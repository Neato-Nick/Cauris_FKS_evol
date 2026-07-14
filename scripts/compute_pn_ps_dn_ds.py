#!/usr/bin/env python3
"""
compute_pn_ps_v2.py

Compute per-site polymorphism (Pn, Ps, Pn_indel, Pn_indel_fs) and optionally
per-site fixed-difference (Dn, Ds, Dn_indel, Dn_indel_fs) counts from a
SnpEff-annotated VCF and a sample-to-clade population map.

Counting modes (--counting-mode):

  classical (default):
    Each (VCF position, variant category) is scored at most once for P and once
    for D, across all clades and all alt alleles.  A nonsyn position segregating
    in three clades still contributes Pn=1.  This is the original McDonald &
    Kreitman (1991) definition of a segregating or fixed site.

    Note on recurrence: this mode does not distinguish a mutation present in one
    clade from the same mutation convergently fixed in five clades.  If you are
    working with a clonal organism or suspect recurrent mutation, analyse mutation
    recurrence separately (e.g. tabulate per-clade allele presence as a binary
    matrix) before interpreting Pn/Ps ratios.

  per-clade (--counting-mode per-clade):
    Each (VCF position, variant category, clade) is scored at most once.  The
    same nonsyn position segregating in three clades contributes Pn=3, reflecting
    three putatively independent mutational events under a clonal model.  This is
    a common adaptation of MK for clonal organisms (e.g. C. auris) where gene
    flow between clades is absent, but it is not the classical MK formulation.
    Multiple distinct alt alleles of the same category at the same position in
    the same clade are still collapsed to a single count.

  per-allele (--counting-mode per-allele):
    Each (VCF position, alt allele, clade) is scored separately.  A multiallelic
    nonsyn site with two distinct alt alleles in one clade contributes Pn=2.
    This reflects the raw number of distinct mutational outcomes but inflates
    counts at hotspot positions and deviates furthest from classical MK.

Polymorphic vs fixed:
  Polymorphic (Pn/Ps):  allele NOT fixed (AF < --fixed-af) in any clade.
  Fixed (Dn/Ds):        allele fixed (AF >= --fixed-af) in at least one clade.
  If fixed in one clade but polymorphic in another, it contributes to Dn
  (for the fixed clade) and is excluded from Pn (matches R cladecons logic).

Variant type classification (by allele length, not SnpEff annotation):
  len(alt) == len(ref)                  -> SNP
    SnpEff impact LOW                   -> Ps / Ds
    SnpEff impact MODERATE or HIGH      -> Pn / Dn
  |len(alt) - len(ref)| % 3 == 0, != 0 -> in-frame indel  (Pn_indel / Dn_indel)
  |len(alt) - len(ref)| % 3 != 0       -> frameshift      (Pn_indel_fs / Dn_indel_fs)

Min-AC filter: applied across ALL clades (so the same allele in 2 isolates
across 2 clades counts as AC=2, not a per-clade singleton).

Usage:
    python scripts/compute_pn_ps.py \\
        --vcf   data/vcfs/Cauris_...snpEff.FKS1.vcf.gz \\
        --popmap data/sample_information/clade_assignments_n22388.tsv \\
        --popmap-clade-col 2 \\
        --out   FKS1_pn_ps.csv \\
        --exclude-clades unknown \\
        [--fixed-af 0.95]          AF threshold to call clade-fixed (default 0.95)
        [--min-ac 1]               min total alt-allele count across all clades (default 1)
        [--no-fixed]               omit Dn/Ds columns
        [--include-invariant]      fill in AA positions with all-zero counts
        [--popmap-sample-col 0]    0-based column for sample ID (default 0)
        [--popmap-clade-col 1]     0-based column for clade label (default 1)
        [--popmap-sep TAB]         column separator in popmap (default: tab)
        [--exclude-clades LIST]    comma-separated clade labels to drop
        [--clade-tree FILE]        Newick tree to correct Dn via Fitch parsimony
        [--clade-tree-outgroup CL] clade name to use as outgroup for rooting
        [--counting-mode MODE]     classical (default) | per-clade | per-allele
                                   (see module docstring for trade-offs)
"""

import argparse
import gzip
import sys
from collections import defaultdict

import pandas as pd


# ── constants ──────────────────────────────────────────────────────────────────

IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

# (kind, category) → output column name
COL_MAP = {
    ("poly",  "syn"):         "Ps",
    ("poly",  "nonsyn"):      "Pn",
    ("poly",  "inframe"):     "Pn_indel",
    ("poly",  "frameshift"):  "Pn_indel_fs",
    ("fixed", "syn"):         "Ds",
    ("fixed", "nonsyn"):      "Dn",
    ("fixed", "inframe"):     "Dn_indel",
    ("fixed", "frameshift"):  "Dn_indel_fs",
}

COUNT_COLS = ["Pn", "Ps", "Pn_indel", "Pn_indel_fs", "Dn", "Ds", "Dn_indel", "Dn_indel_fs"]


# ── helpers ────────────────────────────────────────────────────────────────────

def open_vcf(path):
    """Open a plain or bgzip/gzip VCF for line-by-line reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_ann(ann_str):
    """
    Parse a SnpEff ANN= field value.

    Returns:
        dict: alt_allele_str -> {gene_name: {impact, aa_pos, aa_len}}
        Multiple annotations per (allele, gene) are reduced to the one with
        the highest impact; ties broken by preferring a non-None aa_pos.
    """
    best = defaultdict(dict)  # allele -> gene -> ann_dict

    for entry in ann_str.split(","):
        fields = entry.split("|")
        if len(fields) < 14:
            continue

        allele  = fields[0]
        effect  = fields[1]   # e.g. splice_donor_variant, missense_variant
        impact  = fields[2]
        gene    = fields[4]   # gene ID (e.g. CJI82_02060); field[3] is URL-encoded display name
        hgvsc   = fields[9].strip() if len(fields) > 9 else ""   # e.g. c.53-1G>C

        if not allele or allele == "*" or not gene:
            continue

        aa_pos = aa_len = None

        aa_field = fields[13] if len(fields) > 13 else ""
        if aa_field and "/" in aa_field:
            try:
                a, b   = aa_field.strip().split("/")
                aa_pos = int(a)
                aa_len = int(b)
            except ValueError:
                pass

        cur = best[allele].get(gene)
        if cur is None:
            best[allele][gene] = {"impact": impact, "effect": effect,
                                  "aa_pos": aa_pos, "aa_len": aa_len, "hgvsc": hgvsc}
        elif IMPACT_RANK.get(impact, 0) > IMPACT_RANK.get(cur["impact"], 0):
            best[allele][gene] = {"impact": impact, "effect": effect,
                                  "aa_pos": aa_pos, "aa_len": aa_len, "hgvsc": hgvsc}
        elif (IMPACT_RANK.get(impact, 0) == IMPACT_RANK.get(cur["impact"], 0)
              and aa_pos is not None and cur["aa_pos"] is None):
            cur["aa_pos"] = aa_pos
            cur["aa_len"] = aa_len
            cur["hgvsc"]  = hgvsc

    return dict(best)


def normalize_allele(ref, alt):
    """
    Trim common prefix and suffix to get the minimal allele representation.
    Used as a fallback when the ANN allele doesn't match the full VCF allele.
    """
    # trim common suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # trim common prefix
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    return alt


def allele_type(ref, alt):
    """Return 'snp', 'inframe', or 'frameshift'."""
    diff = abs(len(alt) - len(ref))
    if diff == 0:
        return "snp"
    return "inframe" if diff % 3 == 0 else "frameshift"


def load_popmap(path, sample_col, clade_col, sep, exclude):
    """Return {sample_id: clade} excluding clades in the exclusion set."""
    result  = {}
    exclude = set(exclude)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(sep)
            if max(sample_col, clade_col) >= len(parts):
                continue
            sample = parts[sample_col]
            clade  = parts[clade_col]
            if clade not in exclude:
                result[sample] = clade
    return result


def load_exclude_samples(path):
    """Return a set of sample IDs to exclude (one per line; # lines are comments)."""
    exclude = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                exclude.add(line)
    return exclude


def resolve_gene_symbols(gene_ids, gene_map_path=None, ncbi_lookup=False,
                         ncbi_email=None, ncbi_organism="Candida auris"):
    """
    Return {gene_id: symbol} for each gene_id in gene_ids.
    Resolution order:
      1. User-supplied TSV (gene_id <tab> symbol [<tab> ...])
      2. NCBI Entrez Gene lookup by locus tag (requires BioPython)
    Gene IDs with no resolved symbol are mapped to themselves (pass-through).
    """
    symbol_map = {}

    # 1. User-supplied mapping file
    if gene_map_path:
        with open(gene_map_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    gid, sym = parts[0].strip(), parts[1].strip()
                    if gid and sym:
                        symbol_map[gid] = sym
        print(f"Gene map loaded: {len(symbol_map)} entries from {gene_map_path}", flush=True)

    # 2. NCBI Entrez lookup for gene IDs not resolved by the map file
    unresolved = [g for g in gene_ids if g not in symbol_map]
    if ncbi_lookup and gene_map_path and not unresolved:
        print("Note: all gene IDs resolved via --gene-map; --ncbi-lookup not needed.",
              flush=True)
    if ncbi_lookup and gene_map_path and unresolved:
        print(f"Note: --gene-map resolved {len(symbol_map)} IDs; "
              f"querying NCBI for the remaining {len(unresolved)}: {unresolved}",
              flush=True)
    if ncbi_lookup and unresolved:
        if not ncbi_email:
            print("WARNING: --ncbi-email is required for NCBI lookup; skipping.",
                  file=sys.stderr)
        else:
            try:
                from Bio import Entrez
            except ImportError:
                print("WARNING: BioPython not installed; skipping NCBI lookup.",
                      file=sys.stderr)
                unresolved = []

            if unresolved:
                import re
                import time
                Entrez.email = ncbi_email
                n_before = len(symbol_map)
                print(f"Querying NCBI protein db for {len(unresolved)} locus tags ...",
                      flush=True)
                for gid in unresolved:
                    try:
                        # CJI82_ locus tags are not in the Gene DB; search protein
                        # records by [locus_tag] and parse the /gene= qualifier.
                        handle = Entrez.esearch(db="protein", term=f"{gid}[locus_tag]")
                        record = Entrez.read(handle)
                        handle.close()
                        sym = None
                        if record["IdList"]:
                            pid = record["IdList"][0]
                            handle = Entrez.efetch(
                                db="protein", id=pid, rettype="gp", retmode="text"
                            )
                            text = handle.read()
                            handle.close()
                            for line in text.splitlines():
                                m = re.search(r'/gene="([^"]+)"', line)
                                if m:
                                    sym = m.group(1)
                                    break
                        if sym:
                            symbol_map[gid] = sym
                            print(f"  NCBI: {gid} -> {sym}", flush=True)
                        time.sleep(0.4)   # NCBI: ≤3 requests/sec without API key
                    except Exception as exc:
                        print(f"WARNING: NCBI lookup failed for {gid}: {exc}",
                              file=sys.stderr)
                n_found = len(symbol_map) - n_before
                print(f"NCBI lookup: resolved {n_found}/{len(unresolved)} locus tags",
                      flush=True)

    # Detect duplicate symbols (two gene IDs → same symbol); fall back to ID
    symbol_counts = {}
    for sym in symbol_map.values():
        symbol_counts[sym] = symbol_counts.get(sym, 0) + 1
    for sym, count in symbol_counts.items():
        if count > 1:
            dupes = [gid for gid, s in symbol_map.items() if s == sym]
            print(
                f"WARNING: symbol '{sym}' maps to multiple gene IDs "
                f"({', '.join(sorted(dupes))}); falling back to gene IDs for these genes.",
                file=sys.stderr,
            )
            for gid in dupes:
                del symbol_map[gid]

    return symbol_map


# ── clade-tree helpers ─────────────────────────────────────────────────────────

def build_clade_tree(tree_path, popmap, outgroup_clade=None):
    """
    Load a Newick tree and return a rooted Bio.Phylo tree whose leaf names
    are clade labels.  Accepts two input formats (auto-detected):

    1. Pre-collapsed clade tree  — leaf names are already clade labels
       (e.g. output of collapse_tree_to_clades.py).  Used directly after
       optional re-rooting.

    2. Sample-level tree — leaf names are sample IDs that appear in the
       popmap.  Pruned to one representative per clade and leaves renamed
       to clade labels.

    Detection: if the majority of leaf names match known clade labels from
    the popmap, the tree is treated as pre-collapsed.

    Raises ValueError if neither format can be detected.
    """
    from Bio import Phylo
    import copy

    tree = Phylo.read(tree_path, "newick")

    # Index all leaf names present in the tree
    leaf_names = {t.name for t in tree.get_terminals() if t.name}

    # Clade labels known from the popmap
    known_clades = set(popmap.values())

    # ── auto-detect: pre-collapsed clade tree? ────────────────────────────────
    leaves_matching_clades = leaf_names & known_clades
    if len(leaves_matching_clades) / max(len(leaf_names), 1) > 0.5:
        # Leaf names are already clade labels — use tree as-is
        print(
            f"Clade tree: detected pre-collapsed tree "
            f"({len(leaf_names)} clade leaves: {sorted(leaf_names)})",
            flush=True,
        )
        unrecognised = leaf_names - known_clades
        if unrecognised:
            print(
                f"WARNING: {len(unrecognised)} leaf name(s) not in popmap clades "
                f"and will be ignored in parsimony: {sorted(unrecognised)}",
                file=sys.stderr,
            )
        # Skip to rooting
        clade_to_rep = None   # sentinel: no renaming needed

    else:
        # ── sample-level tree: prune to one rep per clade ─────────────────────
        clade_to_rep = {}
        for sample, clade in sorted(popmap.items()):
            if clade not in clade_to_rep and sample in leaf_names:
                clade_to_rep[clade] = sample

        if not clade_to_rep:
            raise ValueError(
                "No popmap sample IDs match tree leaf names, and leaf names do "
                "not match clade labels either.  Pass either the original "
                "sample-level tree or a pre-collapsed clade tree (e.g. from "
                "scripts/collapse_tree_to_clades.py)."
            )

    if clade_to_rep is not None:
        # Prune tree in a single recursive pass (O(N) vs O(N²) for iterative removal)
        keep         = set(clade_to_rep.values())
        rep_to_clade = {v: k for k, v in clade_to_rep.items()}
        tree         = copy.deepcopy(tree)

        def prune_node(node):
            """Return True if node or any descendant is in `keep`."""
            if node.is_terminal():
                return node.name in keep
            node.clades = [c for c in node.clades if prune_node(c)]
            return len(node.clades) > 0

        prune_node(tree.root)

        # Collapse single-child internal nodes (artifacts of pruning)
        def collapse_singles(node):
            while len(node.clades) == 1:
                child = node.clades[0]
                node.clades = child.clades
                node.branch_length = (node.branch_length or 0) + (child.branch_length or 0)
            for child in node.clades:
                collapse_singles(child)

        collapse_singles(tree.root)

        # Rename leaves: sample ID → clade label
        for leaf in tree.get_terminals():
            leaf.name = rep_to_clade.get(leaf.name, leaf.name)

    # Root the tree
    if outgroup_clade:
        clade_leaves = {t.name for t in tree.get_terminals()}
        if outgroup_clade in clade_leaves:
            tree.root_with_outgroup(outgroup_clade)
        else:
            print(
                f"WARNING: --clade-tree-outgroup '{outgroup_clade}' not found "
                f"in pruned tree ({sorted(clade_leaves)}); falling back to midpoint root.",
                file=sys.stderr,
            )
            tree.root_at_midpoint()
    else:
        tree.root_at_midpoint()

    final_leaves = sorted(t.name for t in tree.get_terminals())
    print(f"Clade tree ({len(final_leaves)} leaves): {final_leaves}", flush=True)

    return tree


def fitch_min_gains(root, fixed_set):
    """
    Minimum number of independent gain events (0→1 transitions) on a rooted tree,
    with state 0 assumed at the root (reference = not fixed).

    Algorithm: two-pass Fitch parsimony biased toward state 0 to minimise
    the number of 0→1 transitions counted.

    Args:
        root:      Bio.Phylo Clade (root node of a pruned clade tree)
        fixed_set: set of leaf names (clade labels) where the allele is fixed

    Returns:
        int ≥ 1: minimum number of independent fixation events
    """
    # Phase 1 — bottom-up: compute Fitch character sets
    def bottom_up(node):
        if node.is_terminal():
            node._fs = frozenset({1} if node.name in fixed_set else {0})
            return
        for child in node.clades:
            bottom_up(child)
        child_sets = [c._fs for c in node.clades]
        inter = child_sets[0]
        for cs in child_sets[1:]:
            inter = inter & cs
        node._fs = inter if inter else frozenset().union(*child_sets)

    bottom_up(root)

    # Phase 2 — top-down: assign states, count 0→1 transitions
    # Prefer state 0 at every ambiguous node to minimise gains.
    gains = [0]

    def top_down(node, parent_state):
        # Inherit parent if possible; otherwise prefer 0 (fewest future gains)
        node_state = parent_state if parent_state in node._fs else min(node._fs)
        if parent_state == 0 and node_state == 1:
            gains[0] += 1
        if not node.is_terminal():
            for child in node.clades:
                top_down(child, node_state)

    root_state = 0 if 0 in root._fs else 1
    if root_state == 1:
        gains[0] += 1   # the root itself is a gain
    for child in root.clades:
        top_down(child, root_state)

    return max(gains[0], 1)  # at least 1 if any clade is fixed


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--vcf",    required=True,
                    help="SnpEff-annotated VCF (plain or .gz / .bgz)")
    ap.add_argument("--popmap", required=True,
                    help="Sample-to-clade map (TSV by default)")
    ap.add_argument("--out",    required=True,
                    help="Output CSV path")
    ap.add_argument("--fixed-af", type=float, default=0.95,
                    help="Clade AF threshold to call an allele fixed (default: 0.95)")
    ap.add_argument("--min-ac",   type=int,   default=1,
                    help="Min total alt-allele count across ALL clades (default: 1)")
    ap.add_argument("--no-fixed", action="store_true",
                    help="Omit Dn/Ds/Dn_indel/Dn_indel_fs columns")
    ap.add_argument("--include-invariant", action="store_true",
                    help="Output rows for every AA position 1..max, filling 0 for invariant sites")
    ap.add_argument("--popmap-sample-col", type=int, default=0,
                    help="0-based column index for sample ID in popmap (default: 0)")
    ap.add_argument("--popmap-clade-col",  type=int, default=1,
                    help="0-based column index for clade label in popmap (default: 1)")
    ap.add_argument("--popmap-sep", default="\t",
                    help="Popmap field separator (default: tab)")
    ap.add_argument("--exclude-clades", default="",
                    help="Comma-separated clade labels to exclude (e.g. 'unknown')")
    ap.add_argument("--exclude-samples", default=None, metavar="FILE",
                    help="File with one sample ID per line to drop before counting "
                         "(equivalent to R omit_misclade_isos)")
    ap.add_argument("--gene-map", default=None, metavar="FILE",
                    help="TSV file mapping gene_id -> gene_symbol (columns: gene_id, symbol). "
                         "Applied after VCF parsing to rename gene IDs in output.")
    ap.add_argument("--ncbi-lookup", action="store_true",
                    help="Query NCBI Entrez Gene by locus tag for any gene IDs not in --gene-map. "
                         "Requires BioPython and --ncbi-email.")
    ap.add_argument("--ncbi-email", default=None, metavar="EMAIL",
                    help="Email address for NCBI Entrez queries (required with --ncbi-lookup)")
    ap.add_argument("--ncbi-organism", default="Candida auris", metavar="ORGANISM",
                    help="Organism name for NCBI locus-tag search (default: 'Candida auris')")
    ap.add_argument("--pn-pops", action="store_true",
                    help="Add a Pn_clades column listing (sorted, comma-separated) the "
                         "clades in which each site has a polymorphic nonsynonymous allele")
    ap.add_argument("--clade-tree", default=None, metavar="FILE",
                    help="Sample-level Newick tree used to correct Dn counting via Fitch "
                         "parsimony. Reduces over-counting when the same allele is fixed "
                         "in multiple closely related clades (shared ancestry vs. parallel "
                         "evolution). Requires BioPython.")
    ap.add_argument("--clade-tree-outgroup", default=None, metavar="CLADE",
                    help="Clade label to use as outgroup for rooting the pruned clade tree "
                         "(e.g. 'Clade_I'). Falls back to midpoint rooting if not specified "
                         "or if the clade is not found in the pruned tree.")
    ap.add_argument("--counting-mode",
                    choices=["classical", "per-clade", "per-allele"],
                    default="classical",
                    help=(
                        "How to count segregating and fixed sites (default: classical). "
                        "classical: each (position, category) scored once across all clades "
                        "and alleles — the original McDonald & Kreitman (1991) definition. "
                        "If the same position is mutated in multiple clades you cannot "
                        "distinguish recurrence from a single shared ancestor under this mode; "
                        "analyse mutation recurrence separately. "
                        "per-clade: each (position, category, clade) scored once; a mutation "
                        "present in N clades contributes N counts, reflecting N independent "
                        "events under a strict clonal model. "
                        "per-allele: each (position, alt allele, clade) scored separately; "
                        "inflates counts at multiallelic/hotspot positions and is furthest "
                        "from classical MK."
                    ))
    args = ap.parse_args()

    # handle literal \t from shell
    if args.popmap_sep in (r"\t", "\\t"):
        args.popmap_sep = "\t"

    exclude_clades = [c for c in args.exclude_clades.split(",") if c]

    # ── load sample exclusion list ─────────────────────────────────────────────
    exclude_sample_ids = set()
    if args.exclude_samples:
        exclude_sample_ids = load_exclude_samples(args.exclude_samples)
        print(f"Excluding {len(exclude_sample_ids)} samples from --exclude-samples",
              flush=True)

    # ── load popmap ────────────────────────────────────────────────────────────
    popmap = load_popmap(
        args.popmap,
        args.popmap_sample_col,
        args.popmap_clade_col,
        args.popmap_sep,
        exclude_clades,
    )
    all_clades = sorted(set(popmap.values()))
    print(f"Popmap: {len(popmap)} samples across {len(all_clades)} clades: {all_clades}",
          flush=True)

    # ── parse VCF ──────────────────────────────────────────────────────────────
    events       = []   # list of dicts, one per (gene, chrom, pos, alt, clade)
    gene_aa_len  = {}   # gene -> max aa_len seen (for --include-invariant)
    vcf_samples  = []
    sample_idx_to_clade = {}   # column offset among GT cols -> clade

    with open_vcf(args.vcf) as fh:
        for line in fh:

            # ── meta-header ──────────────────────────────────────────────────
            if line.startswith("##"):
                continue

            # ── column header ────────────────────────────────────────────────
            if line.startswith("#CHROM"):
                vcf_samples = line.rstrip("\n").split("\t")[9:]
                sample_idx_to_clade = {
                    i: popmap[s]
                    for i, s in enumerate(vcf_samples)
                    if s in popmap and s not in exclude_sample_ids
                }
                print(
                    f"VCF: {len(vcf_samples)} samples total, "
                    f"{len(sample_idx_to_clade)} matched to popmap",
                    flush=True,
                )
                continue

            # ── data line ────────────────────────────────────────────────────
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom = cols[0]
            pos   = int(cols[1])
            ref   = cols[3]
            alt_s = cols[4]

            if alt_s in (".", "*"):
                continue

            alts = alt_s.split(",")

            # FORMAT: find GT field index
            fmt   = cols[8].split(":") if len(cols) > 8 else ["GT"]
            gt_i  = fmt.index("GT") if "GT" in fmt else 0

            # Per-alt counts:
            # clade_ac[alt_1based][clade] = # samples in clade carrying that allele
            # clade_n[clade]             = # samples in clade with a called genotype
            clade_ac = defaultdict(lambda: defaultdict(int))
            clade_n  = defaultdict(int)

            for col_i, clade in sample_idx_to_clade.items():
                raw = cols[9 + col_i] if 9 + col_i < len(cols) else "."
                gt_str = (raw.split(":")[gt_i] if ":" in raw else raw)
                gt_str = gt_str.replace("|", "/")
                calls  = gt_str.split("/")

                # Skip fully missing
                if all(c == "." for c in calls):
                    continue

                clade_n[clade] += 1

                # Count each unique alt allele once per sample (not per copy)
                seen = set()
                for c in calls:
                    if c in (".", ""):
                        continue
                    try:
                        idx = int(c)
                    except ValueError:
                        continue
                    if idx > 0 and idx not in seen:
                        clade_ac[idx][clade] += 1
                        seen.add(idx)

            # ── ANN field ───────────────────────────────────────────────────
            ann_str = ""
            for token in cols[7].split(";"):
                if token.startswith("ANN="):
                    ann_str = token[4:]
                    break
            if not ann_str:
                continue

            ann_data = parse_ann(ann_str)   # allele -> gene -> {impact, aa_pos, aa_len}

            # ── per alt allele ───────────────────────────────────────────────
            for ai0, alt in enumerate(alts):
                ai1 = ai0 + 1
                if alt in (".", "*"):
                    continue

                atype = allele_type(ref, alt)

                ac_by_clade = dict(clade_ac[ai1])   # {clade: count}
                total_ac    = sum(ac_by_clade.values())

                if total_ac < args.min_ac:
                    continue

                # Annotation lookup: try exact match, then normalized form
                allele_anns = ann_data.get(alt)
                if allele_anns is None:
                    norm = normalize_allele(ref, alt)
                    allele_anns = ann_data.get(norm)
                if allele_anns is None:
                    continue

                for gene, ann in allele_anns.items():
                    impact  = ann["impact"]
                    effect  = ann.get("effect", "")
                    aa_pos  = ann["aa_pos"]
                    aa_len  = ann["aa_len"]
                    hgvsc   = ann.get("hgvsc", "")

                    # Skip non-coding (MODIFIER)
                    if impact == "MODIFIER":
                        continue

                    # Splice-site and other HIGH/MODERATE variants without an AA
                    # position (e.g. splice_donor, splice_acceptor) can't be
                    # attributed to a specific codon.  Assign aa_pos = 0 as a
                    # sentinel so they appear as a distinct row in the output.
                    # Downstream: filter aa_pos > 0 to exclude splice-site rows.
                    if aa_pos is None:
                        if impact in ("HIGH", "MODERATE"):
                            aa_pos = 0
                            cds_str = hgvsc if hgvsc else "CDS position unknown"
                            print(
                                f"WARNING: {gene} {chrom}:{pos} "
                                f"({effect}, {impact}) has no AA position; "
                                f"tabulated at aa_pos=0 (CDS: {cds_str})",
                                file=sys.stderr,
                            )
                        else:
                            continue

                    # Variant category
                    if atype == "snp":
                        if impact == "LOW":
                            category = "syn"
                        elif impact in ("MODERATE", "HIGH"):
                            category = "nonsyn"
                        else:
                            continue
                    else:   # indel
                        if impact not in ("MODERATE", "HIGH"):
                            continue
                        category = atype   # "inframe" or "frameshift"

                    # Track protein length for --include-invariant
                    if aa_len:
                        gene_aa_len[gene] = max(gene_aa_len.get(gene, 0), aa_len)

                    # Is this allele fixed in ANY clade?
                    fixed_in_any = any(
                        ac_by_clade.get(cl, 0) / clade_n[cl] >= args.fixed_af
                        for cl in clade_n
                        if clade_n[cl] > 0 and ac_by_clade.get(cl, 0) > 0
                    )

                    # Emit one event per clade that carries this allele
                    for clade, count in ac_by_clade.items():
                        if count == 0:
                            continue
                        n = clade_n.get(clade, 0)
                        if n == 0:
                            continue

                        clade_af   = count / n
                        fixed_here = clade_af >= args.fixed_af

                        if fixed_here:
                            kind = "fixed"
                        elif not fixed_in_any:
                            kind = "poly"
                        else:
                            # Fixed elsewhere, polymorphic here → exclude from Pn
                            # (mirrors R behaviour: cladecons sites removed from Pn)
                            continue

                        events.append({
                            "gene":     gene,
                            "chrom":    chrom,
                            "pos":      pos,
                            "alt":      alt,
                            "clade":    clade,
                            "category": category,
                            "kind":     kind,
                            "aa_pos":   aa_pos,
                        })

    if not events:
        sys.exit("ERROR: No usable annotated variants found. "
                 "Check that the VCF has an ANN= field and that sample IDs "
                 "match the popmap.")

    print(f"Collected {len(events)} allele×clade observations", flush=True)

    # ── Fitch parsimony correction for Dn ─────────────────────────────────────
    # Without a tree, each (pos, alt, clade) with AF≥threshold emits one Dn
    # event regardless of clade relationships — identical to counting once per
    # clade.  The same allele fixed in all 6 clades gives Dn=6 when it is
    # almost certainly one ancestral fixation.  With --clade-tree we replace
    # the per-clade raw counts with the minimum number of independent gain
    # events inferred by Fitch parsimony on the pruned clade topology.
    if args.clade_tree:
        clade_tree_obj = None
        try:
            clade_tree_obj = build_clade_tree(
                args.clade_tree, popmap,
                outgroup_clade=args.clade_tree_outgroup,
            )
        except Exception as exc:
            print(
                f"WARNING: Could not load clade tree ({exc}); "
                f"Dn will be counted per clade (uncorrected).",
                file=sys.stderr,
            )

        if clade_tree_obj is not None:
            # Group fixed events by (gene, chrom, pos, alt, aa_pos, category)
            fixed_groups = defaultdict(set)
            for e in events:
                if e["kind"] == "fixed":
                    key = (e["gene"], e["chrom"], e["pos"], e["alt"],
                           e["aa_pos"], e["category"])
                    fixed_groups[key].add(e["clade"])

            # Warn about clades present in events but absent from the clade tree
            tree_clades = {t.name for t in clade_tree_obj.get_terminals()}
            event_clades = {e["clade"] for e in events if e["kind"] == "fixed"}
            missing = event_clades - tree_clades
            if missing:
                print(
                    f"WARNING: {len(missing)} clade(s) in fixed events not found "
                    f"in pruned clade tree: {sorted(missing)}. "
                    f"These clades will be treated as state 0 (not fixed) in "
                    f"parsimony — Dn may be underestimated.",
                    file=sys.stderr,
                )

            # Replace fixed events with parsimony-corrected ones
            events = [e for e in events if e["kind"] != "fixed"]
            total_raw = sum(len(v) for v in fixed_groups.values())
            total_parsimony = 0

            clade_root = clade_tree_obj.root
            for (gene, chrom, pos, alt, aa_pos, category), fixed_clades in fixed_groups.items():
                n = fitch_min_gains(clade_root, fixed_clades)
                total_parsimony += n
                for i in range(n):
                    events.append({
                        "gene":     gene,
                        "chrom":    chrom,
                        "pos":      pos,
                        "alt":      alt,
                        "clade":    f"_gain_{i}",   # synthetic; unique per gain
                        "category": category,
                        "kind":     "fixed",
                        "aa_pos":   aa_pos,
                    })

            print(
                f"Parsimony Dn correction: {total_raw} raw events → "
                f"{total_parsimony} minimum gains "
                f"(reduction factor {total_raw/max(total_parsimony,1):.1f}×)",
                flush=True,
            )

    df = pd.DataFrame(events)

    # ── resolve gene IDs to symbols ────────────────────────────────────────────
    if args.gene_map or args.ncbi_lookup:
        gene_ids = sorted(df["gene"].unique())
        print(f"Resolving gene symbols for: {gene_ids}", flush=True)
        symbol_map = resolve_gene_symbols(
            gene_ids,
            gene_map_path  = args.gene_map,
            ncbi_lookup    = args.ncbi_lookup,
            ncbi_email     = args.ncbi_email,
            ncbi_organism  = args.ncbi_organism,
        )
        if symbol_map:
            df["gene"] = df["gene"].map(lambda g: symbol_map.get(g, g))
            gene_aa_len = {symbol_map.get(g, g): v for g, v in gene_aa_len.items()}

    # Deduplicate events according to --counting-mode.
    #
    # classical:   each (gene, aa_pos, chrom, pos, category, kind) scored once —
    #              clade and alt allele identity are ignored.  One segregating
    #              nonsynonymous site = Pn+1 regardless of how many clades carry it
    #              or how many alt alleles exist there.  This is the original MK
    #              definition.
    #
    # per-clade:   each (gene, aa_pos, chrom, pos, category, kind, clade) scored
    #              once.  The same position mutated in N clades = N counts.
    #              Multiple alt alleles of the same category in the same clade still
    #              collapse to one count.
    #
    # per-allele:  each (gene, aa_pos, chrom, pos, alt, clade) scored separately.
    #              Two nonsynonymous alleles at the same site in the same clade = 2.
    if args.counting_mode == "per-allele":
        dedup_subset = ["gene", "aa_pos", "chrom", "pos", "alt", "clade"]
    elif args.counting_mode == "per-clade":
        dedup_subset = ["gene", "aa_pos", "chrom", "pos", "category", "kind", "clade"]
    else:  # classical
        dedup_subset = ["gene", "aa_pos", "chrom", "pos", "category", "kind"]
    print(f"Counting mode: {args.counting_mode} (dedup key: {dedup_subset})", flush=True)
    df = df.drop_duplicates(subset=dedup_subset)
    print(f"After deduplication: {len(df)} distinct events", flush=True)

    # Map to output column
    df["col"] = df.apply(lambda r: COL_MAP.get((r["kind"], r["category"])), axis=1)
    df = df[df["col"].notna()].copy()

    # Aggregate: count events per (gene, aa_pos, output_col)
    if df.empty:
        result = pd.DataFrame(columns=["gene", "aa_pos"] + COUNT_COLS)
    else:
        pivoted = (
            df.groupby(["gene", "aa_pos", "col"])
            .size()
            .reset_index(name="n")
            .pivot_table(
                index=["gene", "aa_pos"],
                columns="col",
                values="n",
                fill_value=0,
            )
            .reset_index()
        )
        pivoted.columns.name = None
        # Ensure all columns exist
        for col in COUNT_COLS:
            if col not in pivoted.columns:
                pivoted[col] = 0
        result = pivoted

    # ── include invariant positions ────────────────────────────────────────────
    if args.include_invariant:
        existing = set(zip(result["gene"], result["aa_pos"]))
        inv_rows = []
        for gene, max_len in gene_aa_len.items():
            for ap in range(1, max_len + 1):
                if (gene, ap) not in existing:
                    row = {"gene": gene, "aa_pos": ap}
                    for col in COUNT_COLS:
                        row[col] = 0
                    inv_rows.append(row)
        if inv_rows:
            result = pd.concat([result, pd.DataFrame(inv_rows)], ignore_index=True)
            result = result.sort_values(["gene", "aa_pos"]).reset_index(drop=True)

    # ── Pn_clades: populations carrying each polymorphic nonsyn allele ────────
    if args.pn_pops:
        pn_clades = (
            df[(df["kind"] == "poly") & (df["category"] == "nonsyn")]
            .groupby(["gene", "aa_pos"])["clade"]
            .apply(lambda s: ",".join(sorted(s.unique())))
            .reset_index(name="Pn_clades")
        )
        result = result.merge(pn_clades, on=["gene", "aa_pos"], how="left")

    # ── select and order output columns ───────────────────────────────────────
    out_cols = ["gene", "aa_pos", "Pn", "Ps", "Pn_indel", "Pn_indel_fs"]
    if not args.no_fixed:
        out_cols += ["Dn", "Ds", "Dn_indel", "Dn_indel_fs"]
    if args.pn_pops:
        out_cols += ["Pn_clades"]

    result = result[[c for c in out_cols if c in result.columns]]
    result["aa_pos"] = result["aa_pos"].astype(int)
    result = result.sort_values(["gene", "aa_pos"]).reset_index(drop=True)

    result.to_csv(args.out, index=False)
    print(f"\nWrote {len(result)} sites to {args.out}", flush=True)

    # ── summary ───────────────────────────────────────────────────────────────
    print("\nSummary per gene:")
    for gene in sorted(result["gene"].unique()):
        g = result[result["gene"] == gene]
        msg = (
            f"  {gene}: {len(g)} variant sites, "
            f"Pn={g['Pn'].sum()}  Ps={g['Ps'].sum()}  "
            f"Pn_indel={g['Pn_indel'].sum()}  Pn_indel_fs={g['Pn_indel_fs'].sum()}"
        )
        if not args.no_fixed:
            msg += (
                f"\n         "
                f"Dn={g['Dn'].sum()}  Ds={g['Ds'].sum()}  "
                f"Dn_indel={g['Dn_indel'].sum()}  Dn_indel_fs={g['Dn_indel_fs'].sum()}"
            )
        print(msg, flush=True)


if __name__ == "__main__":
    main()
