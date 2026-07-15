#!/usr/bin/env python3
"""
Compute synonymous (S) and nonsynonymous (N) site counts for every CDS
in a FASTA file using the Nei-Gojobori (1986) method, optionally with
transition-transversion rate ratio correction (kappa).

For each sense codon at each of the three positions:
  - enumerate all 3 possible single-nucleotide changes
  - weight each change by kappa (if a transition) or 1 (if a transversion)
  - synonymous fraction at that position = weighted_syn / (kappa + 2)
Synonymous sites per codon = sum of fractions across positions (range 0-3).
Nonsynonymous sites per codon = 3 - S_codon.

When kappa=1 (default) every change carries equal weight 1/3, which is
identical to the original unweighted Nei-Gojobori method.

Transitions:   A<->G, C<->T  (same biochemical class: purine<->purine or
               pyrimidine<->pyrimidine), weighted by kappa.
Transversions: A<->C, A<->T, G<->C, G<->T, weighted by 1.

Usage:
  python nei_gojobori_sites.py <input.cds.fasta> [output.tsv] \
      [--table N] [--kappa K] [--per-codon]

Arguments:
  input.cds.fasta   CDS sequences (nucleotide, in-frame, may include terminal stop)
  output.tsv        Output TSV (default: <input_stem>_NG_sites.tsv or _NG_per_codon.tsv)
  --table N         NCBI translation table number (default: 12, alternative yeast nuclear)
  --kappa K         Transition/transversion rate ratio (default: 1.0, unweighted NG)
  --per-codon       Output one row per codon instead of per-gene totals
"""

import sys
import argparse
from itertools import product
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

NUCLEOTIDES = ('A', 'T', 'C', 'G')
_STOP = '*'
_PURINES = frozenset('AG')
_PYRIMIDINES = frozenset('CT')


def is_transition(nuc1: str, nuc2: str) -> bool:
    """Return True if nuc1->nuc2 is a transition (purine<->purine or pyrimidine<->pyrimidine)."""
    return (nuc1 in _PURINES and nuc2 in _PURINES) or (nuc1 in _PYRIMIDINES and nuc2 in _PYRIMIDINES)


def build_codon_table(table_id: int) -> dict[str, str]:
    """Return a dict mapping every 3-mer to its amino acid (or '*' for stop)."""
    ct = CodonTable.unambiguous_dna_by_id[table_id]
    codon_to_aa: dict[str, str] = {}
    for nuc_tuple in product(NUCLEOTIDES, repeat=3):
        codon = ''.join(nuc_tuple)
        if codon in ct.stop_codons:
            codon_to_aa[codon] = _STOP
        else:
            codon_to_aa[codon] = ct.forward_table.get(codon, '?')
    return codon_to_aa


def build_sites_table(
    codon_to_aa: dict[str, str],
    kappa: float = 1.0,
) -> dict[str, tuple[float, float]]:
    """
    Precompute (syn_sites, nonsyn_sites) for every unambiguous sense codon.

    Each of the 3 alternative nucleotides at each codon position is weighted
    by kappa if it constitutes a transition, or 1 if a transversion.  The
    synonymous fraction at each position is the sum of synonymous-change
    weights divided by the total weight (kappa + 2).

    When kappa=1 every change carries weight 1/3, reproducing the original
    unweighted Nei-Gojobori formula.  Only sense codons get entries.
    """
    total_weight = kappa + 2.0  # 1 transition + 2 transversions, unnormalised
    sites: dict[str, tuple[float, float]] = {}
    for codon, ref_aa in codon_to_aa.items():
        if ref_aa in (_STOP, '?'):
            continue
        syn_sum = 0.0
        for pos in range(3):
            syn_weight = 0.0
            for nuc in NUCLEOTIDES:
                if nuc == codon[pos]:
                    continue
                mutated = codon[:pos] + nuc + codon[pos + 1:]
                alt_aa = codon_to_aa.get(mutated, '?')
                if alt_aa == ref_aa:
                    syn_weight += kappa if is_transition(codon[pos], nuc) else 1.0
            syn_sum += syn_weight / total_weight
        sites[codon] = (syn_sum, 3.0 - syn_sum)
    return sites


def iter_cds_per_codon(
    seq: str,
    codon_to_aa: dict[str, str],
    sites_table: dict[str, tuple[float, float]],
):
    """
    Yield one dict per codon in seq.

    Stop and ambiguous codons are yielded with empty s_i/n_i so that
    codon numbering stays in sync with other per-site data sources.
    """
    seq = seq.upper()
    for i in range(len(seq) // 3):
        codon = seq[i * 3 : i * 3 + 3]
        codon_num = i + 1   # 1-based
        nt_start = i * 3    # 0-based

        if codon not in codon_to_aa:
            yield {
                "codon_number": codon_num, "nt_start": nt_start,
                "codon": codon, "amino_acid": "?",
                "s_i": "", "n_i": "", "status": "ambiguous",
            }
            continue

        aa = codon_to_aa[codon]
        if aa == _STOP:
            yield {
                "codon_number": codon_num, "nt_start": nt_start,
                "codon": codon, "amino_acid": "*",
                "s_i": "", "n_i": "", "status": "stop",
            }
            continue

        s, n = sites_table[codon]
        yield {
            "codon_number": codon_num, "nt_start": nt_start,
            "codon": codon, "amino_acid": aa,
            "s_i": round(s, 6), "n_i": round(n, 6), "status": "sense",
        }


def compute_cds_sites(
    seq: str,
    codon_to_aa: dict[str, str],
    sites_table: dict[str, tuple[float, float]],
) -> dict:
    """
    Compute S and N for a complete CDS string.

    Terminal and any internal stop codons are detected and excluded.
    Codons with ambiguous bases are tallied but skipped.
    """
    seq = seq.upper()
    n_total_codons = len(seq) // 3
    remainder = len(seq) % 3

    S = 0.0
    N = 0.0
    n_sense = 0
    n_stop = 0
    n_ambiguous = 0

    for i in range(n_total_codons):
        codon = seq[i * 3 : i * 3 + 3]

        if codon not in codon_to_aa:
            # Contains ambiguous nucleotide(s)
            n_ambiguous += 1
            continue

        aa = codon_to_aa[codon]

        if aa == _STOP:
            n_stop += 1
            continue

        s, n = sites_table[codon]
        S += s
        N += n
        n_sense += 1

    return {
        "n_sense_codons": n_sense,
        "n_stop_codons": n_stop,
        "n_ambiguous_codons": n_ambiguous,
        "remainder_bases": remainder,
        "S": round(S, 6),
        "N": round(N, 6),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Nei-Gojobori (1986) synonymous/nonsynonymous site counts per CDS"
    )
    p.add_argument("input", help="CDS FASTA file (nucleotide, in-frame)")
    p.add_argument(
        "output",
        nargs="?",
        default=None,
        help="Output TSV (default: <input_stem>_NG_sites.tsv)",
    )
    p.add_argument(
        "--table",
        type=int,
        default=12,
        help="NCBI translation table number (default: 12, alternative yeast nuclear)",
    )
    p.add_argument(
        "--kappa",
        type=float,
        default=1.0,
        help=(
            "Transition/transversion rate ratio κ for site-count weighting "
            "(default: 1.0 = unweighted Nei-Gojobori).  Transitions are "
            "weighted by κ, transversions by 1; synonymous fraction at each "
            "codon position = weighted_syn / (κ + 2).  Typical coding-sequence "
            "values are 2-4; set to 1 to reproduce the original NG method."
        ),
    )
    p.add_argument(
        "--per-codon",
        action="store_true",
        help="Output one row per codon (sense/stop/ambiguous) instead of per-gene totals",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Error: input file not found: {input_path}")

    if args.kappa <= 0:
        sys.exit("Error: --kappa must be a positive number")

    if args.output:
        output_path = Path(args.output)
    else:
        suffix = "_NG_per_codon.tsv" if args.per_codon else "_NG_sites.tsv"
        output_path = input_path.parent / (input_path.stem + suffix)

    kappa_label = (
        "1.0 (unweighted, original Nei-Gojobori)"
        if args.kappa == 1.0
        else f"{args.kappa} (Ts/Tv-corrected)"
    )
    print(f"Translation table : {args.table}", file=sys.stderr)
    print(f"Kappa (Ts/Tv)     : {kappa_label}", file=sys.stderr)
    print(f"Input             : {input_path}", file=sys.stderr)
    print(f"Output            : {output_path}", file=sys.stderr)

    codon_to_aa = build_codon_table(args.table)
    sites_table = build_sites_table(codon_to_aa, kappa=args.kappa)

    n_written = 0
    n_skipped = 0

    if args.per_codon:
        header = (
            "gene_id\tcodon_number\tnt_start\tcodon\tamino_acid\ts_i\tn_i\tstatus\n"
        )
        with open(output_path, "w") as out:
            out.write(header)
            for record in SeqIO.parse(input_path, "fasta"):
                seq_str = str(record.seq)
                if len(seq_str) < 3:
                    n_skipped += 1
                    continue
                for row in iter_cds_per_codon(seq_str, codon_to_aa, sites_table):
                    out.write(
                        f"{record.id}\t"
                        f"{row['codon_number']}\t"
                        f"{row['nt_start']}\t"
                        f"{row['codon']}\t"
                        f"{row['amino_acid']}\t"
                        f"{row['s_i']}\t"
                        f"{row['n_i']}\t"
                        f"{row['status']}\n"
                    )
                n_written += 1
        print(f"Done: {n_written} CDS written, {n_skipped} skipped.", file=sys.stderr)

    else:
        header = (
            "gene_id\t"
            "seq_length_bp\t"
            "n_sense_codons\t"
            "n_stop_codons\t"
            "n_ambiguous_codons\t"
            "remainder_bases\t"
            "S\t"
            "N\n"
        )
        with open(output_path, "w") as out:
            out.write(header)
            for record in SeqIO.parse(input_path, "fasta"):
                seq_str = str(record.seq)
                if len(seq_str) < 3:
                    n_skipped += 1
                    continue
                stats = compute_cds_sites(seq_str, codon_to_aa, sites_table)
                if stats["n_sense_codons"] == 0:
                    n_skipped += 1
                    continue
                out.write(
                    f"{record.id}\t"
                    f"{len(seq_str)}\t"
                    f"{stats['n_sense_codons']}\t"
                    f"{stats['n_stop_codons']}\t"
                    f"{stats['n_ambiguous_codons']}\t"
                    f"{stats['remainder_bases']}\t"
                    f"{stats['S']}\t"
                    f"{stats['N']}\n"
                )
                n_written += 1
        print(f"Done: {n_written} CDS written, {n_skipped} skipped.", file=sys.stderr)


if __name__ == "__main__":
    main()
