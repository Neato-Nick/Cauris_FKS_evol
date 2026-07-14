# nei_gojobori_sites.py

Counts synonymous (S) and nonsynonymous (N) sites per CDS using the unweighted Nei–Gojobori (1986) method. These counts are used to normalize observed polymorphisms and substitutions when computing pN, pS, dN, and dS.

## Method

For each sense codon (stop codons excluded), at each of the three codon positions:

1. Enumerate all 3 possible single-nucleotide changes at that position
2. Count how many produce a synonymous change (same amino acid)
3. Synonymous fraction at that position = n_syn / 3

Synonymous sites per codon (s) = sum of fractions across the three positions (range 0–3).  
Nonsynonymous sites per codon (n) = 3 − s.

Totals S and N for a CDS are the sums across all sense codons. By construction, S + N = n_sense_codons × 3.

Changes to stop codons count as nonsynonymous. Internal and terminal stop codons are identified and excluded from S/N accumulation. Codons containing ambiguous IUPAC bases are also excluded.

Reference: Nei M, Gojobori T (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Mol Biol Evol* 3(5):418–426.

## Usage

```
python3 scripts/nei_gojobori_sites.py <input.cds.fasta> [output.tsv] [--table N] [--per-codon]
```

| Argument | Description |
|---|---|
| `input.cds.fasta` | CDS sequences in FASTA format (nucleotide, in-frame). May include terminal stop codons. |
| `output.tsv` | Output path (optional). Default: `<input_stem>_NG_sites.tsv` or `<input_stem>_NG_per_codon.tsv`. |
| `--table N` | NCBI translation table number. Default: `12` (alternative yeast nuclear). |
| `--per-codon` | Output one row per codon instead of per-gene totals (see [Per-codon mode](#per-codon-mode)). |

### Example

```bash
# Default table 12, auto-named output (per-gene summary)
python3 scripts/nei_gojobori_sites.py Cauris/B11205/Caur_B11205.cds.fasta

# Explicit output path and table
python3 scripts/nei_gojobori_sites.py \
  Cauris/B11205/Caur_B11205.cds.fasta \
  Cauris/B11205/Caur_B11205_NG_sites.tsv \
  --table 12

# Per-codon output (auto-named _NG_per_codon.tsv)
python3 scripts/nei_gojobori_sites.py \
  Cauris/B11205/Caur_B11205.cds.fasta \
  --per-codon
```

Commonly used translation tables:

| Table | Name |
|---|---|
| 1 | Standard |
| 12 | Alternative yeast nuclear (default) |
| 4 | Mold, protozoan, and coelenterate mitochondrial |

See the [NCBI genetic codes page](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for all table numbers.

## Output

### Per-gene summary (default)

Tab-separated file with one row per CDS (CDSs with zero sense codons are omitted).

| Column | Description |
|---|---|
| `gene_id` | FASTA record ID |
| `seq_length_bp` | Full sequence length in base pairs |
| `n_sense_codons` | Number of sense (non-stop, unambiguous) codons used in calculation |
| `n_stop_codons` | Number of in-frame stop codons found (normally 1 for the terminal stop) |
| `n_ambiguous_codons` | Number of codons skipped due to ambiguous IUPAC bases |
| `remainder_bases` | Trailing bases that do not form a complete codon (normally 0) |
| `S` | Total synonymous sites |
| `N` | Total nonsynonymous sites |

### Example output (per-gene)

```
gene_id         seq_length_bp  n_sense_codons  n_stop_codons  n_ambiguous_codons  remainder_bases  S          N
CJI82_03021.1   2403           800             1              0                   0                477.333333 1922.666667
CJI82_03022.1   1584           527             1              0                   0                343.0      1238.0
```

### Per-codon mode

Activated with `--per-codon`. Outputs one row per codon across all CDS in the file. Stop and ambiguous codons are included as rows (with empty `s_i`/`n_i`) so codon numbering stays in sync with other per-site data sources.

| Column | Description |
|---|---|
| `gene_id` | FASTA record ID |
| `codon_number` | 1-based codon index within the CDS |
| `nt_start` | 0-based nucleotide start position of this codon |
| `codon` | The 3-mer (uppercase) |
| `amino_acid` | Translated residue; `*` for stop; `?` for ambiguous |
| `s_i` | Synonymous sites for this codon (empty for stop/ambiguous) |
| `n_i` | Nonsynonymous sites for this codon (empty for stop/ambiguous) |
| `status` | `sense`, `stop`, or `ambiguous` |

By construction, `s_i + n_i = 3.0` for every `sense` row, and `sum(s_i) = S` / `sum(n_i) = N` for each gene matches the per-gene summary.

### Example output (per-codon)

```
gene_id        codon_number  nt_start  codon  amino_acid  s_i       n_i       status
CJI82_03021.1  1             0         ATG    M           0.0       3.0       sense
CJI82_03021.1  2             3         GAC    D           0.333333  2.666667  sense
CJI82_03021.1  3             6         GCC    A           1.0       2.0       sense
```

## Downstream use

S and N are the site-count denominators needed for the Nei–Gojobori estimates of pN/pS (from polymorphism data) and dN/dS (from divergence data):

```
pS = Ps / S      pN = Pn / N
dS = -3/4 * ln(1 - 4/3 * pS)
dN = -3/4 * ln(1 - 4/3 * pN)
```

where Ps and Pn are the counts of synonymous and nonsynonymous polymorphisms (or substitutions) observed for a given gene.

## Dependencies

- Python ≥ 3.10
- [Biopython](https://biopython.org/) (`pip install biopython`)

## Notes

- All 64 unambiguous codon→AA translations are precomputed once from the chosen table at startup; per-codon processing uses only dictionary lookups, so the script runs in under one second for genomes with ~5000 CDS.
- Genes with multiple in-frame stop codons (e.g. pseudogenes or annotation artifacts) are not excluded; their sense codons are still processed and the stop count is reported in the output for downstream filtering.
