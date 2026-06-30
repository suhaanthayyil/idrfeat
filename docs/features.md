# Feature list

This is the full list of features `idrfeat` computes, one row per IDR segment. For each feature
it gives the definition, the tool or database it comes from, and one line on why it might matter
for whether an IDR is fitness associated or whether a substitution perturbs it. Where the
predictive value is a guess rather than something we have measured, it says so. Nothing here has
been fit to the melanoma fitness labels yet, so treat every "why it matters" as a hypothesis to
test in phase 1, not a result.

## How IDRs are annotated and segments called

A protein sequence is scored per residue for disorder, then contiguous runs at or above a score
threshold and at least a minimum length become IDR segments. Two disorder definitions are
supported so a region can be annotated more than one way:

- **metapredict** (primary). A fast, well calibrated per-residue disorder predictor. It is the
  default segment caller.
- **AIUPred** (secondary). An energy-based per-residue disorder predictor. It also gives the
  fold-on-binding propensity below.

Defaults are a score threshold of 0.5 and a minimum segment length of 10, both configurable.
When metapredict is not installed the segments are called with a transparent built-in heuristic
(hydropathy plus charge plus disorder-promoting composition, window smoothed) so the pipeline
still runs offline. The heuristic is a smoke-test stand-in, not a replacement for metapredict.
Every row records which backend defined it in `disorder_backend`.

Coordinates `idr_start` and `idr_end` are 1-based inclusive into the parent protein. `sequence`
is the segment amino acid string. `seg_length` is its length.

## Missing-value policy

Sequence-based features are always populated. A predictor or database that is not supplied leaves
its columns as missing values (NaN) rather than zero, so a blank means "source absent", not "none
found". A supplied database that simply has no hit for a segment gives a real zero count.

## Composition

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `aa_frac_A` ... `aa_frac_Y` | Fraction of each of the 20 standard amino acids | Sequence | Composition is the baseline description of an IDR and the input to most other features |
| `frac_charged` | Fraction of D, E, K, R | Sequence | Charge density separates linkers from binding and LLPS regions |
| `frac_polar` | Fraction of S, T, N, Q, C, Y, H | Sequence | Polar residues mark solvated, flexible stretches |
| `frac_aromatic` | Fraction of F, W, Y | Sequence | Aromatics drive many binding and phase-separation contacts |
| `frac_proline` | Fraction of P | Sequence | Proline enforces extended, disorder-promoting backbone |
| `frac_glycine` | Fraction of G | Sequence | Glycine gives flexibility and is enriched in spacers |
| `frac_disorder_promoting` | Fraction of A, R, G, Q, S, E, K, P (the positive end of the TOP-IDP scale) | Sequence, Campen et al. 2008 | Directly measures how disorder-promoting the composition is |

## Charge and patterning

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `fcr` | Fraction of charged residues | Sequence | Total charge sets the conformational class of an IDR |
| `ncpr` | Net charge per residue, signed | Sequence | Net charge shifts compaction and binding partners |
| `net_charge` | Count of positive minus negative residues | Sequence | Absolute charge balance |
| `kappa` | Das and Pappu (2013) charge patterning, 0 well mixed to 1 segregated, -1 when undefined | Sequence (native implementation, cross-checks localCIDER) | Charge blockiness changes IDR dimensions and binding even at fixed composition |
| `scd` | Sequence charge decoration, Sawle and Ghosh (2015) | Sequence | A second charge-patterning measure, more negative as opposite charges separate |
| `frac_positive` | Fraction of K, R | Sequence | Positive tracts bind nucleic acids and acidic partners |
| `frac_negative` | Fraction of D, E | Sequence | Acidic tracts are common in transactivation domains |

## Hydropathy

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `hydropathy_mean` | Mean normalized Kyte-Doolittle hydropathy, 0 to 1 | Sequence | Hydrophobicity gates folding and aggregation |
| `hydropathy_patterning` | Blockiness of hydropathy, an SCD-style decoration on mean-centered hydropathy | Sequence | Patchy hydrophobicity can seed folding-on-binding and phase separation |

## Sequence complexity

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `kmer1_distinct`, `kmer2_distinct`, `kmer3_distinct` | Number of distinct k-mers for k of 1, 2, 3 | Sequence | Low k-mer diversity flags repetitive, low-complexity sequence |
| `kmer1_entropy`, `kmer2_entropy`, `kmer3_entropy` | Shannon entropy in bits of the k-mer distribution | Sequence | Entropy summarizes k-mer content in one number per k |
| `kmer1_diversity`, `kmer2_diversity`, `kmer3_diversity` | Distinct k-mers divided by the number of k-mer positions | Sequence | Length-normalized richness, comparable across segment sizes |
| `lowcomplexity_frac` | Fraction of residues in a window whose residue entropy is below threshold | Sequence | Low-complexity regions are common in fitness-associated and LLPS IDRs |
| `longest_single_run` | Length of the longest run of one residue | Sequence | Homopolymer runs (poly-Q, poly-E) are functionally important |

The full k-mer count vectors for k of 2 and 3 are available from `idrfeat.complexity.kmer_spectrum`
for downstream models. The table keeps only the richness and entropy summaries so it stays tidy.

## Disorder and structure

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `disorder_mean`, `disorder_max`, `disorder_frac` | Mean, max, and fraction at or above threshold of the primary backend over the segment | metapredict or heuristic | How confidently the segment is disordered |
| `metapredict_mean`, `metapredict_max`, `metapredict_frac` | Same three statistics from metapredict | metapredict (NaN if not installed) | The reference disorder signal |
| `aiupred_disorder_mean` | Mean AIUPred disorder over the segment | AIUPred (NaN if not installed) | A second, independent disorder call |
| `plddt_mean` | Mean AlphaFold pLDDT over the segment | AlphaFold PDB B-factors (NaN if no model given) | Low pLDDT is a structure-based disorder check; high pLDDT can flag conditional folding |

## Fold on binding

This answers Carter's first open question, how to quantify propensity to fold on binding.

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `aiupred_binding_mean` | Mean AIUPred binding propensity over the segment | AIUPred (NaN if not installed) | High binding propensity marks residues that fold when they bind a partner |
| `aiupred_binding_max` | Max AIUPred binding propensity | AIUPred | Peak fold-on-binding signal in the segment |
| `aiupred_binding_frac` | Fraction of residues at or above the binding threshold | AIUPred | How much of the segment is fold-on-binding |

AIUPred binding is the ANCHOR2-style score for disordered binding regions. It is the most direct
sequence answer to fold on binding. The pLDDT and ELM binding-motif features below give two
independent, partly orthogonal views of the same idea. We have not yet checked which is most
predictive of fitness.

## Motifs and annotation

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `elm_motif_residues` | Residues in the segment covered by any ELM motif instance | ELM (NaN if not given) | Many fitness-associated IDRs carry short linear motifs |
| `elm_motif_density` | `elm_motif_residues` divided by segment length | ELM | Length-normalized motif load |
| `elm_classes_distinct` | Number of distinct ELM classes overlapping the segment | ELM | Motif diversity, a proxy for how many partners a region engages |
| `overlaps_binding_motif` | Whether the segment overlaps a LIG or DOC motif | ELM | Binding and docking motifs are a direct functional handle |
| `ptm_count`, `ptm_density` | Total dbPTM sites in the segment and per-residue density | dbPTM (NaN if not given) | PTM sites are regulatory and often under selection |
| `ptm_phospho_count`, `ptm_phospho_density` | Phosphorylation sites and density | dbPTM | Phosphosites are the most common regulatory PTM in IDRs |
| `ptm_acetyl_count`, `ptm_acetyl_density` | Acetylation sites and density | dbPTM | Acetylation tunes charge and interactions |
| `ptm_ubiq_count`, `ptm_ubiq_density` | Ubiquitination sites and density | dbPTM | Ubiquitination controls turnover, relevant to fitness |
| `ptm_methyl_count`, `ptm_methyl_density` | Methylation sites and density | dbPTM | Methylation modulates binding, notably in RG regions |

## Phase separation

This answers Carter's second open question, how to quantify phase-separation propensity. The short
version is that there is no single trusted number, so `idrfeat` gives three complementary signals.

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `phasepro_overlap` | Whether the segment overlaps a curated PhaSePro LLPS region | PhaSePro (NaN if not given) | Curated, high-confidence LLPS evidence |
| `drllps_member` | Whether the protein is a DrLLPS scaffold or regulator | DrLLPS (NaN if not given) | Broader, lower-confidence LLPS membership |
| `llps_seq_score` | Mean of six LLPS-promoting signals: aromatic, arginine, glycine, glutamine plus asparagine, charge, and low-complexity fractions | Sequence | A transparent sequence guess at LLPS propensity, to be validated, not trusted yet |
| `charge_aromatic_proxy` | Fraction charged times fraction aromatic | Sequence | The charge-times-aromatic proxy Carter suggested, capturing the two main LLPS driving forces |

On the phase-separation question directly: charge and hydrophobicity are part of it, as Carter
guessed, but the field now reads the main sequence grammar as multivalent aromatics (especially
tyrosine), arginine, and glycine or serine spacers, plus charge blockiness for complex
coacervation. `llps_seq_score` folds those compositional signals into one number and
`charge_aromatic_proxy` is the minimal charge-and-aromatic version. Both are heuristics. The
honest move is to compute all of curated membership, the sequence score, and the proxy, then let
phase 1 tell us which actually tracks fitness. Published sequence predictors such as PScore,
catGRANULE, and FuzDrop could be added later if any of these signals looks promising.

## Conservation (optional)

| Feature | Definition | Source | Why it might matter |
| --- | --- | --- | --- |
| `phylop_mean` | Mean phyloP over the segment | UCSC phyloP bigWig (NaN unless provided) | Per-base constraint, expected to track functional importance |
| `phastcons_mean` | Mean phastCons over the segment | UCSC phastCons bigWig (NaN unless provided) | Conserved-element probability, a second constraint view |

Conservation needs a residue to genome coordinate map to read the bigWig at the right positions.
That mapping is not wired up yet, so these columns are present but left as missing values for now.
The bigWig reader is in place for when the mapping is added.

## Sources to fetch

`scripts/fetch_databases.py` downloads ELM, dbPTM, PhaSePro, DrLLPS, the UniProt human reference
proteome, and the UCSC conservation bigWigs, and records versions, dates, and URLs in
`data/raw/SOURCES.txt`. metapredict, localCIDER, and AIUPred are installed as the `predictors` and
`binding` extras, not fetched. The bundled example FASTA is synthetic and exists only so the CLI
runs with no downloads.
