# Project recipe: features that make an IDR fitness associated, and variant effects in IDRs

## Goal

Improve prediction of how missense variants affect intrinsically disordered regions. We start
from IDR-level function: which disordered regions are required for cell fitness, and what features
explain that. We then move to the residue level: whether a substitution perturbs an IDR's
function. The anchor is the Carter lab melanoma cell-line data that labels IDRs as fitness
associated or not.

## Data

- **Labels.** Melanoma cell-line functional data labeling IDRs as required for fitness or not,
  from Kivil. We do not have it yet. We need the IDR definition used (coordinates, predictor, and
  thresholds), how fitness is called (score, cutoff, replicates, screen type), and a sample file.
  No modeling on the labels happens until that arrives. We will not stand in synthetic labels.
- **Sequences and structure.** UniProt human reference proteome. AlphaFold pLDDT where available.
- **Annotation databases.** ELM motifs, dbPTM modification sites, PhaSePro and DrLLPS for phase
  separation. Optional UCSC phyloP and phastCons for conservation.

## Feature set

`idrfeat` computes one feature row per IDR segment covering composition, charge and its
patterning, hydropathy, k-mer content, low complexity, length, disorder, fold on binding, ELM
motifs, PTM sites, and phase-separation signals. Every feature, its source, and its rationale are
in `docs/features.md`. Sequence features always compute. Annotation and predictor features degrade
to missing values when their source is absent. The fold-on-binding signal is AIUPred binding
(ANCHOR2 style). Phase-separation propensity is covered three ways: curated membership, a
transparent sequence score, and the charge-times-aromatic proxy.

## Phase 1: which features make an IDR fitness associated

Join the feature table to the fitness labels at the IDR level. Fit interpretable models first
(regularized logistic regression, then gradient-boosted trees) to predict the fitness label from
features. Read the model, do not just score it: feature importance, SHAP values, and univariate
tests, with attention to confounds like IDR length and overall disorder. The deliverable is a
ranked, honest account of which features carry signal and which do not, including features that
look promising a priori but add nothing here.

## Phase 2: whether a substitution perturbs IDR function

Once we know which features matter, ask whether a missense substitution moves them. For a variant,
recompute the segment features on the mutant sequence and take deltas (change in charge
patterning, motif disruption, PTM-site loss or gain, predicted disorder or binding shift). Add
per-residue variant effect predictor scores and ESM-based features. Model the perturbation against
whatever functional readout the labels support. Frame it as variant effect within IDRs, where
generic predictors are known to be weaker.

## Evaluation

Group-aware cross-validation split by protein so no protein is in both train and test. Report
ROC-AUC and PR-AUC with bootstrap confidence intervals, and calibration. Compare against simple
baselines (length, mean disorder, a single VEP score) so we know the features earn their place.
Hold out a set of proteins untouched until the end. Fix the seed for reproducibility.

## Milestones

1. Feature library and CLI, feature list, this recipe. Done, in draft.
2. Receive the fitness dataset and agree the IDR and label definitions with Kivil.
3. Build the joined IDR feature plus label table. Sanity-check coverage and class balance.
4. Phase 1 models and an interpretation writeup.
5. Phase 2 variant perturbation features and models.
6. Add ESM embeddings and re-evaluate.

## What needs the cluster

ESM and other protein language model embeddings over the proteome, and any large mutant-sequence
rescoring in phase 2, are the compute-heavy steps. They are queued for when cluster access is
ready. Everything in phase 1 and the feature extraction runs on a laptop.
