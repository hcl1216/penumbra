# SUPERSEDED — Fork 3 (cross-modal protein→RNA imputation)

**Status: ARCHIVED, not active.** Superseded on 2026-06-15 by the pivot to Fork 2 (see
`DECISIONS.md`, PIVOT entry). Kept because the blocker is a **real negative finding** worth
preserving: the data substrate Fork 3 needed is inaccessible to a solo researcher.

**Why archived (verbatim from the pivot):** the open paired CosMx breast multiomic substrate is
confirmed inaccessible to a solo researcher (protein lives behind AtoMx, no instrument/tenant
access; the public breast deposit is WTX RNA-only). Fork-3-as-written is blocked at the data layer.

**What carried over to Fork 2 (NOT archived — see CLAUDE.md):** the apparatus, the logging
discipline, the infra/conventions, the checkpoint-and-stop rule, and the **approved 22-protein
Step-0 partition** (repurposed as the protein feature space shared with METABRIC for cross-cohort
validation). The Step-0 script `scripts/phase_0_00_panel_overlap.py` and its outputs
(`results/phase0_step0_*`) remain valid and approved.

**What is dead in Fork 3:** the protein→RNA imputation premise, the CosMx paired-slide dependency,
the single-slide ceiling, the build-vs-borrow imputer, and the imputation Gates A–E (with their
metric and kill criteria). All of that is reproduced verbatim below for the record only.

---

## (archived) Goal — Fork 3

Train a cross-modal imputation (multiplex protein -> whole transcriptome) on same-cell spatial
multiomics, transfer it onto an outcome-linked spatial-proteomic cohort that only measured a thin
protein panel, and discover a TME biomarker invisible to the measured proteins -- validated against
real clinical outcome.

- **Train / bridge:** CosMx SMI Same-Cell Multiomics (paired protein + whole transcriptome on the
  same FFPE cell).
- **Outcome / validation cohort:** METABRIC-IMC (Danenberg 2022) -- ~693 breast tumors, 37-protein
  IMC + genomics + survival.
- **Headline (if positive):** cross-platform imputation reveals a TME biomarker invisible to the
  measured proteins, validated against outcome in METABRIC and additive over standard-of-care.

## (archived) The binding constraint

The only open paired CosMx breast multiomic data is **one demonstration slide** (one patient). So the
only holdout available is held-out **cells from the same patient and batch**. That tests whether
protein -> RNA generalizes to new cells of the same tissue -- necessary, very weak. It **cannot** test
cross-patient or cross-platform transfer, which is what this project actually needs.

Consequence: Phase 0 can only **falsify** (fail within-slide -> dead) or **earn the right to
continue**. A within-slide pass is **not** evidence the project works. Never report or treat a
within-slide result as transfer evidence.

## (archived) Build vs borrow

There is **no drop-in protein -> RNA imputer** for same-cell spatial multiomics. Closest prior art is
RNA -> RNA spatial imputation (SpaFormer, and the scVI / SAVER / Tangram / SpaGE / Sprod family) and
RNA -> protein (scProSpatial) -- wrong direction or wrong modality. So: **borrow the evaluation
protocol and baseline set; build the protein -> RNA map yourself.**

In Phase 0 the candidate map is **per-gene ridge / gradient-boosting** from the shared proteins
(+ spatial features). **No neural imputer in Phase 0.** If a tuned GBM can't clear the floor, a
transformer won't, and the build is saved.

## (archived) Phase 0 plan — the Gate-1 kill test

All within the CosMx slide except Gate E.

**Step 0 -- shared space.** Intersect the CosMx 64-plex IO protein list with METABRIC's 37 markers;
build the protein -> gene map (Ki67 -> MKI67; CD3 -> CD3D/E/G; PanCK -> KRT5/8/18/19; etc.). Partition
the transcriptome into **panel-adjacent** (a shared protein measures the gene) vs **panel-distal**
(none does). This partition is the whole point of Gate D.
*(NOTE: Step 0 was completed and APPROVED; the partition carries over to Fork 2. See CLAUDE.md.)*

**Gate A -- noise / detection floor.** No technical cell replicate exists -> build the floor by
**binomial split-half** of each cell's counts; compute per-gene split-half reliability at this depth.
Drop genes with ~0 reliability (unpredictable in principle). Report median counts/cell, n cells, n
detectable genes. Evaluation runs only on detectable genes.

**Gate B -- baselines + the bar.** On held-out cells, the map must beat, by a margin whose bootstrap
CI over cells excludes zero:
1. global mean profile,
2. library-size-scaled mean,
3. **spatial k-NN average** (neighbours' mean RNA, no protein -- the killer),
4. **cell-type mean** (type from protein or RNA clusters -> per-type mean RNA).
Beating 1-2 is trivial; 3 and 4 are the real bars.

**Gate C -- leakage controls.** Shuffle protein across cells -> must collapse to the no-protein
baseline. Regress out total counts + segmentation area on both sides -> signal must survive removing
the cell-size axis (the most likely artifact). Lateral spillover is partly absorbed by baseline 3.

**Gate D -- headroom (centerpiece + kill criterion).** Decompose (model - best baseline) by
panel-adjacent vs panel-distal genes, against the cell-type mean.
**KILL:** if the map beats baselines only on panel-adjacent genes and collapses to cell-type-mean on
panel-distal genes, the premise is falsified -- protein carries no transferable information about
programs it doesn't directly measure, so any downstream "biomarker invisible to the proteins" is
cell-type in disguise. Write that negative and stop.

**Gate E -- cross-cohort check (only if A-D pass; separate phase, needs the build + harmonization).**
Harmonize the shared proteins CosMx <-> METABRIC IMC (fluorescence Ab-oligo vs Hyperion metal-tag --
run the batch canary on the shared markers first). Impute onto METABRIC, pseudobulk per patient,
correlate vs each patient's real METABRIC bulk RNA (METABRIC patients already have bulk
transcriptomes; Curtis 2012). Floor = cohort-mean bulk / scrambled patient labels. This is the only
cross-patient, cross-platform signal available without new data.

## (archived) Metric — Fork 3

cell-eval / PDS does **not** port (that's perturbation *discrimination*, not imputation accuracy).
Use **per-gene correlation of held-out cells after subtracting the cell-type mean** (does the map
predict within-type residual variation?). Report the **distribution across genes + bootstrap CIs**;
signal = (model - best baseline) per gene. **Never report a single mean correlation alone** -- it is
flattered by the easy panel-adjacent genes.

## (archived) Kill criteria — Fork 3

Project is dead at Phase 0 if either:
- (i) no protein -> RNA map beats spatial-kNN **and** cell-type-mean on held-out cells (Gate B), or
- (ii) the map beats them **only** on panel-adjacent genes (Gate D).

## (archived) Fork-3 "Do not" items

- Do not build a neural imputer (or any non-GBM model) in Phase 0.
- Do not treat or report a within-slide result as transfer evidence.
- Do not reintroduce surface-protein CITE-seq / spatial-CITE-seq as imputation training data --
  wrong protein class; CosMx (intracellular-capable) is the deliberate bridge.
- Do not try to "fix" the single-slide limit by pretending held-out cells test transfer.
