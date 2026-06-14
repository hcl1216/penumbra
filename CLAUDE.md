# Penumbra

Cross-modal imputation for TME biomarker discovery. Personal research project, independent of any
employer. Fresh repo; reuses the *apparatus* from the prior perturbation work (Operon) but none of
its code or conclusions.

## Goal

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
- **Both outcomes are publishable by design.** A clean negative (imputation adds nothing over
  measured protein) is a contrarian, useful result. Set every test up so positive AND negative are
  interesting.

---

## STATUS: Phase 0 -- Gate-1 kill test. Skeleton committed (6fccf97). All 3 panel name-lists IN HAND (authoritative): CosMx IO protein (Bruker), CosMx WTx genes (Bruker), METABRIC-IMC 37-marker (Zenodo 6036188, Danenberg 2022). Step 0 (panel overlap) DONE -- 22 shared protein channels, partition = 21 panel-adjacent / 18,912 panel-distal genes, **PENDING REVIEW** (results/phase0_step0_*). STOP for sign-off before Gate A. CosMx slide still absent (not needed for Step 0; required for Gates A-D).

Do **not** build the imputer until the information is proven present. Phase 0 is a sequence of cheap
kill gates run almost entirely within a **single** CosMx slide. The Phase-0 deliverable is a verdict:
KILL (write the negative) or PROCEED.

**First actions for Claude Code (blocking -- do these before any modeling):**
1. Confirm the CosMx breast multiomic slide is downloaded and record its path (see Data).
2. Obtain the two panel lists (CosMx 64-plex IO protein targets; METABRIC/Ali-2020 37-marker list)
   and run Step 0 (panel overlap) -- everything downstream depends on the panel-adjacent / panel-
   distal partition it produces.
Do not start coding the gates against absent data; if either input is missing, say so and stop.

---

## Logging discipline -- the decision log is not optional

Every session **must** append an entry to `DECISIONS.md` (reverse-chronological, newest on top)
whenever ANY of these happen -- no exceptions:
- a gate produces a verdict (pass/fail) or **any** quantitative result;
- a kill criterion fires;
- we pivot, change scope, or abandon an approach;
- a load-bearing assumption is confirmed or falsified;
- a design decision is made that future sessions must not relitigate.

Rules:
- **Design goes in before running; the result goes in after.** Two entries if needed.
- Every logged RESULT/VERDICT must **also update the STATUS line** at the top of this file.
- Results live in `DECISIONS.md` and `results/`, **never only** in commit messages or chat.
- Keep entries short: date, type (RESULT / VERDICT / PIVOT / DECISION), one-line summary, the
  numbers if any, and the consequence (what it changes / what we do next).

This complements apparatus rule 8 ("record design before running, results after"); `DECISIONS.md`
is the durable, human-readable trail and `results/` holds the artifacts.

---

## The binding constraint -- read before writing any code

The only open paired CosMx breast multiomic data is **one demonstration slide** (one patient). So the
only holdout available is held-out **cells from the same patient and batch**. That tests whether
protein -> RNA generalizes to new cells of the same tissue -- necessary, very weak. It **cannot** test
cross-patient or cross-platform transfer, which is what this project actually needs.

Consequence: Phase 0 can only **falsify** (fail within-slide -> dead) or **earn the right to
continue**. A within-slide pass is **not** evidence the project works. Never report or treat a
within-slide result as transfer evidence.

---

## Apparatus (non-negotiable; inherited)

1. **Floors-first.** Establish a noise floor and an honest non-novel baseline before any model.
   Signal counts only above the floor. Effect ~= floor => inconclusive, not a win.
2. **Beat the baseline or it isn't real.** Baseline = the honest non-novel predictor (mean,
   library-size, spatial neighbours, cell-type mean). Beat it by a margin stable beyond seed noise.
   Never weaken a baseline to flatter a model.
3. **Size the effect before building to exploit it.** Confirm an effect is real (survives controls)
   and has headroom before building anything bigger.
4. **Leakage controls.** Shuffle / regress-out / hold-out to prove a signal isn't a trivial coupling
   or an artifact.
5. **Batch canary** for any cross-study / cross-platform stitch. If batch ~= signal, the stitch is
   confounded.
6. **Both outcomes interesting.** Always.
7. **Don't make an unsolved problem a prerequisite.** Borrow the hard parts; contribute the edge.
8. **Record design before running, results after** -- in this file and in `results/`.

---

## Build vs borrow

There is **no drop-in protein -> RNA imputer** for same-cell spatial multiomics. Closest prior art is
RNA -> RNA spatial imputation (SpaFormer, and the scVI / SAVER / Tangram / SpaGE / Sprod family) and
RNA -> protein (scProSpatial) -- wrong direction or wrong modality. So: **borrow the evaluation
protocol and baseline set; build the protein -> RNA map yourself.**

In Phase 0 the candidate map is **per-gene ridge / gradient-boosting** from the shared proteins
(+ spatial features). **No neural imputer in Phase 0.** If a tuned GBM can't clear the floor, a
transformer won't, and the build is saved.

---

## Phase 0 plan -- the Gate-1 kill test

All within the CosMx slide except Gate E.

**Step 0 -- shared space.** Intersect the CosMx 64-plex IO protein list with METABRIC's 37 markers;
build the protein -> gene map (Ki67 -> MKI67; CD3 -> CD3D/E/G; PanCK -> KRT5/8/18/19; etc.). Partition
the transcriptome into **panel-adjacent** (a shared protein measures the gene) vs **panel-distal**
(none does). This partition is the whole point of Gate D.

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

---

## Metric

cell-eval / PDS does **not** port (that's perturbation *discrimination*, not imputation accuracy).
Use **per-gene correlation of held-out cells after subtracting the cell-type mean** (does the map
predict within-type residual variation?). Report the **distribution across genes + bootstrap CIs**;
signal = (model - best baseline) per gene. **Never report a single mean correlation alone** -- it is
flattered by the easy panel-adjacent genes.

---

## Kill criteria (fixed before running)

Project is dead at Phase 0 if either:
- (i) no protein -> RNA map beats spatial-kNN **and** cell-type-mean on held-out cells (Gate B), or
- (ii) the map beats them **only** on panel-adjacent genes (Gate D).

Gate E is the transfer gate, but A-D come first and cost ~a day. On a KILL, the deliverable is the
written negative result.

---

## Checkpoints

After **each** gate: write the result to `results/`, update this file's STATUS, and **STOP for
review** before proceeding. Do not advance Gate D -> E, and do not build any non-baseline / non-GBM
model, without an explicit go.

---

## Data

- **CosMx breast multiomic (Phase 0):** single FFPE slide, 64-plex IO protein + >18k WTx RNA, same
  cell. Open demonstration dataset from Bruker Spatial Biology.
  **TODO: download, record path here, verify license terms.** Raw lives on Google Drive
  (Colab-mounted).
- **METABRIC-IMC (Gate E / later):** Danenberg 2022, public. Plus METABRIC bulk RNA (Curtis 2012)
  for the Gate-E pseudobulk check.
- Raw data and checkpoints live on Drive (Colab is ephemeral). Code lives in git. Keep only small /
  derived artifacts under `data/`.

---

## Infrastructure / conventions

- **Split:** Claude Code = authoring + CPU work, runs locally (Windows, `C:\penumbra`). Colab = GPU.
  git = the seam. Edit locally -> push -> Colab pulls -> runs -> results back to Drive.
- Colab is ephemeral (~12h, disk wiped) -> data / checkpoints / results on mounted Drive; any
  training must checkpoint + resume.
- Scripts are OS-agnostic Python (must run on both Windows-local and Colab-Linux). Use `pathlib`, no
  hard-coded drive letters in committed code -- read paths from a config / env.
- `pip` on Colab. **License-check every dependency and every dataset before use** (personal project
  -> non-commercial is OK, but verify and record).
- Record design before running, results after -- in this file and `results/`.

---

## Suggested layout

```
penumbra/
  CLAUDE.md
  README.md
  config.(yaml|py)                 # paths (slide, drive, results) -- not hard-coded in scripts
  data/                            # small / derived only; raw on Drive
  scripts/
    phase_0_00_panel_overlap.py    # Step 0: overlap + panel-adjacent/distal partition
    phase_0_01_load_qc.py          # load slide, QC, counts/cell, n cells, n genes
    phase_0_02_noise_floor.py      # Gate A: binomial split-half reliability -> detectable gene set
    phase_0_03_baselines.py        # Gate B: mean / libsize / spatial-kNN / cell-type
    phase_0_04_model.py            # candidate map: per-gene ridge / GBM from shared proteins
    phase_0_05_leakage.py          # Gate C: shuffle protein, regress out size + seg area
    phase_0_06_headroom.py         # Gate D: adjacent vs distal decomposition + verdict
  results/                         # gate outputs, figures, verdicts (mirror to Drive)
```

Phase 0 is **done** when Gates A-D have a recorded verdict and this file's STATUS is updated.

---

## Do not

- Do not build a neural imputer (or any non-GBM model) in Phase 0.
- Do not advance past a gate without the recorded verdict + review.
- Do not treat or report a within-slide result as transfer evidence.
- Do not weaken a baseline to make the map look better.
- Do not report a single mean correlation alone.
- Do not reintroduce surface-protein CITE-seq / spatial-CITE-seq as imputation training data --
  wrong protein class; CosMx (intracellular-capable) is the deliberate bridge.
- Do not try to "fix" the single-slide limit by pretending held-out cells test transfer.

The full master design doc (survey, the falsified Operon arc, and the wider list of what not to
relitigate) lives outside this repo. If context seems missing, ask -- don't reconstruct it.
