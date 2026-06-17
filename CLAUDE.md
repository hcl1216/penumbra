# Penumbra

Single-cohort spatial-proteomic TME biomarker discovery, validated in METABRIC-IMC. Personal
research project, independent of any employer. Fresh repo; reuses the *apparatus* from the prior
perturbation work (Operon) but none of its code or conclusions.

**This is Fork 2.** The project pivoted off Fork 3 (cross-modal protein→RNA imputation) on
2026-06-15 because the paired CosMx substrate Fork 3 needed is inaccessible to a solo researcher
(protein behind AtoMx; public breast deposit is WTX RNA-only). The full Fork-3 spec is archived,
not deleted, in [`docs/fork3_superseded.md`](docs/fork3_superseded.md) — it is a real negative
finding. See `DECISIONS.md` (PIVOT, 2026-06-15).

## Goal (Fork-2 thesis)

Discover a tumor-microenvironment (TME) biomarker in **one** outcome-linked spatial-proteomic
cohort (IMC / CODEX / MIBI), and **validate it in METABRIC-IMC** (Danenberg 2022; ~693 breast
tumors, 37-protein IMC + genomics + survival). No imputation, no paired data, no cross-modal
aligner — so the imputation/circularity threat that sank Fork 3 is gone. Discovery and validation
both live in **protein space**, and cross-cohort validation runs over the **markers shared between
the discovery cohort and METABRIC** (the approved 22-protein partition is the starting shared set).

- **Discovery cohort:** one public spatial-proteomic breast cohort with linked outcome — **TBD**
  (cohort scoping is the current open decision; PENDING approval).
- **Validation cohort:** METABRIC-IMC (Danenberg 2022), public. Authoritative 37-marker panel in
  hand (`data/panels/`).
- **Headline (if positive):** a TME biomarker discovered in one cohort and validated against real
  outcome in METABRIC, additive over standard-of-care.
- **Both outcomes are publishable by design.** A clean negative (the biomarker does not replicate /
  adds nothing over standard markers) is a contrarian, useful result. Set every test up so positive
  AND negative are interesting.

---

## STATUS: Fork-2 cohorts + 20-marker feature space + 88/34 validation subset LOCKED; gate sequence LOCKED (primary = CD8↔tumor proximity; see Phase 0 plan). Building front end. **PART 2 data verify = PASS** (Meyer SCE: coords Pos_X/Y, 19/20 markers, PID, DFS+OS, 215 pts/48 DFS events; METABRIC SingleCells.fst: coords, 20/20 markers, metabric_id join, TNBC subset 88 pts/87 imaged/34 events/172k cells; results/fork2_data_verify.md). Flag: METABRIC has no age → SoC composite = grade+stage. PART 3 front-end: script (scripts/fork2_01_typing_floors.R) + thin Colab runner (notebooks/fork2_01_colab_runner.ipynb) READY. Local execution blocked (bioconductor.org unreachable; no Rtools) → runs on Colab (script paths via env vars). **Front end RAN on Colab (R1–R4 in Drive results/).** R1 markers 20→17 (guard passed). R3 Gate-B canary PASS (CD8 profile corr 0.93, CD8-frac batch margin 4.60). R4 Gate-A floors fit on Meyer DFS (n=152/40 ev): composition CD8-frac p=0.69 (FLAT = low bar, good); SoC = nodes (pN2 HR6.2, pN3 HR8.6). R2 cell-typing PASS after refinement (25 metaclusters -> 7 coherent types: Tumor 33/53%, CD8 5.6/4.2%, +B/Myeloid/Stroma/CD4; doublet quarantined; ~42/31% Other but CD8/tumor masks clean). R3 re-confirmed (divergence 0.074; all real-lineage profile corr >=0.92). **Gate A + Gate B PASS; typing sound.** Gate C BUILT (pre-registered, Meyer only): scripts/fork2_02_gateC_proximity.R + notebook cells 7-8. PRIMARY = CD8 infiltrated-fraction within 20um of CK+ tumour (imcRtools minDistToCells/buildSpatialGraph, per-core); SECONDARY = median CD8->tumour distance (exploratory); type-preserving within-core permutation null (N=499) + density confound check; fit on Meyer DFS vs flat composition floor. HARD STOP after Meyer fit+permutation; METABRIC untouched; feature NOT locked for Gate D. RAN (R5): pre-registered PRIMARY FAILS -- HR 1.18/SD p=0.27 (does NOT beat flat composition floor) AND fails type-preserving permutation (perm p=0.86, real |z|=1.11 < null mean 2.16). KILL criterion (ii) fired. CD8 tumour-distal (median 62um) so 20um near-floor (median infiltr-frac 2.3%). METABRIC NOT touched; radius NOT re-picked (p-hacking). Negative result -> write up. C-index NA was a reporting bug (fixed). Next: Secondary-2 PRE-REGISTERED (Ki67 tumour proliferative architecture; assortativity enrichment, kNN k=6 tumour-only, Otsu Ki67+, vs Ki67-fraction floor + type-preserving permutation) -- the LAST feature; RAN (R6): Ki67 assortativity FAILS -> KILL (n=85/28 ev; floor Ki67-frac flat C=0.50; primary HR 0.92 p=0.71 dC +0.01; perm p=0.75 real|z|=0.37<null 0.93; no density confound). BOTH Gate-C features now dead (R5 CD8-proximity + R6 Ki67-architecture). FORK-2 DISCOVERY = NEGATIVE: no pre-registered spatial feature beats composition or survives the spatial null in Meyer TNBC. Nothing earns Gate D; METABRIC untouched. Deliverable = written negative -- BUT first a pre-specified COHORT-VALIDITY / positive-control check (R7, scripts/fork2_04_cohort_validity.R): can Meyer detect KNOWN prognostic signal (clinical pN/pT/grade; IMC total-TIL+immune) on DFS AND OS, full 215 cohort? Required to interpret the R5/R6 nulls. Awaiting Colab run. Gate before writing the negative; no third spatial feature. Pivoted off Fork 3 (2026-06-15; data-layer blocked, archived in docs/). Carried over: apparatus, APPROVED Step-0 partition (commit a1dec02). **LOCKED: discovery = Meyer 2025 TNBC IMC (USZ/Zurich, n=215, survival); validation = METABRIC TNBC/basal subset.** Chosen for same-platform-as-METABRIC (cleanest batch canary), CD8+ER+HER2 in shared space, population-independent; FOXP3 deliberately NOT required (scope = CD8/cytotoxic-architecture + receptor-context; Treg-subtype thesis out). **Working feature space = Meyer ∩ METABRIC = 20 protein markers** (immune 8: CD3/CD4/CD8/CD11c/CD15/CD20/CD68/HLA-DR; stromal 2: SMA/vWF; epi-tumor 7: CK5/CK8/CK18/panCK/ER/HER2/Ki-67; other 3: HistoneH3/cl-Casp3/cl-PARP). Meyer CKs resolved from deposited SCE clean_target (Zenodo 15304181). **Validation cohort = METABRIC ER−/HER2− subset: 88 patients, 34 BC-death events** (power-limited → few pre-specified hypotheses). Fork-2 Phase-0 gate sequence = TODO (next design pass, on this locked feature set). Single-cell expression (SingleCells 849MB/341MB) range-extractable, NOT pulled yet. WTX RNA-only CosMx slide is NOT a Fork-2 input. Artifacts: docs/fork2_cohort_scoping.md, results/fork2_marker_overlap.md, data/metabric/IMCClinical.fst.

---

## Logging discipline — the decision log is not optional

Every session **must** append an entry to `DECISIONS.md` (reverse-chronological, newest on top)
whenever ANY of these happen — no exceptions:
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
   confounded. (Binding for Fork 2: discovery cohort → METABRIC is exactly such a stitch.)
6. **Both outcomes interesting.** Always.
7. **Don't make an unsolved problem a prerequisite.** Borrow the hard parts; contribute the edge.
8. **Record design before running, results after** — in this file and in `results/`.

---

## The shared protein space (carried over from Step 0; APPROVED, final)

Step 0 intersected the CosMx IO protein panel with METABRIC's 37-marker IMC panel and built the
authoritative protein→gene map. Result (commit a1dec02, APPROVED — do not relitigate):
**22 shared protein channels → 25 cognate genes**; partition of the WTx transcriptome into
**23 panel-adjacent / 18,910 panel-distal** genes (artifacts in `results/phase0_step0_*`;
script `scripts/phase_0_00_panel_overlap.py`).

**Repurposing under Fork 2.** This is no longer an imputation adjacency partition. It now defines
the **protein feature space shared with METABRIC** — the candidate markers along which a
discovery-cohort finding can be cross-validated in METABRIC-IMC. When the discovery cohort is
chosen, intersect ITS panel with METABRIC the same way (vendor/official panel files only, HARD
RULE below) to get the discovery↔METABRIC shared markers; the 22-protein CosMx∩METABRIC set is the
reference for what METABRIC contributes.

> **HARD RULE (panels).** Every panel / marker fact comes from a real file or official source —
> **never reconstructed or approximated from model knowledge.** A wrong panel silently corrupts the
> shared-marker set and every cross-cohort claim built on it. If a panel is behind a wall
> (login / form / unparseable PDF), report the wall and stop — do not substitute a remembered list.

---

## Phase 0 plan (Fork 2) — LOCKED gate sequence

> **DESIGN LOCKED (2026-06-15; do not relitigate). See DECISIONS.md.**

**BINDING CONSTRAINT: 34 validation events** (METABRIC ER−/HER2− subset, 88 patients). Consequences:
- (a) **≤3 pre-specified spatial features, one PRIMARY** (no wide scans).
- (b) The METABRIC TNBC subset (88/34) is **underpowered independent REPLICATION**; Meyer (215)
  carries effect-existence. Honest framing = **"discovered in Meyer, replicated additively in
  independent METABRIC"**, NOT "p<0.001 in METABRIC".
- (c) Every comparison is **1-feature-vs-1-scalar-floor** — floors are **pre-collapsed to single risk
  scores and FIT on Meyer**, then applied as **FIXED scalars** to METABRIC, so the scarce events are
  spent almost entirely on the one feature's coefficient.

**ROLE SPLIT.** Meyer = **discovery** (explore freely; fit + lock features and floor coefficients
here). METABRIC TNBC subset = **validation** (touched once; estimates only the incremental spatial
effect).

**PRE-SPECIFIED FEATURES.**
- **PRIMARY = CD8↔tumor proximity** — per-patient distance from CD8 T cells to nearest CK+/panCK+
  tumor cell, aggregated to one scalar (median, or infiltrated-fraction within a locked radius).
  Evidence-grounded; beats composition in TNBC literature. Headline, multiplicity-protected.
- **Secondary-1 = proximity CONSISTENCY** (uniform vs patchy). Exploratory, multiplicity-flagged.
- **Secondary-2 (receptor-context) = Ki67-tumor spatial structure / CD8 relative to proliferative
  tumor.** Exploratory, multiplicity-flagged.

**GATES.**
- **A — floors + reliability.** Floors = clinical SoC composite (grade/age; ER/HER2 constant in TNBC)
  + composition scalar (e.g. CD8 fraction), **fit on Meyer**. Plus within-patient feature reliability
  (split-half / across-FOV).
- **B — batch canary.** Harmonize the 20 markers, ONE cross-cohort cell-typing; confirm a
  biologically-stable quantity varies **LESS** cross-cohort than the survival signal. **batch ≈ signal
  ⇒ DEAD.**
- **C — lock features on Meyer + position-permutation control.** Permute cell positions within-patient
  → the feature must **collapse to the composition floor**; if it doesn't, it is leaking composition,
  not spatial.
- **D — validate on METABRIC, touch once.** Feature predicts BC-death, **additive over the clinical
  composite, beats the composition floor.** **METRIC = ΔC-index** (feature+floor vs floor) with
  bootstrap CI; **not bare p<0.05**.

**KILL CRITERIA (fixed).** Dead if: (i) batch ≈ signal (B); (ii) primary doesn't beat the composition
floor on Meyer OR collapses under position-permutation (C); (iii) doesn't replicate additively over
both floors on METABRIC (D) → write the negative (*"spatial immune architecture adds no prognostic
value over composition + SoC in independent TNBC validation"*). Both outcomes publishable.

**STATED ASSUMPTION.** Discovery (Meyer) and validation (METABRIC) endpoints differ —
**recurrence vs breast-cancer-specific death**. Carried as an explicit assumption of the replication.

Decided inputs in place: apparatus, locked 20-marker Meyer∩METABRIC feature space, the Meyer/METABRIC
cohorts, the 88/34 validation subset, and the panel HARD RULE.

---

## Checkpoints

After **each** gate/step: write the result to `results/`, update this file's STATUS, and **STOP for
review** before proceeding. Do not finalize the discovery cohort, and do not build any model, without
an explicit go.

---

## Data

- **Validation cohort — METABRIC-IMC (Danenberg 2022):** public. Authoritative 37-marker IMC panel
  acquired (Zenodo 6036188, CC-BY-4.0) — `data/panels/metabric_markers.txt` (+ raw CSVs,
  provenance in `data/panels/PROVENANCE.md`). METABRIC bulk RNA (Curtis 2012) and clinical/survival
  available if needed downstream.
- **Discovery cohort — TBD:** one public spatial-proteomic breast cohort (IMC/CODEX/MIBI) with
  linked outcome; see cohort scoping (PENDING approval). Panel must be obtained from a real
  file / official source (HARD RULE).
- **Reference panels in hand:** CosMx IO protein (Bruker) and CosMx WTx genes (Bruker) — used for
  the Step-0 partition; CosMx itself is no longer a data dependency.
- **NOT an input:** the public CosMx WTX RNA-only breast slide. It carries no protein and is not a
  Fork-2 dependency; any in-progress download of it can be abandoned.
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

## Layout

```
penumbra/
  CLAUDE.md
  README.md
  config.py                        # paths (drive, results, panels) -- not hard-coded in scripts
  data/panels/                     # authoritative panel name-lists (tracked) + provenance
  docs/
    fork3_superseded.md            # archived Fork-3 spec (negative finding; not active)
  scripts/
    phase_0_00_panel_overlap.py    # Step 0 (DONE, APPROVED): shared-protein partition -- carries over
    phase_0_01_load_qc.py ...      # SUPERSEDED Fork-3 imputation-gate stubs; await Fork-2 redesign
  results/                         # step/gate outputs, figures, verdicts (mirror to Drive)
```

The `phase_0_01..06` scripts are empty Fork-3 imputation-gate stubs (superseded). Leave them until
the Fork-2 Phase-0 plan is designed, then repurpose/replace.

---

## Do not

- Do not finalize the discovery cohort, or build any Fork-2 gate/model, without the written +
  approved Fork-2 Phase-0 plan and an explicit go.
- Do not advance past a step/gate without the recorded verdict + review.
- Do not reconstruct or approximate any panel/marker list from memory (HARD RULE) — real file or
  official source only; report walls.
- Do not make a cross-cohort claim before the batch canary on the shared markers passes.
- Do not weaken a baseline (incl. cell-type mean and standard clinical markers) to make a finding
  look better.
- Do not relitigate the approved Step-0 22-protein partition or the pivot off Fork 3.
- Do not treat the CosMx WTX RNA-only breast slide as a Fork-2 input.

The full master design doc (survey, the falsified Operon arc, the Fork-3 negative, and the wider
list of what not to relitigate) lives outside this repo. If context seems missing, ask — don't
reconstruct it.
