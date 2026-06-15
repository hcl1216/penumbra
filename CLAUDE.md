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

## STATUS: Fork 2 — Phase 0 (cohort scoping). Pivoted off Fork 3 (2026-06-15; data-layer blocked, archived in docs/). Carried over: apparatus, the APPROVED Step-0 22-protein shared partition (final; commit a1dec02; repurposed as the protein space shared with METABRIC). Fork-2 Phase-0 gate sequence = TODO (next design pass). Current open item: discovery-cohort scoping DONE (docs/fork2_cohort_scoping.md) — recommend **Jackson/Basel IMC** (same-platform, open, ER+HER2; immune-thin caveat), alternatives Keren/Risom; **PENDING APPROVAL, not finalized**. The WTX RNA-only CosMx breast slide is NOT a Fork-2 input.

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

## Phase 0 plan (Fork 2) — TODO (next design pass)

> **TODO STUB — do not invent the gate sequence yet.** The Fork-2 Phase-0 kill-gate sequence is to
> be designed jointly in the next session. It must instantiate the apparatus above for a
> *discovery-in-one-cohort → validate-in-METABRIC* design (noise/detection floor; honest non-novel
> baselines incl. cell-type and standard clinical markers; leakage controls; **batch canary on the
> shared markers** before any cross-cohort claim; size-the-effect before building). Do NOT start
> building any Fork-2 gate or model until this plan is written and approved.

Decided inputs already in place: the apparatus, the approved 22-protein shared partition, the
METABRIC-IMC validation target, and the panel HARD RULE. Open input: the discovery cohort (scoping
below).

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
