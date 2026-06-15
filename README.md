# Penumbra

Spatial-proteomic tumor-microenvironment (TME) biomarker discovery. **This is Fork 2:** discover a
TME biomarker in **one** outcome-linked spatial-proteomic breast cohort (IMC / CODEX / MIBI) and
**validate it in METABRIC-IMC** (Danenberg 2022) — discovery and validation both in protein space,
over the markers shared between the two cohorts. No imputation, no paired data. Both outcomes are
publishable by design: a clean negative (the biomarker doesn't replicate / adds nothing over
standard markers) is a useful, contrarian result. The project is at **Phase 0** (cohort scoping;
the gate sequence is the next design pass). See [CLAUDE.md](CLAUDE.md) for the full design and
apparatus.

The project pivoted off **Fork 3** (cross-modal protein→RNA imputation) on 2026-06-15 — the paired
CosMx substrate it needed is inaccessible to a solo researcher. That spec is archived as a real
negative finding in [docs/fork3_superseded.md](docs/fork3_superseded.md).

## Where work happens (the Claude Code / Colab / Drive split)

- **Claude Code (local, Windows, `C:\penumbra`)** — authoring + CPU work. Edits code, writes
  configs, runs the cheap CPU gates.
- **Colab (ephemeral, GPU)** — heavier compute. Disk is wiped (~12h sessions), so anything that
  trains must checkpoint and resume.
- **Google Drive (Colab-mounted)** — durable home for **raw data, checkpoints, and results**.
  These never go in git.
- **git** is the seam: edit locally → push → Colab pulls → runs → results land back on Drive.

So the repo holds **code and small/derived artifacts only**. All paths are read from
[`config.py`](config.py) (override via `PENUMBRA_*` environment variables) — no script hard-codes a
drive letter, and the same code runs on Windows-local and Colab-Linux unchanged.

## Layout

```
penumbra/
  CLAUDE.md          design doc, apparatus, kill criteria (read this first)
  config.py          all paths (env-overridable, pathlib, no drive letters)
  requirements.txt   minimal Phase-0 dependencies
  data/              small/derived only (gitignored); data/panels/ tracked
  results/           gate outputs/figures/verdicts (gitignored; mirrored to Drive)
  scripts/           phase_0_00..06 -- one file per step/gate
```

## Status

Phase 0 — scaffolding in place; gate logic not yet written. See `STATUS` in
[CLAUDE.md](CLAUDE.md).
