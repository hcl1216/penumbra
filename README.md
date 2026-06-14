# Penumbra

Cross-modal imputation for tumor-microenvironment (TME) biomarker discovery. We train a
**multiplex protein → whole-transcriptome** map on same-cell spatial multiomics (CosMx SMI),
transfer it onto an outcome-linked spatial-proteomic cohort that measured only a thin protein
panel (METABRIC-IMC), and test whether the imputation reveals a TME biomarker invisible to the
measured proteins — validated against real clinical outcome. Both outcomes are publishable by
design: a clean negative (imputation adds nothing over measured protein) is a useful, contrarian
result. The project is currently at **Phase 0**, a sequence of cheap kill gates run within a
single CosMx slide; see [CLAUDE.md](CLAUDE.md) for the full design, apparatus, and kill criteria.

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
