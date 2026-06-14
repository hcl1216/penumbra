"""
Phase 0 -- Step 1: Load slide + QC.

Implements: slide loading and quality control (pre-gate; feeds Gates A-D).

What it does (design only -- NO logic yet):
  Load the single CosMx breast multiomic slide (paired protein + whole-transcriptome,
  same cell), run basic QC, and report median counts/cell, n cells, n genes. Produces
  the in-memory object (e.g. AnnData) all later gates consume.

Inputs:
  - CosMx slide raw flat files (PATHS["cosmx_slide"])

Outputs (to PATHS["results_dir"] / PATHS["data_dir"]):
  - QC summary: median counts/cell, n cells, n detectable genes pre-floor
  - QC'd cell x gene (RNA) and cell x protein matrices for downstream gates

Reminder: only a within-slide holdout exists (one patient). A within-slide result is
NEVER transfer evidence (see CLAUDE.md "binding constraint").

Status: STUB. No logic until scaffolding/plan is approved.
"""
