"""
Phase 0 -- Candidate map: per-gene ridge / GBM from shared proteins.

Implements: the candidate protein -> RNA map evaluated by Gates B-D.
NOTE: per CLAUDE.md, NO neural imputer in Phase 0 -- per-gene ridge / gradient
boosting (+ spatial features) only. If a tuned GBM can't clear the floor, a
transformer won't, and the build is saved.

What it does (design only -- NO logic yet):
  For each detectable gene, fit a per-gene regressor (ridge and/or LightGBM) from the
  shared proteins (+ spatial features) on training cells; predict held-out cells.

Inputs:
  - shared proteins (Step 0, phase_0_00_panel_overlap)
  - detectable-gene set from Gate A
  - QC'd protein matrix + spatial features (phase_0_01_load_qc)
  - same held-out split as the baselines (phase_0_03_baselines)

Outputs (to PATHS["results_dir"]):
  - per-gene held-out predictions from the candidate map
  - per-gene scores comparable to the baselines (same metric)

Status: STUB. No logic until Gate B is reviewed and a build is explicitly approved.
"""
