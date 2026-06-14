"""
Phase 0 -- Gate B: Baselines + the bar.

Implements: Gate B (the honest non-novel predictors any map must beat).

What it does (design only -- NO logic yet):
  On held-out cells, compute the baselines the candidate map must beat by a margin
  whose bootstrap CI over cells excludes zero:
    1. global mean profile               (trivial)
    2. library-size-scaled mean          (trivial)
    3. spatial k-NN average (neighbours' mean RNA, no protein)  -- the killer
    4. cell-type mean (type from protein or RNA clusters -> per-type mean RNA)
  Beating 1-2 is trivial; 3 and 4 are the real bars.

Metric (CLAUDE.md): per-gene correlation of held-out cells AFTER subtracting the
cell-type mean; report the distribution across genes + bootstrap CIs. Never a single
mean correlation alone.

Inputs:
  - detectable-gene set from Gate A (phase_0_02_noise_floor)
  - QC'd RNA + protein matrices, spatial coordinates (phase_0_01_load_qc)

Outputs (to PATHS["results_dir"]):
  - per-gene baseline predictions/scores for baselines 1-4 on held-out cells
  - recorded Gate-B verdict

Status: STUB. No logic until Gate B is approved.
"""
