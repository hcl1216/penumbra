"""
Phase 0 -- Gate D: Headroom (centerpiece + kill criterion).

Implements: Gate D (the decisive Phase-0 verdict: KILL or PROCEED).

What it does (design only -- NO logic yet):
  Decompose (model - best baseline) by PANEL-ADJACENT vs PANEL-DISTAL genes, against the
  cell-type mean (from Step 0's partition).

  KILL criterion: if the map beats baselines only on panel-adjacent genes and collapses
  to cell-type-mean on panel-distal genes, the premise is falsified -- protein carries no
  transferable information about programs it doesn't directly measure, so any downstream
  "biomarker invisible to the proteins" is cell-type in disguise. Write that negative and
  stop.

Inputs:
  - panel-adjacent / panel-distal partition (Step 0, phase_0_00_panel_overlap)
  - candidate-map scores (phase_0_04_model)
  - baseline scores incl. cell-type mean (phase_0_03_baselines)
  - (passing Gate C is a precondition)

Outputs (to PATHS["results_dir"]):
  - adjacent-vs-distal decomposition of (model - best baseline), with bootstrap CIs
  - recorded Gate-D verdict: KILL (write the negative) or PROCEED (earn Gate E)

Gate E (cross-cohort CosMx <-> METABRIC transfer) is a SEPARATE phase, only if A-D pass,
and is not implemented in this file.

Status: STUB. No logic until Gate D is approved.
"""
