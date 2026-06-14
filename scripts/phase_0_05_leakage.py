"""
Phase 0 -- Gate C: Leakage controls.

Implements: Gate C (prove the signal isn't a trivial coupling or a size artifact).

What it does (design only -- NO logic yet):
  - SHUFFLE protein across cells -> the map must collapse to the no-protein baseline
    (if it doesn't, the "signal" is not coming from protein).
  - REGRESS OUT total counts + segmentation area on both sides -> signal must survive
    removing the cell-size axis (the most likely artifact).
  (Lateral spillover is partly absorbed by the spatial-kNN baseline in Gate B.)

Inputs:
  - candidate-map predictions/scores (phase_0_04_model)
  - baseline scores (phase_0_03_baselines)
  - QC'd matrices incl. total counts + segmentation area (phase_0_01_load_qc)

Outputs (to PATHS["results_dir"]):
  - shuffled-protein control scores (expect collapse to no-protein baseline)
  - regress-out control scores (expect signal to survive)
  - recorded Gate-C verdict

Status: STUB. No logic until Gate C is approved.
"""
