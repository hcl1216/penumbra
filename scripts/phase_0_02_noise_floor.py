"""
Phase 0 -- Gate A: Noise / detection floor.

Implements: Gate A (defines the detectable-gene set all evaluation runs on).

What it does (design only -- NO logic yet):
  No technical cell replicate exists, so build the floor by BINOMIAL SPLIT-HALF of each
  cell's counts and compute per-gene split-half reliability at this depth. Drop genes
  with ~0 reliability (unpredictable in principle).

Inputs:
  - QC'd RNA matrix from phase_0_01_load_qc

Outputs (to PATHS["results_dir"] / PATHS["data_dir"]):
  - per-gene split-half reliability
  - detectable-gene set (genes evaluation is allowed to use)
  - report: median counts/cell, n cells, n detectable genes
  - recorded Gate-A verdict

Status: STUB. No logic until Gate A is approved (STOP after each gate per CLAUDE.md).
"""
