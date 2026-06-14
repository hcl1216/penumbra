"""
Phase 0 -- Step 0: Shared space (panel overlap).

Implements: Step 0 (not a kill gate; it produces the partition every later gate uses).

What it does (design only -- NO logic yet):
  Intersect the CosMx 64-plex IO protein target list with METABRIC's 37-marker
  (Ali/Danenberg) list, build the protein -> gene(s) map (e.g. Ki67 -> MKI67;
  CD3 -> CD3D/E/G; PanCK -> KRT5/8/18/19), and partition the CosMx transcriptome into
  PANEL-ADJACENT (a shared protein measures the gene) vs PANEL-DISTAL (none does).
  This partition is the centerpiece of Gate D.

Inputs:
  - CosMx 64-plex IO protein target list        (PATHS["panels"])
  - METABRIC / Ali-2020 37-marker list          (PATHS["panels"])
  - CosMx transcriptome gene list (>18k WTx)     (from PATHS["cosmx_slide"])

Outputs (to PATHS["results_dir"] / PATHS["data_dir"]):
  - shared protein list (CosMx INTERSECT METABRIC)
  - protein -> gene(s) mapping table
  - panel-adjacent / panel-distal gene partition
  - a short written summary (counts, ambiguous mappings, dropped targets)

Status: STUB. No logic until Step-0 plan is approved.
"""
