"""
Phase 0 -- Step 0: Shared space (panel overlap).

Implements: Step 0 (not a kill gate; it produces the partition every later gate uses).

What it does:
  Intersect the CosMx 64-plex IO protein panel with the METABRIC-IMC 37-marker panel
  (matched at the level of HGNC gene symbol), build the protein -> gene(s) map, and partition
  the CosMx Whole Transcriptome into PANEL-ADJACENT (cognate of a *shared* protein -- a protein
  available on BOTH platforms, hence usable at transfer) vs PANEL-DISTAL.
  This partition is the centerpiece of Gate D.

  STRUCTURAL LOCK (binding rule): the imputer's input feature set == the adjacency-defining set
  == the shared proteins. Polygenic/combined shared channels expand to ALL their cognate genes,
  so a directly-measured molecule can never hide in panel-distal. An all-62-protein run is a
  separate upper bound with its own larger adjacency partition, NOT the binding gate.

Inputs (authoritative; see data/panels/PROVENANCE.md):
  - data/panels/cosmx_io_protein_targets.tsv  (Bruker: ProteinName, Gene -- vendor HGNC symbols)
  - data/panels/cosmx_wtx_genes.tsv           (Bruker: DisplayName, GeneSymbols; 18,934 targets)
  - data/panels/metabric_markers.txt          (Zenodo 6036188: Epitope names; 37 markers + DNA1/2)

Outputs (PATHS["results_dir"]):
  - phase0_step0_protein_gene_map.tsv  full METABRIC-marker -> gene mapping + CosMx match + WTx flag
  - phase0_step0_shared_proteins.tsv   the shared protein set (CosMx IO INTERSECT METABRIC)
  - phase0_step0_panel_partition.tsv   every WTx gene labelled panel-adjacent / panel-distal
  - phase0_step0_SUMMARY.md            counts + every ambiguous/dropped decision; PENDING REVIEW

HARD RULE: panel lists come from real files only (above). The only curation here is normalizing
METABRIC *epitope display names* to HGNC gene symbols via the explicit, reviewable table below
(the approved Step-0 plan calls for this). CosMx protein->gene comes straight from the vendor
'Gene' column. Every non-trivial / ambiguous / dropped mapping is surfaced for review. The
resulting partition is marked PENDING REVIEW -- not final until signed off.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PATHS  # noqa: E402


# ---------------------------------------------------------------------------
# METABRIC-IMC epitope (display name, verbatim) -> HGNC gene symbol(s).
# Basis for each is the standard HGNC-approved symbol for the antibody's target.
# `genes`: tuple of symbols (empty => not gene-mappable / dropped).
# `note`:  surfaced in the summary; flags ambiguity, combined channels, drops.
# Where the same protein also exists in the CosMx IO panel, the symbol is chosen to
# match the CosMx vendor symbol so the gene-level intersection is exact.
# ---------------------------------------------------------------------------
METABRIC_MARKER_MAP: dict[str, dict] = {
    "Histone H3":        {"genes": (),                       "note": "DROP: histone (H3 multigene cluster), structural/housekeeping, not a TME target"},
    "SMA":               {"genes": ("ACTA2",),               "note": "smooth muscle actin"},
    "CK5":               {"genes": ("KRT5",),                "note": "cytokeratin 5"},
    "CD38":              {"genes": ("CD38",),                "note": ""},
    "HLA-DR":            {"genes": ("HLA-DRA", "HLA-DRB1", "HLA-DRB5"), "note": "MHC-II complex (DRA+DRB); SHARED via CosMx HLA-DR; all cognate genes adjacent. In WTx: HLA-DRA, HLA-DRB5 (HLA-DRB1 absent from WTx)"},
    "CK8-18":            {"genes": ("KRT8", "KRT18"),        "note": "combined cytokeratin 8/18; NOT shared (CosMx IO has no keratin) -> KRT8/KRT18 stay panel-distal"},
    "CD15":              {"genes": ("FUT4",),                "note": "SSEA-1/Lewis-x; FUT4"},
    "FSP1":              {"genes": ("S100A4",),              "note": "fibroblast-specific protein 1 = S100A4"},
    "CD163":             {"genes": ("CD163",),               "note": ""},
    "ICOS":              {"genes": ("ICOS",),                "note": "CD278"},
    "OX40":              {"genes": ("TNFRSF4",),             "note": "CD134"},
    "CD68":              {"genes": ("CD68",),                "note": ""},
    "HER2 (3B5)":        {"genes": ("ERBB2",),               "note": "HER2, clone 3B5 (one of two HER2 channels)"},
    "CD3":               {"genes": ("CD3D", "CD3E", "CD3G"), "note": "polyclonal CD3; CD3D/E/G (matches CosMx)"},
    "Podoplanin":        {"genes": ("PDPN",),                "note": ""},
    "CD11c":             {"genes": ("ITGAX",),               "note": ""},
    "PD-1":              {"genes": ("PDCD1",),               "note": "CD279"},
    "GITR":              {"genes": ("TNFRSF18",),            "note": ""},
    "CD16":              {"genes": ("FCGR3A",),              "note": "FCGR3A (matches CosMx); FCGR3B also detected"},
    "HER2 (D8F12)":      {"genes": ("ERBB2",),               "note": "HER2, clone D8F12 (second HER2 channel) -> same gene ERBB2"},
    "CD45RA":            {"genes": ("PTPRC",),               "note": "PTPRC isoform RA"},
    "B2M":               {"genes": ("B2M",),                 "note": ""},
    "CD45RO":            {"genes": ("PTPRC",),               "note": "PTPRC isoform RO (same gene as CD45RA/CD45)"},
    "FOXP3":             {"genes": ("FOXP3",),               "note": ""},
    "CD20":              {"genes": ("MS4A1",),               "note": ""},
    "ER":                {"genes": ("ESR1",),                "note": "estrogen receptor (Ho165 Epitope; AbPanel reagent row mislabels as Rabbit IgG)"},
    "CD8":               {"genes": ("CD8A",),                "note": "CD8A (matches CosMx); CD8B also"},
    "CD57":              {"genes": ("B3GAT1",),              "note": "HNK-1 = B3GAT1"},
    "Ki-67":             {"genes": ("MKI67",),               "note": ""},
    "PDGFRB":            {"genes": ("PDGFRB",),              "note": "CD140b"},
    "Caveolin-1":        {"genes": ("CAV1",),                "note": ""},
    "CD4":               {"genes": ("CD4",),                 "note": ""},
    "CD31-vWF":          {"genes": ("PECAM1", "VWF"),        "note": "combined CD31(PECAM1)+vWF(VWF) endothelial channel; SHARED via CosMx CD31->PECAM1. The METABRIC channel measures both, so both cognate genes are adjacent (conservative vs false headroom). Both in WTx"},
    "CXCL12":            {"genes": ("CXCL12",),              "note": "SDF-1"},
    "HLA-ABC":           {"genes": ("HLA-A", "HLA-B", "HLA-C"), "note": "MHC-I pan class-I; NOT shared (CosMx IO has no MHC-I) -> HLA-A/B/C stay panel-distal"},
    "panCK":             {"genes": ("KRT5", "KRT8", "KRT18", "KRT19"), "note": "pan-cytokeratin; NOT shared (CosMx IO has NO keratin/panCK at all) -> all KRT genes correctly panel-distal (verified against cosmx_io_protein_targets.tsv)"},
    "c-Caspase3c-PARP":  {"genes": ("CASP3", "PARP1"),       "note": "run-together channel split into cleaved Caspase3 (CASP3) + cleaved PARP (PARP1); NOT shared (CosMx IO has no caspase/PARP) -> both stay panel-distal"},
    "DNA1":              {"genes": (),                       "note": "DROP: Ir191 DNA intercalator, not an antibody target"},
    "DNA2":              {"genes": (),                       "note": "DROP: Ir193 DNA intercalator, not an antibody target"},
}


def _read_commented_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a TSV that may have leading '#' comment lines. Returns (header, rows)."""
    header: list[str] = []
    rows: list[list[str]] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or line == "":
                continue
            fields = line.split("\t")
            if not header:
                header = fields
            else:
                rows.append(fields)
    return header, rows


def _read_marker_list(path: Path) -> list[str]:
    out: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or line == "":
                continue
            out.append(line)
    return out


def main() -> None:
    panels = PATHS["panels"]
    results = PATHS["results_dir"]
    results.mkdir(parents=True, exist_ok=True)

    # ---- Load CosMx IO protein -> gene (vendor; may be multi-gene, comma-separated) ----
    _, io_rows = _read_commented_tsv(panels / "cosmx_io_protein_targets.tsv")
    cosmx_protein_genes: dict[str, tuple[str, ...]] = {}
    for r in io_rows:
        name = r[0].strip()
        gene = r[1].strip() if len(r) > 1 else ""
        if gene == "":
            continue  # isotype controls / footer rows
        genes = tuple(g.strip() for g in gene.split(",") if g.strip())
        cosmx_protein_genes[name] = genes
    cosmx_io_gene_set = {g for genes in cosmx_protein_genes.values() for g in genes}

    # ---- Load CosMx WTx transcriptome gene set (authoritative HGNC symbols) ----
    _, wtx_rows = _read_commented_tsv(panels / "cosmx_wtx_genes.tsv")
    wtx_genes: set[str] = set()
    for r in wtx_rows:
        sym = r[1].strip() if len(r) > 1 else ""
        if sym == "":
            continue  # negative-control probes / footer
        for g in sym.split(";"):
            g = g.strip()
            if g:
                wtx_genes.add(g)

    # ---- Load METABRIC markers (verbatim) and normalise via the table above ----
    metabric_markers = _read_marker_list(panels / "metabric_markers.txt")
    unmapped = [m for m in metabric_markers if m not in METABRIC_MARKER_MAP]
    if unmapped:
        raise SystemExit(f"METABRIC markers missing from override table (refusing to guess): {unmapped}")

    # ---- Build the protein->gene map rows + determine shared proteins ----
    # STRUCTURAL LOCK (Gate-D binding rule): the imputer's INPUT FEATURE SET is exactly the
    # shared proteins (everything METABRIC has at transfer), and ADJACENCY IS DEFINED BY THAT
    # SAME SET. A channel is SHARED iff >=1 of its cognate genes is measured by a CosMx IO
    # protein. For a shared channel, ALL of its cognate genes (the polygenic/combined expansion)
    # are panel-adjacent -- so a measured molecule can never sit in panel-distal and let a model
    # post fake headroom. Genes measured only by NON-shared CosMx proteins are NOT adjacent here
    # (an all-62 upper-bound run would carry its own, larger adjacency partition and is NOT the
    # binding gate).
    map_rows = []          # full per-marker mapping
    shared_cognate_genes: set[str] = set()   # union of ALL cognate genes of shared channels
    shared_marker_count = 0
    for marker in metabric_markers:
        entry = METABRIC_MARKER_MAP[marker]
        genes = entry["genes"]
        matched = sorted(set(genes) & cosmx_io_gene_set)   # which cognate genes hit CosMx
        is_shared = len(matched) > 0
        if is_shared:
            shared_marker_count += 1
            shared_cognate_genes.update(genes)             # ALL cognate genes -> adjacency
        cosmx_proteins = sorted({p for p, gs in cosmx_protein_genes.items()
                                 if set(gs) & set(matched)})
        in_wtx = sorted(g for g in genes if g in wtx_genes)
        adjacency_genes = sorted(g for g in genes if g in wtx_genes) if is_shared else []
        map_rows.append({
            "metabric_marker": marker,
            "mapped_genes": ";".join(genes) if genes else "",
            "shared_with_cosmx": "yes" if is_shared else "no",
            "matched_genes": ";".join(matched),
            "cosmx_protein(s)": ";".join(cosmx_proteins),
            "mapped_genes_in_wtx": ";".join(in_wtx),
            "adjacency_genes": ";".join(adjacency_genes),
            "note": entry["note"],
        })

    # ---- Partition the WTx transcriptome ----
    # panel-adjacent = WTx gene that is a cognate of a SHARED protein (transfer-available).
    adjacent_genes = sorted(shared_cognate_genes & wtx_genes)
    shared_genes_missing_from_wtx = sorted(shared_cognate_genes - wtx_genes)
    # broader, for reviewer context: genes measured by ANY CosMx IO protein (not just shared)
    cosmx_any_adjacent = sorted(cosmx_io_gene_set & wtx_genes)

    adjacent_set = set(adjacent_genes)
    distal_genes = sorted(g for g in wtx_genes if g not in adjacent_set)

    # ---- Write outputs ----
    with (results / "phase0_step0_protein_gene_map.tsv").open("w", encoding="utf-8", newline="\n") as f:
        cols = ["metabric_marker", "mapped_genes", "shared_with_cosmx", "matched_genes",
                "cosmx_protein(s)", "mapped_genes_in_wtx", "adjacency_genes", "note"]
        f.write("\t".join(cols) + "\n")
        for row in map_rows:
            f.write("\t".join(row[c] for c in cols) + "\n")

    with (results / "phase0_step0_shared_proteins.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("metabric_marker\tmatched_genes\tcosmx_protein(s)\n")
        for row in map_rows:
            if row["shared_with_cosmx"] == "yes":
                f.write(f"{row['metabric_marker']}\t{row['matched_genes']}\t{row['cosmx_protein(s)']}\n")

    with (results / "phase0_step0_panel_partition.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("gene\tpartition\tmeasured_by_shared_protein\tmeasured_by_any_cosmx_io_protein\n")
        cosmx_any_set = set(cosmx_any_adjacent)
        for g in sorted(wtx_genes):
            part = "panel-adjacent" if g in adjacent_set else "panel-distal"
            f.write(f"{g}\t{part}\t{'yes' if g in adjacent_set else 'no'}"
                    f"\t{'yes' if g in cosmx_any_set else 'no'}\n")

    # ---- Summary ----
    dropped = [r for r in map_rows if r["mapped_genes"] == ""]
    ambiguous = [r for r in map_rows if r["note"] and ("ambiguous" in r["note"].lower()
                 or "combined" in r["note"].lower() or "multi" in r["note"].lower()
                 or "representative" in r["note"].lower() or "run-together" in r["note"].lower())]
    not_shared = [r for r in map_rows if r["shared_with_cosmx"] == "no" and r["mapped_genes"] != ""]

    lines = []
    lines.append("# Phase 0 -- Step 0: panel overlap  (STATUS: PENDING REVIEW -- partition not final)\n")
    lines.append("Generated by scripts/phase_0_00_panel_overlap.py from authoritative panel lists.\n")
    lines.append("## Counts\n")
    lines.append(f"- CosMx IO protein targets (with gene): **{len(cosmx_protein_genes)}** "
                 f"({len(cosmx_io_gene_set)} distinct genes)")
    lines.append(f"- CosMx WTx transcriptome genes: **{len(wtx_genes)}**")
    lines.append(f"- METABRIC-IMC channels: **{len(metabric_markers)}** "
                 f"({sum(1 for m in metabric_markers if METABRIC_MARKER_MAP[m]['genes'])} gene-mappable, "
                 f"{len(dropped)} dropped)")
    lines.append(f"- **Shared proteins (CosMx IO ∩ METABRIC, gene-level): {shared_marker_count} METABRIC channels "
                 f"-> {len(shared_cognate_genes)} distinct cognate genes**")
    lines.append("")
    lines.append("## STRUCTURAL LOCK (binding Gate-D rule -- do not relitigate)\n")
    lines.append("- The imputer's INPUT FEATURE SET == the ADJACENCY-DEFINING SET == the "
                 f"**{shared_marker_count} shared proteins** (all METABRIC has at transfer).")
    lines.append("- A gene is panel-adjacent iff it is a cognate of a shared protein (polygenic/combined")
    lines.append("  channels expand to ALL cognate genes), so a directly-measured molecule can never sit in")
    lines.append("  panel-distal and let the model post fake headroom.")
    lines.append("- Training on all 62 CosMx proteins is NOT the binding gate: the ~40 non-shared IO proteins")
    lines.append("  would directly measure their own genes, which must then move OUT of panel-distal. Any")
    lines.append("  all-62 upper-bound run carries its OWN larger adjacency partition and is reported separately.")
    lines.append("")
    lines.append("## Partition of the WTx transcriptome (PENDING REVIEW)\n")
    lines.append(f"- **panel-adjacent** (cognate of a shared protein): **{len(adjacent_genes)}** genes")
    lines.append(f"- **panel-distal**: **{len(distal_genes)}** genes")
    lines.append(f"- (context only, NOT the gate) genes measured by ANY CosMx IO protein and in WTx: "
                 f"{len(cosmx_any_adjacent)}")
    if shared_genes_missing_from_wtx:
        lines.append(f"- shared-protein cognate genes NOT in WTx (cannot anchor adjacency): "
                     f"{', '.join(shared_genes_missing_from_wtx)}")
    lines.append("")
    lines.append(f"panel-adjacent genes ({len(adjacent_genes)}): {', '.join(adjacent_genes)}\n")
    lines.append("## Shared proteins (the transfer-available set)\n")
    for r in map_rows:
        if r["shared_with_cosmx"] == "yes":
            lines.append(f"- {r['metabric_marker']} -> adjacency {{{r['adjacency_genes']}}}  "
                         f"(CosMx match: {r['matched_genes']} via {r['cosmx_protein(s)']})")
    lines.append("")
    lines.append("## Gate-D composition notes\n")
    lines.append("- The shared/adjacent set is **immune- and stromal-dominated**; **ERBB2 is the only")
    lines.append("  tumor-intrinsic anchor** among the adjacent genes. Gate D's panel-distal signal therefore")
    lines.append("  tests mostly whether protein predicts non-cell-type programs in immune/stromal context.")
    lines.append("- **Priority watch-gene: ESR1 (ER).** Present in METABRIC (Ho165) but NOT in the CosMx IO")
    lines.append("  panel -> it is panel-DISTAL. ER is a key tumor-intrinsic axis in breast cancer, so ESR1")
    lines.append("  behaviour in panel-distal is a headline target to watch at Gate D / transfer.")
    lines.append("")
    lines.append("## METABRIC markers NOT shared with CosMx (mapped, but no CosMx counterpart)\n")
    for r in not_shared:
        lines.append(f"- {r['metabric_marker']} -> {r['mapped_genes']}   {('['+r['note']+']') if r['note'] else ''}")
    lines.append("")
    lines.append("## Dropped channels (not gene-mappable)\n")
    for r in dropped:
        lines.append(f"- {r['metabric_marker']}   [{r['note']}]")
    lines.append("")
    lines.append("## Ambiguous / combined / representative mappings to confirm at review\n")
    for r in ambiguous:
        lines.append(f"- {r['metabric_marker']} -> {r['mapped_genes']}   [{r['note']}]")
    lines.append("")
    lines.append("## Resolved review items (from the two panel files; no memory reconstruction)\n")
    lines.append("- **panCK / CK8-18 (keratins): CosMx IO contains NO keratin/panCK channel at all** "
                 "(verified against cosmx_io_protein_targets.tsv) -> KRT5/KRT8/KRT18/KRT19 are in WTx but")
    lines.append("  CORRECTLY panel-distal (not a dropped mapping; no false-headroom fix needed).")
    lines.append("- **HLA-DR (shared): expanded to all cognate genes** HLA-DRA, HLA-DRB1, HLA-DRB5 -> "
                 "adjacent = {HLA-DRA, HLA-DRB5} (HLA-DRB1 absent from WTx).")
    lines.append("- **HLA-ABC (NOT shared; no MHC-I in CosMx IO):** HLA-A/HLA-B/HLA-C stay panel-distal.")
    lines.append("- **CD31-vWF (shared): expanded to PECAM1 + VWF** (both in WTx) -> both adjacent.")
    lines.append("- **c-Caspase3c-PARP: split into CASP3 + PARP1** (map correctness); NOT shared "
                 "(no caspase/PARP in CosMx IO) -> both stay panel-distal.")
    lines.append("- **FCGR3A (CD16): kept as-is** -- shared protein, but FCGR3A absent from WTx, so it cannot")
    lines.append("  anchor adjacency.")
    lines.append("- **Histone H3 / DNA1 / DNA2: kept dropped** (structural / intercalator, not gene targets).")
    lines.append("")
    lines.append("## Review notes\n")
    lines.append("- CosMx protein->gene is the Bruker vendor 'Gene' column (authoritative).")
    lines.append("- METABRIC epitope->gene is the explicit table in this script; every entry above is reviewable.")
    lines.append("- 'panel-adjacent' uses SHARED proteins only (the set available after transfer to METABRIC),")
    lines.append("  which is what Gate D's kill test requires. The broader 'any CosMx IO protein' column is")
    lines.append("  context only.")
    lines.append("- Partition is PENDING REVIEW; approvable after this resolution pass.")

    (results / "phase0_step0_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # console digest
    print(f"CosMx IO proteins: {len(cosmx_protein_genes)} | WTx genes: {len(wtx_genes)} | "
          f"METABRIC channels: {len(metabric_markers)}")
    print(f"Shared: {shared_marker_count} METABRIC channels -> {len(shared_cognate_genes)} cognate genes")
    print(f"Partition: panel-adjacent={len(adjacent_genes)}  panel-distal={len(distal_genes)}  (PENDING REVIEW)")
    print(f"Outputs written to {results}")


if __name__ == "__main__":
    main()
