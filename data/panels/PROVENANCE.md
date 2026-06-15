# Panel list provenance

Authoritative-source record for every panel/marker list in `data/panels/`.
HARD RULE (CLAUDE.md / task): panel lists come from a real file or official source only —
never reconstructed or approximated from model knowledge. A wrong panel silently corrupts the
Gate-D panel-adjacent / panel-distal partition.

Retrieved 2026-06-15.

---

## 1. CosMx 64-plex Human Immuno-Oncology Protein Panel  ✅ IN HAND (authoritative)

- **Source:** Bruker Spatial Biology, official resource page
  https://www.brukerspatialbiology.com/resources/cosmx-human-protein-immuno-oncology-protein-details/
- **Raw file:** `cosmx_io_protein_panel_RAW.xlsx` (gitignored — large/binary; re-fetchable)
  - Direct URL: https://brukerspatialbiology.com/wp-content/uploads/2023/02/LBL-11122-06_Human-Protein-Immuno-Oncology-BSB.xlsx
  - sha256: `873c2cac96f99cc0e86c748157d92ebab8c905cee69e317feb63e351c5f17410`
- **Extracted (verbatim, raw names):** `cosmx_io_protein_targets.tsv`
  - Sheet `Protein Details`, columns `Protein Name` (col A) and `Gene` (col B).
- **Counts:** 62 biological protein targets (each with a gene) + 2 isotype controls
  (`Ms IgG1`, `Rb IgG`, no gene) = the **64-plex**. The verbatim TSV also retains 5 trailing
  footer/trademark rows from the spreadsheet (Abcam/CST notices, copyright) — these are NOT
  targets and must be dropped at Step 0 cleaning, not here.
- **Note:** raw `Protein Name` values are antibody/display names (e.g. `4-1BB`, `B7-H3`, `SMA`,
  `Tim-3`); the `Gene` column gives the HGNC symbol (e.g. `TNFRSF9`, `CD276`, `ACTA2`). Name→symbol
  normalization happens at Step 0.

## 2. CosMx Human Whole Transcriptome (WTx) Panel  ✅ IN HAND (authoritative)

- **Source:** Bruker Spatial Biology, official resource page
  https://brukerspatialbiology.com/resources/cosmx-human-whole-transcriptome-wtx-genelist/
- **Raw file:** `cosmx_wtx_genelist_RAW.xlsx` (gitignored — large/binary; re-fetchable)
  - Direct URL: https://brukerspatialbiology.com/wp-content/uploads/2025/05/LBL-11220-01_CosMx_Human_Whole_Transcriptome_WTX_GeneList.xlsx
  - sha256: `7642474a4488b06d9adafaeb22b7243227540616dea8c3ccea646f432393434a`
- **Extracted (verbatim, raw names):** `cosmx_wtx_genes.tsv`
  - Sheet `Gene and Probe Details`, columns `Display Name` (col A) and `Gene Symbol(s)` (col B).
- **Counts:** 18,986 raw rows = **18,934 gene targets** (matches the panel's stated
  "18,934 targets") + 50 negative-control probes (`Negative1..50`, no symbol) + 1 footer row.
  Negative controls and footer must be dropped at Step 0 cleaning.

## Fork-2 discovery-cohort panels (for protein-marker overlap) — all from real source files

- **Jackson/Basel IMC** — `basel_zuri_stainingpanel_RAW.csv`, range-extracted from Zenodo **3518284**
  (`SingleCell_and_Metadata.zip` → `Data_publication/Basel_Zuri_StainingPanel.csv`), CC-BY 1.0.
  sha256 `3706c0c7b07e2947…`. NOTE: Gd156 = `Rabbit IgG H L` (control), **no Estrogen Receptor** in
  this panel. Patient metadata (`Basel_PatientMetadata.csv`) confirms `OSmonth` + `DFSmonth` columns.
- **Keren MIBI-TNBC** — `keren_mibi_panel_RAW.txt`, verbatim `cellData.csv` column header (the channels
  ARE the columns), open mirror github.com/jhausserlab/NIPMAP.
- **Engelhardt/Chang CycIF** — `engelhardt_cycif_proteins_RAW.csv`, fetched from authors' repo
  github.com/engjen/cycIF_TMAs (`data/proteins_by_frame.csv`), 42 markers. sha256 `719f2128c89702de…`.
- **Meyer 2025 IMC** — `meyer2025_imc_panel_RAW.tsv`, from deposited raw IMC acquisition `.txt` headers
  (Zenodo **10890543**, `ZTMA174/249_raw.zip`), CC-BY-4.0, read by research agent via HTTP range.
  Cytokeratin channel identities are clone-truncated in the export → flagged provisional.

Overlap computed by `scripts/fork2_marker_overlap.py` → `results/fork2_marker_overlap.md`.

## METABRIC-IMC validation data (Zenodo 6036188) — partial pull
- **Clinical/survival:** `IMCClinical.fst` range-extracted to `data/metabric/` (gitignored,
  patient-level). 709 patients × 11 cols; endpoint = `yearsToStatus` + `DeathBreast`
  (breast-cancer-specific survival); 230 events. sha256 `ba435cabb9e462e0…`.
  (Read with R `fst` package, **AGPL-3** — logged; install.packages('fst').)
- **Single-cell expression:** `SingleCells.csv` (849 MB) / `SingleCells.fst` (341 MB) — INVENTORIED
  only, range-extractable, not pulled (pull when the gate sequence needs it).

## 3. METABRIC-IMC (Danenberg 2022) marker panel  ✅ IN HAND (authoritative)

- **Source:** Zenodo record **6036188** ("Breast tumour microenvironment structures are associated
  with genomic features and clinical outcome", Danenberg et al. 2022, Nat Genet), license
  **CC-BY-4.0**. https://zenodo.org/records/6036188
  - Container: `MBTMEStrIMCPublic.zip` (6.65 GB; Zenodo MD5 `992d04caf3cefcca7bfc5bb64813297f`).
  - **Not fully downloaded.** Extracted only the two small panel members via HTTP **range
    requests** (`remotezip` 0.12.3, MIT) — Zenodo serves `accept-ranges: bytes` (HEAD 403, GET OK).
- **Raw files (verbatim, tracked):**
  - `metabric_markerStackOrder_RAW.csv` — member `MBTMEIMCPublic/markerStackOrder.csv`;
    cols `Isotope,Epitope`. sha256 `af8ed2f3936e36af0749a09b1a3fb56c8e30dd9e6a57e4e13f26e4ec79b25f7a`.
  - `metabric_AbPanel_RAW.csv` — member `MBTMEIMCPublic/AbPanel.csv`; full wet-lab antibody sheet
    (target, clone, metal tag). sha256 `ea86c1da04667bad706e368ef8583d5b9386b17d4048f6b398d5c08aacaee02a`.
- **Extracted list (verbatim, raw):** `metabric_markers.txt` — the `Epitope` column of
  `markerStackOrder.csv`, the canonical image-layer panel the published single-cell data columns
  correspond to. **39 channels = 37 protein markers + DNA1/DNA2** (Ir191/Ir193 intercalator).
  This is the authoritative "37-marker" panel; normalization to gene symbols happens at Step 0.
- **Discrepancies to handle at Step 0 (recorded, not silently cleaned):**
  - Channel `Ho165`: `markerStackOrder.csv` Epitope = **`ER`** (estrogen receptor, the analysis
    label); `AbPanel.csv` target = `Rabbit IgG (H+L)` (the wet-lab reagent row). Use the
    `markerStackOrder` Epitope (`ER` → `ESR1`); the AbPanel label is not a biological target.
  - Channel `Yb176`: Epitope verbatim = **`c-Caspase3c-PARP`** (a run-together of cleaved
    Caspase-3 / cleaved PARP; AbPanel target = `Cleaved Caspase3`). Apoptosis marker; no single
    clean gene — likely panel-distal-irrelevant, decide at Step 0.
  - `DNA1`/`DNA2` are nuclear intercalators, not antibody targets → not gene-mappable (drop).
  - Two HER2 channels (`Eu151` `HER2 (3B5)`, `Tb159` `HER2 (D8F12)`) both → `ERBB2`.
- **Note:** the files originally in `~/Downloads` presumed to be METABRIC-IMC were a
  misidentification — they are **Bassez-2021 BIOKEY scRNA-seq** (RNA gene × cell `dgCMatrix`) and
  a 10x scRNA-seq matrix, containing no IMC panel. See DECISIONS.md (2026-06-15). They are NOT the
  source used here; the authoritative panel above came from the Zenodo deposit.
