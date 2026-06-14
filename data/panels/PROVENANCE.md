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

## 3. METABRIC / Ali-2020 37-marker IMC panel  ❌ NOT IN HAND — BLOCKED

- **Assumption falsified.** The local files in `~/Downloads` previously presumed to be
  Danenberg-2022 METABRIC-IMC are **not** METABRIC-IMC and contain **no IMC protein panel**:
  - `1863-counts_cells_cohort1.rds`, `1864-counts_tcell_cohort1.rds`,
    `1867-counts_cells_cohort2.rds`: sparse `dgCMatrix` of **scRNA-seq counts**, rows =
    25,288/22,889 **RNA gene symbols**, columns = cells named `BIOKEY_##_Pre_…` →
    **Bassez et al. 2021 BIOKEY** anti-PD1 breast scRNA-seq cohort (Nat Med), not Danenberg.
  - `2102-Breastcancer_counts.tar.gz`: 10x `matrix.mtx`/`genes.tsv`/`barcodes.tsv`, 33,694
    **RNA genes** × 44,024 cells (`sc5r…` barcodes) — scRNA-seq, not IMC.
- **Tooling note:** `pyreadr` 0.5.0 (**AGPLv3** — acceptable for personal, non-distributed use;
  recorded) cannot read these (`dgCMatrix` is S4). Read instead via local R 4.3.1
  (`C:\Program Files\R\R-4.3.1`) `readRDS()` → dimnames, which is how the misidentification
  above was established.
- **Consequence:** no authoritative METABRIC 37-marker source is available locally, and the
  HARD RULE forbids reconstructing it from memory. **Step 0 cannot run** (needs all 3 name-lists).
  Need either the real Danenberg-2022 IMC artifact (e.g. the `SingleCells.csv` / IMC `.rds` with
  37 marker channels, Zenodo) or the Ali-2020 marker table from an official source.
