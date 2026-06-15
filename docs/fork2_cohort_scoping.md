# Fork-2 discovery-cohort scoping — PENDING APPROVAL

Status: **report only — no cohort finalized, nothing built.** Recommendation at the bottom is
PENDING the user's go. Compiled 2026-06-15 by parallel source-backed investigation.

**Goal context.** Fork 2 = discover a TME biomarker in ONE public spatial-proteomic breast cohort
with linked outcome, then validate in METABRIC-IMC (Danenberg 2022, 37-marker IMC). Cross-cohort
validation runs over the markers **shared between the discovery cohort and METABRIC**, so marker
overlap + platform comparability (batch canary, apparatus rule 5) matter as much as accessibility.

**HARD RULE applied.** Every panel fact below was taken from a real file / official source (cited);
items that could not be verified from a source are flagged "unverified / behind wall." METABRIC-IMC
37-marker panel used for overlap: see `data/panels/metabric_markers.txt`.

---

## Comparison

| Cohort | Platform | Access | License | n (outcome) | Outcome | Tissue | Overlap w/ METABRIC-37 | Overlap character |
|---|---|---|---|---|---|---|---|---|
| **Jackson/Basel 2020** | **IMC (Hyperion)** — same as METABRIC | **Direct, open** (Zenodo 3518284; 36.8 GB) | **CC-BY 1.0** | 281 Basel (+72 Zurich) | **OS, DFS**, grade, TNM, ER/PR/HER2 | Invasive, all subtypes | **15** | epithelial-biased; incl. **ER + HER2**, keratins, CD31/vWF, Ki-67; immune-thin (CD3/CD20/CD68 only) |
| **Keren 2018** | MIBI-TOF (cross-platform vs METABRIC) | Images **gated** (free reg, mibi-share); single-cell `cellData.csv` **open mirror** | **Unverified** | 41 TNBC (38 in reanalyses) | **OS + recurrence** | TNBC only | **16** (+1 soft H3) | **immune-rich** (CD3/4/8/16/45RO/68/163, FOXP3, HLA-DR/ABC) + SMA/CD31/panCK/Ki-67; **no ER, no HER2** |
| **Risom 2022** | MIBI-TOF (cross-platform) | **Direct, open** (Mendeley d87vg86zd8.2) | **CC-BY 4.0** | 79 tissues (14 progressor / 44 non-prog DCIS) | **DCIS→invasive progression/recurrence** | DCIS (pre-invasive) | **16 confirmed** (panel partly behind wall) | mixed: immune (CD3/4/8/11c/68, PD-1, HLA-DR) + **ER + HER2** + CK5/panCK/SMA/CD31/Ki-67/H3 |
| OHSU CycIF 2025 | **CycIF (fluorescence)** — NOT IMC/MIBI/CODEX | Synapse syn50134757 (light reg wall) | CC-BY 4.0 (article) | 102 (OS/RFS) | OS, RFS | Invasive, all subtypes | **not computed** (42-marker panel behind supplement; HARD-RULE: not guessed) | likely substantial but unverified |

**Excluded (with reason):**
- **Ali 2020 (Nat Cancer, METABRIC IMC):** same METABRIC cohort + same 37-marker IMC panel as the
  validation set's lineage → using it for discovery would **leak into validation**. Not independent.
- **Wang/Hartmann spatial proteomics:** LC-MS region proteomics, no antibody panel — wrong protein
  class, doesn't bridge to METABRIC-IMC.

---

## Per-cohort notes

### Jackson/Basel IMC 2020 (Nature; Bodenmiller)
- Source: Zenodo **10.5281/zenodo.3518284** (https://zenodo.org/records/3518284) — open, CC-BY 1.0;
  includes single-cell data, metadata, and the panel file `Basel_Zuri_StainingPanel.csv`. Also via
  Bioconductor `imcdatasets::JacksonFischer_2020_BreastCancer`. Paper PDF paywalled; **data is not**.
- Panel verified from Supp Table 2 (`41586_2019_1876_MOESM1_ESM.pdf`): 35 measured biomarkers,
  tumor/epithelial-focused (CK5/7/8-18/14/panCK, ER, PR, HER2, EGFR, GATA3, E/P-cadherin, p53,
  c-Myc, Twist, Slug, CAIX, CD44, vimentin, fibronectin, SMA, CD31, vWF, CD3, CD20, CD68, CD45,
  Ki-67, pS6/pmTOR, cleaved-Casp3/PARP, Histone H3).
- **Shared with METABRIC (15):** Histone H3, SMA, CK5, CK8-18, CD68, HER2, CD3, CD20, ER, Ki-67,
  CD31, vWF, panCK, cleaved-Caspase3, cleaved-PARP. (CD45 is pan-CD45 vs METABRIC's CD45RA/RO
  isoforms — excluded as not the same epitope.)

### Keren 2018 MIBI-TNBC (Cell; Angelo)
- Raw images gated behind free registration at mibi-share.ionpath.com (the wall = account creation).
  Processed single-cell `cellData.csv` is openly mirrored (panel read from its header). **License
  unverified** — no explicit data-license text found; confirm before any redistribution.
- **Shared with METABRIC (16):** SMA, HLA-DR, CD68, CD163, CD3, CD16, CD45RO, FOXP3, CD20, CD8,
  Ki-67, CD4, CD31, HLA-ABC, panCK, CD45 (+ soft Histone H3 ≈ H3K9ac/H3K27me3). **No HER2, no ER.**

### Risom 2022 MIBI DCIS (Cell; Angelo)
- Open on Mendeley (CC-BY 4.0). 32 of 37 markers verified verbatim from open full text; the full
  Table S2 / Fig 2B sit behind the Cell WAF → **~5 panel slots unverified** (pull Table S2 from the
  Mendeley deposit to finalize). Outcome is DCIS→invasive **progression**, not survival.

---

## Recommendation (PENDING APPROVAL — reasoning, not a decision)

**Primary recommendation: Jackson/Basel IMC 2020.** Rationale, in priority order:
1. **Same platform as the validation cohort (IMC/Hyperion).** METABRIC-IMC is IMC; Basel is IMC. A
   discovery→validation stitch within one platform is far less batch-confounded than a MIBI→IMC
   cross-platform stitch (apparatus rule 5 — the batch canary is much more likely to pass). Keren
   and Risom are MIBI; OHSU is CycIF — all cross-platform vs METABRIC.
2. **Fully open, no wall** (Zenodo, CC-BY 1.0) — best accessibility; license documented.
3. **Largest outcome-linked invasive cohort** (281 + 72) with OS/DFS — matches METABRIC's invasive,
   survival-linked design (Keren is TNBC-only n=41; Risom is pre-invasive DCIS).
4. **Contains the tumor-intrinsic anchors** ER and HER2 — so **ESR1 (ER), the panel-distal priority
   watch-gene, is actually measured in both cohorts**, and ERBB2 too.

**The honest caveat:** Basel↔METABRIC overlap (15) is **epithelial-biased and immune-thin** (only
CD3/CD20/CD68 on the immune side). If the intended biomarker thesis is **immune/stromal**, the
shared space is richer in **Keren** (16, immune-heavy) or **Risom** (mixed) — but both cost a
cross-platform MIBI→IMC stitch, and Keren adds TNBC-only scope + n=41 + unverified license, Risom
adds DCIS-not-invasive + n=79. That is a **biology-vs-comparability tradeoff for the user to set**:
- tumor-intrinsic / receptor / epithelial biomarker, clean platform match → **Basel**.
- immune-microenvironment biomarker, accept cross-platform → **Keren** (TNBC) or **Risom** (DCIS).

## Verification gaps to close before finalizing (whichever is chosen)
- Confirm **Basel vs METABRIC patient non-overlap** (Basel/Zurich hospitals vs METABRIC UK/Canada —
  expected distinct, but verify, to avoid validation leakage).
- Pull `Basel_Zuri_StainingPanel.csv` from Zenodo and re-derive the discovery↔METABRIC shared set
  byte-exactly (resolve the pan-CD45 vs CD45RA/RO question).
- If Risom is chosen, pull Table S2 (Mendeley) to lock the full 37-marker panel.
- Re-intersect the chosen cohort's panel with METABRIC using the same authoritative method as Step 0.

Sources: Zenodo 3518284; mibi-share.ionpath.com + Keren cellData.csv mirror; Mendeley
10.17632/d87vg86zd8.2; PMC8792442; nature.com/articles/s41586-019-1876-x; cell.com S0092-8674(18)31100-0.
