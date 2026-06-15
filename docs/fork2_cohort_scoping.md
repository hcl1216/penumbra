# Fork-2 discovery-cohort scoping — ELIGIBILITY-FIRST — PENDING APPROVAL

Status: **report only — no cohort finalized, nothing built.** Recommendation at the bottom is
PENDING the user's go. Revised 2026-06-15 to apply two **disqualifying eligibility gates BEFORE any
ranking** (the first pass wrongly ranked on panel overlap and deferred these). A candidate must
**PASS BOTH GATES** to be eligible; panel overlap + platform only break ties among survivors.

**Goal context.** Fork 2 = discover a TME biomarker in ONE public spatial-proteomic breast cohort
with linked outcome, then validate in METABRIC-IMC (Danenberg 2022). METABRIC's validation endpoints
are **overall survival (OS) / disease-free survival (DFS)**; its patients are UK (multiple sites) +
Canada (Vancouver) via the METABRIC consortium (Curtis 2012); the IMC subset is Ali 2020 / Danenberg
2022. HARD RULE applied: every outcome/accrual/panel fact below is from a real file or official
source (cited); unverifiable items are flagged as walls.

---

## The two gates (must pass BOTH)

- **GATE 1 — ENDPOINT MATCH (disqualifying).** The discovery cohort's outcome must be on the SAME
  axis as METABRIC's (OS/DFS survival). A different-axis endpoint (e.g. DCIS→invasive *progression*)
  cannot validate a survival biomarker.
- **GATE 2 — POPULATION INDEPENDENCE (disqualifying).** No shared patients and no shared
  selection/accrual pipeline with METABRIC. (Shared *antibody-panel lineage* is NOT leakage — only
  shared patients or selection pipeline disqualifies.)

## Eligibility table (gates applied to every candidate)

| Cohort | Platform | GATE 1 (endpoint) | GATE 2 (independence) | ELIGIBLE? |
|---|---|---|---|---|
| **Jackson/Basel 2020** | IMC | **PASS** — "patient overall survival"; Cox on "disease-specific overall survival" (Fig 4/4i) | **PASS** — Univ. Hospitals **Basel + Zurich (Switzerland)**, Swiss ethics 2014-397/2012-0553; METABRIC = citation only | **YES** |
| **Keren 2018** | MIBI | **PASS** — "recurrence and overall survival (OS)" | **PASS** — **Stanford Hospital, USA, 2002–2015**, single-institution | **YES** (TNBC-only caveat) |
| **Meyer 2025 (Bodenmiller TNBC)** | **IMC** | **PASS** — phenotype "correlated with rapid disease recurrence" (recurrence/DFS axis) | **likely PASS — WALL** — raw files "ZTMA174/ZTMA249" ⇒ Zurich TMAs (independent), but origin sentence is behind paywalled Methods; not confirmed | **PROVISIONAL** (2 walls) |
| **Engelhardt/Chang CycIF (VUMC/OHSU)** | **CycIF** (not IMC/MIBI/CODEX) | **PASS** — "overall survival (OS) and recurrence-free survival (RFS)" | **PASS** — **Vanderbilt UMC, USA** surgical TMAs | **YES** (platform caveat) |
| **Risom 2022 (DCIS)** | MIBI | **FAIL** — only "ipsilateral breast event (iBE)"/DCIS→IBC progression; no OS/DFS | PASS (US RAHBT/TBCRC, independent) | **NO — Gate 1** |
| **Ali 2020 (METABRIC-IMC)** | IMC | PASS (METABRIC OS/DFS) | **FAIL** — IS the METABRIC cohort (same patients) | **NO — Gate 2** |
| **Wang 2023** | (reuses Basel/METABRIC) | — | **FAIL** — method on existing METABRIC/Basel data, not independent | **NO — Gate 2** |
| Sorin 2023 | IMC | — | — | **NO** — lung, not breast |

**Net: four survivors** — Jackson/Basel, Keren, Meyer-2025 (provisional, walls), Engelhardt/Chang
(platform caveat). Risom is **disqualified** (the deferred gate now decides it). Ali-2020 is the
clean example of a Gate-2 failure (it IS METABRIC).

---

## Ranking the survivors (tie-break only: platform comparability + panel overlap)

Platform comparability matters because discovery→METABRIC is a cross-cohort stitch and the batch
canary (apparatus rule 5) is far more likely to pass IMC→IMC than across platforms.

| Survivor | Platform vs METABRIC (IMC) | n (outcome-linked) | Scope | METABRIC-37 overlap | Notes |
|---|---|---|---|---|---|
| **Jackson/Basel** | **same (IMC)** | 281 Basel (+72 Zurich), OS | all subtypes (invasive) | **15** (verified) — epithelial-biased, incl. **ER + HER2** | fully open, CC-BY 1.0; immune-thin overlap |
| **Meyer 2025** | **same (IMC)** | 215, recurrence | **TNBC only** | **unverified (panel wall)** — same lab as Danenberg ⇒ likely high | open Zenodo CC-BY-4.0; origin + panel both walls |
| **Keren** | cross (MIBI) | 41 TNBC, OS+recurrence | **TNBC only** | **16** (verified) — immune-rich, **no ER/HER2** | images gated (free reg); license unverified |
| **Engelhardt/Chang** | cross (CycIF) | 102, OS+RFS | all subtypes | **unverified (panel wall)** | open (Synapse, light reg), CC-BY-4.0; platform not IMC/MIBI/CODEX |

Overlap detail for the two verified panels (from real panel files):
- **Basel↔METABRIC (15):** Histone H3, SMA, CK5, CK8-18, CD68, HER2, CD3, CD20, ER, Ki-67, CD31,
  vWF, panCK, cleaved-Caspase3, cleaved-PARP. (pan-CD45 vs METABRIC CD45RA/RO excluded as different
  epitopes.)
- **Keren↔METABRIC (16):** SMA, HLA-DR, CD68, CD163, CD3, CD16, CD45RO, FOXP3, CD20, CD8, Ki-67,
  CD4, CD31, HLA-ABC, panCK, CD45 (+ soft Histone H3).

---

## Recommendation (PENDING APPROVAL — reasoning, not a decision)

**Primary: Jackson/Basel IMC 2020.** It is the only survivor that is (a) clean on **both gates with
no walls** (OS confirmed; Swiss, independently accrued), (b) **same platform as METABRIC** (lowest
batch-confound risk), (c) the **largest outcome-linked all-comers invasive cohort** (matching
METABRIC's invasive, survival-linked design), and (d) measures **ER + HER2**, so the panel-distal
watch-gene **ESR1** and the tumor anchor **ERBB2** exist in both. Caveat: Basel↔METABRIC overlap (15)
is epithelial-biased / immune-thin.

**Watch closely — Meyer 2025 (Bodenmiller TNBC IMC).** Same platform, larger TNBC discovery n (215),
recurrence endpoint, openly deposited. If its two walls clear favourably — **(i)** patient origin
confirmed independent of METABRIC (likely Zurich), **(ii)** antibody panel pulled from the deposited
IMC files — it could rival or beat Basel for a TNBC/immune thesis (same-platform AND likely high
overlap). Not recommendable until those walls are closed.

**Immune-thesis alternative — Keren.** If the biomarker thesis is immune/stromal, Keren's overlap is
richer (16, immune-heavy) — at the cost of a cross-platform MIBI→IMC stitch, TNBC-only scope, n=41,
and an unverified license.

**Engelhardt/Chang** is eligible and has good outcome data (OS+RFS, n=102, all subtypes) but is
**CycIF, not IMC/MIBI/CODEX** — a different multiplex class; consider only if that platform is
acceptable, and pull its 42-marker panel first.

### Keren TNBC-scope caveat (calibration, not disqualifying)
Discovering in triple-negative (Keren, Meyer) and validating in all-comers METABRIC is a **subtype
confound**. It is handleable — validate in METABRIC's **basal/TNBC subset** rather than the whole
cohort — so it does not disqualify, but it shrinks the effective validation n and must be declared up
front.

---

## SECONDARY — is the METABRIC-IMC validation data actually in hand? **NO (panels only).**
Local `data/panels/` holds only the panel CSVs (CosMx + METABRIC marker lists). A filesystem check
found **no SingleCells / clinical / survival tables** anywhere local. The full validation substrate —
single-cell expression + per-patient clinical/survival — lives **only inside the Zenodo 6036188 zip
`MBTMEStrIMCPublic.zip` (6.65 GB)**, not yet pulled. **Fork-2 validation cannot run until those cell
+ clinical tables are obtained** (selective extraction via range requests is feasible, as done for
the panel files). Reported, not downloaded.

## Verification gaps / walls to close before finalizing
- **Basel:** confirm whether a DFS column exists in the Zenodo metadata CSV (paper analyzes OS /
  disease-specific OS, not DFS); re-derive Basel↔METABRIC shared set byte-exactly from
  `Basel_Zuri_StainingPanel.csv`; resolve pan-CD45 vs CD45RA/RO.
- **Meyer 2025:** read patient origin from the STAR Methods (confirm METABRIC-independence) and pull
  the antibody panel from the deposited IMC files (Zenodo 10890543) — both currently walls.
- **Engelhardt/Chang:** pull the 42-marker panel from Suppl. Table 1 / Synapse before computing
  overlap.
- **METABRIC-IMC:** pull SingleCells + clinical/survival from Zenodo 6036188 (validation prerequisite).

Sources: Zenodo 3518284 (Basel) · zenodo 10890543 + cell.com/cancer-cell S1535-6108(25)00269-7
(Meyer) · PMC8271023 + cell.com S0092-8674(18)31100-0 (Keren) · PMC9772081 + data.mendeley.com
d87vg86zd8 (Risom) · PMC11948582 + Synapse syn50134757 (Engelhardt/Chang) · repository.cam.ac.uk
Basel accepted manuscript · Zenodo 6036188 (METABRIC-IMC validation).
