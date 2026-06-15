# Fork-2 cohorts — LOCKED (2026-06-15)

> **DECISION LOCKED (do not relitigate; see DECISIONS.md):**
> - **Discovery = Meyer 2025** (TNBC IMC, USZ/Zurich, n=215, survival endpoint).
> - **Validation = METABRIC, TNBC/basal subset** (subtype-matched).
> - **Working feature space = Meyer ∩ METABRIC = 20 shared protein markers** (see §2.1).
> - Chosen for: same-platform-as-METABRIC (cleanest batch canary), shared space carries
>   **CD8 + ER + HER2**, population-independent. **FOXP3 deliberately NOT required** → Treg/
>   exhaustion-subtype thesis is out of scope; scope = **CD8/cytotoxic-architecture + receptor-context**.
>   Keren & Engelhardt set aside.

The scoping that led here is retained below for the record. This pass also (a) resolved Meyer's
clone-walled cytokeratins from the study's own deposited SCE data → final byte-exact overlap, and
(b) pinned the METABRIC validation-subset size (§4). Revised 2026-06-15.

**Protein-marker, not gene.** Fork 2 has no RNA. Overlap = "which antibody targets are measured in
BOTH the discovery cohort and METABRIC-IMC," canonicalized at the marker/epitope level
(CD8a=CD8, ER=estrogen receptor, HER2=c-erbB-2, panCK=Pan-Keratin; combined channels expand:
CK8-18→{CK8,CK18}, CD31-vWF→{CD31,vWF}, c-Caspase3c-PARP→{CASP3,PARP}). The full override +
compartment table is in `scripts/fork2_marker_overlap.py`; the run is `results/fork2_marker_overlap.md`.
Pan-CD45 is NOT matched to METABRIC's CD45RA/CD45RO isoforms (different epitope). HARD RULE: every
panel below is parsed from a real source file (see `data/panels/PROVENANCE.md`); Meyer's
clone-truncated cytokeratins are flagged provisional, not counted strict.

---

## 1. Eligibility gates (must pass BOTH; run BEFORE ranking)

- **GATE 1 — endpoint match** to METABRIC OS/DFS axis (not a different axis like DCIS progression).
- **GATE 2 — population independence** (no shared patients/accrual pipeline; shared antibody lineage
  is fine).

| Cohort | Platform | GATE 1 | GATE 2 | Eligible? |
|---|---|---|---|---|
| **Jackson/Basel 2020** | IMC | **PASS** — `OSmonth` **and** `DFSmonth` columns present in deposited metadata | **PASS** — Univ. Hospitals Basel + Zurich (CH) | **YES** |
| **Keren 2018** | MIBI | **PASS** — OS + recurrence | **PASS** — Stanford, USA | **YES** (TNBC-only) |
| **Meyer 2025 (Bodenmiller)** | IMC | **PASS** — OS/DFS + recurrence (`OS_data`/`DFS_data`) | **PASS** — Univ. Hospital **Zurich** (USZ, accrual 2005–2017); ZTMA174/249 = Zurich TMAs; independent of METABRIC | **YES** (walls now resolved; CK identities clone-walled) |
| **Engelhardt/Chang CycIF** | CycIF | **PASS** — OS + RFS | **PASS** — Vanderbilt UMC, USA | **YES** (platform = CycIF, not IMC/MIBI/CODEX) |
| **Risom 2022 (DCIS)** | MIBI | **FAIL** — only DCIS→invasive progression / iBE; no OS/DFS | (pass) | **NO — Gate 1** |
| **Ali 2020** | IMC | pass | **FAIL** — IS the METABRIC cohort | **NO — Gate 2** |
| Wang 2023 | reuse | — | **FAIL** — method on Basel/METABRIC | **NO** |

**Four eligible survivors:** Jackson/Basel, Keren, Meyer-2025, Engelhardt/Chang. (Meyer's two prior
walls are now closed: origin = USZ Zurich, independent of METABRIC; panel obtained from the deposited
raw IMC headers — only the exact cytokeratin clone identities remain behind the paywalled Methods.)

---

## 2. Compartment-split shared protein markers (the decision substrate)

METABRIC-IMC validation panel = 37 markers → 39 canonical tokens.

| Cohort | Platform vs METABRIC | n (outcome) | Scope | **Shared** | immune | stromal | epi-tumor | CD8 | FOXP3 | ER | HER2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Jackson/Basel** | **same (IMC)** | 281+72, OS+DFS | all subtypes | **15** | 3 | 3 | 6 | – | – | – | **✔** |
| **Meyer 2025** (LOCKED) | **same (IMC)** | 215, OS/DFS+recur | **TNBC** | **20** | 8 | 2 | 7 | **✔** | – | **✔** | **✔** |
| **Keren** | cross (MIBI) | 41, OS+recur | **TNBC** | **18** | 14 | 2 | 2 | **✔** | **✔** | – | – |
| **Engelhardt** | cross (CycIF) | 102, OS+RFS | all subtypes | **15** | 7 | 3 | 5 | **✔** | **✔** | **✔** | **✔** |

Exact shared markers per compartment (from the real panel files):

- **Jackson/Basel (15):** immune = CD3, CD20, CD68 · stromal = SMA, CD31, vWF · epithelial-tumor =
  CK5, CK8, CK18, HER2, Ki-67, panCK · other = Histone H3, cleaved-Casp3, cleaved-PARP.
- **Meyer 2025 (20) — LOCKED feature space:** immune = CD3, CD4, CD8, CD11c, CD15, CD20, CD68,
  HLA-DR · stromal = SMA, vWF · epithelial-tumor = CK5, CK8, CK18, panCK, ER, HER2, Ki-67 ·
  other = Histone H3, cleaved-Casp3, cleaved-PARP. (CK identities resolved from the study's own
  deposited SCE `clean_target`, Zenodo 15304181: Pr141=CK5, Nd144=CK8/18, Lu175=panCK; CK7/CK14 also
  on the panel but not in METABRIC. Vendor clone names remain paywalled; overlap is by target.)

### 2.1 The LOCKED working feature space (Meyer ∩ METABRIC = 20 markers)
- **immune (8):** CD3, CD4, CD8, CD11c, CD15, CD20, CD68, HLA-DR
- **stromal (2):** SMA, vWF   *(note: no CD31 — Meyer carries vWF, not CD31; no FOXP3 by design)*
- **epithelial-tumor (7):** CK5, CK8, CK18, panCK, ER, HER2, Ki-67
- **other/functional (3):** Histone H3, cleaved-Caspase3, cleaved-PARP
This 20-marker set is what every Fork-2 gate and the cross-cohort validation will run on.
- **Keren (18):** immune = CD3, CD4, CD8, CD11c, CD16, CD20, CD45RO, CD68, CD163, FOXP3, HLA-DR,
  HLA-ABC, OX40, PD-1 · stromal = SMA, CD31 · epithelial-tumor = Ki-67, panCK.
- **Engelhardt (15):** immune = CD3, CD4, CD8, CD20, CD68, FOXP3, PD-1 · stromal = SMA, CD31,
  Podoplanin · epithelial-tumor = CK5, CK8, ER, HER2, Ki-67.

### Corrections this pass forced (why byte-exact panels matter)
- **Basel does NOT measure ER.** The deposited `Basel_Zuri_StainingPanel.csv` has **Rabbit IgG
  (control) at Gd156, no Estrogen Receptor** — only Progesterone Receptor + HER2 among receptors.
  The earlier "Basel has ER+HER2" (from a Supp-Table reading) is **wrong per the actual file**. So
  the ER/ESR1 axis is **not available in Basel**.
- **Keren overlap is 18, not 16** — byte-exact canonicalization caught OX40 and CD11c too.
- **Meyer has CD8 + ER + HER2 but NOT FOXP3** (no Treg marker in its deposited panel).

---

## 3. The tradeoff (explicit; NO recommendation — thesis axis is open)

Read the axis off the lists:

- **Immune / T-cell thesis** → **Keren** is richest (14 immune incl. **CD8, FOXP3, PD-1, OX40**,
  CD11c, CD16, CD163, HLA-DR/ABC) — but cross-platform MIBI, TNBC-only, n=41, license unverified.
  **Engelhardt** also has the cytotoxic+Treg pair (CD8+FOXP3+PD-1) and is all-subtype — but CycIF
  (cross-platform). **Meyer** has CD8 (+CD4/CD11c/CD15/CD68) **but no FOXP3**, same-platform IMC.
  **Basel** is unusable for an immune thesis (only CD3/CD20/CD68 lineage markers, no subsets).
- **Tumor-intrinsic / receptor / epithelial thesis** → needs ER and/or HER2 + keratins.
  **Engelhardt** has the full axis (ER, HER2, CK5, CK8, Ki-67) but CycIF. **Meyer** has ER+HER2+Ki-67
  (+ likely CKs) and is same-platform IMC, TNBC. **Basel** has HER2 + CK5/CK8-18/panCK + Ki-67 but
  **no ER**, all-comers, same-platform IMC. **Keren** is weakest here (only Ki-67 + panCK, no ER/HER2).
- **Stromal thesis** → thin everywhere (all have SMA; +CD31 except Meyer which has vWF not CD31;
  Engelhardt adds Podoplanin). Not a strong discriminator.
- **Platform comparability (batch canary)** → **Basel and Meyer are same-platform IMC** (lowest
  batch-confound vs METABRIC IMC); **Keren (MIBI)** and **Engelhardt (CycIF)** are cross-platform.
- **ESR1/ER watch-gene** → only **Engelhardt** and **Meyer** measure ER; Basel does not, Keren is TNBC.

Net shape of the choice: **same-platform + receptor axis → Meyer (TNBC) or Basel (no ER);
same-platform + immune → none strong (Basel/Meyer thin on T-cell subsets); immune-rich → Keren or
Engelhardt at a cross-platform cost.** Subtype confound (Keren/Meyer are TNBC; validate in METABRIC's
basal/TNBC subset) is handleable, not disqualifying.

---

## 4. METABRIC-IMC validation data status (Zenodo 6036188)

- **Clinical/survival table — IN HAND.** `IMCClinical.fst` (10 KB) range-extracted to
  `data/metabric/` (gitignored — patient-level). Schema: **709 patients × 11 cols** —
  `metabric_id, ERStatus, LymphNodesOrdinal, sizeOrdinal, Grade, ERBB2_pos, yearsToStatus,
  DeathBreast, isValidation, PAM50, IntClust`. Endpoint = **breast-cancer-specific survival**
  (`yearsToStatus` years + `DeathBreast` event). Events: **230 BC deaths / 479 censored**; follow-up
  median 8.69 yr (range 0.06–27.68); ER **542 pos / 151 neg / 16 NA**; METABRIC's own
  discovery/validation split `isValidation` = 532 / 177. **Endpoint is usable for survival validation.**
- **Single-cell expression — INVENTORIED, not pulled.** `SingleCells.csv` (849 MB) and
  `SingleCells.fst` (341 MB) inside the same zip; both range-extractable (Zenodo serves byte ranges).
  Pull when the gate sequence calls for it, not now.

### 4.1 Validation subset size + power (the cohort the design runs on)
From `IMCClinical.fst` (709 patients). No PR field → not strict triple-negative; the closest
subtype-match to Meyer's TNBC discovery is the **receptor-defined ER−/HER2− subset**:

| Definition | n | BC-death events | censored | follow-up median |
|---|---|---|---|---|
| **A. ER−/HER2−** (`ERStatus=neg & ERBB2_pos=FALSE`) — **primary** | **88** | **34** | 54 | 6.25 yr |
| B. PAM50 = Basal | 88 | 35 | 53 | — |
| C. ER−/HER2− AND PAM50 Basal (strictest) | 62 | 23 | 39 | — |

**Validation cohort = 88 patients, 34 breast-cancer-death events** (definition A). This is modest —
the Fork-2 gate sequence and any survival test must be designed for ~34 events (power-limited; favors
a small number of pre-specified marker hypotheses over wide scans). Definitions B/C are near-identical
/ stricter alternatives to pre-register. (METABRIC's own internal `isValidation` split within A is
56/32, but Fork 2 uses METABRIC purely as the external validation cohort, so all 88 are in play.)

## 5. Remaining walls / verification to-dos
- **Meyer:** exact cytokeratin identities (CK5 vs CK8-18 vs panCK) are clone-truncated in the deposit;
  confirm from the paywalled Key-Resources-Table before counting CK overlap as strict.
- **Engelhardt:** journal Supp-Table-1 (clones/vendors) behind a PMC anti-bot wall; the 42-marker
  list itself is verified from the authors' GitHub (`proteins_by_frame.csv`).
- **Basel:** pan-CD45 vs METABRIC CD45RA/RO left unmatched by design (epitope mismatch).

Marked **PENDING REVIEW** — thesis axis open; no cohort selected. Sources & sha256 in
`data/panels/PROVENANCE.md`; computation in `scripts/fork2_marker_overlap.py` →
`results/fork2_marker_overlap.md`.
