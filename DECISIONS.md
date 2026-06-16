# Penumbra -- Decision Log

Persistent, reverse-chronological log (newest on top). Append an entry whenever a gate
produces a verdict/result, a kill criterion fires, we pivot/rescope/abandon, a load-bearing
assumption is confirmed or falsified, or a design decision is made that future sessions must
not relitigate. See "Logging discipline" in CLAUDE.md. Every RESULT/VERDICT entry must also
update the STATUS line in CLAUDE.md. Keep entries short.

Entry format:
`### YYYY-MM-DD - TYPE - one-line summary`  (TYPE = RESULT | VERDICT | PIVOT | DECISION)
then: **Numbers:** (if any) - **Consequence:** (what it changes / what we do next).

---

### 2026-06-16 - DECISION - Cell-typing refined after R2 review (priority + doublet quarantine + higher k)
**Change to scripts/fork2_01_typing_floors.R (logic only; design unchanged):** metacluster k 12 -> 25;
annotation rule rewritten to the design priority **Tumor > B > Myeloid > CD8 > CD4 > Stroma > Other**
with a **doublet quarantine** (>=2 of the mutually-exclusive anchors CD3/CD20/panCK/CD68 at >=0.5 -> Other,
fixes the all-high cluster mc12) and CD4 tested LAST among leukocytes (it is promiscuous, was swallowing
B cells mc11/mc12 and myeloid mc8/mc10). Thresholds raised (defining marker >=0.35). Re-run on Colab
(notebook does git pull) and re-review R2 + re-confirm R3 before Gate C. Floors (R4) and guard (R1)
unchanged in logic.

### 2026-06-16 - RESULT - Front end ran on Colab (R1-R4); Gate A/B OK, but cell-typing too coarse -> Gate C NOT ready
**R1 (markers):** 20 -> 19 (Histone H3 absent in Meyer) -> 17 channels; HARD GUARD passed (CD8 + cytokeratins present). Clean.
**R2 (typing) — PROBLEM:** 12 metaclusters collapsed to only 4 types. CD8_T 5.6%/4.2% (Meyer/METABRIC), Tumor_epithelial 24.6%/35.2% (both plausible), CD4_T 24.9%/16.2%, **Other 44.9%/44.4% (~half unclassified)**. The CD4_T bucket is a GARBAGE DRAWER: it absorbs B cells (mc11 CD20=0.88, mc12 CD20=1.0) and myeloid (mc8 CD15=1.0/CD68=0.69, mc10 CD68=0.83/CD11c=0.70/HLA-DR=0.80) because the annotation rule tests CD4>=0.2 before CD20/CD68 and CD4 is promiscuous. **No B_cell / Myeloid / Stroma types emerge.** Tumor(CK+) and CD8(mc9, CD8=0.79) masks look usable, but typing must be fixed (resolution + rule priority) before it underpins the proximity feature.
**R3 (Gate B canary):** PASS (provisional). Per-type cross-cohort profile corr: CD8 0.93, CD4 0.90, Tumor 0.77, Other 0.63; stable-quantity divergence 0.166; CD8-frac cross-patient SD 0.0575 vs cross-cohort shift 0.0125 -> **margin 4.60** (batch below signal axis). Robust on the CD8/tumor axes that matter; re-confirm after typing fix.
**R4 (Gate A floors, Meyer DFS, n=152/40 events; dropped from 215 for missing SoC covariates):** composition floor CD8-fraction coef -1.18, **p=0.69 (FLAT / non-prognostic)** -- this is the LOW BAR the spatial feature must beat (good for the thesis). SoC composite: grade p=0.26 (n.s.), **nodes dominate (pN2 HR 6.2 p=1e-4; pN3 HR 8.6 p<1e-4)**; transferable as grade+stage.
**VERDICT:** Gate A + Gate B informative and provisionally pass; **Gate C BLOCKED on cell-typing refinement** (fix R2: prioritize Tumor>B>Myeloid>CD8>CD4>Stroma, raise metacluster k, handle the all-high doublet cluster, reduce ~45% Other). Artifacts in Drive results/fork2_R1-R4. STOP for review before re-running + Gate C.

### 2026-06-16 - RESULT - Colab runner + data-to-Drive manifest ready (front-end runs on Colab)
**Built:** thin Colab notebook `notebooks/fork2_01_colab_runner.ipynb` (mount Drive → clone github.com/
hcl1216/penumbra → BiocManager install SingleCellExperiment/FlowSOM/imcRtools/cytomapper + fst/survival
→ `Rscript scripts/fork2_01_typing_floors.R` → print R1–R4). Logic stays in the git `.R` script (not
reimplemented in the notebook). Script made Colab-ready (minimal edits, still runs locally): input/
output paths from env vars (MEYER_SCE / METABRIC_SC / PENUMBRA_RESULTS_DIR) defaulting to repo layout.
**Script now emits 4 review artifacts (stdout + Drive results/):** R1 marker-status (resolves locked-20
→ 19 drop Histone H3 absent in Meyer → 17 analysis channels: CK8&CK18 share CK8-18, cl-Casp3&cl-PARP
share one, vWF via CD31-vWF); R2 cluster→type + Meyer% vs METABRIC% side-by-side; R3 Gate-B canary with
numbers + margin (stable per-type profile divergence vs CD8-fraction cross-patient SD / cross-cohort
shift); R4 Gate-A Meyer-fit Cox floors (CD8-fraction composition + grade+stage SoC composite, labeled
"additive over grade+stage" since METABRIC has no age). **HARD GUARD:** halts if CD8 or all cytokeratins
drop from the harmonized set (primary CD8↔CK+ proximity would be impossible). **HARD STOP after Gate A**
(no Gate-C proximity, no hand-rolled clustering). FlowSOM SOM trained on 200k balanced subsample then
mapped to all ~1.21M cells (RAM-safe; High-RAM Colab recommended).
**Manual next steps (user):** (1) move 3 objects to Drive per manifest; (2) push repo to GitHub
(no remote configured locally — needed before the notebook can clone); (3) run notebook; (4) paste
R1–R4 back for review BEFORE Gate C. FLAG: repo public/private unverified → notebook supports both.

### 2026-06-16 - RESULT - PART 3 build BLOCKED locally (Bioconductor unreachable); front-end script ready for Colab
**Blocker (environment):** the mandated IMC cell-typing ecosystem can't be installed on the local
Windows box. **bioconductor.org is UNREACHABLE** (curl HTTP 000 / 25s timeout; R available.packages
hangs) while CRAN/GitHub/Zenodo work → FlowSOM/flowCore/imcRtools/cytomapper (Bioconductor-only)
uninstallable. The design's other allowed clusterer, **Rphenograph (GitHub), needs a C++ compiler and
Rtools is ABSENT** (`has_build_tools=FALSE`). Design forbids hand-rolling clustering → did NOT
substitute.
**Consequence:** cell-typing → Gate B canary → Gate A floors (all depend on the typing) are NOT run
locally. Front-end script authored and committed (`scripts/fork2_01_typing_floors.R`,
harmonize→FlowSOM→canary→floors, stops before Gate C) — **run on Colab** (network-open, installs the
IMC stack, handles ~1.2M cells), per the Claude-Code/CPU vs Colab/heavy split. Harmonization decided:
Meyer `exprs` (arcsinh cofactor 1) + METABRIC arcsinh(cofactor 1) of raw intensity, per-cohort 99th-
pct 0-1 scaling; 17 shared analysis markers (20 minus Histone H3 absent in Meyer, with CD31-vWF and
c-Casp3/PARP each one combined channel). STOP for review; nothing hand-rolled.

### 2026-06-16 - RESULT - Fork-2 single-cell data pulled + VERIFIED; both pass (coords + linkage + survival)
**Meyer (Zenodo 15304181, sce_ALL_sub.rds):** SCE 36 markers × 92,899 cells, 215 patients; coords
Pos_X/Pos_Y ✓; 19/20 shared markers (Histone H3 absent — structural "other", irrelevant); endpoint
status_DFS+DFS_months (48 events/198; recurrence) and status_OS+OS_months (42 events); SoC age/grade/
pTNM. **METABRIC (Zenodo 6036188, SingleCells.fst):** 1,123,466 cells × 55 cols; coords
Location_Center_X/Y ✓; 20/20 markers; metabric_id joins IMCClinical (707/709 imaged); **TNBC subset
ER−/HER2− = 88 pts, 87 imaged, 34 events, 172,074 cells**; SoC Grade/size/nodes — **no age**.
**VERDICT: PASS — build may proceed.** Flags: (1) METABRIC has no age → transferable SoC composite =
grade+stage, not age (adjusts Gate-A floor); (2) endpoint difference Meyer recurrence vs METABRIC
BC-death (stated assumption); (3) METABRIC isTumour flags ~92% cells tumour → unreliable, use own
typing. Raw objects gitignored. Report: results/fork2_data_verify.md.

### 2026-06-15 - DECISION - Fork-2 Phase-0 gate sequence LOCKED; primary feature = CD8↔tumor proximity
**Binding constraint:** 34 validation events → ≤3 pre-specified spatial features (1 primary);
METABRIC TNBC subset (88/34) = underpowered independent REPLICATION (Meyer n=215 carries effect-
existence); framing = "discovered in Meyer, replicated additively in independent METABRIC", not
"p<0.001 in METABRIC"; every test = 1-feature-vs-1-scalar-floor, floors pre-collapsed to single risk
scores FIT on Meyer and applied as FIXED scalars to METABRIC.
**Role split:** Meyer = discovery (fit+lock features & floor coefficients); METABRIC TNBC = validation
(touched once; incremental spatial effect only).
**Features:** PRIMARY = CD8↔tumor proximity (per-patient CD8→nearest CK+/panCK+ tumor distance → one
scalar: median or infiltrated-fraction within locked radius; multiplicity-protected). Secondary-1 =
proximity consistency (uniform vs patchy). Secondary-2 = Ki67-tumor structure / CD8 vs proliferative
tumor. Secondaries exploratory + multiplicity-flagged.
**Gates:** A = floors (SoC composite grade/age + composition scalar, fit on Meyer) + within-patient
feature reliability. B = batch canary (harmonize 20 markers, ONE cross-cohort cell-typing; stable
quantity must vary LESS cross-cohort than the survival signal; batch≈signal⇒DEAD). C = lock features
on Meyer + position-permutation (permute positions within-patient → feature must collapse to
composition floor). D = METABRIC, touch once; ΔC-index (feature+floor vs floor) w/ bootstrap CI, not
bare p<0.05.
**Kill criteria (fixed):** dead if (i) batch≈signal; (ii) primary fails vs composition floor on Meyer
or collapses under permutation; (iii) no additive replication over both floors on METABRIC → write
the negative ("spatial immune architecture adds no prognostic value over composition + SoC in
independent TNBC validation"). Both outcomes publishable.
**Stated assumption:** discovery/validation endpoints differ (Meyer recurrence vs METABRIC BC-specific
death). Recorded in CLAUDE.md (Phase 0 plan, LOCKED). Build the front end (harmonize + cell-typing +
Gate B canary + Gate A floors), STOP before the Gate-C proximity feature.

### 2026-06-15 - RESULT - Meyer CKs resolved (feature space final = 20 markers); METABRIC validation subset pinned
**Task 1 — Meyer cytokeratins resolved** from the study's own deposited SCE `rowData/clean_target`
(Zenodo 15304181, sce_ALL_sub.rds; real source, not memory): Pr141=CK5, Nd144=CK8/18, Sm147=KRT14,
Yb174=CK7, Lu175=panCK (vendor clones still paywalled; overlap by target). Re-ran
scripts/fork2_marker_overlap.py → **Meyer ∩ METABRIC = 20 markers** (was 16+prov): immune 8
{CD3,CD4,CD8,CD11c,CD15,CD20,CD68,HLA-DR}, stromal 2 {SMA,vWF}, epithelial-tumor 7
{CK5,CK8,CK18,panCK,ER,HER2,Ki-67}, other 3 {Histone H3, cleaved-Casp3, cleaved-PARP}. **This is the
LOCKED working feature space.** (No CD31 — Meyer has vWF not CD31; no FOXP3 by design.)
**Task 2 — METABRIC validation subset** (from IMCClinical.fst; no PR field → ER−/HER2− is closest to
Meyer's TNBC): **primary def A ER−/HER2− (ERStatus=neg & ERBB2_pos=FALSE) = 88 patients, 34 BC-death
events**, median follow-up 6.25 yr. Alternatives: PAM50 Basal = 88/35; A∩Basal = 62/23. Power is
limited (~34 events) → gate design must favor few pre-specified hypotheses over wide scans.
**Consequence:** feature space + validation cohort size are now fixed; ready for Fork-2 gate-sequence
design. Single-cell expression still NOT pulled (range-extractable when needed).

### 2026-06-15 - DECISION - Fork-2 cohorts LOCKED: discovery = Meyer 2025 TNBC IMC; validation = METABRIC TNBC/basal subset
**Decision (do not relitigate):**
- **Discovery = Meyer 2025** (TNBC IMC, USZ/Zurich, n=215, survival endpoint).
- **Validation = METABRIC, restricted to the TNBC/basal subset** (subtype-matched to Meyer's TNBC).
- The **Meyer ∩ METABRIC shared marker space is the working feature set** for the whole project.
**Reasoning:** Meyer is same-platform as METABRIC (IMC→IMC, cleanest batch canary), its shared space
carries **CD8 + ER + HER2**, and it is population-independent of METABRIC (USZ Zurich vs UK/Canada).
**FOXP3 deliberately NOT required** → the Treg/exhaustion-subtype thesis is out of scope; the project
scope is **CD8/cytotoxic-architecture + receptor-context**. Keren and Engelhardt set aside.
**Consequence:** next pass = Fork-2 gate-sequence design on this locked feature set. Before that, two
closing items: resolve Meyer's clone-walled cytokeratins (finalize epithelial overlap) and pin the
METABRIC TNBC/basal validation-subset n + event count (statistical power).

### 2026-06-15 - RESULT - Compartment-split protein-marker overlap for all 4 survivors; Meyer walls resolved; METABRIC clinical pulled
**Numbers (shared antibody targets w/ METABRIC-37, byte-exact from real panel files; immune/stromal/
epi-tumor):** Basel 15 (3/3/6) CD8– FOXP3– ER– HER2+ ; Keren 18 (14/2/2) CD8+ FOXP3+ ER– HER2– ;
Engelhardt CycIF 15 (7/3/5) CD8+ FOXP3+ ER+ HER2+ ; Meyer 16 (8/2/3 +CK provisional) CD8+ FOXP3– ER+ HER2+.
**Key findings:**
- **CORRECTION:** Basel does NOT measure ER (deposited panel Gd156 = Rabbit IgG control; only PR+HER2).
  Earlier "Basel ER+HER2" was wrong; real file overrides. Basel Gate 1 confirmed OS+DFS (both columns).
- **Meyer-2025 walls resolved → ELIGIBLE:** origin = Univ. Hospital Zurich (USZ, 2005–2017), independent
  of METABRIC (Gate 2 PASS); panel obtained from deposited raw IMC headers (IMC, ~37 Ab; CK clone
  identities still walled → provisional). Has CD8+ER+HER2 but NO FOXP3.
- Keren overlap is 18 (byte-exact caught OX40, CD11c), immune-richest incl CD8+FOXP3+PD1+OX40.
**METABRIC validation data:** clinical IMCClinical.fst IN HAND (709 pts; endpoint yearsToStatus+
DeathBreast = BC-specific survival; 230 events; ER 542+/151−). Single-cell SingleCells.csv(849MB)/
.fst(341MB) INVENTORIED, range-extractable, NOT pulled.
**Consequence:** tradeoff surfaced in docs/fork2_cohort_scoping.md (no recommendation; thesis axis
open): immune thesis → Keren/Engelhardt (cross-platform) or Meyer (same-platform, no FOXP3);
receptor/epithelial thesis → Engelhardt/Meyer (have ER) or Basel (no ER); same-platform IMC → Basel
or Meyer. **PENDING REVIEW; no cohort picked.** New panels tracked in data/panels/; overlap script
scripts/fork2_marker_overlap.py.

### 2026-06-15 - RESULT - Eligibility gates applied to Fork-2 cohort candidates; Risom disqualified, 4 survivors
**Gates (both disqualifying, run BEFORE ranking):** G1 endpoint must match METABRIC OS/DFS axis;
G2 population independent of METABRIC (no shared patients/accrual pipeline).
**Outcomes (source-verified):**
- Jackson/Basel IMC — G1 PASS (overall survival; disease-specific OS Cox), G2 PASS (Swiss Basel+Zurich,
  independent) → **ELIGIBLE**.
- Keren MIBI — G1 PASS (OS + recurrence), G2 PASS (Stanford USA) → **ELIGIBLE** (TNBC-only caveat).
- Meyer 2025 Bodenmiller TNBC IMC (re-sweep find) — G1 PASS (recurrence), G2 likely PASS (Zurich ZTMA;
  WALL: origin behind paywalled Methods) → **PROVISIONAL** (origin + panel walls).
- Engelhardt/Chang CycIF VUMC/OHSU (re-sweep find) — G1 PASS (OS+RFS), G2 PASS (Vanderbilt) →
  **ELIGIBLE** (platform = CycIF, not IMC/MIBI/CODEX; panel wall).
- **Risom 2022 DCIS — G1 FAIL** (only DCIS→invasive progression / ipsilateral breast event; no OS/DFS)
  → **DISQUALIFIED** (G2 would have passed; the deferred gate decides it).
- Ali 2020 — G2 FAIL (IS the METABRIC cohort). Wang 2023 — G2 FAIL (method on Basel/METABRIC).
**Consequence:** report restructured eligibility-first in docs/fork2_cohort_scoping.md. Survivors
ranked by platform comparability + overlap; recommendation = **Jackson/Basel IMC** (only survivor clean
on both gates with no walls; same platform as METABRIC; largest invasive OS cohort; has ER+HER2 so
ESR1 watch-gene present) — Meyer-2025 a same-platform TNBC contender pending its two walls; Keren the
immune-rich cross-platform alternative. **Pick still PENDING APPROVAL — not finalized.**
**SECONDARY (load-bearing):** METABRIC-IMC validation data is NOT in hand — only panel CSVs are local;
SingleCells + clinical/survival live only in Zenodo 6036188 zip (6.65 GB), not pulled. Validation
cannot run until obtained.

### 2026-06-15 - DECISION - Fork-2 discovery-cohort scoping reported (recommendation PENDING APPROVAL; NOT finalized) (superseded by eligibility-gated RESULT above)
**Numbers:** candidates w/ METABRIC-37 marker overlap (source-verified): Jackson/Basel IMC 15
(epithelial-biased, incl. ER+HER2; same platform as METABRIC; open CC-BY 1.0; n=281+72; OS/DFS),
Keren MIBI-TNBC 16 (immune-rich, no ER/HER2; images gated, license unverified; n=41; OS+recurrence),
Risom MIBI-DCIS 16 (mixed; open CC-BY 4.0; n=79; progression). OHSU CycIF excluded on platform; Ali
2020 excluded (= METABRIC lineage, validation leakage).
**Consequence:** Report in `docs/fork2_cohort_scoping.md`. Recommendation = **Jackson/Basel IMC**
(same-platform IMC->IMC minimizes batch confound, fully open, largest outcome-linked invasive cohort,
has ER/HER2 so ESR1 watch-gene measured in both) — with the explicit caveat that its overlap is
immune-thin, so an immune-biomarker thesis would favor Keren/Risom at a cross-platform cost.
**NOT finalized** — awaiting user go. Verification gaps to close before commit (Basel↔METABRIC patient
non-overlap; re-derive shared set from Basel panel CSV; Risom Table S2). Do NOT build until cohort
approved + Fork-2 Phase-0 plan written.

### 2026-06-15 - PIVOT - Off Fork 3, onto Fork 2 (single-cohort biomarker discovery)
**Reason (verbatim):** the open paired CosMx breast multiomic substrate is confirmed inaccessible to
a solo researcher (protein lives behind AtoMx, no instrument/tenant access; the public breast deposit
is WTX RNA-only). Fork-3-as-written is blocked at the data layer. Fork 2 = discover a TME biomarker in
ONE spatial-proteomic cohort (IMC/CODEX/MIBI), validate in METABRIC-IMC. No imputation, no paired
data, no aligner — so the imputation/circularity threat is gone. The apparatus and the 22-protein
partition carry over unchanged (still protein space, METABRIC still the validation target).
**Consequence:** Fork-3 imputation machinery (protein->RNA premise, CosMx paired-slide dependency,
one-slide ceiling, build-vs-borrow imputer, Gates A-E) is superseded/archived, not deleted (a real
negative finding). Apparatus, logging discipline, infra/conventions, checkpoint-and-stop, and the
approved 22-protein partition carry over. New Fork-2 Phase-0 gate sequence is a TODO for the next
design pass (not invented here). **Note:** the WTX RNA-only public breast slide is NOT a Fork-2 input;
any in-progress download of it is NOT a Phase-0 dependency and can be abandoned.

### 2026-06-15 - DECISION - Step-0 panel partition APPROVED (final)
**Numbers:** 22 shared proteins; 23 panel-adjacent / 18,910 panel-distal; structural lock recorded
(commit a1dec02).
**Consequence:** The Step-0 partition is FINAL — do not relitigate. Under Fork 2 it is repurposed:
it defines the **protein feature space shared with METABRIC** for cross-cohort validation (no longer
an imputation adjacency partition). Carries over unchanged.

### 2026-06-15 - RESULT - Step 0 partition RESOLVED (review pass); structural lock recorded; approvable
**Numbers:** 22 shared protein channels -> 25 distinct cognate genes; **panel-adjacent = 23**,
**panel-distal = 18,910** (WTx 18,933). Changes vs first pass (+2 adjacent): HLA-DR expanded to
{HLA-DRA, HLA-DRB5} and CD31-vWF to {PECAM1, VWF}.
**Consequence:** Review items resolved from the two panel files (no memory reconstruction):
(1) **STRUCTURAL LOCK** recorded in CLAUDE.md + Step-0 output -- input feature set == adjacency-
defining set == the 22 shared proteins; polygenic/combined shared channels expand to ALL cognate
genes; an all-62 run is a separate upper bound with its own larger partition, NOT the binding gate.
(2) **panCK/CK8-18 keratins:** CosMx IO has NO keratin channel at all -> KRT5/8/18/19 in WTx but
CORRECTLY panel-distal (not a dropped mapping). HLA-ABC likewise not shared (no MHC-I in CosMx) ->
HLA-A/B/C distal. (3) c-Caspase3c-PARP split into CASP3+PARP1 (not shared -> distal). (4) FCGR3A/CD16
kept (shared but absent from WTx); Histone H3/DNA1/DNA2 kept dropped. Gate-D notes: adjacent set is
immune-/stromal-dominated, **ERBB2 the only tumor-intrinsic anchor**; **ESR1 (ER) flagged as panel-
distal priority watch-gene** (in METABRIC, not CosMx IO). Partition **approvable**; Gate A remains
blocked on the absent CosMx slide.

### 2026-06-15 - RESULT - Step 0 panel overlap produced; partition PENDING REVIEW (superseded by resolution above)
**Numbers:** CosMx IO 62 proteins (65 distinct genes) ∩ METABRIC 39 channels (36 gene-mappable,
3 dropped: Histone H3, DNA1, DNA2) = **22 shared channels -> 22 genes**. WTx = 18,933 genes.
Partition: **panel-adjacent = 21**, **panel-distal = 18,912**. (One shared gene, FCGR3A/CD16, is
absent from the WTx panel so cannot anchor adjacency -> 22 shared genes, 21 adjacent.)
**Consequence:** Gate-D substrate exists: 21 adjacent vs 18,912 distal genes. Panel-adjacent set =
{ACTA2, CD163, CD38, CD3D/E/G, CD4, CD68, CD8A, ERBB2, FOXP3, FUT4, HLA-DRA, ICOS, ITGAX, MKI67,
MS4A1, PDCD1, PECAM1, PTPRC, TNFRSF18}. Outputs in results/ (gitignored): protein_gene_map.tsv,
shared_proteins.tsv, panel_partition.tsv, SUMMARY.md. METABRIC epitope->gene used an explicit,
reviewable override table (CosMx side is vendor 'Gene' column). **PENDING REVIEW** -- partition not
final; ambiguous/combined channels flagged (HLA-DR/HLA-ABC multi-gene, CK8-18, CD31-vWF,
c-Caspase3c-PARP, panCK). STOP for sign-off before Gate A. Gate A also needs the CosMx slide (absent).

### 2026-06-15 - RESULT - METABRIC 37-marker panel ACQUIRED from Zenodo; Step 0 UNBLOCKED
**Numbers:** Zenodo 6036188 (Danenberg 2022, CC-BY-4.0), members `markerStackOrder.csv` +
`AbPanel.csv` extracted via HTTP range (remotezip 0.12.3, MIT) from a 6.65 GB zip without full
download. 39 channels = 37 protein markers + DNA1/DNA2. sha256 recorded in PROVENANCE.md.
**Consequence:** All 3 Phase-0 name-lists now in hand from authoritative sources (CosMx IO,
CosMx WTx, METABRIC-IMC). Supersedes the earlier BLOCKED entry below: the block was on the *local*
files (which were misidentified BIOKEY scRNA-seq); the panel itself was obtained from the official
deposit. Recorded discrepancies for Step-0 handling: Ho165 Epitope=`ER` (use; AbPanel's
`Rabbit IgG (H+L)` is the reagent row), Yb176=`c-Caspase3c-PARP` (run-together), DNA1/DNA2 not
gene-mappable, two HER2 channels -> ERBB2. **Step 0 may now run** (CosMx slide not required for it).

### 2026-06-15 - RESULT - METABRIC 37-marker list BLOCKED; "local METABRIC-IMC" assumption falsified
**Numbers:** 3 local `.rds` files = sparse dgCMatrix, 25,288 / 22,889 RNA genes x cells named
`BIOKEY_##_Pre_…`; tarball = 10x matrix, 33,694 RNA genes x 44,024 `sc5r…` cells.
**Consequence:** LOAD-BEARING ASSUMPTION FALSIFIED. The files in ~/Downloads presumed to be
Danenberg-2022 METABRIC-IMC are actually **Bassez-2021 BIOKEY anti-PD1 breast scRNA-seq** + a 10x
scRNA-seq matrix -- they contain RNA genes, **no IMC protein panel**. (pyreadr 0.5.0, AGPLv3, can't
read dgCMatrix S4; established via local R 4.3.1 readRDS dimnames.) No authoritative METABRIC
37-marker source is in hand; HARD RULE forbids reconstructing it from memory. **Step 0 cannot run**
(needs all 3 name-lists). Need the real Danenberg/Ali-2020 IMC marker list from an official source.

### 2026-06-15 - RESULT - CosMx IO protein + WTx gene lists acquired (authoritative, Bruker)
**Numbers:** IO protein = 62 biological targets (+2 isotype controls = 64-plex); WTx = 18,934 gene
targets (+50 negative controls). Raw xlsx sha256 recorded in data/panels/PROVENANCE.md.
**Consequence:** 2 of 3 Phase-0 name-lists now in hand from the official Bruker Spatial Biology
resource pages (verbatim extracts in data/panels/*.tsv; raw xlsx gitignored, re-fetchable). Footer
/control rows retained raw, to be dropped at Step-0 cleaning. Does not unblock Step 0 (METABRIC
still missing).

### 2026-06-15 - DECISION - Repo skeleton committed; all three Phase-0 inputs inventoried MISSING
**Numbers:** skeleton commit `6fccf97` (config.py, .gitignore, README, requirements.txt,
scripts/phase_0_00..06 docstring-only stubs).
**Consequence:** Scaffolding done, no gate logic written. Inventory of the Phase-0 blockers:
(1) CosMx breast multiomic slide -- MISSING; (2) CosMx 64-plex IO protein target list -- MISSING;
(3) METABRIC/Ali-2020 37-marker list -- MISSING as a usable list. Step 0 cannot run until the
slide + both panel lists are in hand. Danenberg-2022 METABRIC-IMC `.rds` files exist in the
user's Downloads but are Gate-E-stage material, not a Phase-0 input. Waiting on go + panel lists.

### 2026-06-15 - DECISION - Keep Fork 3 (Penumbra); Gate 1 is the kill test
**Consequence:** Proceed with the cross-modal imputation fork. Phase 0 = the Gate-1 kill test run
within a single CosMx slide (Gates A-D), deliverable = a KILL/PROCEED verdict. Do not relitigate
the fork choice. Do not build the imputer (or any non-GBM model) before the information is proven
present.

### 2026-06-15 - RESULT - Step-0 scoping readout: Check #1 fails on open data, Check #2 qualified pass
**Numbers:** Only open paired CosMx breast multiomic data = one demonstration slide (one patient).
**Consequence:** Check #1 (a holdout that tests cross-patient / cross-platform transfer) is
impossible on open data -- the only available holdout is held-out cells from the same patient/batch.
Check #2 (within-slide protein->RNA generalization is testable) is a qualified pass. Binding
constraint confirmed: Phase 0 can only falsify (fail within-slide -> dead) or earn the right to
continue; a within-slide pass is NEVER transfer evidence and must never be reported as such.

### 2026-06-15 - DECISION - Project = Penumbra (Fork 3, cross-modal imputation)
**Consequence:** Scope set: train protein->whole-transcriptome imputation on same-cell CosMx SMI
multiomics, transfer to outcome-linked METABRIC-IMC, hunt a TME biomarker invisible to the measured
proteins, validate against clinical outcome. Reuses the Operon *apparatus* (floors-first, beat the
baseline, leakage/batch controls, both-outcomes-interesting) but none of its code or conclusions.
Both a positive and a clean negative are publishable by design.
