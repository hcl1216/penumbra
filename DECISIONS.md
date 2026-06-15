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

### 2026-06-15 - DECISION - Fork-2 discovery-cohort scoping reported (recommendation PENDING APPROVAL; NOT finalized)
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
