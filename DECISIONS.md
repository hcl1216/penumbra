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
