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
