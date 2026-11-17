# ubuntu-develop — Development Log

**Range**: 2025-05-22 → 2026-10-06 (HEAD)  
**Commits**: 438 (389 non-merge) · +40,633 / −37,002 lines since 2025-05-23  
**Contributors**: SantosJGND (365), Joao Santos (46), dsobral (18), lore-elias (9)  
**Activity peaks**: 2026-04 (58), 2025-07 (52), 2026-02 (46), 2026-08 (42); quiet stretch late 2025 (Sep–Dec nearly idle), then a sustained push through 2026.

This log covers the whole `ubuntu-develop` line in one continuous narrative. Work concentrates in `pathogen_identification` (the TELEVIR path).

---

## 1. TELEVIR pipeline and tree orchestration

- **PK-based tree navigation (2026-04)**: deprecated use of `SoftwareTreeNode` ordered indices in favour of direct PK lookups in pipeline trees; fixed mapping of software-tree nodes to local pipeline leaves. Node registration for deployment no longer relies on ordered indices (local sort + match to software nodes).
- **Pipeline type registration (2026-02)**: introduced `SoftwareTree` pipeline-type registration and available-node querying to avoid unnecessary tree-build steps; pipeline types registered on availability modification. Migration `0056` renamed `Model`→`pipeline_type`.
- **Tree deployment rework**: branching tree support, node-to-pipeline mapping, improved logging, and safe per-leaf error handling.
- **Pipeline layout (2026-05)**: assembly + contig classification moved before the read-class step; pipeline layout refactored out of utilities for visibility.
- **`RunBatching` (2026-10)**: introduced to freeze ensemble runs at point of deployment, replacing deprecated commands; sample batch registration (`SampleBatchRegister`, migration `0008`) and kill-sample-batch-jobs added.
- **Leaf parameters**: leaf-parameter model created and node leaf parameters registered in DB, loaded on demand to accelerate load (2026-06-08).
- **Validation steps**: validation QC and depletion pipeline steps created and integrated (2026-09-28); validation pre-processing / preparatory steps and extra-QC marked not DB-dependent (2026-09-03).

## 2. Reporting: compound and combined reports

- **Aggregate models (2026-04-27)**: `FinalReportCompound` / `ReportGroup` / `ReportAggregate` tables host aggregate and grouped final-report data; compounds store runs instead of raw references (migrations `0005`–`0009`, incl. `runs` on `ReportAggregate` and `shared_proportion_threshold`).
- **Heatmaps**: JSON clade heatmap store/retrieve (2026-04-28); taxon heatmap setup with table/card lines removed from combined reports (2026-08-27), debugged (2026-08-28).
- **QCrun stack (2026-06)**: `QCrun` / `TelevirRunQCStack` introduced for reporting simplification, with a `performed` field (migration `0005`) wired into frontend requests.
- **Taxon articulation**: `taxon_species` added to `ReferenceTaxid` (2026-06-08); `main_species` on `ReportGroup` (migration `0003`); report-group view toggled off when a single group is present.
- **Ordering and aggregation**: reports sorted by private reads (`private_propotion`→`private_proportion`, migration `0007`), report groups by highest coverage; aggregate sample job counts shown in view tables (2026-09-03).
- **Display and controls**: controls-collection and better reporting of references in controls (2026-08-26); accordion fixes, missing references added to sample detail (2026-10-03), assembly/mapped-reads files unblocked in compound pages, group-by-species toggle fixed, extra column display defaults to true.
- **Roots (2025)**: report colours and sample-detail compound layout (2025-09-08); error rate on sample/project report pages (2025-06-27).

## 3. Reference and taxonomy management

- **`ReferenceTaxid` and taxonomy (2026-05)**: taxonomy constants + `Taxon` table (2026-05-12), basic lineage on `ReferenceTaxid` (2026-05-13), `tax_species` (migration `0004`); tax fields / `related_names` fixed alongside a migration reset.
- **Entrez wrapper**: moved and refactored (2026-04-24); output processing made specific to `EntrezMethod` (2026-04-15); refactored for batch Entrez requests and applied to reference registration (2026-06-08); optimised taxid-library checks (2026-06-04); prevented from deletions (2026-09-22).
- **NCBI recovery**: recovery of taxids not present in the database (2026-04-24); cached reference source files (`ReferenceSourceFile.is_cache`, migration `0011`, 2026-04-29).
- **Download and matching**: better handling of reference sequence downloads (2026-08-03); object reference creation post-download fixed; taxids mapped limited to a user-defined maximum; more than one reference-source accid accounted for (2026-08-11); accid checks added (2026-10-02).
- **Raw references**: references without a DB match (acc or taxids) saved as a txt file to media; `RawReference.status` tracked (migration `0009`, 2026-10-02); soft-clip remap recall threshold keeps remaining references raw (2026-09-04).
- **Roots (2025)**: `register_references` with `--force` and process/monitor info (2025-07); improved reference-table file filters (2025-07-04); virosaurus accid exclude-genes fix (2025-05-26); MPXV Rivers added to default references (2025-07-02).

## 4. Metagenomics and classification

- **Classifier standardization (2025)**: Metaphlan conda deployment + paired-end support + output/tree parsing (2025-05-28→05-30); Voyager paired-end method, ONT-technology config, JSON taxonomy parse and empty-output handling (2025-05-27→06-03); Kraken2/Centrifuge processor (2025-08-02); Centrifuge phage filter (2025-08-04); unique-taxon read count on the Kraken2 report (2025-08-03); `min-hits` defaults and Kraken2 `nreads` min of 5 (2025-08-25); accids used from classifier reports when available (2025-06-12).
- **Classifier output tracking**: `ClassifierOutput` / `ClassifierOutputFile` tables (2026-03, migration `0002`); classifier reports moved to media and displayed in intermediate reports (2026-04-23).
- **Composition modelling (2026-07)**: `ml_api_client` utility and constants; taxid taxonomy + predict-cutoff setup; `model`, `recall` and tax-level added to constants, with tax-level properties on `ReferenceTaxid` (2026-07-10). Recall-model API calls implemented on metagenomics map planning (2026-07-14); BAM→distance-matrix calculation and composition clustering modelling integrated (2026-07-15); matrix parse + API communication debugged, recall-predict column added, taxid type checks (2026-07-16); reference-taxid order query and model methods (2026-07-16).
- **Prediction plumbing**: read and contig count sums added to compound reference models (2026-07-10); empty distance matrices worked around (2026-08-06); node-statistics calculation before the composition API request fixed to match training (2026-08-11); metagenomics tree queue refactored to prevent sequential waits (2026-07-10).

## 5. Software, pipeline steps and validation

- **GATK (2026-06-03)**: mark-and-remove-duplicates added to the remap processing module; GATK registered as software and linked to remap process steps, with tests and search options; redundant GATK run removed (2026-08-04).
- **Validation pipeline**: validation QC + depletion steps (2026-09-28) and preparatory pre-processing steps (2026-09-03) as above.
- **Parameter/step wiring**: pipeline-step parsing fixed for gate decisions on parameter setting (2026-10-03); pipeline step specified in software parameter retrieval (fixes the enmesh bug, 2026-10-06).
- **Process typing (2026-05-13)**: `PROCESS` type definition moved from the `Profile` class to `Constants`; job-specific memory/CPU params set and implemented for TELEVIR jobs; SGE references removed from `process_SGE`.
- **Settings**: assembly can run without contig classification; soft-clip remap recall threshold (2026-09-04); strategy changed to specify an extra-QC BWA filter (2025-08-27/28); `metagenomics_settings` pipeline step removed (2026-02-17).
- Focus collection updated to use the new software name table column (2026-09-08).

## 6. Reference, mapping and assembly

- **Host depletion**: BWA filter registered as polyvalent software and pipeline steps added to the polyvalent dictionary (2025-05-23→05-26); extra-QC BWA filter.
- **Reference builds**: MPXV All-Clades Nextstrain build (2025-06-25) and MPXV Rivers (2025-07-02).
- **Primers / viruses**: RSV paired-primers table (2025-06-25); primer-set management functionality (2025-05-23); IRMA and iVAR default settings (2025-06-12).
- **Consensus / calling**: Medaka model updates, consensus generation and variant calling/annotation fixes (2026-02-10→02-26); Snippy configured to use the TELEVIR binary (2025-06-05).
- **Abricate**: `db_projects` update (2025-06-16); different-column-count fix (2025-10-10).
- **Mapping/remap**: remapping added to the workflow monitor (2025-07-24); remap targets filtered (2025-06-16); BWA host mapping removed from reference map (2025-07-21); mapping target generation refactored (2026-10-01); map controls integrated (2026-10-06); teleflu focus sample select fixed (2026-10-02).

## 7. Infrastructure

- **SLURM migration (2025-07-31)**: job tracking and wait queues adapted for SLURM (`submit_job`, `sbatch` config); a BASH-SCRIPT file type added for temp files; SGE references removed from `process_SGE` (2026-05-13).
- **Django upgrades**: `is_ajax()` replaced, lazy-translation updates for Django >3.2, and breadcrumbs migrated from `bootstrap_breadcrumbs` to `view_breadcrumbs` across apps — all landed around 2025-07-28; URL patterns moved to `re_path` (2025-07-31); `python-decouple`/dotenv replaced by Django's native `.env` (2026-04-24).
- **Config**: `NCBI_EMAIL` placeholder added to the default `.env` (2026-05-04).
- **Nextstrain**: called by full path instead of relying on PATH, temp-dir handling and tree output adjusted (2026-07-22); improved logs (2025-07-22); Nextclade multiple-dataset link (2025-07-18).
- **PKs and migrations**: `DEFAULT_AUTO_FIELD` set to `BigAutoField` with PKs migrated (2026-07-22), then reverted to `AutoField` with the setting kept explicit; migration resets (2026-02-20, 2026-05-04, 2026-05-13) and PK-alter migrations across `datasets`, `extend_user`, `log_login`, `manage_virus`, `managing_files`, `settings`.
- **Path handling**: `pathlib` paths preferred over `str` concatenation (2026-05-22); `Path`/`PosixPath` concatenation fixes in utils.

## 8. UI and frontend

- **Report pages (2025)**: report colours and sample-detail compound layout (2025-09-08); error rate on sample/project report pages (2025-06-27); project/main frontend edits (2025-08-26/27).
- **Project tags (2026-04-24)**: `ProjectTag` / `ProjectTagAssignment` tables (migration `0003`), `tag.name`, and AJAX fetch-all-user-tags.
- **Navigation**: Nextclade multiple-dataset links (2025-07-18); project-settings breadcrumbs (2026-04-09); datasets/projects/references breadcrumb URLs fixed (2026-03-11, 2026-05-05).
- **Report UX**: report-wrapper run index and extra-column default display (2026-09-22); main-page paginate number reset (2026-09-22); map controls (2026-10-06).

## 9. Refactoring, bug fixes and operations

- **Utilities**: `general_utils` split into `pathutils` with explicit staticmethods, and `Classification._get_report` given a safety wrapper (2026-05-21); TELEVIR utilities refactor and description handling (2026-08-13).
- **`run_main`**: deprecated attributes converted to properties; attribution dropped its reliance on username in favour of the user object (2026-10-01).
- **Code quality**: repeated `lynt` passes across `pathogen_identification/` throughout the range.
- **Dead code**: unused `pathogen_identification` scripts deleted/deprecated (2026-02-12); deprecated commands removed with the introduction of `RunBatching` (2026-10-01); `metagenomics_settings` step removed (2026-02-17); `LockedAtomicTransaction` replaced by `transaction.atomic` (2025-07-29).
- **Safety guards**: empty report aggregates (2026-08-28), empty taxids in metadata (2026-01-05), empty Voyager output (2025-05-27), empty distance matrices (2026-08-06); `delete_media` removed from failed pipeline-node registration (2025-06-26); remap made safe (2025-06-26).
- **Tests / deps**: new TELEVIR BFS-persistence and Entrez-optimization tests (2026-06); `scikit-learn` added to requirements (2025-06-05).
- **Operations**: parallel preprocessing with BWA filter (2025-05-23); sample re-upload bug fixed by removing the dependency on sample name (2025-06-30); reads recovery in preprocessing for the multi-method structure; processes killed on project deletion.
- **Migrations**: `pathogen_identification` `0002`–`0011` (ReportGroup.main_species, ReferenceTaxid.tax_species, TelevirRunQCStack.performed, RawReferenceCompoundModel counts, RunBatch, SampleBatchRegister, RawReference.status, RunRemapMain.missing_log) plus `0054`–`0056` (PK alters, `SoftwareTree`→`pipeline_type`).
- **Integration branches**: PRs #186 (2025-07-31), #187 (2025-09-18), #195 (2025-12-17), #197 (2026-01-07), #199–#200 (2026-02), #202/#204/#205/#206/#207/#208 (2026-02), #3 from INSaFLU (2026-01-05), and #5/#6/#7 from `lore-elias/start` (2026-05→08); merges from `filtermito`, `dev-filter` and `develop`.

---

**Key architectural shifts**: job scheduling migrated from SGE to SLURM, and TELEVIR's tree navigation was rebuilt around direct PK lookups (plus pipeline-type registration) instead of fragile index-based traversal. The compound/combined report then became the main TELEVIR output surface — `QCrun`/`TelevirRunQCStack`, `ReportGroup`/`ReportAggregate`, taxon heatmaps, controls and coverage/read ordering — while `RunBatching` froze ensemble runs at deployment and a new `ml_api_client` opened a composition-modelling / recall-prediction path built on BAM-derived distance matrices.
