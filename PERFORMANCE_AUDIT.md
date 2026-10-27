# Workflow performance audit — 2026-10-27

Scope: active Snakemake rules, six workflow Python scripts, helpers, configs, CLI/build/Slurm wiring and development tests. All changes were made in source. No production FASTQs, completed biological results, installed SIF or remote cluster were modified. Existing unrelated checkout changes were retained.

The reported approximately 265 GB prepared FASTQs and five-hour single-core QC are user-supplied observations, not measurements from this development session. Exact parity and failure behavior were tested with deterministic fixtures; production throughput and peak RAM remain unmeasured.

## Findings and implemented changes

| Component | Original behavior | Bottleneck | Implemented optimization | Risk / limitation | Priority |
|---|---|---|---|---|---|
| FASTQ QC | Parse both raw mates and both prepared mates after whitelist construction | Repeated gzip inflation and serial parsing | Routine QC verifies a sample/file-bound receipt and whitelist hash; no FASTQ traversal | Independent prepared-R2 integrity checking is opt-in | P0 |
| Whitelist | Full prepared R1 traversal; sequence lines only | Duplicated by QC | Validate complete records and 27+8 geometry during the existing whitelist scan; keep the same sorted observed panel-valid barcode set | Memory scales with distinct observed valid barcodes | P0 |
| Unsampled lanes | Concatenate raw gzip files, run one Cutadapt job, later reread original lanes for pairing/counts | Raw-sized copy plus another raw scan | Independent paired Cutadapt units, strict names checked in its existing iterator, sum exact counts, assemble prepared gzip members in configured order | Assembly still reads/writes prepared compressed bytes; downsampling keeps the original path | P0 |
| Cropped R1 | Write cropped gzip, decompress it, cut columns, recompress | Extra compressed intermediate and gzip cycle | Stream Cutadapt R1 stdout through the exact existing `cut` columns into pigz | Adapter tied to pinned Cutadapt 4.9/Linux fork; pipeline failure propagation tested | P0 |
| Chemistry | Classify early, reject non-enhanced/unknown only in late QC | Expensive work after a known incompatible classification | Raise in bead validation before Cutadapt; retain diagnostic log | Same sampled prefix and classifier; not an exhaustive chemistry assay | P0 |
| Reference | Check contents only when quantification starts; record fixed GRCh38/v49 labels | Late failure and potentially incorrect provenance | Shared reference preflight before preparation, optional hash-bound manifest, explicit unknown identity when absent | Manifest is declared identity, not biological authentication; large index files are stat-checked, not fully hashed | P0 |
| Downsampling | Concatenate sample lanes, run seqtk separately per mate with a shared seed | Raw copy, two sequential sampling commands, full original-input integrity scan | Preserve current selection and original-input denominator; stream subsequent preparation | Restarting RNG per lane would change reads; intentionally not done | Keep |
| Native counts/H5AD | Sparse export, serialized reopen, independent native comparison | Simultaneous sparse allocations and repeated matrix materialization | Preserve checks and count definitions | Chunked validation/memory reduction remains proposed | P2 |
| Resources | Preprocessing had no explicit memory reservation | Ambiguous concurrency budgeting | Preprocessing 8192 MB per unit; whitelist/QC 4096 MB; assembly/concat/reference preflight 512 MB; downsample 1024 MB; bead scan 256 MB | Starting reservations, not measured peak memory | P1 |
| Intermediate retention | Prepared and combined raw files retained; cropped R1 temporary | Large disk footprint | Prepared mates persist; lane mates are temporary after assembly; no cropped R1; raw concat only for downsampling | Multi-lane assembly transiently holds unit and final prepared files | P1 |

## Four newly implemented upgrades

1. **Independent unsampled lanes.** `prepare_fastq_unit` consumes one configured raw pair. Its report records sample/unit identity, input stat identities, strict mate-validation mode and exact Cutadapt input/output counts. `standardize_cb_umis_cutadapt` checks these reports against configured sources and assembles units in lane-list order. Single-unit output uses hard links; multiple units use bounded binary copying of gzip members, without inflation. Successful unit files survive a peer failure and are cleaned after assembly. Unit reports/logs remain available. A change to a shared dependency, such as chemistry validation or source code, can correctly invalidate multiple units; reuse is not guaranteed after arbitrary input/config changes.
2. **Streamed R1 preparation.** The adapter, overlap, no-indel requirement, error tolerance, crop action, untrimmed discard and pair-filter mode remain unchanged. Cropped R1 is sent through `cut -c1-9,14-22,27- | pigz -p 1`. This retains the original transformation of every line, including long headers. Logs are isolated from FASTQ stdout and Bash `pipefail` rejects a failure anywhere in the pipeline. Outputs are staged on the same filesystem and renamed after success. Renaming multiple files is not a transactional commit; declared dependencies and receipts remain the completion gate. Hard kills can leave hidden staging directories, although ordinary exceptions clean them.
3. **Early chemistry rejection.** The unchanged classifier scans the same sample-wide prefix of at most `bead_scan_n_reads` R1 sequences. For unsampled lanes it streams them in configured order without raw concatenation; parity tests compare against the original concatenated scan across a lane boundary. Non-enhanced and unknown classifications fail before any Cutadapt launch, with diagnostics retained in the log. Downsampled inputs retain classification of the existing downsampled R1. The late gate remains as a defensive check.
4. **Reference preflight and honest provenance.** The shared preflight verifies piscem metadata, expected index files and the three-column S/U map before chemistry/preparation, including raw concatenation/downsampling. `genome` and `GENCODE_release` are no longer constants. An optional manifest supplies their declared values and must match SHA256 of the t2g map and index JSON. No manifest means `unknown` plus an explicit unverified status, preserving compatibility with existing configs. Quantification rechecks preflight input stat identities before using its recorded metadata; it does not reparse the map or hash huge index files. Preflight uses ordinary dependencies: reference updates can invalidate preparation as well as quantification. An ordering-only `ancient()` variant was rejected because the pinned scheduler skipped missing downstream targets in a small recovery test; it is not enabled.

Implementation is in `workflow/scripts/prepare_reads.py`, `workflow/scripts/simpleaf_quant.py`, `workflow/src/prepare_fastq.snmk`, `workflow/src/required_helpers.py` and `workflow/Snakefile`.

## Exact counts, pairing and scientific behavior

The installed Cutadapt 4.9 report implementation emits exact pair counts in `read_counts.input` and `output`. Its paired reader parses and pairs records before filtering. The existing QC's name check was stricter than Cutadapt's default: for example, Cutadapt accepts `read1`/`read2` while legacy QC rejects them. The new wrapper therefore adds the legacy first-token comparison and asymmetric `/1`/`/2` normalization to `dnaio.pairedend.TwoFilePairedEndReader.__iter__` inside the Cutadapt process. Fork workers inherit it. It is restored after processing and guarded against unsupported Cutadapt/start-method combinations. This is a version-bound parser hook, not a change to trimming or filtering. Tests exercise serial and multiple-worker cases, including these formerly over-permissive names.

Statistics are not blindly assumed to equal final files. Each R2 unit is written directly by Cutadapt; R1 column removal preserves record order/count. Ordered gzip-member assembly preserves those records. The final R1 parser checks complete four-line records, nonempty sequence, matching quality length, exactly 35 DNA bases and agreement with the sum of output-pair counts. Ordinary input edits during preparation/validation are rejected by before/after stat identities. Aggregated reports check sample, unit order, source paths/sizes/mtimes and schema/strict-validation mode before avoiding raw scans.

The whitelist is still the sorted set of observed 27nt barcodes with all three segments in the original panels. Lowercase/ambiguous barcode segments remain excluded; N-containing UMIs remain allowed, and no extra reads are filtered by whitelist generation. Raw and prepared counts include records with barcodes absent from the whitelist. For downsampling, raw counts still come from original-lane scans, so `retained_percent = 100 * final_prepared / all_original_pairs`. Counts are exact, not estimates. `corrected_barcodes` remains `unavailable` in the pre-quantification QC table. Legacy direct-report helper calls retain their compatibility path and may require original-lane scans.

The receipt records sample, schema, count sources, exact counts, ordered file stat identities and whitelist SHA256. It is atomically published after checks; Snakemake touches output timestamps, so whitelist association deliberately uses content hashing. `fastq_qc.tsv` retains its five columns and is consumed only as the quantification gate.

**Routine guarantees:** strict raw pairing established during each unsampled unit's existing parse; validated final R1 structure/count/geometry; prepared-pair consistency by the inspected record-preserving construction; exact whitelist/counts and configured sample/file association. Routine mode does not independently decompress final R2 after generation. Path/size/mtime identity is not a content checksum and cannot detect edits deliberately preserving those attributes. Sample labels and manifest labels are not biological authentication.

**Exhaustive mode:** YAML `fastq_qc_exhaustive: true` rereads all original and prepared mates, checks FASTQ/gzip integrity, strict transformed mate-header consistency and exact counts, and compares them with the receipt. This can take hours and is never applied automatically to external sequencing files by development tests. Changed receipt input identity fails before scanning; an intentional replacement requires reviewed regeneration, not trusting stale metadata.

Unchanged scientific contracts: original BD panel order and 1-based BD4 barcode conversion; 27nt barcode/8nt UMI; enhanced linkers/error/crop settings; observed whitelist; SimpleAF spliceu/USA commands and `cr-like` UMI resolution; both orientations; `--unfiltered-pl --min-reads 1`; gene IDs/order; no new cell/gene filtering; H5AD X/spliced = S+A and unspliced = U. No original BD expression matrix is imported.

## Expected performance impact

Let R be compressed original mates, B1/B2 compressed final prepared mates, and C the old cropped R1 gzip size.

- Original whitelist + QC read `R + 2*B1 + B2`. Current routine unsampled whitelist + QC reads `B1`, whether input has one lane or many. Downsampled samples still read `R + B1` to preserve original-input validation/denominator.
- Unsampled multi-lane raw concatenation no longer reads and writes R bytes. Prepared assembly instead reads/writes B1+B2 compressed bytes, with unit files removed after success. Single-unit assembly uses hard links without copying data.
- Cropped-R1 gzip creation, its C-byte reread and corresponding inflate/compress cycle are removed. Final R1 remains gzip-compressed and persistent.
- Quantification still reads prepared FASTQs once. Strict name checking adds comparisons within the existing parse, not another decompression/FASTQ parsing pass.

For the reported approximately 265 GB prepared pair, routine QC avoids that full additional prepared-pair traversal. Unsampled multi-lane inputs now also avoid the extra original-lane QC scan. These are pass-count/I/O deductions, not measured production speedups. Storage contention, compression and barcode diversity still determine runtime; no numerical speedup or production memory reduction is claimed.

## Reference manifest, resources and resume

Existing sample configs remain valid. Optional `simpleaf.reference_manifest` points to a JSON manifest; otherwise `reference_manifest.json` is discovered in the resolved index directory. The manifest requires `version: 1`, `index_type: spliceu`, nonempty `genome` and `GENCODE_release`, and exact `reference_t2g_sha256`/`reference_index_meta_sha256`. Build records must supply identities; do not derive them from directory names. The README shows its schema. The preflight records manifest path/hash and matching status, alongside the map and index-metadata hashes. It does not write to the reference.

Each streaming unit requests `nthreads`, capped by allocation, and needs at least two cores: `threads-1` Cutadapt workers plus one pigz thread, with helper/I/O overhead. Lane jobs may overlap when CPU/memory budgets allow. `preprocess_mem_mb` defaults to 8192 per unit; `fastq_qc_mem_mb` defaults to 4096. Quantification keeps `max_mem_mb: 64000`; export/final validation keep `export_mem_mb: 16000`. Whitelist RAM depends on distinct observed valid barcodes, and sparse export memory depends on counts/nonzeros. Slurm must provide the total scheduled resources; these are unprofiled starting reservations.

Final prepared FASTQs remain non-temp for retrying quantification. Prepared lane files are temp after successful assembly, and multiple lane outputs coexist with assembled outputs during copying. Combined/downsampled raw files persist only on the downsampling path. No automatic deletion of existing sequencing data was added. Normal failed jobs retain completed steps and can be retried with `--rerun-incomplete`; the regression test verifies a completed lane report is unchanged after a peer fails and its R2 is repaired. If code, shared chemistry inputs or other dependencies change, Snakemake may invalidate more work. Failed SimpleAF work is still restarted at the quantification step, not an internal mapping checkpoint.

New dependencies and report schemas deliberately prevent silently trusting older preparation outputs. Do not run this migration against completed results: inspect a dry run and select an isolated `working_dir`. A previously prepared dataset without current unit/preflight metadata requires regeneration or a separate reviewed migration. No production migration or container replacement was performed.

## Validation and rebuild

```bash
export XDG_CACHE_HOME="$PWD/.cache"
export VELOCITY_SNAKEMAKE="$PWD/.pixi/envs/default/bin/snakemake"
export VELOCITY_SIMPLEAF="$PWD/.pixi/envs/default/bin/simpleaf"
export VELOCITY_CUTADAPT_PYTHON="$PWD/.pixi/envs/default/bin/python"
export VELOCITY_PREPROCESS_BIN="$PWD/.pixi/envs/default/bin"
.pixi/envs/default/bin/python -W error::DeprecationWarning -W error::FutureWarning tests/velocity/test_velocity.py
.pixi/envs/default/bin/python -m unittest discover -s tests -p test_cli.py
PATH="$PWD/.pixi/envs/default/bin:$PATH" .pixi/envs/default/bin/python tests/velocity/synthetic_e2e.py
```

Coverage: upstream prepared-stream parity over all diversity offsets and long headers; same QC table/whitelist; strict malformed/truncated/mismatched mates including Cutadapt-permitted terminal-digit differences; serial/multiple-worker processing; explicit cut/pigz failure cleanup; no redundant raw traversal for unit reports; wrong-sample/reordered receipts; early bad-reference/V1/unknown-chemistry failure before Cutadapt; missing/matched/mismatched manifests and stale preflight; failed-lane reuse with unchanged shared inputs; unchanged sampling commands and original-pair denominator; serialized sparse/fractional S+A/U output. End-to-end checks use 40 generated pairs presented as single, split-lane and sample-wide-downsampled inputs, with matching single/split counts, temporary-unit cleanup and unchanged-output completed-run resume.

Final local checks passed: **24 velocity regression tests**, **two CLI tests**, Python compilation, documentation/link/syntax checks, and the three synthetic input layouts described above, including failed-lane recovery and completed-run reuse. Logs are retained in `.cache/validation/provenance-upgrade.2i2vlfvi/{tests,synthetic}.log`. No production throughput benchmark or rebuilt-SIF validation is represented.

Rebuild a new artifact when ready; the existing recipe already copies the new source and tests:

```bash
./container_build.sh "$PWD/rhapsodist-simpleaf.performance.sif"
```

The build refuses to overwrite its destination, tests fixtures/synthetic workflows and writes a checksum. No installed SIF was modified; rebuilding/deployment and production/HPC benchmarking remain outstanding.

GitNexus impact was LOW for resolved reference/driver functions and UNKNOWN for dynamic/new Snakemake/preparation symbols; source/text inspection identified and updated their callers. The old index does not cover all current files. Graph change analysis cannot be treated as a complete all-clear, and no commit is made as part of this task.

## Additional implemented upgrade: execution provenance

Native quantification provenance and H5AD now retain aligned workflow source paths/SHA256 values, the lockfile SHA256, canonical effective-config SHA256 and preparation-receipt SHA256. Hashes cover the actual runtime source bytes and small metadata; they do not reread FASTQs. Source/lock/receipt changes during quantification fail before publication. Snakemake declares the source, lock and receipt dependencies and tracks the config digest as a parameter. Direct library calls without config/receipt information explicitly record `unavailable` for those fields. Regression checks cover stable and changed digests, incomplete input rejection, failed publication after a mid-run change, and round-trip preservation in H5AD.

## Remaining proposals

The four original upgrades and the execution-provenance upgrade are implemented. Four proposals remain: reduce sparse export/validation peak memory; parallelize the unchanged mate-downsampling commands within a CPU budget; extend hard-interruption recovery tests and cleanup; benchmark compression levels against storage cost. Staging normal preparation writes was necessary for streamed failure handling, but does not claim transactional multi-file publication or hard-kill cleanup. None of these deferred items was silently enabled.
