# Validation performed on 29 October 2026

- Rust 1.99.1, Linux x86_64; optimized build completed.
- `cargo test --locked`: 31 library tests and 4 CLI integration tests passed (35 total).
- `cargo fmt --check` passed.
- Original 24 library tests preserved and passing.
- CLI checks validate exported sequence content, reference-order sorting, preservation of unresolved sequences, merged reference coverage (exactly 6,000 of 12,000 synthetic reference bases), insertion tolerance, mixed orientation rejection, circular-origin labeling, duplicate-ID rejection and existing-output protection.
- `example_output/` contains a reproducible synthetic demonstration from the included `test_data/demo_*` FASTAs.
- Real-data smoke test: NCBI chromosome AP006716.1 (JCSC1435) versus CP025031.1 (SGAir0252), k=31, circular mode. Downloaded via NCBI EFetch by accession on 29 October 2026.
- Reproduced 2,092,571 exact anchor matches. Reference: 2,207,379/2,685,015 bases = 82.211% exact-anchor coverage. Query: 2,231,218/2,632,932 bases = 84.743%.
- Query classified ambiguous_blocks and left unchanged. No organism-level completeness or confirmed rearrangement claim follows from this result.
- Initial optimized real-data run: about 11.3 seconds and 532 MiB peak resident memory on this execution environment. Not a cross-machine performance guarantee.
- `validation_real/` includes the report and tabular results, without the downloaded genomes or exported FASTAs. Regenerate the complete outputs with the documented command.

Remaining validation: comparison against an independent aligner and broader curated benchmarks, calibrated thresholds, and macOS execution. Circular-origin chains are flagged conservatively rather than automatically merged. This is an integrated research prototype, not production validation.

Report checks: both HTML files parsed successfully, with populated tables and SVG plots. Browser visual inspection was unavailable because the browser download failed. Both legacy CLI commands were also run successfully on the synthetic example.
