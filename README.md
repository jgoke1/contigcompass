# ContigCompass

I am learning Rust through this guided, AI-assisted microbial genomics tutorial. I use bacterial genome assemblies to practise reading FASTA records, indexing k-mers with a hash map, and testing how a program handles contigs that match the reference in either direction.

The three commands build on each other. `contigcompass` summarises each sequence, `anchors` finds exact matches, and `place` uses those matches to suggest contig placements and draw a dot plot. The reported coverage measures the bases touched by exact anchor matches. Assembly completeness and biological rearrangements require independent evidence.

## Run the tutorial

With Rust and Cargo installed, run these commands from the repository root.

```bash
cargo test --locked
mkdir -p results
cargo run --locked --release --bin place -- \
  test_data/demo_reference.fasta test_data/demo_queries.fasta 31 results/demo
```

Open `results/demo/report.html` in a browser. On macOS, run `open results/demo/report.html`. Choose a fresh output directory for each run. A completed run writes a `SUCCESS` file.

To try the same comparison on your own circular bacterial chromosome, replace the two example FASTA paths, choose a new output directory and append `--circular` to the command.

## Outputs

- `placements.tsv`: every input query, exact-anchor coverage, status, exported strand and reference anchor position.
- `ordered_oriented.fasta`: all queries, supported placements first in reference-anchor order, unresolved queries last in original input order. Supported reverse placements are reverse-complemented. Unresolved sequences retain original orientation. IDs are preserved; header descriptions are omitted. No joining, gap filling, polishing or sequence correction occurs.
- `unresolved.fasta`: an additional copy of unresolved queries, unchanged. Do not concatenate it with ordered_oriented.fasta, which already contains those queries.
- `blocks.tsv`: chains of non-overlapping anchors. Block span endpoints are NOT an alignment or a claim that intervening bases match.
- `coverage.tsv`: union of all exact matched k-mer intervals on the reference and on each query. Percentages are not identity, sequencing depth or genome completeness.
- `reference_uncovered.bed`: reference intervals with no exact unique-anchor coverage. They are not confirmed deletions or missing assembly sequence.
- `candidate_events.tsv`: transitions between substantial blocks, labeled orientation_switch, discontinuous_blocks, overlapping_query_blocks or possible_circular_boundary. These can reflect biological differences, assembly errors or matching limitations.
- `report.html`: offline, dependency-free report with sampled anchor dot plots. No network services required.
- `run_parameters.txt` and `SUCCESS`: settings and completion marker.

All intervals are zero-based and end-exclusive. Placement position is the first reference coordinate covered by the selected sampled chain, not an extrapolated whole-contig start. Dot plots show original query coordinates even when exported contigs are reversed.

## Algorithm and thresholds

1. Build a single-reference index of canonical exact k-mers, excluding repeated canonical keys and strand-ambiguous palindromes. Skip windows with ambiguous DNA bases.
2. Stream query k-mers through the index. Keep all hits for coverage.
3. Thin hits to non-overlapping query windows. Extend a compatible recent chain on the same strand with monotonic reference positions. Search the 32 most recently created chains; choose the lowest local coordinate discrepancy. This is a bounded greedy heuristic, not optimal dynamic-programming chaining or base alignment.
4. Maximum inter-anchor step: max(10,000 bp, k); maximum local reference/query step discrepancy: 500 bp.
5. A substantial chain requires >=3 anchors and >=max(500,3k) covered bases.
6. Place only if the dominant chain covers >=50% of query bases and accounts for >=90% of sampled-anchor covered bases, with no substantial secondary chain >=10% of its support. These are development heuristics, not calibrated confidence probabilities.
7. Otherwise report insufficient_support, partial_support or ambiguous_blocks. No anchors yields no_matching_anchors.

`--circular` flags possible crossings of the reference origin when both substantial blocks approach opposite ends (within min(10,000 bp,10% reference length)). It does not rotate chromosomes or merge crossing chains. Such placements remain conservative. Reference topology must be chosen by the user.

## What coverage means

Reference anchor coverage = union length of matched reference k-mer intervals / reference length.
Query anchor coverage = sum of per-query union lengths / sum of query lengths.
Overlapping k-mers never double-count a base within the same coordinate system. Identical query contigs count separately in the query total, but not in reference union coverage. Exact matching is sensitive to k, divergence, ambiguous bases and repeat exclusion. This tool alone cannot estimate organism-level completeness or contamination.

## Limits and validation

Single reference FASTA record only; query can contain many contigs. Duplicate query IDs and invalid DNA are rejected. All query sequences are held in memory, and the reference index still uses allocated string keys. Very large genomes or many thousands of contigs need further engineering. Report plot sampling may omit rare events; inspect blocks and candidate events as well.

Synthetic tests cover forward/reverse placement, union coverage, absent segments, unknown contigs, indel tolerance, mixed orientations, circular boundaries, duplicate IDs and refusing to overwrite output. See VALIDATION.md for the executed checks. Biological structural events require independent alignment and, where available, read evidence. No independent-aligner accuracy benchmark has been performed.

## Existing commands

```bash
cargo run --bin contigcompass -- assembly.fasta
cargo run --bin anchors -- reference.fasta assembly.fasta 31 NEW_summary.tsv > anchors.tsv
```

The original strict anchor summary remains available for learning and diagnostics.
