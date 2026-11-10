use bio::io::fasta;
use contigcompass::index::ReferenceIndex;
use contigcompass::matching::{Strand, match_indexed_query};
use contigcompass::placement::{PlacementStatus, assess_consistency};

use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};

fn validate_dna(sequence: &[u8], id: &str) -> Result<(), Box<dyn Error>> {
    if sequence.is_empty() {
        return Err(format!("Empty sequence: {id}").into());
    }

    for base in sequence {
        if !matches!(
            base.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'N'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
        ) {
            return Err(format!("Invalid DNA character in sequence: {id}").into());
        }
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 5 {
        eprintln!("Usage: anchors <reference.fasta> <query.fasta> <k> <summary.tsv>");
        std::process::exit(1);
    }

    let k: usize = args[3].parse()?;

    if k == 0 {
        return Err("k must be greater than zero".into());
    }

    // Read exactly one reference record.
    let reader = fasta::Reader::from_file(&args[1])?;
    let mut records = reader.records();

    let reference = records
        .next()
        .ok_or("Reference FASTA contains no records")??;

    if let Some(extra) = records.next() {
        extra?;
        return Err("Reference FASTA must contain exactly one record".into());
    }

    validate_dna(reference.seq(), reference.id())?;

    if reference.seq().len() < k {
        return Err("Reference sequence is shorter than k".into());
    }

    // Build once, then reuse for every query.
    let reference_index = ReferenceIndex::new(reference.seq(), k)?;
    let queries = fasta::Reader::from_file(&args[2])?;

    // Create both writers outside the query loop.
    let mut output = io::BufWriter::new(io::stdout().lock());

    let summary_file = File::create_new(&args[4])?;
    let mut summary = io::BufWriter::new(summary_file);

    writeln!(
        output,
        "query_id\tquery_start_0based\treference_id\treference_start_0based\tstrand\tk"
    )?;

    writeln!(
        summary,
        "query_id\tlength_bp\tanchor_count\tforward_anchors\treverse_anchors\tstatus"
    )?;

    let mut query_count = 0;

    for result in queries.records() {
        // Everything using `query` stays inside this loop.
        let query = result?;
        validate_dna(query.seq(), query.id())?;
        query_count += 1;

        let matches = match_indexed_query(&reference_index, query.seq())?;

        let forward_count = matches
            .iter()
            .filter(|hit| hit.strand == Strand::Forward)
            .count();

        let reverse_count = matches.len() - forward_count;

        let status = if matches.is_empty() {
            "no_matching_anchors"
        } else {
            match assess_consistency(&matches) {
                PlacementStatus::InsufficientAnchors => "insufficient_anchors",
                PlacementStatus::ConsistentForward => "consistent_forward",
                PlacementStatus::ConsistentReverse => "consistent_reverse",
                PlacementStatus::ConflictingAnchors => "conflicting_anchors",
            }
        };

        writeln!(
            summary,
            "{}\t{}\t{}\t{}\t{}\t{}",
            query.id(),
            query.seq().len(),
            matches.len(),
            forward_count,
            reverse_count,
            status
        )?;

        eprintln!("{}: {} anchor matches", query.id(), matches.len());

        for hit in matches {
            let strand = match hit.strand {
                Strand::Forward => "+",
                Strand::Reverse => "-",
            };

            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}",
                query.id(),
                hit.query_position,
                reference.id(),
                hit.reference_position,
                strand,
                k
            )?;
        }
    }

    if query_count == 0 {
        return Err("Query FASTA contains no records".into());
    }

    summary.flush()?;
    output.flush()?;

    Ok(())
}
