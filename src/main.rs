use bio::io::fasta;
use std::env;
use std::error::Error;
use std::io::{self, Write};

fn main() -> Result<(), Box<dyn Error>> {
    // The first argument is the program name; the second is the FASTA path.
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: contigcompass <assembly.fasta>");
        std::process::exit(1);
    }

    let reader = fasta::Reader::from_file(&args[1])?;
    let mut output = io::BufWriter::new(io::stdout().lock());

    writeln!(output, "contig_id\tlength_bp\tgc_percent\tn_count")?;

    for result in reader.records() {
        let record = result?;
        let sequence = record.seq();

        let mut gc_count: usize = 0;
        let mut acgt_count: usize = 0;
        let mut n_count: usize = 0;

        for base in sequence {
            match base.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc_count += 1;
                    acgt_count += 1;
                }
                b'A' | b'T' => {
                    acgt_count += 1;
                }
                b'N' => {
                    n_count += 1;
                }
                // Accept other standard IUPAC ambiguity symbols.
                b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {}
                other => {
                    return Err(format!(
                        "Invalid DNA character '{}' in contig '{}'",
                        other as char,
                        record.id()
                    )
                    .into());
                }
            }
        }

        // GC percentage uses only unambiguous A, C, G and T bases.
        let gc_percent = if acgt_count > 0 {
            format!("{:.2}", 100.0 * gc_count as f64 / acgt_count as f64)
        } else {
            String::from("NA")
        };

        writeln!(
            output,
            "{}\t{}\t{}\t{}",
            record.id(),
            sequence.len(),
            gc_percent,
            n_count
        )?;
    }

    output.flush()?;
    Ok(())
}
