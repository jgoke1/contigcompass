use bio::io::fasta;
use contigcompass::{
    analysis::{Chain, chains, covered, decide, gaps},
    dna::reverse_complement,
    index::ReferenceIndex,
    matching::{Strand, match_indexed_query},
};
use std::{
    collections::HashSet,
    env,
    error::Error,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};
type Result<T> = std::result::Result<T, Box<dyn Error>>;
fn strand(s: Strand) -> &'static str {
    if s == Strand::Forward { "+" } else { "-" }
}
fn valid(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|b| b"ACGTRYSWKMBDHVN".contains(&b.to_ascii_uppercase()))
}
fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn writer(dir: &Path, name: &str) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(dir.join(name))?))
}
fn fasta_write(w: &mut impl Write, id: &str, seq: &[u8]) -> Result<()> {
    writeln!(w, ">{id}")?;
    for chunk in seq.chunks(80) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}
fn boundary(a: &Chain, b: &Chain, k: usize, len: usize) -> bool {
    if a.strand() != b.strand() {
        return false;
    }
    let x = a.hits.last().unwrap().reference_position;
    let y = b.hits[0].reference_position;
    let window = 10000.min(len / 10);
    match a.strand() {
        Strand::Forward => x + k >= len.saturating_sub(window) && y < window,
        Strand::Reverse => x < window && y + k >= len.saturating_sub(window),
    }
}
fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 || args.iter().any(|a| a == "--help") {
        eprintln!(
            "Usage: place <reference.fasta> <query.fasta> <k> <NEW_output_directory> [--circular]\nSingle reference record. Conservative exact-anchor placement; not a completeness estimator."
        );
        return if args.iter().any(|a| a == "--help") {
            Ok(())
        } else {
            Err("Missing arguments".into())
        };
    }
    if args.len() > 6 || (args.len() == 6 && args[5] != "--circular") {
        return Err("Unknown option".into());
    }
    let circular = args.len() == 6;
    let k: usize = args[3].parse()?;
    if k < 3 {
        return Err("k must be at least 3".into());
    }
    let rr: Vec<_> = fasta::Reader::from_file(&args[1])?
        .records()
        .collect::<std::result::Result<_, _>>()?;
    if rr.len() != 1 || !valid(rr[0].seq()) || rr[0].seq().len() < k {
        return Err("Provide one nonempty valid reference record, length >= k".into());
    }
    let reference = &rr[0];
    let rlen = reference.seq().len();
    let queries: Vec<_> = fasta::Reader::from_file(&args[2])?
        .records()
        .collect::<std::result::Result<_, _>>()?;
    if queries.is_empty() {
        return Err("No query records".into());
    }
    let mut ids = HashSet::new();
    for q in &queries {
        if !valid(q.seq()) || !ids.insert(q.id()) {
            return Err(format!("Invalid DNA or duplicate query ID: {}", q.id()).into());
        }
    }
    let dir = Path::new(&args[4]);
    fs::create_dir(dir)?;
    let mut summary = writer(dir, "placements.tsv")?;
    let mut blocks = writer(dir, "blocks.tsv")?;
    let mut events = writer(dir, "candidate_events.tsv")?;
    writeln!(
        summary,
        "query_id\tlength_bp\tanchor_count\texact_anchor_covered_bp\texact_anchor_coverage_percent\tstatus\texport_strand\treference_anchor_start_0based\tdominant_chain_support_fraction"
    )?;
    writeln!(
        blocks,
        "query_id\tblock\tquery_start_0based\tquery_end_exclusive\treference_start_0based\treference_end_exclusive\tstrand\tnonoverlapping_anchors\tanchor_covered_bp"
    )?;
    writeln!(
        events,
        "query_id\tleft_block\tright_block\tquery_left_end_exclusive\tquery_right_start_0based\ttype"
    )?;
    eprintln!("Building reference index...");
    let index = ReferenceIndex::new(reference.seq(), k)?;
    let mut ref_intervals = Vec::new();
    let mut total_covered = 0;
    let mut export = Vec::new();
    let mut rows = String::new();
    let mut plots = String::new();
    let mut placed = 0;
    let mut event_count = 0;
    for (qi, q) in queries.iter().enumerate() {
        let hits = match_indexed_query(&index, q.seq())?;
        let qintervals: Vec<_> = hits
            .iter()
            .map(|h| (h.query_position, h.query_position + k))
            .collect();
        let qcovered = covered(&qintervals);
        total_covered += qcovered;
        // Merge per contig before accumulation to bound memory on large assemblies.
        let ri: Vec<_> = hits
            .iter()
            .map(|h| (h.reference_position, h.reference_position + k))
            .collect();
        ref_intervals.extend(contigcompass::analysis::merged(&ri));
        let cs = chains(&hits, k, 10000.max(k), 500);
        let p = decide(&cs, q.seq().len(), k, 500.max(3 * k), 0.5);
        let status = if hits.is_empty() {
            "no_matching_anchors"
        } else {
            p.status
        };
        if p.strand.is_some() {
            placed += 1;
        }
        let orientation = p.strand.map(strand).unwrap_or("unchanged");
        let pos = p
            .position
            .map(|x| x.to_string())
            .unwrap_or_else(|| "NA".into());
        writeln!(
            summary,
            "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.4}",
            q.id(),
            q.seq().len(),
            hits.len(),
            qcovered,
            100.0 * qcovered as f64 / q.seq().len() as f64,
            status,
            orientation,
            pos,
            p.support
        )?;
        export.push((p.position.unwrap_or(usize::MAX), qi, p.strand));
        rows.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{:.2}%</td><td>{}</td><td>{}</td></tr>",
            esc(q.id()),
            q.seq().len(),
            100.0 * qcovered as f64 / q.seq().len() as f64,
            status,
            orientation
        ));
        for (i, c) in cs.iter().enumerate() {
            writeln!(
                blocks,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                q.id(),
                i + 1,
                c.query_start(),
                c.query_end(k),
                c.ref_start(),
                c.ref_end(k),
                strand(c.strand()),
                c.hits.len(),
                c.covered(k)
            )?;
        }
        let strong: Vec<_> = cs
            .iter()
            .enumerate()
            .filter(|(_, c)| c.hits.len() >= 3 && c.covered(k) >= 500.max(3 * k))
            .collect();
        for pair in strong.windows(2) {
            let ((ai, a), (bi, b)) = (pair[0], pair[1]);
            let typ = if a.query_end(k) > b.query_start() {
                "overlapping_query_blocks"
            } else if circular && boundary(a, b, k, rlen) {
                "possible_circular_boundary"
            } else if a.strand() != b.strand() {
                "orientation_switch"
            } else {
                "discontinuous_blocks"
            };
            writeln!(
                events,
                "{}\t{}\t{}\t{}\t{}\t{}",
                q.id(),
                ai + 1,
                bi + 1,
                a.query_end(k),
                b.query_start(),
                typ
            )?;
            event_count += 1;
        }
        plots.push_str(&format!("<details {}><summary>{} — anchor dot plot</summary><p>X: reference position (0–{} bp); Y: original query position (0–{} bp). Blue +; orange −. Up to 5,000 sampled anchors.</p><svg viewBox='0 0 900 400' role='img' aria-label='Anchor dot plot'><rect x='30' y='15' width='850' height='355' fill='#f1f5f9'/>",if qi==0{"open"}else{""},esc(q.id()),rlen,q.seq().len()));
        let step = hits.len().div_ceil(5000).max(1);
        for h in hits.iter().step_by(step) {
            let x = 30.0 + 850.0 * h.reference_position as f64 / rlen as f64;
            let y = 370.0 - 355.0 * h.query_position as f64 / q.seq().len() as f64;
            plots.push_str(&format!(
                "<circle cx='{x:.2}' cy='{y:.2}' r='1.3' fill='{}'/>",
                if h.strand == Strand::Forward {
                    "#2563eb"
                } else {
                    "#ea580c"
                }
            ));
        }
        plots.push_str("</svg></details>");
        eprintln!("{}: {} anchors; {}", q.id(), hits.len(), status);
    }
    let rc = covered(&ref_intervals);
    let total: usize = queries.iter().map(|q| q.seq().len()).sum();
    let mut cov = writer(dir, "coverage.tsv")?;
    writeln!(
        cov,
        "reference_id\treference_length_bp\treference_anchor_covered_bp\treference_anchor_coverage_percent\tquery_total_bp\tquery_anchor_covered_bp\tquery_anchor_coverage_percent"
    )?;
    writeln!(
        cov,
        "{}\t{}\t{}\t{:.3}\t{}\t{}\t{:.3}",
        reference.id(),
        rlen,
        rc,
        100.0 * rc as f64 / rlen as f64,
        total,
        total_covered,
        100.0 * total_covered as f64 / total as f64
    )?;
    let mut bed = writer(dir, "reference_uncovered.bed")?;
    for (s, e) in gaps(&ref_intervals, rlen) {
        writeln!(bed, "{}\t{}\t{}", reference.id(), s, e)?;
    }
    export.sort_by_key(|&(pos, qi, _)| (pos, qi));
    let mut fa = writer(dir, "ordered_oriented.fasta")?;
    let mut unresolved = writer(dir, "unresolved.fasta")?;
    for (_, qi, s) in export {
        let q = &queries[qi];
        let seq = if s == Some(Strand::Reverse) {
            reverse_complement(q.seq())
        } else {
            q.seq().to_vec()
        };
        fasta_write(&mut fa, q.id(), &seq)?;
        if s.is_none() {
            fasta_write(&mut unresolved, q.id(), q.seq())?;
        }
    }
    let mut report = writer(dir, "report.html")?;
    writeln!(
        report,
        "<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width'><title>ContigCompass report</title><style>body{{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px;color:#163047}}h1{{color:#086f76}}table{{border-collapse:collapse;width:100%}}th,td{{padding:12px;border-bottom:1px solid #cbd5e1;text-align:left}}details{{margin:24px 0}}svg{{width:100%}}.note{{background:#fff4d8;padding:18px;border-radius:8px}}</style><h1>ContigCompass</h1><p>Reference: {} · k={} · {} contigs · {} placed · {} candidate block transitions</p><h2>{:.2}% exact-anchor reference coverage</h2><p class='note'>This is not genome completeness, sequence identity, or read depth. Coverage counts the union of exact matched k-mer intervals. Repeats and sequence divergence reduce detectable coverage. Candidate events require independent investigation. Dot plots are sampled; metrics use all matches.</p><table><tr><th>Contig</th><th>Length</th><th>Query anchor coverage</th><th>Status</th><th>Export direction</th></tr>{}</table>{}<h2>Method</h2><p>Non-overlapping query anchors; bounded greedy chaining (32 recent chains), maximum step 10,000 bp or k if larger; local drift 500 bp. Placement requires at least 3 anchors, 500 covered bases (or 3k if larger), 50% query coverage, 90% dominant-chain support, and no substantial secondary chain. These are heuristic thresholds, not calibrated confidence. Circular mode labels possible origin crossings; it does not rotate or force placement of crossing contigs.</p></html>",
        esc(reference.id()),
        k,
        queries.len(),
        placed,
        event_count,
        100.0 * rc as f64 / rlen as f64,
        rows,
        plots
    )?;
    let mut run = writer(dir, "run_parameters.txt")?;
    writeln!(
        run,
        "ContigCompass {}\nk={}\ncircular={}\nmax_gap={}\ndrift=500\nmin_block_bp={}\nmin_query_fraction=0.5\nmin_chain_support=0.9\nreference={}\nquery={}\nCoordinates: 0-based, half-open intervals. Reference placement coordinate is an anchor start, not an extrapolated full-contig start.\nSUCCESS",
        env!("CARGO_PKG_VERSION"),
        k,
        circular,
        10000.max(k),
        500.max(3 * k),
        args[1],
        args[2]
    )?;
    for w in [
        &mut summary,
        &mut blocks,
        &mut events,
        &mut cov,
        &mut bed,
        &mut fa,
        &mut unresolved,
        &mut report,
        &mut run,
    ] {
        w.flush()?;
    }
    fs::write(dir.join("SUCCESS"), b"Completed\n")?;
    Ok(())
}
