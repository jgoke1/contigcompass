use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
fn workspace() -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "compass-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    p
}
fn genome(n: usize) -> String {
    let mut x = 123456789u64;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x % 4) as usize] as char
        })
        .collect()
}
fn rc(s: &str) -> String {
    String::from_utf8(contigcompass::dna::reverse_complement(s.as_bytes())).unwrap()
}
fn run(p: &PathBuf, extra: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_place"))
        .args([
            p.join("ref.fa").to_str().unwrap(),
            p.join("query.fa").to_str().unwrap(),
            "31",
            p.join("output").to_str().unwrap(),
        ])
        .args(extra)
        .output()
        .unwrap()
}
#[test]
fn exported_sequences_coverage_and_unresolved_are_correct() {
    let p = workspace();
    let r = genome(12000);
    fs::write(p.join("ref.fa"), format!(">reference\n{r}\n")).unwrap();
    let reverse = rc(&r[6000..9000]);
    let forward = &r[1000..4000];
    fs::write(
        p.join("query.fa"),
        format!(
            ">reverse\n{reverse}\n>forward\n{forward}\n>unknown\n{}\n",
            "N".repeat(700)
        ),
    )
    .unwrap();
    let result = run(&p, &[]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let t = fs::read_to_string(p.join("output/placements.tsv")).unwrap();
    assert!(t.contains("\tplaced\t-\t6024\t"));
    assert!(t.contains("\tplaced\t+\t1000\t"));
    assert!(t.contains("no_matching_anchors"));
    let fasta = fs::read_to_string(p.join("output/ordered_oriented.fasta")).unwrap();
    let flat = fasta.lines().collect::<String>();
    assert!(flat.starts_with(&format!(">forward{forward}>reverse{}", &r[6000..9000])));
    let coverage = fs::read_to_string(p.join("output/coverage.tsv")).unwrap();
    assert!(coverage.contains("reference\t12000\t6000\t50.000\t6700\t6000"));
    assert!(p.join("output/SUCCESS").exists());
    assert!(!run(&p, &[]).status.success());
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn indels_can_place_but_inversions_stay_unresolved() {
    let p = workspace();
    let r = genome(20000);
    fs::write(p.join("ref.fa"), format!(">r\n{r}\n")).unwrap();
    let insertion = format!("{}{}{}", &r[1000..2500], "A".repeat(50), &r[2500..4000]);
    let mixed = format!("{}{}", &r[5000..7000], rc(&r[9000..11000]));
    fs::write(
        p.join("query.fa"),
        format!(">indel\n{insertion}\n>mixed\n{mixed}\n"),
    )
    .unwrap();
    let result = run(&p, &[]);
    assert!(result.status.success());
    let t = fs::read_to_string(p.join("output/placements.tsv")).unwrap();
    assert!(
        t.lines()
            .find(|s| s.starts_with("indel\t"))
            .unwrap()
            .contains("\tplaced\t")
    );
    assert!(
        t.lines()
            .find(|s| s.starts_with("mixed\t"))
            .unwrap()
            .contains("ambiguous_blocks")
    );
    assert!(
        fs::read_to_string(p.join("output/candidate_events.tsv"))
            .unwrap()
            .contains("orientation_switch")
    );
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn circular_crossing_is_flagged_and_not_forced() {
    let p = workspace();
    let r = genome(12000);
    fs::write(p.join("ref.fa"), format!(">r\n{r}\n")).unwrap();
    fs::write(
        p.join("query.fa"),
        format!(">wrap\n{}{}\n", &r[10000..], &r[..2000]),
    )
    .unwrap();
    assert!(run(&p, &["--circular"]).status.success());
    assert!(
        fs::read_to_string(p.join("output/candidate_events.tsv"))
            .unwrap()
            .contains("possible_circular_boundary")
    );
    assert!(
        fs::read_to_string(p.join("output/placements.tsv"))
            .unwrap()
            .contains("ambiguous_blocks")
    );
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn duplicate_ids_rejected_before_output() {
    let p = workspace();
    fs::write(p.join("ref.fa"), format!(">r\n{}\n", genome(100))).unwrap();
    fs::write(p.join("query.fa"), ">x\nACGT\n>x\nACGT\n").unwrap();
    assert!(!run(&p, &[]).status.success());
    assert!(!p.join("output").exists());
    fs::remove_dir_all(p).unwrap();
}
