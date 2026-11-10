/// Extract uppercase A/C/G/T k-mers and their zero-based positions.
///
/// Ambiguous windows are skipped. Positions still refer to the
/// original sequence. Returns an error when k is zero.
pub fn extract_kmers(sequence: &[u8], k: usize) -> Result<Vec<(usize, Vec<u8>)>, &'static str> {
    if k == 0 {
        return Err("k must be greater than zero");
    }

    let uppercase = sequence.to_ascii_uppercase();

    let kmers = uppercase
        .windows(k)
        .enumerate()
        .filter(|(_, window)| {
            window
                .iter()
                .all(|base| matches!(*base, b'A' | b'C' | b'G' | b'T'))
        })
        .map(|(position, window)| (position, window.to_vec()))
        .collect();

    Ok(kmers)
}

#[cfg(test)]
mod tests {
    use super::extract_kmers;

    #[test]
    fn extracts_overlapping_kmers() {
        let result = extract_kmers(b"ATGCA", 3).unwrap();

        assert_eq!(
            result,
            vec![
                (0, b"ATG".to_vec()),
                (1, b"TGC".to_vec()),
                (2, b"GCA".to_vec()),
            ]
        );
    }

    #[test]
    fn skips_ambiguous_windows_preserving_positions() {
        let result = extract_kmers(b"ATGNACGT", 3).unwrap();

        assert_eq!(
            result,
            vec![
                (0, b"ATG".to_vec()),
                (4, b"ACG".to_vec()),
                (5, b"CGT".to_vec()),
            ]
        );
    }

    #[test]
    fn converts_lowercase_to_uppercase() {
        let result = extract_kmers(b"atgc", 3).unwrap();

        assert_eq!(result, vec![(0, b"ATG".to_vec()), (1, b"TGC".to_vec())]);
    }

    #[test]
    fn returns_no_kmers_when_sequence_is_too_short() {
        assert!(extract_kmers(b"AT", 3).unwrap().is_empty());
    }

    #[test]
    fn rejects_zero_k() {
        assert!(extract_kmers(b"ATGC", 0).is_err());
    }
}

/// Streaming counterpart used by genome-scale indexing and matching.
pub fn iter_kmers(
    sequence: &[u8],
    k: usize,
) -> Result<impl Iterator<Item = (usize, Vec<u8>)> + '_, &'static str> {
    if k == 0 {
        return Err("k must be greater than zero");
    }
    Ok(sequence
        .windows(k)
        .enumerate()
        .filter(|(_, w)| {
            w.iter()
                .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        })
        .map(|(i, w)| (i, w.to_ascii_uppercase())))
}
