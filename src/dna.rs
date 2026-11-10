/// Return the reverse complement of an IUPAC DNA sequence.
///
/// Input is expected to contain validated DNA symbols.
/// Output is always uppercase.
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            other => panic!("Invalid DNA byte: {}", other),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::reverse_complement;

    #[test]
    fn complements_and_reverses_dna() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
    }

    #[test]
    fn accepts_lowercase() {
        assert_eq!(reverse_complement(b"atgc"), b"GCAT");
    }

    #[test]
    fn handles_ambiguous_bases() {
        assert_eq!(reverse_complement(b"ARYKMBVDHSWN"), b"NWSDHBVKMRYT");
    }

    #[test]
    fn handles_empty_sequence() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn reversing_twice_restores_original() {
        let original = b"ATGCRYN";
        let reversed = reverse_complement(original);

        assert_eq!(reverse_complement(&reversed), original);
    }
}
