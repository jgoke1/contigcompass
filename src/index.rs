use std::collections::HashMap;

use crate::dna::reverse_complement;
use crate::kmers::iter_kmers;

/// A unique anchor in one reference sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchor {
    pub position: usize,
    pub canonical_reversed: bool,
}

/// Return the canonical sequence and whether reverse complementation
/// was needed to obtain it.
///
/// Palindromic k-mers return None because their strand is ambiguous.
pub fn canonical_kmer(kmer: &[u8]) -> Option<(Vec<u8>, bool)> {
    let forward = kmer.to_ascii_uppercase();
    let reverse = reverse_complement(&forward);

    if forward == reverse {
        return None;
    }

    if forward < reverse {
        Some((forward, false))
    } else {
        Some((reverse, true))
    }
}

/// Index k-mers that occur exactly once across both orientations.
///
/// This initial version accepts ONE reference sequence.
/// Positions are zero-based within that sequence.
pub fn build_reference_index(
    sequence: &[u8],
    k: usize,
) -> Result<HashMap<Vec<u8>, Anchor>, &'static str> {
    // Some(anchor) means seen once; None means seen multiple times.
    let mut candidates: HashMap<Vec<u8>, Option<Anchor>> = HashMap::new();

    for (position, kmer) in iter_kmers(sequence, k)? {
        let Some((key, canonical_reversed)) = canonical_kmer(&kmer) else {
            continue;
        };

        let anchor = Anchor {
            position,
            canonical_reversed,
        };

        candidates
            .entry(key)
            .and_modify(|existing| *existing = None)
            .or_insert(Some(anchor));
    }

    // Keep only keys that were observed exactly once.
    let index = candidates
        .into_iter()
        .filter_map(|(key, anchor)| anchor.map(|value| (key, value)))
        .collect();

    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::{Anchor, build_reference_index, canonical_kmer};

    #[test]
    fn reverse_complements_share_a_key() {
        let (forward_key, forward_reversed) = canonical_kmer(b"ATG").unwrap();

        let (reverse_key, reverse_reversed) = canonical_kmer(b"CAT").unwrap();

        assert_eq!(forward_key, reverse_key);
        assert!(!forward_reversed);
        assert!(reverse_reversed);
    }

    #[test]
    fn stores_unique_positions_and_orientations() {
        let index = build_reference_index(b"ATGC", 3).unwrap();

        assert_eq!(index.len(), 2);

        assert_eq!(
            index.get(b"ATG".as_slice()),
            Some(&Anchor {
                position: 0,
                canonical_reversed: false,
            })
        );

        // TGC at position 1 becomes the canonical key GCA.
        assert_eq!(
            index.get(b"GCA".as_slice()),
            Some(&Anchor {
                position: 1,
                canonical_reversed: true,
            })
        );
    }

    #[test]
    fn excludes_repeats_across_both_strands() {
        // ATG/CAT share one key; TGC/GCA share another.
        let index = build_reference_index(b"ATGCAT", 3).unwrap();
        assert!(index.is_empty());

        // A third occurrence must not restore a repeated key.
        let repeated = build_reference_index(b"AAAAA", 3).unwrap();
        assert!(repeated.is_empty());
    }

    #[test]
    fn excludes_palindromes_and_rejects_zero_k() {
        assert!(canonical_kmer(b"ATAT").is_none());
        assert!(build_reference_index(b"ATAT", 4).unwrap().is_empty());
        assert!(build_reference_index(b"ATGC", 0).is_err());
    }
}

/// A reference index that keeps its k-mer size with its anchors.
pub struct ReferenceIndex {
    k: usize,
    anchors: HashMap<Vec<u8>, Anchor>,
}

impl ReferenceIndex {
    pub fn new(sequence: &[u8], k: usize) -> Result<Self, &'static str> {
        Ok(Self {
            k,
            anchors: build_reference_index(sequence, k)?,
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn get(&self, key: &[u8]) -> Option<&Anchor> {
        self.anchors.get(key)
    }
}
