use crate::index::{ReferenceIndex, canonical_kmer};
use crate::kmers::iter_kmers;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// One exact k-mer match to a unique reference anchor.
///
/// Both positions are zero-based starts in the original sequences.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnchorMatch {
    pub query_position: usize,
    pub reference_position: usize,
    pub strand: Strand,
}

/// Find query k-mers that match unique reference anchors.
///
/// Builds the reference index using the same k as the query.
/// Results follow query-position order.
/// Convenience wrapper used by our existing tests.
#[cfg(test)]
pub fn match_query(
    reference: &[u8],
    query: &[u8],
    k: usize,
) -> Result<Vec<AnchorMatch>, &'static str> {
    let index = ReferenceIndex::new(reference, k)?;
    match_indexed_query(&index, query)
}

/// Match a query using an already-built reference index.
pub fn match_indexed_query(
    index: &ReferenceIndex,
    query: &[u8],
) -> Result<Vec<AnchorMatch>, &'static str> {
    let mut matches = Vec::new();

    for (query_position, kmer) in iter_kmers(query, index.k())? {
        let Some((key, query_reversed)) = canonical_kmer(&kmer) else {
            continue;
        };

        if let Some(anchor) = index.get(&key) {
            let strand = if query_reversed == anchor.canonical_reversed {
                Strand::Forward
            } else {
                Strand::Reverse
            };

            matches.push(AnchorMatch {
                query_position,
                reference_position: anchor.position,
                strand,
            });
        }
    }

    Ok(matches)
}
#[cfg(test)]
mod tests {
    use super::{AnchorMatch, Strand, match_query};

    #[test]
    fn reuses_one_index_for_multiple_queries() {
        let index = crate::index::ReferenceIndex::new(b"ATGC", 3).unwrap();

        let forward = super::match_indexed_query(&index, b"ATGC").unwrap();

        let reverse = super::match_indexed_query(&index, b"GCAT").unwrap();

        let unmatched = super::match_indexed_query(&index, b"AAAA").unwrap();

        assert_eq!(forward.len(), 2);
        assert_eq!(reverse.len(), 2);
        assert!(unmatched.is_empty());

        assert!(forward.iter().all(|hit| hit.strand == Strand::Forward));
        assert!(reverse.iter().all(|hit| hit.strand == Strand::Reverse));

        // Reusing the index must not change subsequent results.
        let forward_again = super::match_indexed_query(&index, b"ATGC").unwrap();

        assert_eq!(forward, forward_again);
    }

    #[test]
    fn matches_forward_contig() {
        let matches = match_query(b"ATGC", b"ATGC", 3).unwrap();

        assert_eq!(
            matches,
            vec![
                AnchorMatch {
                    query_position: 0,
                    reference_position: 0,
                    strand: Strand::Forward,
                },
                AnchorMatch {
                    query_position: 1,
                    reference_position: 1,
                    strand: Strand::Forward,
                },
            ]
        );
    }

    #[test]
    fn matches_reverse_complement_contig() {
        let matches = match_query(b"ATGC", b"GCAT", 3).unwrap();

        assert_eq!(
            matches,
            vec![
                AnchorMatch {
                    query_position: 0,
                    reference_position: 1,
                    strand: Strand::Reverse,
                },
                AnchorMatch {
                    query_position: 1,
                    reference_position: 0,
                    strand: Strand::Reverse,
                },
            ]
        );
    }

    #[test]
    fn preserves_query_positions_after_ambiguous_bases() {
        let matches = match_query(b"ATGC", b"NNATG", 3).unwrap();

        assert_eq!(
            matches,
            vec![AnchorMatch {
                query_position: 2,
                reference_position: 0,
                strand: Strand::Forward,
            }]
        );
    }

    #[test]
    fn handles_unmatched_queries_and_invalid_k() {
        assert!(match_query(b"ATGC", b"AAA", 3).unwrap().is_empty());
        assert!(match_query(b"ATGC", b"AT", 3).unwrap().is_empty());
        assert!(match_query(b"ATGC", b"ATG", 0).is_err());
    }
}
