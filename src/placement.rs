use crate::matching::{AnchorMatch, Strand};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlacementStatus {
    InsufficientAnchors,
    ConsistentForward,
    ConsistentReverse,
    ConflictingAnchors,
}

/// Check whether anchors support one exact coordinate relationship.
///
/// Requires at least two distinct query positions.
/// This strict version does not accommodate indels or circular wrapping.
pub fn assess_consistency(matches: &[AnchorMatch]) -> PlacementStatus {
    let Some(first) = matches.first() else {
        return PlacementStatus::InsufficientAnchors;
    };

    // Repeated copies of the same query position do not add evidence.
    if !matches
        .iter()
        .any(|hit| hit.query_position != first.query_position)
    {
        return PlacementStatus::InsufficientAnchors;
    }

    let expected = coordinate_key(first);

    for hit in matches {
        if hit.strand != first.strand || coordinate_key(hit) != expected {
            return PlacementStatus::ConflictingAnchors;
        }
    }

    match first.strand {
        Strand::Forward => PlacementStatus::ConsistentForward,
        Strand::Reverse => PlacementStatus::ConsistentReverse,
    }
}

fn coordinate_key(hit: &AnchorMatch) -> i128 {
    // Convert before arithmetic to allow negative forward offsets.
    let reference = hit.reference_position as i128;
    let query = hit.query_position as i128;

    match hit.strand {
        Strand::Forward => reference - query,
        Strand::Reverse => reference + query,
    }
}

#[cfg(test)]
mod tests {
    use super::{PlacementStatus, assess_consistency};
    use crate::matching::{AnchorMatch, Strand};

    fn anchor(query: usize, reference: usize, strand: Strand) -> AnchorMatch {
        AnchorMatch {
            query_position: query,
            reference_position: reference,
            strand,
        }
    }

    #[test]
    fn recognises_consistent_forward_anchors() {
        let hits = vec![
            anchor(0, 100, Strand::Forward),
            anchor(5, 105, Strand::Forward),
        ];

        assert_eq!(
            assess_consistency(&hits),
            PlacementStatus::ConsistentForward
        );
    }

    #[test]
    fn recognises_consistent_reverse_anchors() {
        let hits = vec![
            anchor(0, 105, Strand::Reverse),
            anchor(5, 100, Strand::Reverse),
        ];

        assert_eq!(
            assess_consistency(&hits),
            PlacementStatus::ConsistentReverse
        );
    }

    #[test]
    fn detects_incompatible_positions_on_same_strand() {
        let hits = vec![
            anchor(0, 100, Strand::Forward),
            anchor(5, 200, Strand::Forward),
        ];

        assert_eq!(
            assess_consistency(&hits),
            PlacementStatus::ConflictingAnchors
        );
    }

    #[test]
    fn detects_mixed_strands() {
        let hits = vec![
            anchor(0, 100, Strand::Forward),
            anchor(5, 95, Strand::Reverse),
        ];

        assert_eq!(
            assess_consistency(&hits),
            PlacementStatus::ConflictingAnchors
        );
    }

    #[test]
    fn requires_distinct_query_positions() {
        let hit = anchor(0, 100, Strand::Forward);

        for hits in [vec![], vec![hit], vec![hit, hit]] {
            assert_eq!(
                assess_consistency(&hits),
                PlacementStatus::InsufficientAnchors
            );
        }
    }
}
