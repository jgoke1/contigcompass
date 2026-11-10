//! Conservative, bounded greedy chaining. These are anchor chains, not alignments.
use crate::matching::{AnchorMatch, Strand};

#[derive(Debug, Clone)]
pub struct Chain {
    pub hits: Vec<AnchorMatch>,
}
impl Chain {
    pub fn strand(&self) -> Strand {
        self.hits[0].strand
    }
    pub fn query_start(&self) -> usize {
        self.hits[0].query_position
    }
    pub fn query_end(&self, k: usize) -> usize {
        self.hits.last().unwrap().query_position + k
    }
    pub fn ref_start(&self) -> usize {
        self.hits
            .iter()
            .map(|h| h.reference_position)
            .min()
            .unwrap()
    }
    pub fn ref_end(&self, k: usize) -> usize {
        self.hits
            .iter()
            .map(|h| h.reference_position)
            .max()
            .unwrap()
            + k
    }
    pub fn covered(&self, k: usize) -> usize {
        covered(
            &self
                .hits
                .iter()
                .map(|h| (h.query_position, h.query_position + k))
                .collect::<Vec<_>>(),
        )
    }
}

pub fn merged(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (s, e) in sorted {
        if s >= e {
            continue;
        }
        if let Some(last) = out.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
        }
        out.push((s, e));
    }
    out
}
pub fn covered(intervals: &[(usize, usize)]) -> usize {
    merged(intervals).iter().map(|(s, e)| e - s).sum()
}
pub fn gaps(intervals: &[(usize, usize)], length: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut end = 0;
    for (s, e) in merged(intervals) {
        if s > end {
            out.push((end, s));
        }
        end = end.max(e);
    }
    if end < length {
        out.push((end, length));
    }
    out
}

/// Thin to non-overlapping query anchors, then extend the best compatible recent chain.
/// Bounded to the 32 most recently created chains; not a global optimal alignment.
pub fn chains(hits: &[AnchorMatch], k: usize, max_gap: usize, drift: usize) -> Vec<Chain> {
    let mut out: Vec<Chain> = Vec::new();
    let mut next = 0;
    for &hit in hits {
        if hit.query_position < next {
            continue;
        }
        next = hit.query_position + k;
        let mut best = None;
        for i in (out.len().saturating_sub(32)..out.len()).rev() {
            let last = out[i].hits.last().unwrap();
            if hit.strand != last.strand {
                continue;
            }
            let dq = hit.query_position - last.query_position;
            let dr = match hit.strand {
                Strand::Forward => hit.reference_position.checked_sub(last.reference_position),
                Strand::Reverse => last.reference_position.checked_sub(hit.reference_position),
            };
            if let Some(dr) = dr {
                if dr >= k && dq <= max_gap && dr <= max_gap && dq.abs_diff(dr) <= drift {
                    let score = dq.abs_diff(dr);
                    if best.is_none_or(|(_, s)| score < s) {
                        best = Some((i, score));
                    }
                }
            }
        }
        if let Some((i, _)) = best {
            out[i].hits.push(hit);
        } else {
            out.push(Chain { hits: vec![hit] });
        }
    }
    out.sort_by_key(|c| c.query_start());
    out
}

#[derive(Debug)]
pub struct Placement {
    pub status: &'static str,
    pub strand: Option<Strand>,
    pub position: Option<usize>,
    pub support: f64,
}
/// Only one dominant chain may drive export. Substantial secondary blocks prevent placement.
pub fn decide(
    chains: &[Chain],
    query_len: usize,
    k: usize,
    min_bp: usize,
    min_fraction: f64,
) -> Placement {
    let mut ranked: Vec<_> = chains
        .iter()
        .filter(|c| c.hits.len() >= 3 && c.covered(k) >= min_bp)
        .collect();
    ranked.sort_by_key(|c| std::cmp::Reverse(c.covered(k)));
    let Some(best) = ranked.first() else {
        return Placement {
            status: "insufficient_support",
            strand: None,
            position: None,
            support: 0.0,
        };
    };
    let total: usize = chains.iter().map(|c| c.covered(k)).sum();
    let support = best.covered(k) as f64 / total.max(1) as f64;
    let fraction = best.covered(k) as f64 / query_len.max(1) as f64;
    let competing = ranked
        .iter()
        .skip(1)
        .any(|c| c.covered(k) as f64 >= best.covered(k) as f64 * 0.1);
    let status = if competing || support < 0.9 {
        "ambiguous_blocks"
    } else if fraction < min_fraction {
        "partial_support"
    } else {
        "placed"
    };
    Placement {
        status,
        strand: if status == "placed" {
            Some(best.strand())
        } else {
            None
        },
        position: if status == "placed" {
            Some(best.ref_start())
        } else {
            None
        },
        support,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn h(q: usize, r: usize, s: Strand) -> AnchorMatch {
        AnchorMatch {
            query_position: q,
            reference_position: r,
            strand: s,
        }
    }
    #[test]
    fn intervals_do_not_double_count() {
        assert_eq!(covered(&[(0, 31), (1, 32), (60, 70)]), 42);
        assert_eq!(gaps(&[(10, 20), (15, 30)], 40), vec![(0, 10), (30, 40)]);
    }
    #[test]
    fn forward_with_indel() {
        let v = vec![
            h(0, 100, Strand::Forward),
            h(31, 131, Strand::Forward),
            h(62, 198, Strand::Forward),
        ];
        let c = chains(&v, 31, 1000, 50);
        assert_eq!(c.len(), 1);
        assert_eq!(decide(&c, 93, 31, 60, 0.5).status, "placed");
    }
    #[test]
    fn reverse_placement() {
        let v = vec![
            h(0, 200, Strand::Reverse),
            h(31, 169, Strand::Reverse),
            h(62, 138, Strand::Reverse),
        ];
        let c = chains(&v, 31, 1000, 50);
        assert_eq!(decide(&c, 93, 31, 60, 0.5).strand, Some(Strand::Reverse));
    }
    #[test]
    fn sparse_matches_not_whole_contig() {
        let v = vec![
            h(0, 100, Strand::Forward),
            h(1000, 1100, Strand::Forward),
            h(2000, 2100, Strand::Forward),
        ];
        let c = chains(&v, 31, 10000, 50);
        assert_eq!(decide(&c, 5000, 31, 60, 0.5).status, "partial_support");
    }
    #[test]
    fn mixed_blocks_stay_unresolved() {
        let v = vec![
            h(0, 100, Strand::Forward),
            h(31, 131, Strand::Forward),
            h(62, 162, Strand::Forward),
            h(93, 900, Strand::Reverse),
            h(124, 869, Strand::Reverse),
            h(155, 838, Strand::Reverse),
        ];
        let c = chains(&v, 31, 1000, 50);
        assert_eq!(decide(&c, 186, 31, 60, 0.5).status, "ambiguous_blocks");
    }
    #[test]
    fn circular_boundary_splits_conservatively() {
        let v = vec![
            h(0, 938, Strand::Forward),
            h(31, 969, Strand::Forward),
            h(62, 0, Strand::Forward),
        ];
        assert_eq!(chains(&v, 31, 1000, 50).len(), 2);
    }
    #[test]
    fn no_matches() {
        assert_eq!(decide(&[], 100, 31, 60, 0.5).status, "insufficient_support");
    }
}
