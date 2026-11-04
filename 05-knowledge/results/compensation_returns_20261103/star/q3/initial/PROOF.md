# Completing the other q3 endcap at order13

2026-11-03. **FINITE-EXACT complete initial-cap star classification.**
For D=C3->TT9, no strong one-apex extension with H1409 or1427 has a vertex
of start defect1 or2. This remains true when any vertex of end defect0 may
serve as the prospective outer root. The exact target counts are attained;
the failure is the positive-entry condition.

## 1. The last carrier and its precise scope

The closest mechanism is the parent [deleted-joint formula](../PROOF.md).
The terminal-cap hostile attains both counts but lacks positive entries;
the corrected near miss is the relocated tail, which supplies entries but
forces a larger path count. The remaining carrier reverses their positions:

    C3(201->401->400->201) -> TT(99,10,11,2,1,3,4,12,200).

The physical vertex set is exactly that of the original thirteen-vertex
C2 children, with retained inner mark0. Moving the cap before the chain
changes27 specified non-star pairs; no vertex is deleted or identified.
This is a fixed-carrier classification, not a native conversion program.

The exact formula

    H(D+w)=sum_i(1-x_i)E_i + sum_j x_jS_j
                +sum_(i,j)(1-x_i)x_jJ_ij

classifies all4096 actual star cuts x_i=[w->i]. The coefficients are freshly
counted on the literal deleted graph D, whose Hamilton count is3. Independent
directed reachability gives1792 strong extensions. Only strong candidates
whose count is1409 or1427 are sent to the full literal subset path profiler;
the two unchanged original children serve as additional positive controls.

## 2. Result, with both port types retained

| H | Strong target-count graphs | With some end-zero port | With some start1/2 port | With both |
|---:|---:|---:|---:|---:|
| 1409 | 6 | 6 | 0 | 0 |
| 1427 | 3 | 3 | 0 | 0 |

Every target graph has a possible root-end-zero port. None has a positive
entry anywhere, so choosing a different root or relabelling physical99
cannot meet the joint requirement. All nine are literally counted, with
the full E1,S1,E2,S2,J matrix and fresh named field observers. The observers
are not used as the acceptance predicate.

`census.json` saves all4096 words, the physical bit order, deleted graph,
deleted coefficient profile, strongness, and exact H. `witnesses.json`
saves every target's actual rows/names, cut, every start/end defect, all
candidate root names, and full joint/phase sidecars. `receipt.json` records
the totals and input hashes. No claim is made about arbitrary source
transformations, other deletion counts, or other tournament orders by this
finite packet alone.

Together with the terminal and interior carrier arguments, this is the
last endcap needed for a separately justified classification of q3 deleted
SCC shapes at this order. That combination also needs its SCC theorem and
the name-preserving order bound; neither is inferred from the nine graphs.

## 3. Reproduction

    python 05-knowledge/results/compensation_returns_20261103/star/q3/initial/reproduce.py

Normal and `-O` outputs must be identical. The parent q3 packet is pinned
and remains byte-for-byte unchanged. The complete inverse universe is a
4096-cut family, not an arbitrary graph-order or word-length cutoff.
