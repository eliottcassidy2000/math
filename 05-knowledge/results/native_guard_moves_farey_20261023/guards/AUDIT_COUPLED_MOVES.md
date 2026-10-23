# Independent audit of native pinned cyclic-cell moves

23 October 2026. Reviewer: original-guard/attachment lane.

**PASS after a receipt-format correction.** Reviewed all of
[`COUPLED_PHASE_MOVES.md`](../moves/COUPLED_PHASE_MOVES.md) and its checker.
The proof establishes the stated all-size regular-host family results;
electrical conclusions beyond the supplied actual order-eight example
remain explicitly conditional on guard factorization.

## Proof audit

1. The host is a genuine tournament for every odd `q,t>=3`, including
   composite `q`; its displayed score is `(qt-1)/2`. The chosen native
   contiguous-cell paths have exactly binary consecutive phase differences.
   Fixed endpoint owners determine the actual endpoint phases.
2. A one-cell change preserves the sum of the two adjacent binary
   differences modulo odd `q>=3`; their possible integer sums `0,1,2`
   are distinct residues. The only nontrivial replacement is `01<->10`.
   Thus integer weight, not merely endpoint residue, exactly indexes the
   one-cell components, and adjacent swaps connect each layer.
3. A cumulative difference walk with increments `-1,0,1` and displacement
   `gq` must visit each intermediate integer. At least `g(q-1)` visited
   levels are nonzero modulo `q`, at distinct actual internal clocks.
   This proves the lower bound even for arbitrary nonlocal coupled moves.
   Canonical seam flips attain it. There is no assumption of legal
   one-cell intermediate states during an atomic block event.
4. Consequently gaps between accepted entire weight layers give the exact
   scalar connectivity budget. For internal pins, the interval supports
   are disjoint, so exact costs add as
   `(q-1)/q * sum_j |k_j-l_j|`. The accepted weight-vector quotient is
   equivalent to connectivity of the actual state graph for budget at
   least one, because the within-layer graphs are connected. Singleton
   and empty boundaries are handled separately in the text.
5. These quotient statements require acceptance of ENTIRE weight layers.
   The proof does not silently replace the actual electrical character by
   an integer weight. The order-27 correlated guard is plainly declared
   synthetic. The order-eight example instead rebuilds actual arcs,
   source/target owners, original harmonic values, and real endpoint pins.
6. Converse transports the actual paths and indexed pins and preserves
   support size. The common-phase electrical example stays inside one
   actual automorphism orbit; the larger seam moves need not do so.
   Neither operation is asserted to be an elementary metagraph edge.

## Computation audit and repair

The first normal and optimized replays reached the final receipt equality
and failed because in-memory tuples were compared with JSON lists. The
author normalized the computed result through JSON before comparing it.
No mathematical predicate or frozen numerical value needed repair.
Both repaired frozen replays pass **35,950** checks, comprising thirteen
actual host families, 103 weight filters, the correlated-pin hostile, and
the independently reconstructed original order-eight electrical example.

An independent phase-first control is saved as
[`audit_coupled_controls.py`](audit_coupled_controls.py). It enumerates
cyclic phase assignments directly, filters binary differences and real
pins, then compares every cross-layer minimum and within-layer graph.
The four cases are `(q,r)=(3,6),(3,8),(3,8),(5,6)`, with two different
nonzero/zero pin systems at `(3,8)`. Normal and optimized runs both pass
**33** comparisons on respectively `4,25,25,7` retained states. This
independent finite check supports the proof; it is not its replacement.

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/moves/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/moves/check.py
python -B 05-knowledge/results/native_guard_moves_farey_20261023/guards/audit_coupled_controls.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/guards/audit_coupled_controls.py
```

The remaining mathematical obligation is authenticating stronger original
electrical guard factorizations on the larger regular native hosts. The
unfiltered all-size connectivity result does not discharge that obligation.
