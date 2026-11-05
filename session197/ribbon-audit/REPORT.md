# Independent ribbon-transfer audit

2026-11-05. **PASS for the stated theorem and bounded certificates.**

The audited packet is `session197/ribbon/PROOF.md`, `transfer.py`,
`explore.py`, revised `verify.py`, the three `closure-*.json` files, and
`VALIDATION.json`. Source copies and hashes are retained in this audit folder.
No edits were made to the author's source directory. The audit used only this
named mathematical source, not earlier session claims.

## Independent verification

- Reconstructed every transition by lifting endpoint states to concrete
  directed graphs, adding the candidate incident edges, recomputing graph
  components, rejecting cycles, and projecting the endpoints back to the
  frontier. All 10,210 transition rows agree through widths 1, 2, 3, and 4.
  The independently obtained state counts are 7, 27, 119, and 571.
- Enumerated all 501 directed linear forests on five labeled vertices by
  permutations and cuts. Every coordinate of each path-cover vector matches
  the transfer, not merely its final Hamiltonian-path sum. The comparison
  covers all 16 width-1, 128 width-2, 512 width-3, and 1,024 width-4 order-5
  tournament words with irrelevant padding bits fixed to zero.
- Rebuilt the bounded closures using the independent graph transitions and
  a deliberately different coordinate ordering. All queues exhaust, with
  processed count equal to reached count:
  - Width 1, bound 101: 66 vectors.
  - Width 2, bound 47: 1,488 vectors, odd gaps exactly 7, 21, 47.
  - Width 2, bound 1001: 192,274 vectors, odd gaps exactly 7, 21, 47.
  - Width 3, bound 101: 100,168 vectors, odd gaps exactly 7, 21.
- All 568 stored representatives across the three certificates were checked
  with a separately implemented predecessor-removal subset recurrence. All
  values and shortest representative lengths agree with independent BFS.
- Every width-1 word through order 12 (4,095 words in total) agrees with the
  ordered-join product formula and the claimed strong-connectivity criterion.
  Literal permutations of `R_n` through order 8 and direct composition
  enumeration through order 17 verify the tribonacci sequence independently.
  Its multiplicative monoid through 101 agrees with the width-1 closure.
- 3,500 randomized full-vector tests through length 25 verify saturation
  against exact integer vectors. This supports, but does not replace, the
  arbitrary-length algebraic proof below.
- 1,200 separate distance-based graph constructions agree with the stated
  mask convention, including initial padding and all forward long edges.
- Literal permutation enumeration of the explicit five-reversal order-7
  tournament gives exactly 47 paths, and reachability gives one SCC. Thus
  the number 47 is realized at width 3, and the all-order width-2 exclusion
  makes its minimum realizing backward span exactly 3.
- The author's revised verifier was run from a copied source directory,
  writing only audit-local output, and passed in full.

## Proof audit

The frontier state records exactly the free endpoint data of a path cover.
Forgotten endpoints may collapse as labels, but different components do not
collapse: each has its own retained start, except for the unique permissible
forgotten start. A component with both endpoints forgotten is still live:
its end may point to the next vertex. It must not be discarded as a closed
component. The audited code keeps it. For example, from state `((-1,-1),)`
at width 2 and mask 3, the transition to `((-1,1),)` exists with weight one.

A forgotten start cannot receive any future incoming edge, because all
future-to-forgotten pairs have forward orientation. Two such starts cannot
occur in a restriction of a future Hamiltonian path. Internal vertices need
no state information because their selected indegree and outdegree are both
already occupied. Distinct histories may merge into one state, and the code
correctly retains multiplicities; coefficients greater than one are present
and independently verified.

For cap `C=M+1`, `q(x)=min(C,x)` obeys
`q(x+y)=q(q(x)+q(y))` and `q(cx)=q(c*q(x))` for all nonnegative integers.
Consequently coordinatewise saturation is exact after every finite word,
by induction, not only through a tested length. The terminal sum is exact
whenever it is at most M, and is greater than M otherwise. The deterministic
insertion proof of induced-subtournament Hamiltonian-path monotonicity is
valid and injective, including one-vertex edge cases. Therefore pruning a
prefix above M cannot remove a later output at most M. Completed finite
reachability proves the claimed bounds at every order.

The width-1 decreasing-run composition bijection is valid. Runs are consecutive
integer intervals; an attempted transition to a lower interval is either a
forbidden longer descent or would continue the same maximal run. The only
forbidden forward transition between successive intervals occurs when both
parts have size one. The generating function and its initial coefficients
give exactly the stated tribonacci recurrence. Every forward adjacent edge
is a genuine all-forward join cut at width 1, so the multiplicative-monoid
classification is all-order. This join argument is not claimed at width 2.

All coverage conclusions are explicitly for arbitrary tournaments. Negative
conclusions inherit to the strong subclass; positive coverage is not silently
upgraded to strong coverage. The sole claimed strong 47 example is verified.
There is no finite-to-infinite extrapolation of the *target-value* interval:
the global coverage questions remain explicitly open.

## Concrete warning against an interval-SCC shortcut

At width 2, word `0033` gives backward pairs
`(1,3),(2,3),(2,4),(3,4)`, with all other pairs forward.
Its SCCs are `{1,3,4}` and `{2}`. The first is the directed cycle
`1 -> 4 -> 3 -> 1`, and all its vertices point to 2. Every natural-order
cut is crossed by a backward edge, yet the tournament is not strong, and
one SCC is not an interval. No interval-SCC test appears in the audited code.

## Resolved and limited issues

The original verifier only rechecked regenerated representatives while the
prose said stored representatives. This was reported to the author. The
revised `verify.py` now loads existing certificate files, compares their
completeness, state counts, and exact value-to-word maps with regeneration,
and rechecks their saved words. The revised file passed independently.

Future-use limitation, not an error in any present claim: `explore.py` writes
words by concatenating decimal mask strings. This is unambiguous for the
claimed widths 1–3, but would need a delimiter or list encoding for width 4
or above when masks reach 10. No width-4 bounded closure is claimed.

## Reproduction

From this audit folder, run `python audit.py` and `python edge_checks.py`.
The original source verifier was also run as
`PYTHONPATH=source python source/verify.py`, after creating the audit-local
`session197/ribbon/` output directory. Results are in `AUDIT.json`,
`EDGE_CHECKS.json`, `independent-closure-*.json`, and the associated logs.
