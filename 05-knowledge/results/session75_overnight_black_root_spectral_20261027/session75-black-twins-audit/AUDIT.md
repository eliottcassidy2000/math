# Bounded final addendum audit: the two order10 witnesses

**Verdict: ACCEPTED in the scopes below.** 2026-10-27.
The earlier Phase2 producer milestone and independent audit remain unchanged.
No repository files or remote state were modified.

## What was independently checked

Exactly the two named matrices in the producer's frozen
`/tmp/session75-black-twins` packet were checked. Its manifest SHA256 is
`dcb9b8d6d65122e7e47dba0f9a3984c145bdcf67fbf0dd059db4892b54796d2b`;
all eight listed files match. Incoming theorem/recipe/output bytes match git
commit `87d61f67673a6b10cbafff1bf9e1377f3c41e56c` exactly.

The independent checker imports no producer module. Its adjacency-matrix
primitives are an exact copy of the independently audited Phase2 helper,
bound in `SOURCE_BINDINGS.json`. The new cycle computation uses a different
carrier: all permutations p of the ten named vertices satisfying either
p(v)=v or the arc v->p(v). Fixed points represent unused vertices. Rejecting
any nontrivial even cycle leaves exactly the unordered collections of
disjoint directed odd cycles. Each permutation determines its collection
uniquely; each collection determines one permutation. Rotation and component
ordering therefore introduce no overcounting, while reversal is not merged.

The enumeration gives514 allowed directed permutation covers for the old
witness, of which224 have only odd nontrivial cycles; for E10 it gives517
covers, of which262 meet that condition. It independently reproduces the full
multivariate polynomials:

```
F_old = 1 + 14x3 + 30x5 + 52x7 + 30x9
          + 35x3^2 + 48x3*x5 + 6x3*x7 + 4x5^2 + 4x3^3
F_E10 = 1 + 14x3 + 54x5 + 90x7 + 36x9
          + 22x3^2 + 36x3*x5 + 9x5^2.
```

Every retained producer literal odd-cycle word is checked against the
independent permutation-cycle set:126 distinct words for the old witness,
194 for E10. This is an exact list comparison, not just agreement of totals.
Consequently the component-count polynomials are

```
Q_old(t) = 1 + 126t + 93t^2 + 4t^3
Q_E10(t) = 1 + 194t + 67t^2
Q_E10(t)-Q_old(t) = 2t(2-t)(2t+17).
```

Both evaluate to657 at t=2, independently agreeing with exhaustive Hamilton
word counts. This verifies the numerical evaluation identity for these two
objects; it is not an audit of a general cycle-polynomial theorem.

## Named cuts, colors, and the triangle warning

All1,022 named outgoing cuts agree, as do the named degree vectors. The
independent Hamilton replay gives33 blue/624 black words and six black
one-step exits for the old witness, versus99 blue/558 black words and no
black one-step exits for E10. Both literal two-black core routes are replayed
from their actual named sources, with native K endpoint matrices checked
against the producer. Minimum cuts and semidegrees are independently1 and2.
Thus black nonstrong distances are1 and2 respectively; both core distances
are2 by the endpoint lower bound and exact routes. The root independently
verified the automorphism counts1 and9; this addendum does not need those
counts to establish the polynomial or distance distinction.

**The labelled cyclic triples are NOT equal.** Each tournament has scalar
c3=14, but only two cyclic vertex triples are shared; twelve occur only in
each source. The receipt lists both differences. Equality of all named cuts
and of the scalar triangle count must not be stated as equality of the
labelled cyclic-triple sets.

## Incoming E10 reconciliation and credit

Independently transcribing the pinned `stack`, `C3sub`, and `doubled` recipe
only for M=C3 yields rows

```
(974,1012,1016,1010,993,513,288,96,160,448).
```

The map from the Phase2 E10 labels to incoming labels is
`(4,1,3,2,0,9,6,8,7,5)`. All100 ordered-pair equations are verified directly.
The incoming E10 and Phase2 E10 are therefore the same ordinary isomorphism
type; the old two-R5 witness is not that object. Record concurrent independent
construction/computation; no first-discovery claim is established here.
Only this named recipe and isomorphism are certified. In particular the
incoming full order10 census is NOT certified by this addendum.

## Hypothesis correction and scope

`C3_CORRECTION_NOTE.md` binds the current HYP-9421 text at87d61 to the
already-audited order11 source. Clause(1)'s E-only classification and
no-order11-exception prediction are refuted. Clause(2)'s E-only extra-move
exclusivity is refuted, while its universal lambda+1 upper bound remains
OPEN. The fixed Phase2 witness is re-used; no new search range is introduced.
Global black-core distance=delta remains OPEN.

No statement that F generally determines black distance is made. In fact
these two different F values merely distinguish this pair; their common
value at2 loses that distinction.

## Reproduction and final-byte acceptance

Run `bash run_verification.sh` here. It checks the local inventory and reruns
the independent addendum under normal Python and Python -O, requiring byte
identity with both receipts. Every essential rejection is an explicit error.
`ACCEPTANCE.json` records exact receipt and producer-manifest hashes.
This completes the bounded addendum audit; earlier frozen packets remain
untouched.
