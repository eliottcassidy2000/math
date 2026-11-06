# Independent readiness and descendant-macro audit

Status: PASS for the mathematical classification and macro argument, subject to the explicit existential wording correction below.

Final revision signoff: PASS. The revised proof now states the existential representative quantifier explicitly for each sign, identifies the edge-seeded scope and inherited complete-source-grammar boundary, and records the excluded singleton response. The stale sentence suggesting longer descendant macros were unaddressed has also been corrected. The requested changes are satisfied; no mathematical blocker remains.

## Scope

The states begin at the rooted edge and hence describe nontrivial O/Z/J sources. The classification is about counts admitting at least one suitable representative, not every representative of that count. The corollary must say: **There exists a reachable state of count h having an O/Z/J descendant of count 2h±1 if and only if h lies in the corresponding ready set**, for odd h>1. Appending O to a ready state supplies an immediate counterexample to an unintended universal-over-representatives reading.

If interpreted as a theorem about all nontrivial physical unit-start sources, the complete-source-grammar theorem is an inherited dependency and must be cited. It is not proved just by the four-coordinate induction. The singleton is deliberately outside the state set: its response is (1,1,1,1), its J child has count1, and its Z child does not retain a unit start. Adding singleton sources would add h=1 to the minus-ready count set. This changes no h>1 statement.

## Algebra and necessity

I checked all simultaneous invariant updates. In particular R-e is even, h and d remain odd, and C=R-h+e-d is nonnegative, preserved by Z and reset by J. The auxiliary E=h-2e is positive odd. The identities E_Z=e+2(h-d-e)+(R-h) and E_J=e+(R-h) have the advertised signs.

The last-letter classification of E=1 is exhaustive. An e=0 predecessor can only be initial or end in O; the latter cannot have R=h+1. If e=1 and R=h, the d=1 reduction is valid, including the J-predecessor parity argument. In the resulting e=d=1,R=h section, O makes C positive, Z cannot undo this, and a J from h>1 makes d>1. The listed powers-of-two section and its J images are therefore complete, not just examples.

The C=d-1 and C=d+1 branches in the Z readiness argument are correct. At fixed d between J resets, the first O takes C from0 to at least d+1, later O strictly increases it, and Z preserves it. Equality d+1 requires exactly one O applied at E=1. The closed formula after O Z^t is correct for t>=1, including its first-step boundary e=0.

The J-output readiness defect is -E, so only minus-readiness can occur. All displayed families are physically attained by the stated recipes. No reversal or vertex-symmetry identification is used.

## Arbitrary descendant macros

For every reachable state, e<=(h-1)/2 and R>=h. Every non-O operation therefore sends h to at least (3h+1)/2. Two non-O operations give at least (9h+5)/4, strictly greater than 2h+1. This bound remains valid when O's occur between them because O preserves h and the invariants.

After any O, the readiness defect R-e-h is at least3. Further O only increases it. Hence a sole non-O preceded by an O produces at least2h+3. A successful macro for h>1 must therefore be one immediate Z or J followed by terminal O's. An all-O or empty macro preserves h and cannot meet either target for h>1. For h=1 the minus target is1 and this trivial macro is possible; plus-readiness is already present at the initial rooted edge. The argument remains valid at every stated boundary.

## Both signs and sparse scope

The two displayed direct-ready count sets intersect exactly at {5,9,17}. This permits different representatives for the two signs. To check the intersection, the minus family 3*2^m+3 is0 modulo3 and can intersect only the plus family5*2^t-1; its unique solution is9 by powers-of-two valuation. For the minus family2^m+1, comparing with the three plus families gives5,9,17 respectively, again by the 2-adic valuations of2^m+2. A single fixed state cannot be directly ready for both signs. Including the excluded singleton adds1 to the existential both-sign set.

Thus the no-universal-descendant-binary-step consequence is stronger than just a density observation. Neither this consequence nor the finite evidence proves that the scalar support is or is not cofinite.

## Independent computation

`readiness.cpp` modifies the separately audited positive-coordinate C++ enumeration, independently of the author's h,R,e,d implementation. It exhausts all4,082,826 level-states at scalar bound50,000 and yields26 minus-ready counts and99 plus-ready counts below25,000, exactly matching the displayed formulas. `READINESS.out` records these sets. Prior independent literal DP checks in `literal.py` verify all364 recipes through order7 and include the small boundary constructions used here.

The scalar bound captures all ready parents below25,000 because their required child has count at most49,999; the inherited pruning cannot remove such a continuation. Finite enumeration checks the formulas but does not replace their all-order proof.
