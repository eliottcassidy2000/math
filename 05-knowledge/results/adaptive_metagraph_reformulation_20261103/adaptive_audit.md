# Independent rooted-profile audit

**PASS for ADAPTIVE.md §§2–4 and the stated proof boundaries.** Exact R,N
updates, the reset and defect criterion are all-order proofs. The displayed
same-profile obstruction is independently verified. Its minimal order-six
claim relies on the producer's complete ordinary-class augmentation through
six; this audit does not independently duplicate that census. Section 5's
full A/B/G transfer is outside this receipt's independent computational scope.

For the new-root map N, Hamilton paths starting at the newborn are exactly
old paths not starting at the chosen old root, so s'=h−s. Deleting the newborn
gives q'=h. A two-slot cover ending at it either has singleton first slot
(contribution h), or appends it after the old root in an old ending-root cover
(contribution F). Thus F'=h+F. Together with the inherited birth count,
N(h,s,q,F)=(h+F−s,h−s,h,h+F). Substitution proves
R∘N=(h+2F,h,h,h+2F).

The factor trap is correctly scoped: after resetting an actual balanced state
with count b>1 and b≠7, subsequent keep/new counts equal b times an actual C3
keep/new count. That latter factor is at least three and is never seven,
since its first counts are 3,5,9 and then strictly grow. The product is neither
a prime nor 7p. The b=1 exception and b≠7 guard matter. This is no obstruction
to later choosing another older root.

For δ=F−h−q+s, splitting Hamilton words after the root accounts for h covers.
Among the remaining broken covers, singleton first slots account for q−s.
Removing the root from every other first slot gives exactly the two nonempty
paths ending and starting at two different in-neighbors. Hence the stated
δ=sum W_ab identity is bijective. The cited Session124 pair-deletion bound
b_ab≥H(U−{a,b})≥1 has the same nonempty-slot normalization and yields
δ≥d(d−1), so δ=0 iff d≤1. No stronger degree formula is asserted.

In a strong source, if the root has unique predecessor u, then q−s counts
Hamilton paths of its deletion beginning at u. Every old vertex reaches u
in that deletion, so u lies in its terminal SCC. Such a path exists iff the
deletion is strong. The longest-fixed-start insertion proof in §3 is valid:
an outside vertex that cannot be inserted anywhere after the fixed start or
at the end must dominate the entire maximal path, contradicting strongness
if any outside vertices remain. Thus the profile's legal-inverse criterion
uses both δ=0 and q>s correctly.

## Independent literal universe

`adaptive_independent_audit.py` imports no producer or repository modules.
It enumerates literal vertex words. A word with no broken adjacency contributes
to H and to all F_v by splitting after v. A word with exactly one broken
adjacency contributes only to the F_v at that break; two or more breaks
contribute to none. Deletion counts are separate permutation counts.

All 285 rooted labelled sources at orders one through four satisfy R,N,RN,
the defect laws, the sharp indegree-zero detector, and δ≥d(d−1).
The two six-vertex named sources in ADAPTIVE.md independently have profile
(37,5,13,71). Their first marked children and **all fourteen actual eight-vertex
second children** were counted by full vertex permutations. The next-root
lists are exactly

    A: 287,321,219,329,275,233,179
    B: 289,307,339,219,271,227,179.

Thus equal profiles cannot authorize equal arbitrary-root continuation, while
the keep-root and newborn values agree as required. This is an exact object
witness, not numerical coincidence or a proposed unrealized profile.

Reproduction:

    python3 adaptive_independent_audit.py > adaptive_independent_audit.json

Normal and optimized Python runs produce identical JSON.
