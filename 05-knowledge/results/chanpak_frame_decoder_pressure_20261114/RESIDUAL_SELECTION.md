# Canonical choice needs more than the Hamilton count

## The minimum-order parent rule already fails at 63

An attractive rule is to choose, for every smaller attained count, a
tournament of minimum order. Could every new count arise by one vertex
insertion into a chosen smaller representative, perhaps with a clever
tie-break? No: 63 defeats **every** tie-break.

Use the ascending spine word `100110`, followed by the two named hubs
$x\to y$. Its ordinary count is 63. The eight single-vertex deletion
counts (spine first, then the hubs) are
$$
  (11,37,37,33,33,25,25,1).
$$
The positive counts 11,25,33,37 already have witnesses of order at most six:
the two-hub words `110`, `1121`, `1120`, `1200`, respectively.
The unit count has a singleton witness. Thus all eight order-seven parents
are non-minimal for their own count.

The inherited [order-eight census THM-344](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-344-n8-h63-complete-omega.md)
gives exactly two order-eight classes, of opposite scores and converses;
the score sequence of this word is $(1,2,2,3,3,5,6,6)$.
There are none of smaller order. The
[Busch strong order bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
in its $H(T)\ge5^{(n-1)/3}$ form (see the
[Session 324 application](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session324_busch_inverse_and_gap_structure_20261113/GAP_STRUCTURE.md))
rules out strong cores of order at least nine. Any witness of 63 has exactly
one nonsingleton strongly connected component, since every proper product
of 63 into nonunit factors requires 7 or 21, which are forbidden.
All other vertices are singleton padding. Inserting a padding vertex
does not increase the count; inserting a core vertex when other
padding is present leaves the parent padded and therefore nonminimal.
Converse preserves each deletion count. Hence 63 cannot arise by one
strict insertion from a minimum-order representative in any way.

A shortest-parent **two**-vertex lift does work: delete the two consecutive
zeros at positions 1 and 2 (zero-based) to recover `1110` of count 23
at order six. The replay exhausts all labelled tournaments through five
and verifies that 23 has no smaller witness. Inflating that seam gives
$$
 H(1\,0^k\,110)=11+13\,2^k \quad(k\ge1),\qquad H(1110)=23.
$$
The exact block tensor is $(37,13,0)$ in the basis
$1,2^k-2,3^k-3\,2^k+3$; it certifies all lengths. The real but
nonminimal intermediate $k=1$ has count 37. The minimum-order jump
for this target is exactly two, not an absence of a recursive construction.

More generally deletion of a vertex preserves $H$ exactly when it lies
in no directed triangle: one triangle gives a positive new term in the
odd-cycle packing sum; no triangle means a singleton strong component,
whose presence multiplies the count by one. Thus any minimum-order
representative of $H>1$ strictly descends under **every** deletion, but
the parent need not be minimal at the new scalar. This distinction is
the missing hypothesis in a scalar induction.

## A correct canonical section for appends: 17 suffix probes

Let $f_w(z)=H(wz)$ on normalized two-hub words. Two words lie in the
same right residual fibre if their $f_w$ agree on **all** suffixes.
Replacing a prefix by the first shortest-lexicographic word in this
fibre is safe for arbitrary append continuation; replacing it by a
shortest word with the same $H$ only is not. For instance
$$
 H(10011)=H(1121)=25,\qquad H(100110)=63,\quad H(11210)=79.
$$
The shorter word for 25 is not a valid replacement for this continuation.

In fact a minimal linear right profile has **17** coordinates. Take the
probes
$$
 \varepsilon, 0, 1, 2, 00, 10, 20, 30, 01, 11, 21, 31, 000,100, 200, 010, 210 .
$$
The inherited 36-state exact transfer gives reachable row dimension 28
and observable column dimension 20. The rank of their pairing (the
Hankel matrix $H(wz)$) is 17. The probes above achieve this rank
against the 17 explicit prefixes in the verifier. Therefore they
characterize the full right residual, and 17 is minimal for a
rational linear realization on the state span. This does not assert that
fewer probes could never separate the discrete set of actual word states. This is a finite closure certificate, not a finite-word
heuristic: if each appended letter sends every basis row into the span,
all reachable rows lie in that span; dually for prepended letters
and suffix-testing columns. Exact rational elimination constructs and
checks these closures.

Each letter induces a rational linear map of the 17 profile coordinates.
All four maps have rank 13 and eigenspace dimensions $4,8,4,1$
at eigenvalues $0,1,2,3$ respectively, certified by exact elimination.
The three positive eigenvalues are the same occupancy channels
$1^k,2^k,3^k$ as in the homogeneous block formula. More generally that
formula applied to any $H(u d^k v)$, $k\ge1$, implies the polynomial
$t(t-1)(t-2)(t-3)$ annihilates the operator of a letter on the minimal
linear quotient. The zero channel is the information lost when the
first letter is appended, not an obstruction to finding $k=0$ by
treating it as a different skeleton.

As a consequence, for the shortest and then lexicographically least
two-hub word of a realised target $N$, **every prefix is itself the
shortest-lexicographic word in its 17-coordinate residual fibre**.
Otherwise replace it by a shorter or lexicographically earlier equivalent
prefix and keep the same suffix, contradicting the choice. This is a
genuine compatible section for ordinary count and suffix appends.
The 17 coordinates do not by themselves preserve named anchors, grades,
strongness, reorderings or arbitrary insertions; those require keeping
their respective extra guards.

The remaining existence obligation is clear: prove that the cone of
physically reachable profiles meets each hyperplane whose empty-suffix
coordinate is a desired count (or prove enough of those intersections
and use additional constructions for the rest). The linear profile
and right-fibre choice certify safe replacement, not that any arbitrary
integer vector is physically reachable. See `python residual_machine.py`
for the exact realization and [FRAME_DECODER.md](FRAME_DECODER.md) for absolute target
selection once a marked skeleton is known.
