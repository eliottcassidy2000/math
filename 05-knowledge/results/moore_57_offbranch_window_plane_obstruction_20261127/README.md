# Off-branch window 160–166 and a projective-plane obstruction

2026-11-27. Continue the [rooted degree-57 tournament encoding](../moore_57_slotwise_balance_20261127/README.md).
Every statement about a degree-57 matching table is conditional on
the table obeying exact Moore transport. The Paley-59 benchmark
uses the 57-vertex quotient obtained by deleting two vertices.
We narrow its off-branch optimum to four possible values. The
most direct attempt to attain 162 by making matched branch
interiors exact converses is impossible on the largest Paley
defect clique: it would construct a projective plane of order 14.

## 1. An upper bound valid for every matching table

Use a semiregular tournament for each 56-slot interior, meaning
that signed row sums are \(\pm1\). For example, delete one vertex
from the regular 57-vertex carousel. For a fixed pair of branches,
put \(g=(S^2)_{ij}\), align its direct matching to the identity,
and multiply the block by \(\mathrm{sgn}(g)\). Choose the direction
\(i\to j\) for which \(S_{ij}g>0\). Then the entries are
\[
 W_{ab}=52|g|+4q_{ab}+x_a-y_b-2(C_i[a,b]+C_j[a,b]),
\]
where \(x,y\) are the signed interior row sums, \(q_{aa}=0\) and
\(q_{ab}=\pm1\) off the direct matching. Every row and column of \(q\)
sums to \(|g|\). The exact Moore route partition is what makes
the route term bounded by four.

For a two-vertex deletion of any doubly regular 59-vertex
tournament, the prior note gives \(|g|\le3\), so the same
semiregular choice works for every matching system:
\(|W_{ab}|\le156+4+2+4=166\). The
[five-step defect certificate](../moore_57_slotwise_balance_20261127/README.md)
gives a lower bound of 160 for these quotients, even allowing
arbitrary interiors. Entries are even. Thus the optimum belongs
to **160, 162, 164, 166** for each such fixed matching table and
quotient. The upper bound constructs interiors, not matchings.

## 2. Exact tests for the remaining thresholds

For a candidate cap \(B\in\{160,162,164,166\}\) first fix branch
score vectors in \(\{-3,-1,1,3\}\). For the Paley benchmark
these are required at cap 160 by the prior score bound.
Pairs \(|g|=1\) are automatically safe with these scores.
Pairs \(|g|=3\) have positive scaled entries (at least 142);
only their upper bounds matter. At matched slots test
\(156+x_a-y_a\le B\). At each remaining slot pair, the
exact constraint on the two interior arc signs is
\[
 C_i[a,b]+C_j[a,b]\ge (156+4q_{ab}+x_a-y_b-B)/2.
\]
In original frames one sign uses \(P_{ij}^{-1}b\) in the
first branch and the other uses \(P_{ij}a\) in the second.
The sum is in \(\{-2,0,2\}\); the constraint is either
vacuous, a two-literal OR, two unit clauses, or inconsistent.
Add the exact degree condition in each row that defines
its proposed score. This is a complete feasibility formulation
at fixed scores and matching table, with no remaining path
enumeration. Dropping the degree conditions gives a sound
2-SAT rejection filter, not a sufficient completion test.

For semiregular scores, the \(q\le0\) cases are all harmless,
and the positive route table becomes especially small:

| Cap | \(x_a-y_b=-2\) | \(0\) | \(2\) |
|---|---|---|---|
| 160 | At least one positive | At least one positive | Both positive |
| 162 | Free | At least one positive | At least one positive |
| 164 | Free | Free | At least one positive |
| 166 | Free | Free | Free |

Two local results make the matching dependence explicit.
With equal aligned semiregular interiors, a row of score
\(1\) has 29 positive routes but just 28 outgoing arcs.
At a positive route where its arc is negative, both interior
arcs are negative and \(x_a-y_b\ge0\), forcing at least 164.
With opposite aligned semiregular interiors the arc terms
cancel and the maximum is **exactly 162**: for any positive
score row of the first branch, there are 29 positive routes
among 55 slots but only 28 negative-score slots in that
branch. Some positive route therefore joins two positive
scores in the first frame, contributing \(160+2\).
All other entries are bounded by 162.

## 3. Why the exact-converse proposal fails on the Paley clique

The preceding exact 162 block suggests requiring that every
defect matching transport one interior to the converse of the
other. The Paley-59 quotient has four defect cliques of sizes
14, 14, 14, 15, the sign patterns to the deleted vertices:
within one pattern \(g_{ij}=1+1+1=3\). In fact **no clique of
15 branches of a Moore-perfect table admits opposite
interiors across all its matchings**, even without
semiregularity. Here is the obstruction.

On a fixed clique, use one reference branch; its transports
along closed walks form a permutation group. If every edge
transports an interior to its converse, the sign of a walk is
\((-1)^{\rm length}\). This is well-defined on its permutation:
the identity cannot send a tournament to its own converse
while fixing all labels. Thus there is a sign homomorphism,
surjective by the closed walks of length three. Its positive
kernel is a subgroup of the automorphism group of the reference
interior, so has odd order: an involution would interchange
some two vertices and reverse their arc. The negative coset
flips the sign of every signed row sum; all these sums are
nonzero because the order is 56.
On any slot orbit the positive kernel is transitive on each
score-sign half, and the halves have equal odd size. Indeed
a positive element is necessary to take a slot to one of
the same sign in its orbit, and a negative element gives a
bijection between the halves.

Transport this orbit partition to each clique branch. It is
well-defined since any two paths differ by a closed walk at
the reference; every clique matching respects the transported
orbits and flips intrinsic score signs. For distinct branches
\(i,j\) and any slot in an orbit, the 13 two-step routes
through the other clique branches land at distinct slots
in the same orbit and of the same score sign at \(j\).
Hence each orbit half has odd size at least 13. They sum
to 28, so the orbit sizes are exactly 26 and 30.

Now restrict to the 26-slot orbit. Call its positive slots
points and its negative slots lines in each of the 15
branches. Adjacent point and line mean matched by the Moore
table. Each line meets one point in each of the other
14 branches, none in its own. Any two points from different
branches share one line: the 13 clique routes cover exactly
all 13 same-sign slots in the target. Points in the same
branch share no line, by the matchings. The dual statements
hold as well.

Adjoin a point \(p_i\) for each branch and add it to all old
lines of that branch. Adjoin one more point \(\infty\).
For each \(i\) add a line through its 13 old points,
\(p_i\), and \(\infty\). Finally add the line through all
15 points \(p_i\). The augmented structure has 211 points
and 211 lines, each line containing 15 points. Every pair
of points lies on exactly one line: old pairs use the stated
routes or the new line of their branch, an old point and
\(p_i\) use its unique matched line in branch \(i\) or the
new line if it is their own branch, and the other pairs
use the displayed new lines. All lines meet exactly once
as well, by counting the 210 other lines incident with the
15 points on any given line. This is a projective plane of
order 14. [Bruck–Ryser](https://doi.org/10.4153/CJM-1949-009-2)
rules it out: 14 is \(2\bmod4\) and not a sum of two squares.

So a strategy limited to **equal or opposite** aligned
semiregular interiors throughout the Paley defect clique
must have an equal pair somewhere, forcing at least 164.
This does **not** exclude 160 or 162 for general interiors.
What it excludes is the uniform converse construction that
would have made 162 automatic without inspecting individual
route signs.

Run `node replay.js` for the threshold table, signed-orbit
controls, an exact synthetic routed-block test, and the plane
completion checked on existing projective planes of orders
three and five. The degree-57 table is not manufactured by
these tests.
