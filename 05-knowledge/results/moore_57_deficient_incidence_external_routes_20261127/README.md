# Moore transport rules out a 14-branch converse design before completion

2026-11-27. Follow-up to the [off-branch window and 15-branch plane obstruction](../moore_57_offbranch_window_plane_obstruction_20261127/README.md).
The 14-branch missing-partner design has a clean completion criterion,
but a stronger obstruction comes from counting routes through the other
43 branches. In fact **at most eight branches can be pairwise
score-sign reversing**, with no need to assume exact converse arcs.
A [spectral follow-up](../moore_57_mixed_alignment_spectral_20261127/README.md)
now rules out the entire equal-or-converse template, even with mixed
signs, so the conditional template bound below is superseded.

## Rooted transport and the deficient design

A putative degree-57 diameter-two Moore graph, rooted at any vertex,
has 57 branches of 56 second-shell slots apiece. Each pair of
branches has a perfect matching between its slots. For two slots in
different branches, exactly one of their direct match or their
two-step routes through a third branch occurs. No two distinct
slots in the same branch have a common neighbor among the slots.
These are the only global Moore properties used below.

Suppose every pair in a set of 14 branches transports its interior
tournament to the converse. The signed-holonomy argument in the
preceding note applies with 12 distinct two-step routes in place
of 13. Briefly, closed transports are signed by the parity of
their walk length, with a negative element from any triangle;
an identity of negative sign would identify an interior with
its own converse on fixed labels. The positive kernel consists of
tournament automorphisms and has odd order. On each orbit the
positive and negative row-score halves have equal odd size
(the order of each interior is even, so there is no zero score).
Every matching in the chosen branches preserves the transported
orbit partition and interchanges the signs. The 12 routes force
each half to have at least 12 slots, hence at least 13. The
half-sizes sum to 28, so they must be 13 and 15.

In the 26-slot orbit, take positive slots as points and negative
slots as lines, divided into 14 groups of 13 of each. The incidence
matrix N has zero diagonal blocks and permutation matrices
off-diagonal. Two points in different groups have one common
line unless they are the unique missing partners in each
other's groups. Two points of the same group share no line.
For the 13-regular graph A of missing point partners,
with D block diagonal with fourteen all-ones 13-by-13 blocks,
the exact identity is

    NN^T = 13 I + J - D - A.

Likewise N^T N uses the dual missing-line graph B. In particular
A and B are both covers of K14 by perfect matchings between groups.

## Completion and the information in the Gram matrix

On the space of group-constant vectors, NN^T acts by 169 on
the all-ones direction and by 1 on its 13-dimensional complement
within that space. On the space W of vectors summing to zero in
each group it is 13I - A. Each connected component of A is
regular and meets every group, so the multiplicity of 13 as an
A-eigenvalue is the number c of components. Exactly one of these
eigenvectors is in the group-constant space. Consequently

    nullity(N) = c - 1,      components(A) = components(B).

The second equality uses the square matrix N. Positive
semidefiniteness alone is automatic for any abstract matching cover:
on W it just says the eigenvalues of a 13-regular graph are at most
13. For example start with identity matchings except for one
13-cycle matching; the resulting A is connected. This is a Gram
PSD example, not a claim that it has a binary incidence factor.
The exact determinant condition is

    det(N)^2 = 169 det((13I - A)|_W).

In particular a connected leave must pass a nonzero square
determinant test; for disconnected leaves the determinant already
vanishes. We do not use this test to infer completion.

The design extends by a fifteenth group of 13 points and 13 lines
if and only if A is thirteen copies of K14. Necessity: each new
line meets one point in each old group and its pairs must be
missing, so the 13 new lines give disjoint 14-cliques of A.
For sufficiency, if A is such a union then c=13 and the nullity
identity gives 13 components of B. Each has at least 14 vertices,
so B too is thirteen copies of K14. Add a new line for each
A-clique, incident with precisely its 14 points, and a new point
for each B-clique, incident with precisely its 14 lines. Use no
point-line incidences within the new group. An old line contains
13 points with distinct A-cliques (two on the same line cannot
be missing partners), hence one in each; dually each old point
meets one line of every B-clique. This checks all mixed pairs,
and the old pairs and within-group pairs have the required
incidence by construction. The resulting 15-group design augments
to a projective plane as in the preceding note: add a corner
point to each group's old lines, a line through the group's
points and its corner and a common infinity point, and the line
of all corners. Each pair of points is on one line. [Bruck–Ryser](https://doi.org/10.4153/CJM-1949-009-2)
rules out order 14. Thus even a bare local order-14 deficient
design, if it existed, could not have the cliquewise leave.
We need not decide bare local existence.

## Outside branches overfill the leave

Here is a stronger, elementary obstruction that uses all Moore
transport. For a chosen set of k branches, let there be h marked
slots in each such that there are no marked-to-marked direct
matches, and such that the k-2 internal two-step routes from any
marked slot to another chosen branch always land on marked slots.
This holds for a fixed positive orbit-half as above, or for all
positive slots of a clique of score-sign reversing matchings.

For every pair of chosen branches exactly h(h-k+2) of their
marked pairs still need a route. Each outside branch has 56 slots;
write t_v for the number of marked neighbors of its slot v.
The sum of t_v in that branch is kh, because every chosen
branch contributes its h matched slots. Every pair of marked
neighbors of v is from different branches and must be an
uncovered pair. Conversely each uncovered pair needs exactly
one such v. Thus the following is an exact equality with a
necessary integer-convexity lower bound:

    choose(k,2) h(h-k+2)
      = sum_outside v choose(t_v,2)
      >= (57-k) F(kh),
    F(e) = 56 choose(q,2) + r q,   e=56q+r,   0<=r<56.

The minimum per branch comes by balancing loads to q or q+1.
For k=14, h=13 there are just 1,183 pairs in the leave. Each
outside branch forces at least 210 pairs: its load is 182, so
42 slots with load 3 and 14 with load 4 minimize the pair count.
Over 43 branches this is 9,030. The 15-half independently fails
(4,095 pairs available versus 12,642 forced). In particular no
component structure of the 13-half leave can save the embedding.

There is a stronger cutoff without using the orbit decomposition.
Take nine branches whose matchings all flip row-score signs.
In a triangle, those bijections force 28 positive and 28
negative slots in each branch, since the positive counts of
each two branches sum to 56. Mark every positive slot. Now
k=9, h=28, and the route budget is 21,168 versus 21,504
forced outside. Therefore **no nine branches can be pairwise
score-sign reversing**. This also rules out any larger such
set by restriction. The score signs can even be replaced
by arbitrary binary labels that every chosen matching flips.
For eight exact-converse branches, applying the same test to
each odd orbit-half removes the possible halves 9 and 11; the
remaining partitions of 28 are 7+21, 13+15 and 7+7+7+7.
Passing the count is not a construction.

## Quantitative consequence for the off-branch template

For k chosen branches with 28 positive slots apiece, do not
require any sign flips. Let e be the total number of direct
positive-positive edges between them. There are choose(k,2) 28²
cross-group positive pairs. Exactly e are direct; the rest
are accounted for by counting pairs of positive neighbors at
each slot, both inside and outside. The outside contribution
is at least (57-k)F(28k) regardless of e. The negative slots
inside have total positive-neighbor load 28k(k-1)-2e: each
positive slot has k-1 direct contacts inside, and positive-positive
contacts are counted twice. At a negative slot of load u, use

    choose(u,2) >= choose(k-1,2) - (k-2)(k-1-u).

The positive slots' own route contributions are nonnegative.
Including the e direct pairs gives the necessary inequality

    (2k-5)e >= (57-k) F(28k) - choose(k,2) 28(30-k).

For k=14 the minimum forced by the inequality is e>=427;
for k=15 it is e>=541. In the equal-or-converse semiregular
template, each equal-aligned pair contributes 28 to e, and
each converse pair zero. Therefore each 14-clique needs at
least 16 equal-aligned pairs and the 15-clique needs at least 20.
Across the four Paley defect cliques that is at least 68
equal-aligned blocks. In the oriented defect direction of any
such block, each of the 28 positive-score rows has 29 positive
routes and only 28 outgoing arcs in its aligned interior.
Some positive route therefore hits a negative arc in both
interiors, and its score difference is nonnegative, giving
a scaled entry at least 164 by the original block formula.
This forces at least 1,904 entries at or above 164 in the
chosen directed blocks, conditional on using that template
throughout the four cliques.

The maximum lower bound in that restricted template stays
164 (already known from the 15-clique), but its violations
can no longer be isolated to one block or even one defect
clique. **The general 160/162/164/166 window is unchanged.**
The matching-table constraints are conditional on Moore
transport and restrict only the simultaneous signings of
its rooted encoding; they do not disprove existence of
a degree-57 Moore graph.

Run node replay.js for the arithmetic, the eight-branch
orbit filter, and deficient point and line controls obtained
by deleting one group from PG(2,q), q=3,5,7. The connected
abstract leave control is only a test of the graph/PSD
constraints, not a Moore transport example.
