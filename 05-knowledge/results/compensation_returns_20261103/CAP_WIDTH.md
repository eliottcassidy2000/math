# Constant endpoint mass can require unbounded feedback width

**PROVED:** the zero-channel cap of [structure/PROOF.md](structure/PROOF.md)
has feedback vertex number exactly `ceil(j/2)` at depth j. The separate
endpoint law is `H(D_j)=3^j`, `end_(x0)(D_j)=1`. **FINITE-EXACT:** independent
literal checks through j=6. These are statements about actual graphs, not
a numerical resemblance between sequences.

## Inheritance and exact map

The inherited [fixed-feedback automaton](../feedback_word_hamilton_spectrum_20261020/proofs.md#5-exact-finite-state-grammar-at-arbitrary-fixed-feedback-vertex-number)
is exact at every fixed feedback vertex bound. The later
[active-port construction, Section 6](../adaptive_metagraph_reformulation_20261103/ACTIVE_PORTS.md#6-temporal-boundedness-is-different-from-static-feedback-boundedness)
already distinguishes bounded temporal memory from unbounded static feedback
cost, using disjoint triangles. This note applies that mechanism to the
specific endpoint-preserving cap, with an exact minimum rather than merely
an unbounded lower bound. No novelty in the general distinction is claimed.

Start at physical anchor x0. At step i add x_i beating only x_(i-1), then
y_i beating only x_i; all other incidences point from old to new. The next
auxiliary is x_i, while the endpoint anchor stays x0. Relist the same names
in the order

    x0, y1, x1, y2, x2, ..., y_j, x_j.

All arcs point forward in this order except exactly

    x_i -> x_(i-1),  1 <= i <= j.

This is a relisting, not a change of graph or identification of names.
It retains every joint response. Forgetting the labels or identifying all
x_i with one abstract auxiliary would destroy the future attachment data.

## Exact minimum

For each odd i, the triple `{x_(i-1), y_i, x_i}` is a directed triangle.
These `ceil(j/2)` triangles are vertex-disjoint, so every set whose deletion
is transitive must meet each one. Conversely delete the odd-indexed x_i.
Every backward arc then loses an endpoint; the remaining displayed order
is transitive. Consequently

    fvs(D_j) = ceil(j/2).

The empty depth has fvs=0. The argument allows arbitrary feedback deletions,
not merely subsets of the x_i. The odd triangles prove that deleting y_i
instead cannot reduce the minimum.

For the full reader wrapper T_j with marked w and deletion `U -> D_j`,
where U is transitive in the two current q0=1 sources,

    ceil(j/2) <= fvs(T_j) <= ceil(j/2)+1.

The lower bound is inherited from the induced cap; the upper bound deletes
w and the cap feedback set. Equality for either side is not asserted in
general. At the first actual children1409/1427, fvs is exactly2: the cap
gives a deletion set of size2, while the all-order q1 exclusion forbids
any one-vertex feedback set.

## Consequence for the next construction

The first children fit the inherited two-apex word carrier: delete physical
0 and201, leaving the transitive order
`99,10,11,2,1,3,4,12,200,401,400`.
That observation cannot be carried to arbitrary zero-channel depth with a
fixed two-apex bound. The cap itself has fvs>2 once j>=5, even though one
endpoint count remains1 at every depth.

Thus the endpoint coordinate is an exact useful invariant, but it is not a
bound on graph complexity. A return grammar based on it must keep variable
depth (or a faithful recursive representation). The outstanding arithmetic
obligation is still to pay the preparation charge in equation(3) of the
structure note while preserving the actual named root and joint guards.
Neither small endpoint mass nor an exact finite automaton supplies that
payment or proves global Hamilton-spectrum coverage.

## Reproduction

    python 05-knowledge/results/compensation_returns_20261103/check_cap_width.py
    python -O 05-knowledge/results/compensation_returns_20261103/check_cap_width.py

`cap_width.json` retains all seven actual matrices, names, exact Hamilton
and endpoint counts, and minimum feedback sizes. The checker constructs
each matrix both by literal births and by the displayed ordered rule,
checks the minimum by every vertex subset, and counts paths independently
by endpoint subset DP. Its predicates remain active under optimization.
