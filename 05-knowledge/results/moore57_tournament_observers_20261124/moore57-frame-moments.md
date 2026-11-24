# Multi-branch completion certificates for a Moore(57,2) search

Fix a hypothetical Moore graph and its root neighborhood arrays:
57 branches of 56 sheets with perfect matching blocks P_ij between
distinct branches. Restrict each branch to its 55-dimensional
real zero-sum sheet space and retain all matching operators.
On the resulting direct sum the second-layer adjacency B satisfies
B²+B=56I. These tests are necessary for a partial system of blocks;
passing them is not a completion result.

## Residual Gram matrix and complementary spectrum

For any set of k chosen branches decompose the centered B as

    B = [[S,C],[Cᵗ,T]].

The quadratic relation gives three exact equations:

    C Cᵗ = R_S := 56I−S−S²,
    Cᵗ C = R_T := 56I−T−T²,
    S C+C T+C = 0.

Therefore R_S is positive semidefinite with rank at most 55(57−k);
equivalently the spectrum of S lies in [−8,7] and at most
55(57−k) eigenvalues, with multiplicity, lie in the open interval.
Each interior S-eigenspace of eigenvalue λ is carried by Cᵗ
isomorphically to the T-eigenspace of eigenvalue −1−λ, since the
scale squares to 56−λ−λ²>0 by the first two equations.
As global B has multiplicities 1672 and 1463 at 7 and −8,
interlacing further forces S to have at least 55k−1463 eigenvalues
equal to 7 and 55k−1672 equal to −8 whenever these bounds are positive.
An endpoint vector x must satisfy Cᵗx=0 separately in every still
unfilled branch: its transports into any such branch sum to zero.
This last conclusion is linear propagation on the remaining
matching entries, not merely a rank test after completion.

For a purely exact search, one can use rational zero-sum bases
(e.g., differences of sheet coordinate vectors), with their known
Gram matrix, to test positivity via rational symmetric congruence
and rank via finite-field minors. The finite-field rank of a
well-defined integral matrix gives a rigorous lower bound on its
rational rank. The endpoint kernel computed exactly gives the
linear propagation rules before assigning any matching outside.

## Component sheet moments with integral propagation

Suppose the chosen matchings preserve a common transported
partition into sheet classes of sizes s_1,...,s_r in every
chosen branch, with ∑s_u=56. In particular a connected component
of the partial matching graph gives one such class: all its
inter-branch matchings are bijections between its fibers, hence
the component meets each branch in the same number of sheets.

For each external branch j and sheet a let n_(j,a,u) count how
many of its k matches into the selected branches land in part u.
For each j, its vector literally satisfies ∑_a n_(j,a,u)=ks_u
and each row satisfies ∑_u n_(j,a,u)=k.
Let w_u be the class characteristic vector minus s_u/56 in
each selected branch. Then S w_u=(k−1)w_u, and

    (Cᵗw_u)_(j,a) = n_(j,a,u)−ks_u/56.

The residual Gram identity gives the exact uncentered second
moments, summed over every external sheet (j,a):

    Σ n_(j,a,u)n_(j,a,v)
      = k(k−1)s_u s_v + δ_uv·k(56−k(k−1))s_u.

Indeed the centered source inner product is
k(δ_uv s_u−s_u s_v/56), and the residual on these source vectors
has scalar 56−k(k−1). Adding back the means gives the formula.

In particular for k=8 the residual vanishes on every class vector:
every external sheet must have exactly s_u/7 matches into part u.
Thus each s_u must be divisible by seven, independently of any
candidate filling. The scalar spectral envelope alone cannot
detect this: for k≤8 the matching block row norms automatically
give ||S||≤k−1≤7. For k>8 a nontrivial transported partition
already gives a forbidden centered eigenvalue k−1>7.

## Two locally legal partial systems rejected for different reasons

The companion exact checker lists two actual partial systems.
To read them, partition each branch's sheets by the specified
component sizes and identify each part of size s with Z/s.
In lexicographic order of i<j the listed offsets prescribe
translation P_ij by that offset modulo s in each part, and
P_ji uses its inverse. For each ordered branch pair, within
each sheet part the direct shift and all shifts via one other
chosen branch are distinct. This is exactly the complete
triangle/four-cycle exclusion within a matching-block system:
a collision of direct/two-step words forms a triangle, and
a collision of two distinct two-step words forms a square.
A triangle or square cannot project to a repeated branch at
distance two because every matching is a bijection.

* The eight-branch system uses sheet parts 27 and 29. It passes
  the full local triangle/square tests and the spectral envelope.
  It cannot extend, because an external sheet would need exactly
  27/7 matches into the first part.
* The nine-branch system uses two parts of 28 sheets. It passes
  all those local triangle/square tests but its centered block
  has eigenvalue 8. Thus it fails positivity of R_S. Nine is the
  first possible branch count for rejection by the raw spectral
  interval alone, by the general norm bound at k≤8.

## A five-step identity for a fixed three-branch start

This extra local calibration follows already from the full pair-Latin
identities; it is not an additional independent obstruction. For
three distinct branches i,j,k let H=P_ij P_jk P_ki and sum over
ordered distinct l,m outside those three. On full sheet coordinates,

    Σ P_ij P_jk P_kl P_lm P_mi = 52(J−H)+2I.

Indeed for each l, the sum over m of P_lm P_mi is
J−P_li−P_lj P_ji−P_lk P_ki. After prefix multiplication,
summing the first of the three subtracted terms over l yields
J−H−I; summing the second yields the same; the third yields 54H.
Combining with the 54 copies of J proves the formula. Since H is
a derangement, each starting sheet has 54 pentagonal returns
over this fixed three-branch start. Equivalently for any fixed
four-branch initial path there is exactly one allowed fifth
branch closing it: adjacent endpoints would make a square,
and the unique return via two matches cannot reuse an internal
branch without creating a triangle or reusing a matching edge.

## Applying the tests before assigning the next branch

Construct a k-partite compatibility graph on the chosen sheet sets.
An interpart pair is compatible iff it is not already linked by
a direct edge or a two-step matching path among the chosen branches.
A possible new external sheet has exactly one neighbor in every
chosen branch; these k neighbors must form a clique in the
compatibility graph. Conversely, to insert one entire branch
while preserving the triangle/square exclusions is exactly to
take 56 vertex-disjoint such cliques covering all chosen sheets
and make a new sheet for each clique. The disjointness ensures
the added edges are matchings and precludes a square with two new
sheets; the pair compatibility precludes a triangle or square with
one new sheet. This is a k-partite perfect hypermatching problem.

When seeking a full Moore completion, filter those cliques before
hypermatching: all endpoint vectors x above must have zero sum
across the k selected neighbors of each hypothetical new sheet.
Thus the nullspace of the partial residual produces exact
per-clique equalities, useful as early search constraints even
when the local girth test alone offers abundant choices.
