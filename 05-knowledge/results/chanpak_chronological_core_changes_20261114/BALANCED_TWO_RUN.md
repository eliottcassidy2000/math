# Exact chronological-path cost of a balanced one-apex tournament

Let B_r be the strong tournament on one apex v and two transitive blocks X,Y of order r, with v→X→Y→v. The orders within X,Y are from first to last, with arcs forward. Let κ(T) denote the smallest number of Hamiltonian paths whose union of arcs can be deleted to leave a DAG, the prior session's [chronological path budget](../chanpak_collatz_exact_addresses_20261114/READOUT_63.md). Equivalently there is a common rank order and every rank-backward arc occurs in at least one of those paths.

**Theorem.** κ(B_1)=1 and for all r≥2,
 
    κ(B_r)=r−floor(sqrt(r)).
 
Moreover H(B_r)=(2^r−1)^2+2. Thus a tournament can be one vertex from transitive and nevertheless need roughly half its order in full chronological paths to explain all deviations from one rank. This does not assert a large minimum κ over *other* tournaments with the same H.

## Lower bound: the two extreme edges of a double contact

Fix a common forward rank order and k witnessing paths. Each of the X vertices before v forces a backward arc v→x, and each of the Y vertices after v forces a backward arc y→v. There are at most k vertices of each type, since each path can spend just one outgoing and one incoming port at v. Put A=X after v and D=Y before v; then |A|,|D|≥r−k. Every arc A→D is backward and must be read.

Every full path has at most two X→Y transitions, one on each side of its unique visit to v. If it has two, it consists of an increasing X run, an increasing Y run, v, a second increasing X run and a second increasing Y run. The tails of the two transitions include the last vertex of X; the heads include the first vertex of Y. Indeed these are the extrema of the two increasing runs which together exhaust the corresponding block.

Assign each required A→D arc to one path that reads it. Every path receiving two assignments uses one distinct assigned arc incident to last(X) and one distinct assigned arc incident to first(Y). If either extreme is absent in the rectangle, there are no double assignments. Otherwise there are at most min(|A|,|D|) double assignments. So |A||D|≤k+min(|A|,|D|). For k<r and l=r−k≥1 the right inequality implies k≥l(l−1), hence r≥l². This proves the lower bound; if k≥r it is automatic.

## Upper bound: pack a rectangle into full paths

Take l=floor(sqrt(r)), k=r−l. Rank the vertices in this order:

    first k of X, first l of Y, v, last l of X, last k of Y.

Apart from the k incident backward arcs in each direction at v, the only backward arcs are the rectangle from the last l of X to the first l of Y. If l=1, partition this rectangle into one singleton group. If l≥2, pair the l−1 noncorner entries in its last row with the l−1 noncorner entries in its first column, one from each. Pair the shared corner with one still unused interior entry. All pairs have distinct tails and heads and collectively touch last(X) and first(Y) *within each pair*. The other cells stay single. There are l²−l groups for l≥2, which is at most k. The case l=1 also fits in k for r≥2.

Assign at most one group to each of k paths. Assign the k distinct X-prefix vertices as their required outgoing ports u after v, and the k distinct Y-suffix vertices as their required incoming ports z before v, in any bijective pairing. Now extend each group:

* For no assigned rectangle edge, take X_{<u}, Y_{≤z}, v, X_{≥u}, Y_{>z}, omitting empty runs.
* For one assigned edge x→y, use it before v. The pre-v X run can contain exactly X_{<u} and x, and the pre-v Y run can contain exactly y,z. The remaining X run starts at u and the remaining Y run, if nonempty, follows it.
* For two assigned edges x→y and x'→y', use the first before and the second after v. Partition X into increasing runs with the first ending at x and the second starting at u and ending at x'. This is possible: u precedes both tails and the intervals of allowed entries (all vertices up through x, or the interval from u through x') cover X since one tail is its last vertex. Assign the two tails and u to their required runs before allocating the others. Partition Y into runs with the first starting at y and ending at z, and the second starting at y'. The allowed intervals [y,z] and [y',last(Y)] cover Y because one head is its first vertex and z follows both. Assign the endpoints first. The two concatenated X,Y runs, separated by v, meet both required edges and both ports.

Thus the paths cover all backward arcs at the stated k; deleting them leaves the common order. For r=1, the triangle is not acyclic but any of its Hamiltonian paths already suffices.

## Exact count and interpretation

The cyclic triangle through v has one nonempty increasing run of X and one of Y. For any pair of endpoint choices x in X and y in Y, its odd-cycle packing contribution is 2 when they are consecutive in the spine and 2^(j−i−1) otherwise, where i,j are their spine positions: choose even cardinality among the intermediate positions. Adding over x,y gives (2^r−1)^2 from the geometric weights and one extra from the consecutive pair at the seam. Add the empty packing. This recovers the inherited block-word count (see [one-apex cut law](../session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md) and [Irving–Omar](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20)).

The forcing mechanism is an endpoint resource not visible to feedback vertex number: a path can spend two rectangle contacts only by spending global beginning/end resources of the two monotone blocks. The near-square H is a natural fixed-shape arithmetic output, but its representation may be expensive in chronological paths. This is one reason to avoid replacing the scarce absolute-target arithmetic by a claim that all simple cores can be read with two paths.

## Rectangular generalization and its exact count inverse

For positive a,b let B(a,b)=C3[1,TT_a,TT_b] with the same block order as above. Define a rectangle cost

    g(A,D)=0                         if A=0 or D=0,
           A·D                       if min(A,D)=1,
           A·D−min(A,D)              if A,D≥2.

Then

    κ(B(a,b))=min { k≥0 : g(max(a−k,0),max(b−k,0)) ≤ k }.

For the lower bound, coverage of the rank-backward rectangle needs at least g(|A|,|D|) paths. When one dimension is 1, a directed path cannot use two of its edges (common head or tail); otherwise use the double-contact bound above. The function g is nondecreasing in each coordinate, including the change from 1 to 2. Thus the displayed test is necessary.

For sufficiency first take k<min(a,b), choose the first k of X and last k of Y as hub ports, and use the same rank interleaving. If either dimension of the remaining rectangle is 1, its cells must be singleton groups. Otherwise pair distinct noncorner last-row and first-column cells until the smaller of the two lists is exhausted. Pair the center with any interior cell. This makes min(A,D) valid disjoint two-cell groups. The same extension of a group plus its designated hub ports into a full path works at unequal block sizes. For k≥min(a,b), it is enough by converse to take b=min(a,b), rank v,X,Y, and cover the only backward arcs y→v by b individual paths. The path for y is Y_{≤y},v,X,Y_{>y}. This also shows the minimization has a solution.

The same geometric sum gives H(B(a,b))=(2^a−1)(2^b−1)+2. The unordered pair {a,b} is uniquely read back from H. Specifically H=3 comes from a=b=1. For any other candidate set t=v₂(H−3). For unequal a<b, t=a, and solve 2^b=(H−2)/(2^t−1)+1 with the power and inequality checks. In the equal case t=a+1, and directly test a=t−1≥2. These two branches are disjoint: the minimum unequal count with a=t exceeds the equal count with a=t−1. This is a prescribed-target inverse for this physical family.

For each fixed k≥1, the number of distinct values H≤X from this rectangle family with budget at most k is

    k log₂ X + O_k(1).

Indeed if a,b>k and κ≤k, the test implies (a−k)(b−k)≤2k, giving finitely many such pairs. For each a=1,...,k, the choices b≥a with H≤X number log₂ X+O_k(1). The exact count inverse prevents collisions between these branches. Subtracting adjacent bounds shows even each exact positive budget contributes log₂ X+O_k(1) distinct counts. Thus allowing the shape to change within this two-run construction gives a real hierarchy but no fixed-budget cofiniteness proof.

In particular, the at-most-two-path values in this family are exactly

    {2^b+1, 3·2^b−1 : b≥1} ∪ {51,107,227}.

To see this, min(a,b)≤2 gives the two branches. Otherwise the cost inequality at k=2 requires both residual sides in {1,2}, giving the three unordered pairs from {3,4}. Moreover no value from this family, at any budget, is divisible by 7: each Mersenne factor is 0, 1 or 3 modulo 7, and the products omit 5=−2 mod 7. This obstruction is arithmetic, distinct from the path-budget inequality.
