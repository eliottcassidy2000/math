# Residual scores: cycle/path correction and complete realizability

Status: PROVED (self-contained proofs below); bounded enumerations FINITE-EXACT.
No priority claim. No claim to solve the integer square-magic problem, count
unlabelled tournaments, or determine metagraph wiring. 25 October 2026.

## Inheritance and boundary

The exact inherited native object is the fixed directed Hamiltonian path
0 -> ... -> n-1 in a tournament. Session52 `sources/session52-balanced-metagraph.md`
Sections 2 and 5 prove endpoint transport and the internal-port bijection to
Eulerian orientations of K_(n+1) minus a fixed cycle when n is even. These
statements are recovered, not new discoveries. The new extraction here is an
exact residual-score feasibility criterion and its prescribed-high-set corollary.
Concept board: conservation; endpoint marks; cut feasibility; Eulerian lift;
integer lift versus modular shadow. Hostiles: K4 minus C4 is a matching;
a path deletes one fewer edge than a cycle; score data forget arc wiring.

## R1. The missing unit and parity

Let r_v count outgoing arcs remaining after deleting specified directed arcs.
Deleting a Hamilton cycle from an n-vertex tournament (n >= 3) gives

    sum r_v = n(n-3)/2.

Deleting a Hamilton path gives instead

    sum r_v = (n-1)(n-2)/2 = n(n-3)/2 + 1.

These follow by subtracting n or n-1 from binom(n,2). If n=2p >=4 and every
r_v lies in {p-2,p-1}, writing h for the number with r_v=p-1 gives
sum r_v = n(p-2)+h. Thus h=p for a cycle and h=p+1 for a path.
The latter is an obstruction to applying the former half rule, not a rounding
convention. At n=4 the path residual has three edges, so degrees 0/1 require
three high vertices, whereas cycle deletion requires two.
For odd n the floor and ceiling of (n-3)/2 coincide. A path residual cannot
have every degree equal to this value, because its total is one larger.

## R2. Exact residual-score feasibility, including all endpoint information

Let F be any fixed set of underlying removed edges in K_n, with each edge
of F given a specified orientation. Let G=K_n minus F. For integral r_v,
there is a residual orientation with precisely these outgoing degrees iff

    r_v >= 0,
    sum_v r_v = |E(G)|,
    sum_(v in S) r_v >= binom(|S|,2) - |E(F[S])|  for every S subset V. (1)

The right side is |E(G[S])|. Upper bounds r(S) <= |E(G[S])|+|delta_G(S)|
follow by applying the lower inequality to the complement. In particular
0 <= r_v <= deg_G(v). Thus (1) is both necessary and sufficient, unlike
the total sum alone. For paths F=P, the term |E(P[S])| retains adjacency
and endpoint position; sorted degree multisets cannot replace labelled r.

Proof. Necessity counts the tails of internal edges. For sufficiency, form a
bipartite graph with one left vertex for every edge e of G and r_v right
slots for every vertex v; join e to the slots of its endpoints. For any
collection A of left edges let S be the union of their endpoints. Its
neighborhood has r(S) >= |E(G[S])| >= |A| slots. The finite matching theorem
therefore gives a matching saturating all edges. Equality of total sizes
saturates all slots too. Orient each edge away from the endpoint whose slot
it occupies. This realizes r. The matching theorem itself follows by the
usual augmenting-path proof: a maximum matching with unmatched left vertex
and no augmenting path has its reachable left set larger than its reachable
right neighborhood, contradicting the inequalities. No atlas is used.

Equivalently set b_v=2r_v-deg_G(v). Then the criterion is

    sum b_v=0,  b_v congruent deg_G(v) mod 2,
    |sum_(v in S) b_v| <= |delta_G(S)| for every S,                 (2)

with r=(b+deg)/2. Singleton inequalities include nonnegativity and upper bounds.
The equivalence follows from sum_(S) deg=2|E(G[S])|+|delta_G(S)|.

## R3. Every prescribed high set is realizable (sharp small exception)

(a) Let n=2p >=6. For every chosen subset H of size p, K_n minus C_n has
an orientation with r_v=p-1 on H and p-2 elsewhere. Consequently adjoining
the fixed directed cycle gives a balanced tournament whose high set is H.
At n=4 this is true exactly when H contains one vertex of each edge of the
residual matching {0,2},{1,3}; cardinality alone is insufficient.

Proof. Here deg_G=n-3 and b_v is +1 on H, -1 elsewhere. For S of size
k<=n/2, at most 2k deleted cycle edges cross its cut, so

    |delta_G(S)| >= k(n-k-2) >= k >= |b(S)|.

Every other cut is a complement of one of these and b(V)=0. Apply R2.
For n=4 the two disjoint edges each supply exactly one outgoing unit,
which proves both directions of the stated exception.

(b) Let n=2p >=4. For every chosen subset H of size p+1, K_n minus P_n
has an orientation with r_v=p-1 on H and p-2 elsewhere. There are no
additional endpoint restrictions for residual feasibility alone.

Proof for n>=6. Let epsilon(S) count how many of the path's two endpoints
belong to S. The residual degree is n-3 at internal vertices and n-2 at
endpoints. Thus b_v=2[ v in H ]-1-[v is endpoint], b(V)=0, and
|b(S)| <= k+epsilon(S). Since at most 2k-epsilon(S) path edges cross,

    |delta_G(S)| >= k(n-k-2)+epsilon(S) >= k+epsilon(S)

for k<=n/2. Again use complementary cuts and R2. For n=4 the residual
is the tree with edges 02,03,13; orient every edge toward the unique low
vertex. This yields exactly the required score 0 at the low and 1 elsewhere.

These are existence statements for labelled score patterns, not equality
of their realization counts. Different H can have different numbers of
orientations. At n=4 the cycle exception is a genuine failure of total-only
reasoning even when all degrees are in the two nominated adjacent values.

## R3-odd. The odd-order correction is one excess unit, not a half split

For odd n>=5, write m=(n-3)/2. Every prescribed single vertex w can be the
unique vertex of residual degree m+1 in K_n minus P_n, all other vertices
having residual degree m. At n=3 this is possible exactly when w is a path
endpoint. In K_n minus C_n all degrees can equal m, including n=3.

Proof. For the path prescription, b_v=2[v=w]-[v is endpoint], so b(V)=0
and |b(S)|<=2. For n>=5 and 1<=k<=n/2, the preceding cut lower bound
k(n-k-2)+epsilon(S) is at least 2: for k=1 use n-3>=2, and for k>=2
use n-k-2>=1. Apply R2. For n=3 the residual consists of the single edge
joining the endpoints and an isolated middle vertex, proving the exception.
For the odd cycle prescription b=0 everywhere, and R2 applies immediately.
This describes floor/ceiling balancing about the TRUE residual average,
not the impossible constant assignment about the path's incorrect cycle total.

## R4. How the native internal ports sit inside the path score slice

In the fixed path gauge, full tournament degree is d_v=r_v+1 for v<n-1
and d_(n-1)=r_(n-1). Thus a residual with r in {p-2,p-1} gives a balanced
source exactly when the path's end belongs to H. Its start can be low or
high. Native complement C reverses only the residual edges, so

    d_C(v)=n-1-d_v+[v=0]-[v=n-1].

A balanced source has balanced partner iff d_0=p and d_(n-1)=p-1.
In this residual slice this is exactly {0,n-1} subset H. Conversely every
internal balanced word has this residual form. Therefore for every even
n>=4 all binom(n-2,p-1) candidate internal labelled score patterns occur,
though the number of words realizing a pattern need not be constant.
For n=2 the unique path is a separate degenerate case, with residual scores 0.

Under native C, residual degrees become deg_G-r. Interior membership in H
is interchanged, while the two endpoint residual degrees p-1 remain fixed.
Hence the high set transforms by

    H -> {0,n-1} union ( {1,...,n-2} minus H ).

This is a literal endpoint-preserving score action, not an unmarked quotient.

Add the unique regularizing vertex z, with z -> each source high and each
source low -> z. The marked path closes through z exactly for the internal
HL endpoints. The parent order is N=n+1 (odd), full degree p, and after
removing this marked cycle every residual degree is p-1=(N-3)/2.
Thus the Session52 native cycle carrier is Eulerian on K_N minus C_N.
It is NOT the even-order two-value residual cycle carrier in R3(a).
The original path residual sum and this larger cycle residual sum differ
because the new vertex also adds edges. Identifying their n symbols would
erase a necessary type distinction.

## Loss ledger and the square-magic boundary

Source: an oriented graph with fixed path/cycle edges and vertex labels.
Map: record residual outdegrees. Preserved: sum, degree parities, all cut
inequalities, endpoint transport on the indicated domain. Lost: individual
arc orientations, cycles, automorphisms and native metagraph wiring. Sidecar:
retain the labelled orientation (or actual word), not merely a score histogram.
Cheapest hostile: the n=4 cycle matching and a path's extra edge unit.

The sum constraint is an integer conservation law, not a transfer theorem to
3x3 magic squares of distinct squares. No map here sends an orientation to
nine integer squares while preserving the eight magic-line equations and
distinctness. The incidence/matching feasibility proof does not enforce
quadratic representability or no-wrap carries. Any claimed exact arithmetic
bridge needs those additional data and a proved map; none is asserted here.

## Reproduction

Run `python3 check.py > normal.json` and `python3 -O check.py > optimized.json`;
the JSON outputs must agree. Exhaustive universes and independent tests are
reported in README.md. No unlabelled census, solver, external dependency,
network source, or literature-priority assertion is involved.
