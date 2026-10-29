# Recovered Cayley aggregate: exact components, and why paths do not lift

Session102, 2026-10-29. **INHERITED PROVED, rederived and FINITE-EXACT reverified.**
Source pin `3b6e3ef82382da222fc59aa02914eea72d5a105b`.

## Provenance first

The complete result requested in this lane already appears, with the additional
bipartiteness and injectivity conclusions, as **Theorem3, §8** of
`05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`.
The fixed-path switching transversal is older still: **THM-474**, file
`01-canon/theorems/THM-474-tilings-are-switching-classes.md`, and the spanning-tree
extension in **THM-1467**, file
`01-canon/theorems/THM-1467-oriented-spanning-tree-switching-sum.md`.
Only their actual tree-gauge statements are used; unrelated historical invariant
claims in those files are not dependencies. Session90's `trijection/TREE_COMPRESSION.md`
in `05-knowledge/results/session90_even_graphs_odd_descent_20261029/` led to this
recovery. This note must not be presented as a newly discovered connectivity
classification. The fresh contributions are a compact audit, independent code,
and an explicit representative-level lift obstruction.

This recovery also corrects the novelty boundary in `SIGN_GRAMMAR.md`: its
n>=8 two-step switching obstruction follows from inherited generator
injectivity. Its full n7 legal skeleton/colour classification and general
word-length obstruction are subsequent refinements.

## 1. The exact carrier

Fix labelled vertices [n], n>=3, and one reference total order. Tournament
coordinates x belong to X=F2^(E(K_n)), recording disagreement with that order.
Let B be the cut space, of dimension n−1. The labelled switching-class space
Q=X/B is an affine space with chosen coordinate origin, of dimension

    beta = binom(n−1,2).

This does not quotient by vertex permutations. Relabelling acts affinely on
state coordinates; it acts linearly on displacement vectors. In particular
this is not an unlabelled tournament count or a natural identification with
unlabelled even graphs.

For each ordered spanning path P, exactly one tournament in any class q has
all arcs of P forward: fix the first vertex's switching bit, then solve the
successive switching bits along the path. The global complement of all bits
is the same cut, so the tournament representative is unique.

The native move on that representative toggles J\E(P). Its endpoint class is

    q -> q + g_P,    g_P=[J+E(P)] in Q.                (1)

Consequently the **existential support quotient** is the Cayley graph on Q
with generators g_P, where P can be any undirected Hamilton path. More
precisely, aggregating *all* original labelled native incidences gives

    C[q,r] = 2 #{undirected P : r=q+g_P}.              (2)

Every row sums to n!, but this is not an equitable quotient: different
tournaments inside q can have different Hamilton-path counts. The factor2
in (2) comes from the two directed orientations, each realized uniquely.
A native source/target pair determines its agreement path and its direction,
so these incidences are not duplicate witnesses of the same pair.
No g_P is zero, since J\E(P) cannot be a cut (a consecutive triple has one
toggled edge, whereas cuts have zero or two).

## 2. Generator-span proof

Compare vertex words (a,b,c,rest) and (b,a,c,rest). Their path-edge vectors
differ by e_ac+e_bc. Thus differences of path generators include every wedge
(two distinct edges sharing a vertex). Such wedges span the entire even-edge
subspace X_even: the line graph of K_n is connected, and telescoping along a
line-graph walk expresses any pair of edges as a sum of wedges. Conversely,
every pair of Hamilton-path vectors has even total weight. Therefore

    span{g_P+g_R} = image(X_even in Q).               (3)

Before quotient, span{J+E(P)} is X if beta is odd, and X_even if beta is even:
all generators have weight beta, and their differences already span X_even.

For even n, a singleton cut has odd weight n−1, so image(X_even)=Q.
For odd n, every cut has even weight s(n−s), so edge parity descends to Q
and image(X_even) has codimension one. When n is odd, generator parity beta
is odd exactly at n=3 mod4. It follows that the Cayley support has

* **one component** if n is even or n=3 mod4;
* **exactly two components** if n=1 mod4, distinguished by orientation parity.

The origin naming the two parity sheets depends on the reference order;
the partition itself does not. An odd relabelling can exchange the sheets.
Components are cosets of the generator span. This proves the claimed
classification for every n>=3, not just the tested orders.

The stronger inherited bipartiteness conclusion follows from (3): a Cayley
component is bipartite iff a linear functional equals1 on every generator.
It must vanish on all differences. Thus the connected n=3 mod4 case is
bipartite, using edge parity; for even n and for either n=1 mod4 component,
the differences span the entire component, ruling out such a functional.

## 3. Generator collisions and the new refinement's proper place

Equality g_P=g_R means E(P) symmetric-difference E(R) is a cut. Its degrees
are at most4. For n>=9 a nontrivial cut has a vertex degree at least5.
At n=8 only K4,4 could survive, but a path endpoint has degree at most3.
Hence for n>=8, different undirected paths have different generators;
there are n!/2 generators and aggregate weight2. This is the inherited
proof, not new work of Session102.

At n7 the inherited census gives1260 generators from2520 undirected paths.
`SIGN_GRAMMAR.md` supplies the explicit extremal skeleton proof and tracks
legal representatives, showing why its two native moves must have the same
colour. That is information absent from the colour-forgetting Cayley support.

## 4. Concrete failure to lift a quotient edge

On vertices0,1,2,3 let T be the transitive tournament in increasing order.
Its only directed Hamilton path is P=(0,1,2,3), so it has only one native
neighbour. But the quotient support at n4 has degree6. In particular take
Q=(0,1,3,2). Then g_Q differs from g_P, although Q becomes directed in a
different switching representative of T. Thus the quotient edge labelled Q
cannot be lifted from this specified T, even by choosing a different native
path: the only path P reaches a different quotient vertex.

The checker encodes Q's gauge representative by backward mask42 in pair order
(01,02,03,12,13,23), switching vertex2. Fundamental-triangle coordinates give
source class0, target g_P=1, and target g_Q=2. After an actual native step,
the next desired quotient edge may likewise require a switch to another
representative. Such switches are additional operations, not native moves.
Therefore neither Cayley connectivity nor its component count proves native,
black-only, unlabelled, or prescribed-score connectivity.

### Native disconnection hidden by the quotient

At n3 the entire labelled native graph is **two disjoint copies of K1,3**.
Each of the two directed3-cycles has three Hamilton paths; its three native
partners are transitive. Each transitive tournament has a unique Hamilton
path and hence a unique neighbour. Thus the six transitive tournaments form
the leaves, with no other edges. The switching quotient has just two
vertices and one support edge, so is connected. This gives a literal
connectivity failure, not merely a path-label failure.

The fresh finite census is:

|n|labelled vertices|native edges|native component sizes|quotient components|
|---|---:|---:|---|---:|
|3|8|6|4,4|1|
|4|64|96|64|1|
|5|1024|3840|512,512|2|
|6|32768|368640|16384,16384|1|

The order6 census was computed by directed-path backtracking and union-find,
not an all-permutations scan. Its two components are exactly the even and odd
backward-edge parity sets, independently measured from every component. Each
native move toggles beta=10 pairs, whereas an odd-size switch toggles an odd
number. Thus the connected quotient erases a genuine native obstruction.
`SIGN_GRAMMAR.md` equation(5b) gives the all-order necessary parity guard; the
assertion that exactly two components suffice here is **FINITE-EXACT at n6**,
not an all-order conclusion. No full native census at n7 was attempted.

These are exhaustive labelled counts, not new unlabelled classification
claims. No all-order native connectivity assertion is extrapolated from them.

For n>=8 there is an all-order local obstruction to a covering-map analogy.
Take the transitive T in increasing order and the vertex word
P=(0,2,1,3,...,n−1). Neither P nor reverse(P) is directed in T. But the
quotient still has its P-edge from [T]. Generator injectivity says that any
native edge projecting to that edge would have the same undirected path P,
and thus one of those two orientations. Neither is legal. Hence that quotient
edge has **no native lift from T**. This is stronger than saying a preselected
path label fails; it excludes every possible lifting witness.

## 5. Fresh reproducibility

`check_cayley_quotient.py` uses no inherited implementation. It encodes the
quotient by independent fundamental triangles (0,i,j), computes binary ranks,
and compares them to the proof. Every undirected Hamilton path is enumerated
for n3..8. At n9 and10 it uses an explicit complete ordered-triple wedge
family, **not** an exhaustive path census. That generating family is enough
to certify the full lower rank; the parity proof supplies the upper rank.

| n | dim Q | generator rank | difference rank | components |
|---|---:|---:|---:|---:|
|3|1|1|0|1|
|4|3|3|3|1|
|5|6|5|5|2|
|6|10|10|10|1|
|7|15|15|14|1|
|8|21|21|21|1|
|9|28|27|27|2|
|10|36|36|36|1|

Every labelled tournament and every directed vertex word are tested for
normalization at n3,4; direct native-incidence aggregation verifies (2).
The order4 quotient has48 nonzero directed entries and192 underlying
ordered native incidences. Normal and `python3 -O` produce identical JSON;
all checks remain active under optimization. Reproduce with

    python3 check_cayley_quotient.py > cayley_checks.json
    python3 -O check_cayley_quotient.py > cayley_checks_optimized.json
    cmp cayley_checks.json cayley_checks_optimized.json

The code provides750,053 explicit checks, including exhaustive native connected-component
enumeration at orders3..6 using directed-path backtracking and union-find.
Every enumerated native incidence is also checked against its beta edge-parity
change. These certify the stated finite
universes; the all-order result is the inherited elementary proof above.
