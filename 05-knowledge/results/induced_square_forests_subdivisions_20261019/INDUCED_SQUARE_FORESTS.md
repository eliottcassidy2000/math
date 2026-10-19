# Induced square-sum forests at near-quadratic height

Status: **PROVED**, independently [reviewed](review/INDEPENDENT_REVIEW.md), 2026-10-19. No priority claim.

Let S_N have vertices 1,...,N, distinct a,b adjacent exactly when a+b is a square. This definition differs from literature called square-sum labeling where edge weights are sums of squares.

## Inheritance and the missing predicate

The repository's [elementary square-sum embeddings](../square_sum_elementary_sparse_path_cycle_bounds_20261019.md) embeds every tree with O(n^2) labels, but explicitly permits extra edges. Here the retained predicate is the entire adjacency relation, including nonedges. The additional tool is bounded pairwise common-neighborhood size, supplied by factorization of a difference of squares.

## A general greedy lemma

In a finite simple graph G of minimum degree delta and maximum codegree lambda, every n-vertex tree has an induced embedding if delta > (n-2)lambda+(n-1). Indeed, grow a rooted tree in parent-before-child order. When k vertices are embedded, candidates for the next child are neighbors of its parent. Discard the k occupied labels and the common neighbors of the parent with each of the k-1 earlier nonparents. At most k+(k-1)lambda candidates disappear. This preserves every edge and nonedge.

This is an elementary sufficient condition, not an equivalence. It is a version of a standard greedy sparse-codegree argument; no claim that this abstract lemma is new.

## Exact arithmetic codegree

Write D(M)=max_{1<=d<=M} tau(d). In the induced graph on integer labels [M,2M], two distinct labels b,v have at most tau(|v-b|) common neighbors. For a common neighbor x, write x+b=r^2 and x+v=s^2 with positive roots. Then (s-r)(s+r)=v-b. A divisor s-r determines s+r and hence r,s,x uniquely. Signs and parity only reduce the count. For negative difference reverse b and v.

Each vertex b has at least (2-sqrt(3))sqrt(M)-2 neighbors inside [M,2M]. The interval of admissible roots is [sqrt(M+b),sqrt(2M+b)] and has length at least (2-sqrt(3))sqrt(M), attained at b=2M. Counting integers loses at most one; excluding a possible loop loses at most one more. In particular the degree is at least sqrt(M)/4-2.

Therefore, if

    sqrt(M)/4-2 > (n-2)D(M)+(n-1),

then every n-vertex tree occurs induced in S_(2M), entirely inside [M,2M]. Every n-vertex forest also occurs: connect its components arbitrarily into a tree on the same n vertices, embed that tree, then this does NOT give an induced forest. Instead use a new auxiliary root joined to one root of each component, embed the resulting (n+1)-vertex tree, and delete the auxiliary root. This yields the same asymptotic bound and the displayed condition with n replaced by n+1.

## Near-quadratic upper bound, optimal exponent

For every epsilon>0, tau(d)<=C_epsilon d^epsilon for all positive d. Elementary proof: for primes p>=2^(1/epsilon), a+1<=2^a<=p^(epsilon a); for the finitely many smaller primes take the finite supremum of (a+1)/p^(epsilon a), and multiply these constants. Thus D(M)=M^{o(1)}.

It follows that, for every eta>0, there is C_eta such that every n-vertex forest has an induced square-sum embedding with maximum label <=C_eta n^(2+eta). Equivalently the optimal worst-case label height is n^(2+o(1)): the matching exponent lower bound comes from the n-vertex star. In S_N its center has at most floor(sqrt(2N)) neighbors, so N>=(n-1)^2/2. This lower bound does not imply a sharp constant or remove the subpower loss.

For the upper bound choose epsilon<eta/[2(2+eta)] and M a sufficiently large constant times n^(2+eta), so sqrt(M) dominates n M^epsilon. No unproved prime-distribution estimate is used.

## An explicit elementary quartic bound

For all d>=1, tau(d)<=9d^(1/4). For p>=17 the prime-power ratio (a+1)/p^(a/4)<=1. For p=2,3,5,7,11,13 the maximum occurs at a=5,3,2,1,1,1 respectively. These locations follow by comparing consecutive ratios ((a+2)/(a+1))^4 with p. The fourth power of the product of maxima is

    127401984/25025 < 9^4.

Set M=(64(n+1))^4. The displayed degree-codegree criterion for (n+1)-vertex trees follows since sqrt(M)/4-2=1024(n+1)^2-2, whereas (n-1)D(M)+n <=576(n-1)(n+1)+n. Hence every n-vertex forest embeds induced with maximum label at most 2(64(n+1))^4. This constant is deliberately loose.

## Infinite extension and exact seeds

Every countable forest embeds induced in the positive square-sum graph. Root each component. Start with any enumeration of the countable vertex set; at stage j append every still-unlisted vertex on the finite root-to-jth-vertex path, in ancestor order. This gives a parent-before-child enumeration even when a level is infinite, adding a new component root when needed. At a child step, choose x=r^2-b. For each earlier nonparent v, the difference-of-squares factorization forbids only finitely many r, and collisions forbid finitely many more. At a new-root step, exclude the union of finitely many sets {r^2-v}; each has density zero, so infinitely many positive choices remain. Every existing finite induced seed is preserved exactly if future vertices attach to at most one existing vertex at their insertion. This is an extension statement for forests over the seed, not arbitrary prescribed two-neighbor attachments.

## Boundary and next target

This does not prove graceful labeling or an induced embedding using every consecutive integer. Ordinary interval square-sum Hamiltonian existence already has independent constructions. The induced theorem answers the extra-edge issue in the earlier sparse construction and isolates the next obstruction: simultaneous positive adjacency to two prior vertices needs actual common neighbors, whereas the tree argument only uses upper bounds on them. For example labels 1 and 2 have no common square-sum neighbor because two positive squares cannot differ by 1. A two-parent analogue cannot follow from a generic greedy step.

Open quantitative target: can the n^(o(1)) loss be removed, proving an absolute C such that every n-vertex tree embeds induced in S_(Cn^2)? The star shows that the exponent 2 is unavoidable. The maximum-divisor estimate used here cannot by itself supply this constant-scale theorem.

## Exact-seed extension (additional proof candidate)

Let an arbitrary finite induced seed have distinct positive labels at most H. Add vertices successively with exactly one specified earlier neighbor. The same proof works with new labels in [M,2M] once M>=H: earlier parents may lie below M, which only increases the root-interval length; label differences are now <=2M, so replace D(M) by D(2M). Occupied labels and earlier nonparents are counted as before. For every eta>0, a seed plus additions totaling n vertices therefore has an exact-seed induced realization of maximum label O_eta(max(H,n^(2+eta))). The hidden constant is uniform in the seed. A new isolated root can instead be handled by an auxiliary vertex and subsequent deletion, or by the density-zero argument in the infinite case.

In particular the exact seed 482,3362,6242 (all three pair sums are squares) admits arbitrary attached rooted forests while retaining precisely its triangle and attachment edges. The conclusion does not erase the triangle or imply interval coverage.

## Induced cycle closure with two vertices (additional proof candidate)

The earlier sparse cycle formula also becomes induced after choosing its parameter sufficiently large. Let A,...,B be an induced square-sum path of at least two vertices, with largest label H. Choose an integer r>2H+4 with r congruent to A+B+1 modulo 2. Put

    v=r^2-B, s=(A+B+1-r^2)/2, w=s^2-A.

Then B+v=r^2, v+w=(s-1)^2, w+A=s^2. Both new labels exceed H and are distinct (w>v for this large choice). For an old vertex y different from B, v+y=r^2+(y-B), where 0<|y-B|<H<2r-1. This lies strictly between the nearest squares below and above r², so is nonsquare. Likewise for y different from A, w+y=s²+(y-A); since |s|>H, its nonzero deviation is smaller than 2|s|-1, and is nonsquare. Thus no unwanted new edges occur, and the old induced path closes to an induced cycle with exactly two additional vertices.

Every finite cycle therefore occurs induced: use this construction for lengths at least four, and the explicit seed triangle for length three. Combining with exact-seed leaf attachments realizes every finite connected unicyclic graph induced, with polynomial label bounds. This argument alone does not give near-quadratic height for unicyclic graphs: r=O(H) gives w=O(H^4), and the induced-path bound supplies only a larger polynomial bound. Optimizing this induced closing step is a separate target.

## Universal induced two-subdivisions

The two-vertex closure proof never used the old path, only two distinct old labels A,B. It therefore inserts an induced length-three path A--w--v--B between any two seed vertices, with no other new adjacencies. Any existing A--B edge remains unchanged. Starting with independent original labels 4i+1 and doing this for every edge of an arbitrary graph G realizes its exact two-subdivision as an induced square-sum graph. In particular this supplies literal induced subdivisions of K5, K3,3, and the Petersen graph; it does not transport chromatic number, unit-distance geometry, or tournament orientations.

The labels need not grow by successive fourth powers. A common parameter window gives a polynomial bound. Let G have n vertices and m edges, and K=n+2m. Original labels have maximum H=4n+1. Choose integer R>=4sqrt(H+1). For every edge with endpoints A,B choose r in [R,2R] of parity A+B+1, and use the same v,s,w formula. All constructed labels lie in [1,5R^4]. Each earlier label y forbids at most D(5R^4) parameters through the unwanted v+y square condition (except y=B), and at most D(5R^4) through the unwanted w+y square condition (except y=A). For the second assertion, the factorization determines |s|, and |s| is strictly increasing with r throughout this window. Each of the two possible label collisions forbids at most one parameter per earlier label. Finally w>v>H throughout the window, so the new pair itself has no collision with original labels or each other.

Consequently fewer than 2K(D(5R^4)+1) parameters are forbidden at any step. There are at least R/2-1 candidates of the designated parity. Choose R so

    R/2-1 > 2K(D(5R^4)+1).

Since D(5R^4)=R^{o(1)}, this is possible with R=K^{1+o(1)}. More precisely, for every eta>0 a uniform constant gives maximum label O_eta(K^(4+eta)). Thus exact two-subdivisions of ALL finite simple graphs have induced square-sum realizations at polynomial numerical height. The only arithmetic estimate is the elementary divisor bound proved above. This leaves open a near-quadratic height bound for these subdivisions.

For countable G, enumerate vertices and edges fairly; do not apply finite avoidance to an already infinite seed. Proceed by inserting each new original vertex after avoiding square sums with the finitely many previous labels, then inserting each edge bridge only once its endpoints exist. The finite-avoidance proof then applies. Every countable graph has its two-subdivision as an induced subgraph of the infinite square-sum graph.

Taking A=B in the same formula gives an attached triangle, with r odd. More generally, attach a fresh path at a seed vertex and then use the bridge to return to that vertex. It follows that every finite or countable cactus graph is induced-realizable by finite-stage attachment. A polynomial height bound uniform in the number of cactus cycles does not follow from the simple sequential construction and is not asserted here.
