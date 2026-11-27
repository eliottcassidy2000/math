# A 63-cycle as a spanning scaffold, not a cyclic blow-up

2026-11-27. This follows the [7 → 21 → 63 detour construction](../moore_57_heptagon_63_detours_20261127/README.md), the [closure calculations](../moore_heptagon_petersen_deficit_20261127/README.md), the earlier [tiling and curvature experiments](../../../04-computation/general_pq_tessellations_s18p.py), the [four-algebra comparison](../s92_prime_powers_quaternions_bases_20261111/README.md), and the [Cullen–Woodall idempotents](../cullen_woodall_split_complex_trace_ladder_20261127/README.md). All degree-57 assertions are conditional; [existence remains open](https://www.mdpi.com/2227-7390/14/23/4318).

## 1. An exact sense in which a 63-cycle reaches all 3250 vertices

Let G be a degree-57 Moore graph. Mark any ordinary 63-cycle (chords allowed). Starting with its vertices, repeatedly add the common neighbor of any two marked nonadjacent vertices. Recall the general two-neighbor process in [Ibrahim–LaFayette–McCall](https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p060.pdf). Here is a short spectral proof of spanning, also for any starting set of at least 63 vertices.

Write n=3250. The Moore identity A²+A−56I=J gives eigenvalues 7 and −8 on the space orthogonal to constants. A closed proper set S of size s has at most one edge from each outside vertex into S, hence cut size ≤n−s. The cut lower bound from the eigenvalue 7 is 50s(n−s)/n, by applying 57I−A to the centered indicator of S. Thus s≤65. The induced graph on S has diameter at most two and girth at least five, since it contains every common neighbor required for its nonedges.

For clarity, any such induced graph on more than one vertex is either a star or regular. If B is its adjacency matrix and d_i its internal degrees, then B²+B=J+diag(d_i−1). Commuting with B gives (1−B_ij)(d_i−d_j)=0. Distinct degree classes are therefore fully joined; three classes would give a triangle, and if there are two then one is a singleton, since otherwise there is a square. The singleton is universal so the others form an independent set. A star here has ≤58 vertices. A regular closed graph of degree k has k²+1 vertices by the two-step tree count. Among 63,64,65 only 65 has this form, for k=8. Its nonprincipal eigenvalues would be the conjugate roots of x²+x−7, with equal multiplicities 32 by integrality, contradicting trace zero (8−32≠0). This eliminates every proper closed set of size ≥63.

Thus a **pointwise marked 63-cycle is a graph-labeling seed**, though it is not a standalone description of the graph. A fixed scan order on pairs of labeled vertices gives unique labels as new common neighbors are exposed. Any graph automorphism fixing the cycle pointwise is trivial; the stabilizer of its chosen cycle edges acts faithfully on the dihedral 63-gon, even when the cycle has chords.

## 2. What the other 3187 vertices cannot be

Give each of the 63 cyclic positions a bag, anchored by its marked cycle vertex. A balanced assignment has 37 bags of size 52 and 26 of size 51. The whole graph cannot consist of these bags connected solely to their cyclic neighbors, since even the support on the bags would then have diameter ≥31.

**Quantitative result.** In *every* balanced cyclic assignment, at least **66,578 of the 92,625 edges** join nonconsecutive, distinct bags. To see this, retain only the local edges (inside a bag or between consecutive bags), and let l_v be the remaining degree of v. For v in bag i, all endpoints reached locally in one or two steps lie in bags i−2,…,i+2. By no triangles or squares, neighbors and nonbacktracking two-step endpoints are distinct and none equals v. Sum over v. The first-step total is Σl_v and the nonbacktracking two-step total is Σl_v(l_v−1), so

    Σ l_v² ≤ Σ_i s_i(s_(i−2)+s_(i−1)+s_i+s_(i+1)+s_(i+2)−1)
            ≤ 5 Σ_i s_i² − n = 835120.

The last bound uses Cauchy–Schwarz on each cyclic shift. For integer l_v, (l_v−16)(l_v−17)≥0, so 33Σl_v−272n≤835120. The number of local edges is at most floor((835120+272n)/66)=26047. Thus more than 71% of edges must be supplied by the nonlocal arrangement, not the ring.

**No balanced equitable model at all.** In an equitable partition each vertex of bag i has the same count q_ij into bag j. Along a pair of sizes 51 and 52, 51q_ij=52q_ji and coprimality force either no edges or the complete bipartite graph, which contains squares. Thus an equitable balanced partition would disconnect G.

**Unbalanced packages must vary internally too.** For any equitable partition into k nonempty bags, the quotient is similar to a symmetric compression and its eigenvalues are 57 once and otherwise 7 or −8. Its trace, the sum of internal degrees q_ii, satisfies

    Σ_i q_ii ≡ 50+7k (mod 15).

For k a multiple of three the right side is 2 mod 3. In particular neither 21 nor 63 positions admit a quotient with the same internal degree at each position, even if nonlocal edges are permitted and bag sizes are unrestricted. Seven positions instead require a common internal degree congruent to 12 mod 15. The same calculation rules out a 63-fold circulant equitable quotient without any divisibility argument about sizes.

## 3. Positive, zero, and negative curvature: what 5,6,7 measures

In a regular map with p-sided faces and q corners per vertex, the Euler characteristic per vertex is 1−q/2+q/p. At q=3 this is (6−p)/(2p): the pentagon, hexagon and heptagon give the positive, flat and negative cases, as in [Coxeter's tessellation framework](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/D0219653B540B3A893C4E3D2FDA5CD84/S000843950000165Xa.pdf/the-trigonometry-of-hyperbolic-tessellations.pdf). The cutoff is a statement about face and vertex incidences, not an intrinsic curvature assigned to a bare cycle.

One can realize the 63-cycle as the boundary of a single face in an orientable cellular embedding of G: at each cycle vertex, put the previous neighbor immediately before the next neighbor in the local rotation. The resulting face walks are nonbacktracking and each has length at least 5. With E=92625, the chosen 63-face gives F≤1+floor((2E−63)/5)=37038, hence genus ≥26170 in *every* such embedding. Without designating it as a face the corresponding girth bound is 26164. These are lower bounds, not a genus determination. The picture is a disk rim in a cellular surface of very high genus, not a spherical or flat tiling. This does not identify a metric embedding of the graph.

## 4. Three algebras in the actual adjacency operator

On the real space orthogonal to constants put K=(2A+I)/15. The Moore identity gives K²=I and K is symmetric. As in the Cullen–Woodall split-complex construction, (I±K)/2 are orthogonal idempotents; here they select the 7 and −8 eigenspaces rather than the two signs of an integer.

For the endpoints x,y of a marked edge let u=(e_x−e_y)/√2. Then c=〈u,Ku〉=−1/15. Put s=√(1−c²)=√224/15 and w=(Ku−cu)/s. The plane with basis (u,w) is K-invariant and has matrix [[c,s],[s,−c]]. If R is reflection in the u-axis then KR acts as the rotation [[c,−s],[s,c]]. Hence I_c=(KR−cI)/s squares to −I. The two involutions give split-complex directions separately; their product on this plane gives a genuine complex structure. For a nonadjacent marked pair, the same calculation gives c=+1/15, so the sign of the cosine reads adjacency.

There are always two edges along the 63-cycle whose four endpoints induce 2K₂. Indeed the internal degrees of its 63 vertices have mean ≤√62, by counting disjoint paths of length one and two as in §2. Some cycle edge has sum of internal endpoint degrees at most 15, since the mean sum over cycle edges is at most 2√62<16. Deleting the union of the two internal neighborhoods removes ≤15 vertices, including the two endpoints; the remaining ≥48 vertices must contain a consecutive pair. The two resulting planes are orthogonal: for their difference vectors u,v, both 〈u,v〉 and 〈u,Kv〉 vanish, and K preserves inner products.

Identify the planes by their ordered bases. On their direct sum let I_H act by I_c on both and put J_H(p,q)=(−Rq,Rp). Since R anticommutes with I_c, I_H²=J_H²=−I and I_HJ_H=−J_HI_H, giving an explicit quaternion action on this four-dimensional space. It preserves the Euclidean form for unit quaternions. The four-dimensional space and its structures are built using actual cycle edges and A, but the rotations act on vectors, not necessarily on graph vertices.

There is also a clean flat limit in the separate nonbacktracking transfer picture. For a degree-d regular graph, the walk transfer on an adjacency mode λ has determinant-one normalized form [[λ/√(d−1),−1],[1,0]]. This follows by excluding immediate reversals (the recurrence after the first two steps is P_(m+1)=λP_m−(d−1)P_(m−1)). For a nonscalar transfer T the element J_T=T−Tr(T)I/2 satisfies J_T²=(Tr(T)²/4−1)I. Negative, zero, and positive give complex, dual, and split-complex algebras, respectively. In the degree-57 candidate the two nonprincipal modes are complex (49 and 64 are below 224), while the constant mode is split (57²>224). The constant mode on the bare cycle, d=2, is exactly dual (2²=4). This trace threshold and the 5/6/7 Euler threshold are different invariants; the analogy is the sign trichotomy, not an equality between them.

Run **node replay.js** for degree-2, Petersen, and Hoffman–Singleton controls, using the old explicit 42-cycle in the last. It checks the Moore identity, exact integer versions of the reflection-plane Gram relations, rotation-system face counts, and the conditional numerical bounds. The constructive target left here is a non-equitable placement of the 3187 extra vertices whose nonlocal edges satisfy the unique two-step route for every nonedge.
