# Admitting richer cores: factorial ports, negative borrow, and the multi-apex slot budget

Three complementary tests of what a new physical core can buy. The first gives a sharp finer reader for genuinely new absolute profiles; the second effects a literal signed carry; the third explains why the old alphabet capacity argument does not exclude an expanded architecture. None implies Chan–Pak coverage.

## 1. Homogeneous readers have an exact factorial tariff

Let a host vertex v be replaced by an arbitrary nonempty tournament U, every vertex of U inheriting the same outside contacts as v. Write L_j(U) for the covers by exactly j nonempty **ordered** paths, and Q_j for j ordered path slots with empty slots allowed. Thus Q_1=H. Compress each maximal run in U in a full Hamilton path to a single occurrence of v; its symbolic word is a skeleton, and all exterior vertices appear exactly once.

If c_j counts valid skeletons of j occurrences, the readout is
```
H(host[v←U]) = Σ_j c_j L_j(U),     (j−1)! divides c_j.
```
Indeed every one of the j−1 internal exterior fragments (between two successive occurrences of v) is nonempty, individually starts at an out-neighbor of v and ends at an in-neighbor of v. Permuting these fragments fixes all contacts and the two exterior end fragments. The distinct exterior labels make that action free. The coefficient of Q_j in the same identity is a multiple of (j−1)!, by L_i=Σ_j (−1)^(i−j) binomial(i,j)Q_j and divisibility for every i≥j.

**Top visit coefficient attained.** Give v exactly k out-neighbor portals q_i and k in-neighbors p_i. Orient q_i→p_i, and p_j→q_i when j≠i, with arbitrary orientations within the q and p sets. An internal exterior fragment needs at least two vertices: an out-neighbor at the start and an in-neighbor at the end. At the maximal k+1 visits, every fragment must therefore be one matched pair q_i→p_i, with no exterior end fragments, and there are exactly k! skeletons. This is a strong host after substitution: every q_i reaches U via p_i, while every p_i is reached through q_i.

The previous [slot spectral result](../chanpak_crossing_scales_20261115/SLOT_SPECTRAL.md) says that equal-order cores with equal Q_1,…,Q_k have divided slot polynomial difference vanishing at 1²,…,k². In detail, for order n=2D+r with r∈{1,2}, the odd-cycle packing expansion gives
```
Q_j = j^n + 2 j^r G(j²),   G in Z[z].
ΔQ_j = 2(2k+1)! j^(r−1) binomial(j+k,2k+1) B(j²),   j≥k+1,
```
where the integer polynomial B is the quotient of ΔG by the monic product ∏_(l=1)^k(z−l²). In the last equality use j·∏(j²−l²)=(j+k)!/(j−k−1)!. Combining this with the skeleton tariff proves that **any homogeneous vertex reader** has absolute difference divisible by
```
2(2k+1)!·k!,       odd core order;
2(2k+1)!·(k+1)!,   even core order.
```
For odd order each active j≥k+1 contributes (j−1)!; for even order it contributes j·(j−1)!, giving the stronger factorial. This holds even when the reader can revisit U more than k+1 times. It does not apply to mixed outside contacts on U or to changing the lower slot observations.

**Physical admission and sharp instance.** For the previous two order-eight rows (92,141,32,228,206,19,38,101) and (52,229,208,231,10,148,49,81), the shared profile (H,Q_2)=(345,5640) and ΔQ_3=720 now give exactly 1440 paths of separation, using the matched two-pair portal host. This is the minimum possible nonzero separation between these two cores in any homogeneous host, improving the old 2880 cyclic reader. Also this is not merely a new geometry of an old profile: **(345,5640) is not in the old ordinal one-apex profile alphabet**. To check that finite claim at all allowed orders, the core bound Q_2≥2·3^(n−3) restricts each strong one-apex SCC order to at most 10 when Q_2≤5640. Enumerate its strong opening/closing words using the [exact gap polynomial](../chanpak_heisenberg_target_20261115/GAP_HEISENBERG.md), and close their (H,Q_2) under multiplication, including singleton factor (1,2). This exhausts every ordinal sum at Q_2≤5640: each nontrivial factor increases Q_2, so no excluded branch can return. The replay implements the finite bank and verifies the absence, and independently counts both literal new portal reads.

## 2. A real negative-low borrow with a richer core

For two nonempty cores U,W and singleton q with quotient U→q→W→U, the cyclic run count, or the [single-outlet identity](../chanpak_chronological_core_changes_20261114/ONE_OUTLET.md), is
```
H = (Q_2(U)−h_U)(Q_2(W)−h_W)+2h_Uh_W.
```
To see it through the actual runs: a path whose q seam is interior permits a nonempty ordered two-piece cover of each side when expanded with the endpoint cases; classifying the length 3,4,5 cyclic run words gives coefficients 3,1,1,1 on the products with one run on each side, just U split, just W split, both split. Replacing the nonempty two-cover by Q_2−2H gives the formula above.

Append TT_s ordinally to U. The two-slot count scales by 2^s and H stays fixed, so the signed expression is
```
H = Q_2(U)(Q_2(W)−h_W)·2^s + h_U(3h_W−Q_2(W)).
```
Choose the strong order-six core U with rows (24,53,57,50,32,1) and (H,Q_2)=(17,224), and W=TT_2, to obtain 672·2^s−17. U contains two vertex-disjoint triangles, so it cannot have a single-vertex feedback set. Another branch takes singleton U and W with rows (4,1,2,7), profile (3,24), giving 42·2^s−15. After shifting the widths these yield
```
21·2^j−17,    21·2^j−15,      for every j≥5.
```
At the nominal high coefficient 21, each low part is genuinely negative: both perform one borrow, giving actual quotients 20 and remainders 2^j−17, 2^j−15. The target-first check is exact: test whether N+17 or N+15 is 21 times an allowed power of two. The cores in this example need not supply new (H,Q_2) profiles; the richer geometry of the first is a weaker statement than the profile exclusion in §1. Neither branch has arbitrary placement of its borrow.

## 3. Exact slot cost of simultaneously active hubs

Fix k≥1 apex vertices, with a transitive spine of n vertices; orientations among hubs and from hubs to the spine can otherwise be arbitrary. For j≥1 fixed, every such tournament obeys
```
Q_j ≤ [2j(k+1)]^k (j+k)^n.
```
For the proof use the [odd-cycle packing expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20). An odd-cycle packing has weight 2 per cycle; choose one of j slots for every cycle and every unused vertex, obtaining the Q_j weight. Every cycle visits a hub, so first encode which hubs are used and each used hub's next hub (possibly itself). There are at most (k+1)^k successor codes. Then each spine vertex can be assigned to an uncovered slot or a hub-to-successor gap, at most j+k labels. Within each gap the spine order forces the word. Process the free factor 2 per cycle and the slots for cycles and unused hubs at multiplicative cost ≤(2j)^k. Inconsistent codes are harmless for an upper bound.

There is a large uniform family attaining the exponential rate. Put k opening spine vertices first, all beaten by every hub, and k closing vertices last, all beating every hub. The middle n−2k contacts can be chosen arbitrarily, and fix a hub orientation. For n≥3k, reserve k parity-patch spine vertices in the middle. Use one private opening and closing per hub and assign each remaining middle vertex freely to any of the k cycles or j uncovered slots. If a cycle would have even order, include its private patch; otherwise assign the unused patch to slot 1. This makes all cycles odd and disjoint; they have independent slot choices and two weight choices. Consequently,
```
Q_j ≥ (2j)^k (j+k)^(n−3k).
```
There are k disjoint triangles. Deleting the hubs makes the core transitive, so the feedback vertex number is exactly k. The core is strong: through the last closing and the first opening every hub and spine vertex can reach every other one.

For a fixed spine/hub label convention there are 2^(k(n−2k)) distinct cores. An isomorphism class is represented at most binomial(n+k,k)·k! times: choose images of hubs and their labeling, then the image of the spine has its unique increasing transitive order. Thus at fixed k≥3 the sum of reciprocal Q_2 over these **isomorphism types** grows exponentially through each order, up to polynomial factors: the branching base is 2^k, the exact cost base is k+2, and 2^k>k+2. This is not a statement about the number of distinct (H,Q_2) profiles: many distinct cores can collide, so the [square-root tariff](../chanpak_crossing_scales_20261115/HEIGHT_TARIFF.md) does not reverse just because there are enough input types. Strong guarded words are a plausible supply beyond the old alphabet only after quantifying that collision problem or giving the reader access to more of the retained state.

Replay: `python replay_rich_admission.py`, also run with `-O`, checks the old bank cutoff, exact skeleton profiles through four portal pairs, independent literal reads, and guarded configurations. The universal statements use the proofs, not a profile census.
