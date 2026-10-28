---
id: THM-4835
title: "Marked excess -1 carriers in tournaments: the Euler-characteristic criterion. Let D be an oriented subgraph of K_n whose underlying graph G has n-1 edges (a tree, or any disconnected graph with cycle rank = components - 1), with automorphism group A = Aut(D) acting on V and on E without inversions. The complement of G always has C(n-1,2) pairs. If V ≅ E ⊔ {pt} as A-sets, equivalently chi(G^K) = 1 for every subgroup K <= A, then the tournaments containing D are A-equivariantly and F2-affinely the tournaments on the (n-1)-set E (a 'linear trijection'). This holds for every orientation of every tree (explicit end-edge map) and every functional orientation (explicit root-swap map). It recovers the S3 trijection (line, cycle+solo, star) and the vertex-deleted deck. Conversely, linear trijectivity forces chi(G^K) = 1 for every 2-subgroup and every subgroup of odd prime order. The 2-part of the symmetry is invisible to counts and to the A-set of tournaments; the odd part governs them. Explicit infinite families separate the tiers."
status: >
  PROVED: parts 1-6 and 8. FINITE-EXACT: part 7 (full equivalence for every oriented carrier with n <= 7). Not yet
  audited (pre-audit checkpoint).
  - Sufficiency (part 2) is proved by Burnside marks for A x Z2-sets.
  - The explicit maps (part 3) are elementary.
  - Necessity (part 4) is proved for 2-subgroups (index-2 chain argument) and for subgroups of odd prime order (full
    automorphism group, a transposition argument plus relative homology of graphs). Necessity at the remaining
    subgroups (odd composite order, or mixed order) is OPEN: HYP-9435(1).
  - Census (part 7), carrier_census.py, all 1250 oriented carriers with 3 <= n <= 7. The two uniform stars at n = 7
    (|A| = 720) were not run through the group search; they are covered by part 3:
    - linear and the criterion coincide in every case;
    - the cyclic-subgroup version coincides with the all-subgroup version in every case.
  - Pre-registered confirmations (04-computation/tournament_carriers/PREREGISTRATION.md, P1, pushed in efdf0b0a5f
    before the runs):
    - bigons and alternating cycles up to n = 11;
    - 600 random trees and functional orientations at n = 8, 9.
  - The criterion was formulated from n <= 6 data, and its n = 7 check is exploratory.
  - Novelty, as far as checked (S3 note; session80, session81 packets; RPGFD 2023 abstract), all searched for the
    notion: no prior classification of equivariant trijections for marked carriers was found.
    - The Burnside-ring Euler characteristic [V] - [E] is standard (tom Dieck).
    - The end-edge bijection between non-adjacent pairs and pairs of edges of a tree is elementary and plausibly
      folklore.
session: kind-pasteur-2026-10-27-S8
source: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md
scripts:
  - 04-computation/tournament_carriers/carrier_census.py (+ carrier_census_n7.out, .json): census n <= 7
  - 04-computation/tournament_carriers/carrier_core.py, families_and_controls.py (+ .out, .json): P1 families and controls
  - 04-computation/tournament_carriers/tiers_and_sporadics.py (+ .out, .json): tiers, theta_p, group dependence, identity
related: [THM-4836, HYP-9435, LEM-003, THM-002, trijection_four_nodes_three_edges_20261020, session80_tournament_symmetry_20261027, session81_counterexample_families_20261028, Royle-Praeger-Glasby-Freedman-Devillers arXiv:2204.01947, tom Dieck Transformation Groups]
---

# THM-4835 — The Euler-characteristic criterion for marked excess -1 carriers

## Setting

- **The carrier.** An *oriented carrier* is a loopless oriented graph D on V = [n] with no 2-cycles, whose
  underlying graph G has exactly `n-1` edges.
  - Such G are the trees and the disconnected graphs with cycle rank `β = c - 1`, where c is the number of
    components.
  - Its *excess* is `|E| - |V| = -1`, i.e. Euler characteristic `χ(G) = 1`.
- **The symmetry.** Write `A = Aut(D)`, the permutations of V that map arcs to arcs. A acts on V and on E. No element
  reverses an arc, so G is an A-complex without inversions, and for every subgroup `K <= A` the fixed subgraph
  `G^K = G[V^K]` satisfies `|E^K| = e(G[V^K])`.
- **The defect** is `δ(K) = χ(G^K) - 1 = |V^K| - |E^K| - 1`.
- **The complement pairs.** Π is the set of unordered non-edges of G, `|Π| = C(n,2) - (n-1) = C(n-1,2)`. A tournament
  containing D is exactly D plus a choice of orientation on each pair of Π. A acts on Π, and an element may swap the
  two ends of a pair; this makes Π a *signed A-set*, i.e. an `(A x Z2)`-set of ordered pairs on which `Z2 = ⟨swap⟩` acts
  freely.
- **Linear trijectivity.** D is *linearly trijective* if, for some A-set X with `|X| = n-1`, the signed A-set Π is
  isomorphic to the signed pair set of X. Equivalently, there is an A-equivariant bijection
  `{T ⊇ D} -> Tour(X)` induced by a pair-level bijection, i.e. an F2-affine isomorphism.
- **Weak trijectivity.** D is *weakly trijective* if `{T ⊇ D}` and `Tour(X)` are merely isomorphic A-sets.
- **Marked types.** The number of marked types `|{T ⊇ D}/A|` is the number of isomorphism types of pairs `(T, D')`
  with `D' ≅ D` in T.

## Statement

1. **(Complement count.)** For every carrier, `|Π| = C(n-1,2) = |E(K_{n-1})|`. Hence `|{T ⊇ D}| = |Tour(n-1)|`.

2. **(Sufficiency.)** Suppose `χ(G^K) = 1` for every subgroup `K <= A`. In Burnside-ring language, `[V] = [E] + 1`
   in `B(A)`, i.e. `V ≅ E ⊔ {pt}` as A-sets. Then:
   - A fixes a vertex r;
   - D is linearly trijective with `X = V \ {r} ≅ E`;
   - consequently `{T ⊇ D} ≅ Tour(E)` A-equivariantly, and the number of marked types is `|Tour(n-1)/A|` for A
     acting on the `n-1` edges.

3. **(Explicit trijections.)**
   - **(a) Trees, any orientation: the end-edge map.** Send a non-adjacent pair `{u, w}` to the two end edges of the
     u-w path, and an orientation `u -> w` to `e_u -> e_w`. This is an `Aut(G)`-equivariant bijection
     `Π -> C(E, 2)`. Hence `Tour(V) ≅ Or(G) x Tour(E(G))`, `Aut(G)`-equivariantly, by `T -> (T|_G, T_G)`.
   - **(b) Functional orientations.** Here every vertex except one root r has out-degree exactly 1, so G is one tree
     component plus unicyclic components with directed cycles. The root-swap map does the following:
     - it sends a pair `{v, f(v)}` with `f(v) != r` to `{r, v}`;
     - it is the identity on the other pairs of `V \ {r}`;
     - it transports orientations by `(r -> v) <-> (v -> f(v))`.

     This is an `Aut(D)`-equivariant signed bijection `pairs(V \ {r}) -> Π`. The same holds for the converses
     (cofunctional orientations).
   - **(c) Special cases.** These include the S3 trijection:
     - the directed path (A = 1);
     - the directed `(n-1)`-cycle plus solo (`A = Z_{n-1}`);
     - the uniform star (`A = S_{n-1}`);
     - the deck: an in-star at v gives `T -> T - v`.

4. **(Necessity, for the full automorphism group.)** If D is linearly trijective, then `χ(G^K) = 1` for:
   - (i) every 2-subgroup `K <= A`;
   - (ii) every subgroup of odd prime order.

5. **(The odd-order rule splits the tiers.)**
   - A subgroup `K <= A` of even order fixes no tournament containing D, because an involution of K has a 2-cycle
     `(a b)`, `{a,b}` is a non-edge, and it is reversed.
   - Hence the A-set `{T ⊇ D}`, and every count derived from it, only sees odd-order subgroups.
   - **Lemma.** If A fixes a vertex r and `δ(K) = 0` for every odd-order K, then D is weakly trijective with
     `X = V \ {r}`.

6. **(Explicit families.)**
   - **(a) Bigons `B_k`.** Two directed `s -> t` paths of length k plus an isolated vertex, `n = 2k+1`, `k >= 2`.
     `B_k` is not linearly trijective (the path swap fixes `{s, t, solo}`, so `δ = 2`). It is weakly trijective, and
     has exactly `2^(C(2k,2)-1)` marked types.
   - **(b) Alternating cycles.** `C_{2k} + K_1` with alternating orientation, `k >= 2`, is not linear (a reflection
     fixes two opposite vertices and the solo, `δ = 2`), but it is weak.
   - **(c) Dense cores.** A core H with `|E(H)| = |V(H)| + j` (`j >= 1`) together with `j+1` isolated vertices is not
     linear: a transposition of two isolates has `δ = -2`.
   - **(d) Odd thetas.** `θ_p^k + (p-1)K_1` consists of p directed `u -> w` paths of length k plus `p-1` isolated
     vertices, `n = pk+1`, with p an odd prime and `k >= 2`. It is not even weakly trijective. Its `Z_p` rotation of
     the paths fixes `2^m` tournaments containing D for some m. On any A-set X of size `n-1`, a `Z_p` mark of the
     same size would force `x = |X^{Z_p}| = p+1`, while `x ≡ n-1 ≡ 0 (mod p)`.

7. **(Census, FINITE-EXACT.)** For `3 <= n <= 7` there are 1250 oriented carriers up to isomorphism (3, 10, 41, 187,
   1009). In every case, linear trijectivity holds iff `χ(G^K) = 1` for all `K <= A`.

   | n | carriers | linear | weak only | neither |
   |---|---|---|---|---|
   | 3 | 3 | 3 | 0 | 0 |
   | 4 | 10 | 10 | 0 | 0 |
   | 5 | 41 | 39 | 2 | 0 |
   | 6 | 187 | 167 | 20 | 0 |
   | 7 | 1009 | 816 | 188 | 5 |

   - Every non-linear carrier has an involution defect (215/215).
   - The 5 "neither" carriers at n = 7 are:
     - `K_4` (C3+sink) + 3K_1 and `K_4` (C3+source) + 3K_1;
     - `θ_3^2 + 2K_1`;
     - `K_{3,2}` (bipartite, either direction) + 2K_1. These two orientations are the only carriers with n <= 7 whose
       A fixes no vertex.
   - Each of the 5 has an odd-order defect. Their marked counts (1888, 1888, 2816, 1408, 1408) equal no
     `|Tour(X)/A|` for any A-set X of size 6.

8. **(Group dependence.)** `K_4` (transitive) + 3K_1 is linearly trijective for the subgroup `Z_3 <= A = S_3`, with X
   two free 3-orbits. It is not linear for `Z_2` or for the full `S_3`. The criterion is a property of the full
   symmetry group; the witness that kills it is a transposition the subgroup forgot.

## Proofs

**1.** `C(n,2) - (n-1) = C(n-1,2)`.

**2.**
- *Reduction to marks.* Let `Ω` be the ordered non-edges and `Ω_X` the ordered pairs of distinct points of X, both
  `(A x Z2)`-sets with `Z2` free. A linear trijection is exactly an `(A x Z2)`-isomorphism `Ω ≅ Ω_X`; transport the
  section `T \ D` of `Ω` to a section of `Ω_X`. By Burnside's theorem, two finite G-sets are isomorphic iff
  `|Ω^L| = |Ω_X^L|` for every subgroup `L <= A x Z2`.
- *Which subgroups matter.* Subgroups containing the swap have no fixed points on either side. Every other L is the
  graph `{(k, ψ(k))}` of a homomorphism `ψ: K -> Z2`.
  - If `ψ = 0`, then `|Ω^L| = 2[C(|V^K|,2) - |E^K|]` and `|Ω_X^L| = 2 C(|X^K|, 2)`.
  - If `ψ != 0`, put `K0 = ker ψ`. Then `Ω^L` consists of the ordered pairs (a, b) fixed pointwise by K0 and swapped
    by `K \ K0`. Such a pair is never an edge, since an automorphism would reverse its arc. So
    `|Ω^L| = |V^{K0}| - |V^K|`, and likewise `|Ω_X^L| = |X^{K0}| - |X^K|`.
- *The case of a fixed vertex.* Hypothesis 2 at `K = A` gives `|V^A| >= 1`; choose `r ∈ V^A` and `X = V \ {r}`. Then
  `|X^K| = |V^K| - 1` for all K, so the `ψ != 0` marks agree. The `ψ = 0` marks agree iff
  `C(v,2) - e = C(v-1,2)`, i.e. `e = v - 1`, which is hypothesis 2.
- *The A-set X.* Finally `V ≅ X ⊔ {r}` and `V ≅ E ⊔ {pt}` give `X ≅ E`, by marks. `□`

**3(a).**
- *The map is a bijection.* In a tree, a non-adjacent pair `{u, w}` determines its path `u = x_0 ... x_k = w`, with
  `k >= 2`, and so two distinct end edges `{x_0 x_1}`, `{x_{k-1} x_k}`. Conversely, two distinct edges e, e' lie on a
  unique minimal path. Its outer endpoints are at distance `>= 2`, and its end edges are e and e'. These maps are
  mutually inverse.
- *Equivariance.* The map is natural, hence `Aut(G)`-equivariant. Reversing `u -> w` reverses `e_u -> e_w`.
- *The orientation set.* The orientations of G form `Or(G) = {±}^E`, and `T` corresponds to `(T|_G, T_G)`.
- *Consistency with part 2.* For a finite group acting on a finite tree without inversions, every fixed set is a
  nonempty subtree, so `χ = 1`.

**3(b).**
- *Edge count.* Every edge of G is `{v, f(v)}` for a unique tail v. The edges inside `V \ {r}` are those with
  `f(v) != r`, and they are exactly as many as the non-edges `{r, v}`.
- *The map.* The map is a bijection and respects the swap.
- *Equivariance.* `Aut(D)` fixes r, the unique vertex of out-degree 0, and commutes with f, so the map is
  equivariant.
- *Consistency with part 2.* Fixed sets are f-closed and contain r, so `|E^K| = |V^K| - 1`.

**3(c).** The three S3 carriers and the in-star are functional or cofunctional (path, cycle+solo, star) or trees, so
part 3 applies.

**4(i).**
- *An invariant constant along chains.* For `K0 < K` of index 2, the `ψ != 0` marks give
  `|V^{K0}| - |V^K| = |X^{K0}| - |X^K|`. So `d(K) = |V^K| - |X^K|` is constant along index-2 chains.
- *2-subgroups.* Every 2-group has such a chain from 1, and `d(1) = n - (n-1) = 1`. Hence `|X^K| = |V^K| - 1`, and the
  `ψ = 0` mark gives `|E^K| = |V^K| - 1`.

**4(ii).** Let g have odd prime order p, write `v = |V^g|` and `x = |X^g|`, and put `δ = δ(⟨g⟩)`.
- *Moved components have χ = 0.* Suppose a component C of G lies in a nontrivial g-orbit; the orbit has `p >= 3`
  isomorphic copies of C. The transposition of two copies is an automorphism of D with fixed-set Euler characteristic
  `1 - 2χ(C)`, so 4(i) forces `χ(C) = 0`.
- *Fixed components have nonnegative defect.* For a component C fixed setwise, either `C^g` is nonempty and
  `χ(C) - χ(C^g) = -rank H_1(C, C^g) <= 0`, or `C^g` is empty and C is not a tree (trees have fixed vertices), so
  `-χ(C) >= 0`.
- *Hence `δ >= 0`.*
- *The marks.* The `ψ = 0` mark gives `C(x,2) = C(v,2) - e`, with `e = v - 1 - δ >= 0`. Counting p-cycles,
  `x ≡ n-1` and `v ≡ n (mod p)`, so `x = v - 1 + pa` for an integer a. Then `δ = pa(2v - 3 + pa)/2`.
  - If `a > 0`, then `δ >= 3v > v - 1`, contradicting `e >= 0`.
  - If `a < 0`, then `x >= 0` gives `v >= 1 + p|a|`, so `2v - 3 - p|a| > 0` and `δ < 0`, contradicting `δ >= 0`.
  - Hence `a = 0` and `δ = 0`. `□`

**5.**
- *Even-order subgroups.* The first claim is the odd-order rule: an even cycle of an automorphism has a power that
  reverses a chord.
- *The lemma, odd K.* For K of odd order, Burnside on unordered pairs gives `#K-orbits = (1/|K|) Σ Fix(k)`.
  - On Π, `Fix(k) = C(v_k,2) - e_k`.
  - On `X = V \ {r}`, `Fix(k) = C(v_k - 1, 2)`.
  - These agree when every `δ(⟨k⟩) = 0`.
- *The lemma, even K.* Here K fixes nothing on either side: its involutions have 2-cycles in `V \ {r}`.

**6.**
- *(a)–(c).* The defects are as stated, so part 4(i) applies. Weakness follows from part 5: for `B_k`, `A = Z_2`
  fixes the solo. For the alternating cycles, the only odd-order elements are rotations, whose fixed set is the solo,
  so they have no defect.
- *(a), the count.* There `Z_2` acts freely on `{T ⊇ D}`, giving `2^{C(2k,2)}/2`.
- *(d).*
  - *The mark on Π.* For `ρ` of order p, the Burnside count is `#orbits = (C(n-1,2) + (p-1)·Fix)/p`. On Π,
    `Fix = C(p+1, 2)`: the fixed vertices are u, w and the `p-1` isolates, with no edges among them.
  - *The mark on X.* On X, `Fix = C(x,2)`.
  - *The contradiction.* Equal marks force `x = p+1`, while `x ≡ n-1 = pk ≡ 0 (mod p)`. `□`

**7, 8.** These are finite computations (scripts above), with independent re-checks pending audit.

## Corollary (score-product averaging, PROVED)

- **Counting choice functions.** A functional sub-carrier of T rooted at r is a choice of an out-neighbour `f(v)` for
  every `v != r`. So T contains exactly `∏_{v≠r} s_v(T)` of them, where s is the score.
- **Double counting.** Each labelled functional carrier lies in `2^{C(n-1,2)}` tournaments, so

      Σ_T Σ_r ∏_{v≠r} s_v(T) = F_n 2^{C(n-1,2)},

  with `F_n = n Σ_k (-1)^k (n-1)! (n-1)^{n-1-2k} / (k! 2^k (n-1-2k)!)`. These are the rooted functional digraphs
  without loops or 2-cycles; the first values are `F_n = 9, 72, 815, 11700` for `n = 3..6` (verified by brute force).

## Reading

- **The criterion.** Marking an excess -1 carrier costs exactly one vertex of tournament information, and the cost
  is paid *equivariantly* precisely when every symmetry fixes a sub-carrier that is again of Euler characteristic 1.
- **The obstructions** are the symmetric bigon, the alternating even cycle, the dense core with isolates, and the
  odd theta. In each, a symmetry's fixed set has the wrong Euler characteristic.
  - For involutions the error is always even (`χ(G^τ) ≡ χ(G) mod 2`).
  - For elements of odd prime order the error is a multiple of p (Smith).
- **The tiers.** The 2-part of the symmetry is visible only linearly (signed structure). The odd part is visible to
  counts.
- **The S3 trijection** is the case where both tiers hold.
