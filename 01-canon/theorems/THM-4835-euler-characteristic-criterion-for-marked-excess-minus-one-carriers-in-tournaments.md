---
id: THM-4835
title: "Marked excess -1 carriers in tournaments: the Euler-characteristic criterion. Let D be an oriented subgraph of K_n whose underlying graph G has n-1 edges (a tree, or any disconnected graph with cycle rank = components - 1), with A = Aut(D) acting on V and on E without inversions. The complement of G always has C(n-1,2) pairs. The tournaments containing D are, A-equivariantly and by a pair-level (monomial F2-affine) bijection, the tournaments on an (n-1)-set (a 'linear trijection') IF AND ONLY IF chi(G^K) = 1 for every subgroup K <= A, iff chi(G^g) = 1 for every element g, iff V ≅ E ⊔ {pt} as A-sets. The (n-1)-set can then be taken to be E. The criterion holds for every orientation of every tree (explicit end-edge map) and every functional orientation (explicit root-swap map), recovering the S3 trijection and the vertex-deleted deck. Involution defects are invisible to counts. An odd-order defect is necessary for failure of the weaker A-set tier (when A fixes a vertex) but not sufficient. Explicit infinite families separate the tiers."
status: >
  PROVED: parts 1-6 and 8, and the corollary. Part 4 is now a full equivalence. FINITE-EXACT: part 7.
  Independently audited (audit_20261027: MAJOR 1, MINOR 5, NOTE 8; all repairs applied; MISTAKE-689).
  - Part 2 (sufficiency): Burnside marks for A x Z2-sets.
  - Part 3 (explicit maps): elementary. The end-edge bijection is implicit in the standard proof of Buckley's identity
    W(L(T)) = W(T) - C(n,2) (cf. arXiv:1510.00800, Thm 7.1).
  - Part 4 (full necessity, and the cyclic version): proved by the S8 independent referee (audit_20261027/REPORT.md,
    N1). The session re-checked lemmas L1-L4 before adopting them. This closes HYP-9435(1); the session's own proof
    covered only 2-subgroups and odd prime order.
  - Part 5 (tier law): RESTATED after audit M1. The earlier claim "odd-order defects break counts" was contradicted by
    two n = 7 carriers in the session's own output.
  - Part 7 (census):
    - author: all 1250 oriented carriers with 3 <= n <= 7;
    - referee: the same independently (1250/1250 identical flags), plus the full n = 8 census out of sample (5781
      carriers): linear iff criterion in every case.
  - Part 8: PROVED (finite verification, hand-checkable).
  - Credit and novelty.
    - Session81 (GRAPH_CARRIERS.md §5 "General mark law", commit 134770e3f6, pushed before this session's first
      commit) already states two things. Marked classes are the A-orbits of the affine cube of completions, each element
      fixing 0 or 2^{o} completions. And an actual correspondence with tournaments on n-1 vertices "requires a
      compatible pair-coordinate action".
    - THM-4835 answers exactly when that action exists.
    - The Burnside-ring Euler characteristic [V] - [E] is standard (tom Dieck). Fixed subtrees of finite group actions
      on trees are standard (Serre).
    - No prior classification of equivariant trijections of marked carriers was found.
session: kind-pasteur-2026-10-27-S8
source: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md
scripts:
  - 04-computation/tournament_carriers/carrier_census.py (+ carrier_census_n7.out, .json): census n <= 7
  - 04-computation/tournament_carriers/carrier_core.py, families_and_controls.py (+ .out, .json): P1 families and controls
  - 04-computation/tournament_carriers/tiers_and_sporadics.py (+ .out, .json): tiers, theta_p, group dependence, identity
  - 04-computation/tournament_carriers/audit_20261027/ (independent referee: census n <= 8, necessity lemmas, families, affine example, HYP-9435(3) counterexample)
related: [THM-4836, HYP-9435, LEM-003, THM-002, trijection_four_nodes_three_edges_20261020, session80_tournament_symmetry_20261027, session81_counterexample_families_20261028, Royle-Praeger-Glasby-Freedman-Devillers arXiv:2204.01947, Buckley 1981 (line-graph Wiener index), tom Dieck Transformation Groups]
---

# THM-4835 — The Euler-characteristic criterion for marked excess -1 carriers

## Setting

- **The carrier.** An *oriented carrier* is a loopless oriented graph D on V = [n] with no 2-cycles, whose
  underlying graph G has exactly `n-1` edges. These G are the trees and the disconnected graphs with cycle rank
  `β = c - 1` (c the number of components). Equivalently, G has Euler characteristic `χ(G) = 1`.
- **The symmetry.** `A = Aut(D)` acts on V and on E without inversions: no element reverses an arc. So every fixed
  subgraph `G^K` (K ≤ A) is the induced subgraph `G[V^K]`.
- **The defect** of K is `δ(K) = χ(G^K) - 1`.
- **Completions.** Π is the set of unordered non-edges, with `|Π| = C(n-1,2)`. A tournament containing D is D plus an
  orientation of every pair in Π. The ordered non-edges form an `(A x Z2)`-set Ω, with `Z2 = ⟨swap⟩` acting freely.
- **Linear trijectivity.** D is *linearly trijective* if Ω is isomorphic, as an `(A x Z2)`-set, to the ordered pairs
  of distinct points of some A-set X with `|X| = n-1`. Equivalently, there is an A-equivariant bijection
  `{T ⊇ D} -> Tour(X)` induced by a bijection of pairs that respects the swap: a coordinate-permuting (monomial)
  F2-affine isomorphism.
- **Affine trijectivity is weaker.** It only asks for an A-equivariant F2-affine bijection, and is strictly weaker
  than linear trijectivity. Audit example: a rigid oriented diamond plus 2K_1 (n = 6) is affinely but not linearly
  trijective (audit_20261027 m1, `affine.py`).
- **Weak trijectivity** asks only for an isomorphism of A-sets `{T ⊇ D} ≅ Tour(X)`.

## Statement

1. **(Complement count.)** `|Π| = C(n-1,2) = |E(K_{n-1})|`.

2. **(Sufficiency.)** If `χ(G^K) = 1` for every `K ≤ A`, equivalently `[V] = [E] + 1` in the Burnside ring
   `B(A)`, then A fixes a vertex r, and D is linearly trijective with `X = V \ {r} ≅ E`.

3. **(Explicit trijections.)**
   - **(a) Trees, any orientation: the end-edge map.** Send a non-adjacent pair `{u, w}` to the end edges
     `{e_u, e_w}` of the u-w path, with `u -> w` mapping to `e_u -> e_w`. This is an `Aut(G)`-equivariant signed
     bijection, so `Tour(V) ≅ Or(G) x Tour(E(G))` `Aut(G)`-equivariantly.
   - **(b) Functional orientations: the root-swap map.** Here every vertex except one root r has out-degree 1; then
     G is one tree component plus unicyclic components with directed cycles.
     - The map sends `{v, f(v)}` with `f(v) != r` to `{r, v}`, and is the identity on the other pairs of
       `V \ {r}`.
     - It transports orientations by `(r -> v) <-> (v -> f(v))`.
     - It is an `Aut(D)`-equivariant signed bijection onto Π. The same holds for converses (cofunctional
       orientations).
   - **(c) Special cases.** These include the S3 trijection:
     - the directed path;
     - the directed `(n-1)`-cycle plus solo;
     - the uniform star;
     - the deck: an in-star at v gives `T -> T - v`.

4. **(Full equivalence; referee's theorem, audit_20261027 N1.)** The following are equivalent:
   - (a) D is linearly trijective;
   - (b) `χ(G^K) = 1` for every subgroup `K ≤ A`;
   - (c) `χ(G^{⟨g⟩}) = 1` for every element `g ∈ A`.

   Consequences:
   - in a linearly trijective carrier, components with nonzero Euler characteristic are pairwise non-isomorphic
     (in particular there is at most one isolated vertex);
   - A fixes a vertex.

   The hypothesis that A is the *full* automorphism group is used (lemma L1), and part 8 shows it cannot be dropped.

5. **(Tier law, restated after audit M1.)**
   - **(a) Only odd-order subgroups count, PROVED.** A subgroup of even order fixes no completion: an involution has a
     2-cycle `(a b)`, `{a,b}` is a non-edge, and it is reversed. So the A-set `{T ⊇ D}`, and every count derived from
     it, depends only on marks at odd-order subgroups. Involution defects are invisible to it.
   - **(b) A sufficient condition for weakness, PROVED.** If A fixes a vertex r and `δ(K) = 0` for every odd-order
     K, then D is weakly trijective with `X = V \ {r}`. Hence a non-weak carrier with an A-fixed vertex has an
     odd-order defect.
   - **(c) The converse of (b) is false, FINITE-EXACT.** The two n = 7 carriers `TT_4 + 3K_1` and
     `(strong T_4) + 3K_1` have `δ(Z_3) = -3` but are weakly trijective. Their marked count, 5472, equals
     `|Tour(X)/S_3|` for X the regular `S_3`-set. The mechanism is an alternative fixed count `x = 0 != v - 1`.

6. **(Explicit families.)**
   - **(a) Bigons `B_k`** (two directed `s -> t` paths of length k, plus a solo; `n = 2k+1`, `k >= 2`). Not linear,
     with `δ = 2` at the path swap. Weak, with `2^{C(2k,2)-1}` marked types.
   - **(b) Alternating cycles** `C_{2k} + K_1`, `k >= 2`. Not linear, with `δ = 2` at a vertex reflection. Weak.
   - **(c) Dense cores.** A core H with `|E(H)| = |V(H)| + j` (`j >= 1`) plus `j+1` isolates. Not linear: swapping two
     isolates gives `δ = -2`.
   - **(d) Odd thetas.** `θ_p^k + (p-1)K_1`: p directed `u -> w` paths of length k, p an odd prime, `k >= 2`. Not
     weakly trijective: equal `Z_p` marks would force a fixed count `x = p+1`, while `x ≡ n-1 ≡ 0 (mod p)`.
   - **Counts versus A-sets for the thetas (audit m3).**
     - For `p = 3` (all `k >= 2`) the *count* also matches no `|Tour(X)/A|`: an involution acting trivially on X
       overshoots the Burnside sum, and otherwise equal counts force equal `Z_3` marks.
     - For `(p, k) = (5, 2)` this was checked.
     - For other p the count tier is OPEN.

7. **(Census, FINITE-EXACT.)** Oriented carriers up to isomorphism, with tier counts:

   | n | carriers | linear | weak only | neither |
   |---|---|---|---|---|
   | 3 | 3 | 3 | 0 | 0 |
   | 4 | 10 | 10 | 0 | 0 |
   | 5 | 41 | 39 | 2 | 0 |
   | 6 | 187 | 167 | 20 | 0 |
   | 7 | 1009 | 816 | 188 | 5 |
   | 8 (referee) | 5781 | 4155 | 1469 | 157 |

   - "Linear iff (b) iff (c)" holds in all 7031 cases.
   - Every non-linear carrier with n ≤ 8 has an involution defect (1841/1841). This is NOT a law: the referee's
     n = 10 carrier `K_4(C_3 + source) + K_1 + K_2 + vec P_3`, with `A = Z_3` and no involutions, is non-linear and
     non-weak (HYP-9435(3), refuted).
   - The 5 "neither" carriers at n = 7 are:
     - `K_4(C_3 ± apex) + 3K_1`;
     - `θ_3^2 + 2K_1`;
     - `K_{3,2} + 2K_1` (both directions; the only n ≤ 7 carriers whose A fixes no vertex).

     Their counts match no `|Tour(X)/A|`.

8. **(Group dependence; PROVED by finite verification.)** `TT_4 + 3K_1` is linearly trijective for the subgroup
   `Z_3 ≤ A = S_3`, with X two free orbits: `C(x,2) = C(4,2) - 6` forces `x ∈ {0,1}`, and `x ≡ 0 (mod 3)` gives
   `x = 0`. It is not linear for `Z_2`, where `C(x,2) = 4` has no solution, nor for `S_3`. The transposition that the
   subgroup forgot is exactly the witness of lemma L1.

## Proofs

**1.** `C(n,2) - (n-1) = C(n-1,2)`.

**2.**
- *Reduction to marks.* A linear trijection is an `(A x Z2)`-isomorphism `Ω ≅ Ω_X`. By Burnside's theorem it exists
  iff the marks agree at every subgroup `L ≤ A x Z2`. Subgroups containing the swap fix nothing on either side. Every
  other L is the graph of some `ψ: K -> Z2`.
  - For `ψ = 0`: the marks are `2[C(v_K,2) - e_K]` and `2 C(x_K,2)`.
  - For `ψ != 0`, with `K0 = ker ψ`: the marks are `v_{K0} - v_K` and `x_{K0} - x_K`. A swapped pair is never an edge.
- *Choice of X.* The hypothesis at `K = A` gives a fixed vertex r. With `X = V \ {r}` we have `x_K = v_K - 1`, so the
  `ψ != 0` marks agree. The `ψ = 0` marks agree iff `e_K = v_K - 1`, which is the hypothesis.
- *Identifying X.* `X ≅ E` follows from marks. `□`

**3(a).**
- *Bijection.* A non-adjacent pair determines its path `x_0 ... x_k` (`k >= 2`) and two distinct end edges.
  Conversely, two distinct edges lie on a unique minimal path whose outer ends are at distance `>= 2`. These maps are
  inverse to each other, natural, and swap-compatible.
- *Consistency with part 2.* Fixed sets of finite groups acting on trees without inversion are nonempty subtrees.

**3(b).**
- *Bijection.* Every edge is `{v, f(v)}` for a unique tail v. The edges inside `V \ {r}` are exactly as many as the
  non-edges `{r, v}`.
- *Equivariance.* `Aut(D)` fixes r and commutes with f.

**3(c).** The three S3 carriers and the in-star are trees or (co)functional, so part 3(a) or 3(b) applies.

**4 (referee's proof, re-checked).** Write `v_K, e_K, x_K` and `δ_K = v_K - e_K - 1`. For a K-invariant component C
put `c_K(C) = χ(C ∩ G^K) - χ(C)`.

- **L1 (twins).**
  - Suppose components `C != C'` are isomorphic via `φ`. Then `τ = φ ⊔ φ^{-1} ⊔ id ∈ A`, because components share no
    arcs, and `δ(⟨τ⟩) = -2χ(C)`.
  - Under (a): the `ψ != 0` marks make `v_K - x_K` constant along index-2 chains, so for 2-subgroups
    `x_K = v_K - 1` and `δ = 0`. Hence `χ(C) = 0`.
  - Under (c): `χ(C) = 0` directly.
- **L2 (sign).** If twinned components have `χ = 0`, then for every K, `δ_K = Σ_{C K-invariant} c_K(C) >= 0`:
  - non-invariant components are twinned and contribute 0;
  - for an invariant C with `C^K` nonempty, `c_K(C) = rank H_1(C, C^K) >= 0`, since `H_0(C, C^K) = 0`;
  - if `C^K` is empty, C is not a tree (trees have fixed vertices), so `c_K(C) = -χ(C) >= 0`.

  In particular `v_K >= 1`.
- **L3 (solvable K, assuming (a)).**
  - The `ψ = 0` mark gives `C(x_K,2) = C(v_K-1,2) + δ_K`, with `0 <= δ_K <= v_K - 1`.
  - So `φ_K = x_K - v_K + 1 ∈ {-1, 0, 1}`, and `φ_K = 0` gives `δ_K = 0`.
  - For `H ◁ H'` of prime index p, `φ_H ≡ φ_{H'} (mod p)`, by counting fixed points of `H'/H` on `X^H` and `V^H`.
    For `p >= 3` this forces equality; for `p = 2`, the `ψ != 0` mark gives equality.
  - Since `φ_1 = 0`, every solvable K has `δ_K = 0`. This proves (a) ⇒ (c).
- **L4 (cyclic ⇒ all).** Assume L2's hypothesis and `δ(⟨h⟩) = 0` for every h. Fix K and a K-invariant component C.
  By L2, `c_{⟨h⟩}(C) = 0` for all `h ∈ K`.
  - If C is a tree, `c_K(C) = 0`.
  - If some `C^h` is empty, then `χ(C) = 0` and `C^K` is empty, so `c_K(C) = 0`.
  - Otherwise each `C^h` is connected (`H_1(C, C^h) = 0`) and contains every cycle, hence the 2-core.
  - Then `C^K = ∩_h C^h` contains the 2-core. It is connected, because the unique path from a fixed vertex to the
    2-core is fixed. So `χ(C^K) = χ(C)`.
  - Hence `δ_K = 0`, proving (c) ⇒ (b). Part 2 gives (b) ⇒ (a). `□`

**5.**
- (a) is the odd-order rule.
- (b), odd K: Burnside on unordered pairs gives equal orbit counts on Π and on `pairs(V \ {r})`, because
  `C(v_k,2) - e_k = C(v_k - 1, 2)` for every `k ∈ K`.
- (b), even K: K fixes nothing on either side, since its involutions have 2-cycles inside `V \ {r}`.
- (c) is a finite verification (`tiers_and_sporadics.out` (a); audit M1).

**6.**
- The defects in (a)–(c) are as stated, and part 4 applies.
- Weakness in (a)–(b) follows from 5(b): `B_k` has `A = Z_2` fixing the solo. For the alternating cycles, the odd-order
  elements are rotations whose fixed set is the solo.
- The `B_k` count: `Z_2` acts freely, so there are `2^{C(2k,2)}/2` marked types.
- (d):
  - The marks are `#orbits = (C(n-1,2) + (p-1)·Fix)/p`, with `Fix = C(p+1,2)` on Π and `C(x,2)` on X.
  - So `x = p+1`, contradicting `x ≡ 0 (mod p)`.
  - The count for `p = 3` is in audit m3. `□`

**7, 8.** Finite computations, author and referee; part 8 is also hand-checked above.

## Corollary (score-product averaging, PROVED)

- **Counting.** T contains exactly `∏_{v≠r} s_v(T)` functional sub-carriers rooted at r: choose an out-neighbour for
  every `v != r`.
- **Double counting.** Each labelled functional carrier lies in `2^{C(n-1,2)}` tournaments, so
  `Σ_T Σ_r ∏_{v≠r} s_v(T) = F_n 2^{C(n-1,2)}`.
- **The constant.** `F_n = n Σ_k (-1)^k (n-1)! (n-1)^{n-1-2k} / (k! 2^k (n-1-2k)!)` counts the rooted functional
  digraphs without loops or 2-cycles. Its first values are 9, 72, 815, 11700, 201747 for `n = 3..7`; the referee
  checked n = 7.

## Reading

- **The criterion.** Marking an excess -1 carrier costs exactly one vertex. The cost is paid equivariantly, through
  the carrier's edge set, iff every symmetry fixes a sub-carrier of Euler characteristic 1. Equivalently, the carrier
  is a point in the Burnside ring: `[V] - [E] = 1`.
- **The obstructions.**
  - Symmetric bigons and alternating even cycles have involutions with `δ = 2`.
  - Dense cores with twin isolates have `δ = -2` (L1).
  - Odd thetas have `δ = p` (Smith: `δ ≡ 0 (mod p)`).
- **The tiers.**
  - Involution defects destroy only the signed (linear) structure.
  - An odd-order defect is needed to destroy the A-set structure when A fixes a vertex, but need not destroy it.
- **The S3 trijection** is the case where every tier holds.
