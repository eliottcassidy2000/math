---
id: THM-4846
title: "Run observers of tournaments, the partial structures of the metagraph's base path. M_k(T) counts ordered spanning covers by k directed paths with prescribed starts and ends (overnight_rank_return THREE_RUN T1). (A) One arc, k = n-1: M_{n-1}(T) is an orientation of the star-transposition Cayley graph Cay(S_n, {(i n)}) (the Nauru graph at n = 4), each arc of T orienting (n-1)! of its edges. Its run-symmetric block is A^T, its sign block is (-1)^(|x-y|-1) A[x][y], and det M_{n-1} = prod_lambda det(M_lambda)^(f^lambda). The strong 4-tournament gives a unimodular orientation; each non-strong one loses exactly one free S_3-orbit. (B) Two arcs, k = n-2: the unordered observer is the second permanental compound, P[V-alpha, V-beta] = per A[beta, alpha]. (C) T5's four kernel directions are the vertex stars, the Hodge duals of the tetrahedron boundaries. Each lives on a copy T5 - w of the strong 4-tournament S4 (4 vertices, 6 arcs = a Hamiltonian 4-cycle plus two diagonals) as +1 on the ends of one diagonal and -1 on the ends of the other. The proof is the alternating in-neighbourhood identity, and T5 is the only class of order 5, 6 or 7 whose plus sector kills the boundaries. (D) Ghost lemma: with sectors P (plus), S (minus), B = (S-P)/2, the number of spurious 2-factors created by forgetting the S_3 (or S_2) glue equals rank(B: ker P -> coker P) over F_2. For T5, ker P over F_2 is the cut space of K5 and B maps it isomorphically onto coker P, which gives the four spurious Z/2. With the four spurious Z/3 of the mod-3 glue (R3) these are the four spurious Z/6. (E) Census thresholds as n grows: order 4 has no torsion and no glue loss; order 5 has the first torsion, and 4 of 12 classes lose to the glue; at order 6 the star law fails (12 of 56) and boundary killers vanish; at order 7 unique-Hamiltonian-cycle classes first carry torsion (5 of 21). Through order 8 every strong class with at least two Hamiltonian cycles has torsion in M_2."
status: >
  PROVED (elementary; mac-mini-2026-10-29-growth; pre-audit): (A) the orientation statement, the trivial and sign
  blocks and the Schur factorization of det M_{n-1}; (B) the compound law; (C) the T5 identity, hence the kernel,
  together with its identification with the inherited tetrahedron-boundary kernel; (D) the ghost lemma at p = 2.
  FINITE-EXACT (scripts in 04-computation/experiments/growth_gf_runs_fermat_20261029, every class from gentourng):
  - (A) at n = 4, 5, 6, with the Nauru isomorphism at n = 4 and the determinant table;
  - (B) against direct enumeration for n <= 7;
  - the T5 Smith forms, reproducing THREE_RUN_GLUE (R7) exactly;
  - the ghost lemma on all 144 (class, k) pairs with n <= 6 and k = 2, 3, and the mod-3 counts;
  - the run-observer censuses n = 3..6 for every k, and n = 7 for k = 2, 3;
  - the star-law census n = 4..8;
  - the Hamiltonian-cycle/torsion tables: k = 2 through n = 8, k = 3 through n = 7, k = 4, 5 at n = 5, 6.
  INHERITED: THREE_RUN.md T4 (the T5 kernel is spanned by the lifted tetrahedron boundaries; the order-4 ranks
  18, 18, 24, 18) and THREE_RUN_GLUE.md (R2)-(R8) (the index-54 chart, the T5 Smith forms, the gcd(m,6)^4 image law).
  CITED: Douglas 1970 (tournaments with exactly one Hamiltonian cycle; the Fibonacci count F(2n-6), recalled and
  attributed to Garey 1972, reproduced here for n = 4..8). The conjectural parts are in HYP-9445.
session: mac-mini-2026-10-29-growth
source: 05-knowledge/results/tournament_growth_gaussian_gf_runs_fermat_20261029.md
depends_on:
  - 05-knowledge/results/overnight_rank_return_20261029/THREE_RUN.md
  - 05-knowledge/results/overnight_research_20261029/THREE_RUN_GLUE.md
related:
  - 01-canon/theorems/THM-4823-the-four-4-vertex-tournaments-census-by-scores-and-one-trace-pattern-free-classes-and-w4-free-equals-locally-transitive-plus-a-sink-tail.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 05-knowledge/results/overnight_rank_return_20261029/COMPLEMENT_KERNEL.md
  - 05-knowledge/hypotheses/HYP-9445-run-observer-torsion-and-hamiltonian-cycles-star-law-fermat-routes-and-gaussian-gf-questions.md
scripts:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/runs.py (shared: observers, S_3 chart, PARI interface)
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/t5_check.py, census_runs.py, summarize.py
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/star_orientation.py, block_checks.py
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/compound_kernel.py, ghost_lemma.py
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/fast_runs.py, torsion_vs_hc.py, torsion_vs_hc2.py
outputs:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/out/ (census_n*_k*.out, compound_n*.out, fast_m2_n7.out, fast_m2_n8.out)
---

# THM-4846 — run observers, the Nauru orientation, the T5 stars and the ghost lemma

## 0. Setting: the base path cut into runs

In the tiling model of the metagraph (01-canon/definitions.md, THM-643), a tiling is a tournament `T` together with
its base Hamiltonian path. A metagraph line replaces `T` by `T^op ⊕ E(P)`, the Hamilton-path complement.

Cutting the base path into `k` runs gives a spanning cover by `k` directed paths. These covers are the partial
structures studied here. Their counting matrix is the overnight packets' observer: for ordered `k`-tuples `s, t` of
distinct vertices, `M_k(T)[s,t]` is the number of ordered covers whose `i`-th path runs from `s_i` to `t_i`
(singletons allowed).

`M_k` commutes with simultaneous permutation of the run positions (right action of `S_k`). Each `k`-set of ports
carries a free `S_k`-orbit. For `k = 3`, the integral Fourier chart `F` of THREE_RUN_GLUE (R2) splits `M_3` over `Q` into:
- the plus sector `P`: unordered covers, start set × end set;
- the minus sector `S`: signed sums;
- the standard sector, twice.

The chart's image has index 54, cut out by the parity and ternary glue (R3). The *unglued* matrix `N = F M F⁻¹`
is the direct sum of the sectors, and "spurious factors" are the invariant factors that `N` has and `M` lacks.

## 1. Statements

**(A) One arc: an orientation of the star graph.**
- *The orientation.* Identify an ordered `(n−1)`-tuple `s` with the permutation `(s_1, …, s_{n−1}, w)`, where `w`
  is the missing vertex. A cover by `n−1` paths has exactly one arc, `s_j → w`. Hence `M_{n−1}[s,t] = 1` iff
  `t = s·(j n)` and `s_j → w`. So `M_{n−1}(T) + M_{n−1}(T)ᵀ` is the adjacency matrix of
  `ST_n = Cay(S_n, {(1 n), …, (n−1 n)})`, and `M_{n−1}(T)` orients each edge `{π, π(j n)}` by the arc between
  `π(j)` and `π(n)`.
- *Edge count.* Every unordered pair of vertices labels `(n−1)!` edges of `ST_n`. At `n = 4` this is 6 per pair, and
  `ST_4` is the Nauru graph: 24 vertices, 36 edges.
- *The run-symmetric and sign blocks.* In complement coordinates the run-symmetric block is `Aᵀ`: the cover with arc
  `x → y` has start set `V − y` and end set `V − x`. The sign block is `(−1)^(|x−y|−1) A[x][y]` in the natural vertex
  order.
- *Determinant.* Since `M_{n−1}` commutes with `S_{n−1}`, Schur gives `det M_{n−1}(T) = Π_{λ ⊢ n−1} det(M_λ)^(f^λ)`,
  where `f^λ` is the dimension of the irreducible representation `λ`.
- *Data.* At `n = 4` the strong tournament `S4` (scores 1122) has `det = 1`. Each of the three non-strong ones has
  Smith form `1^18 0^6`, and its kernel is a copy of the regular representation of `S_3`: sector Smith forms plus
  `1^3 0`, minus `1^3 0`, standard `1^6 0^2` in each copy. A source kills the six rows indexed by the orderings of
  the other three vertices, and a sink kills the corresponding six columns.
  - At `n = 5` the nonzero determinants are `1, 1, 1, 2^8` (T5) and `3^6`.
  - At `n = 6`, 33 of the 56 classes are singular. The nonzero determinants are `1` (eight classes, seven of them
    strong), `2^8`, `2^24 3^28 5^10`, `2^8 3^12 5^10 7^10`, `7^10` (twice), `23^10`, `3^30` (twice), `2^42` (twice),
    `2^36 11^10` and `2^40 13^10` (four classes).
  - Every exponent is a combination of the `S_5`-dimensions 1, 4, 5, 6.

**(B) Two arcs: the second permanental compound.** For `k = n − 2` a cover has two arcs, which match the two non-end
vertices bijectively onto the two non-start vertices. Two arcs cannot close a cycle in a tournament. Hence

    P_{n−2}[V − α, V − β] = per A[β, α] = A[β1,α1] A[β2,α2] + A[β1,α2] A[β2,α1],

with `α` the non-starts and `β` the non-ends.

**(C) The T5 stars.** Let `T5` have arcs `i → i+1, i+2` on `Z/5`.
- *The vectors.* For `w ∈ Z/5` let `c_w` be the alternating sign `(−1)^pos` on `Q = Z/5 − {w}` in increasing order,
  and let `v_w = Σ_{y ≠ w} c_w(y) [{w,y}]`. This is the star of `w`, in the complement coordinates of (B).
- *The theorem.* `P_3(T5) v_w = 0` for every `w`, and `v_0, …, v_4` span the whole rational kernel, which has
  dimension 4.
- *Comparison.* Under end set ↔ complement, `v_w` is the run-symmetric lift of the tetrahedron boundary `∂[Q]`.
  This is THREE_RUN T4's kernel, and the single relation among the five is `∂∂[01234] = 0`.
- *Where the stars live.* `T5 − w` is the strong 4-tournament `S4`. Its six arcs are its unique Hamiltonian cycle
  `w+1 → w+2 → w+3 → w+4 → w+1` plus the diagonals `w+1 → w+3` and `w+2 → w+4`. The sign `c_w` alternates along that
  cycle: it is `+1` on the ends of one diagonal and `−1` on the ends of the other.
- *Uniqueness.* Among all classes of orders 5, 6 and 7, `T5` is the only one whose three-run plus sector kills the
  boundaries (censuses).

**(D) The ghost lemma.** Let `M` be an integral matrix commuting with `S_k`, `k ∈ {2, 3}`. Its 2-local glue is the
single parity condition `u ≡ v (mod 2)` between the plus and minus sectors; for `k = 3` the standard sector splits
off 2-locally. Put `B = (S − P)/2`, which is integral because `S ≡ P (mod 2)`. Then

    spur_2 := rank_F2(M) − rank_F2(N) = rank_F2 [[P, 0], [B, P]] − 2 rank_F2 P = rank( B̄ : ker P̄ → coker P̄ ),

and this is the excess of the number of invariant factors of `N` divisible by 2 over that of `M`.

*For T5.*
- `P` has Smith form `1^6 0^4`.
- The five star indicators reduce mod 2 to the cut space of `K5`, which has dimension 4 and equals `ker P̄`.
- `B̄` maps it isomorphically onto `coker P̄`, so `spur_2 = 4`.
- At 3 the glue ties the trivial and sign characters to the invariant flag of the standard sector mod 3 (R3), and the
  count is `spur_3 = 4` (computed).

So dropping the glue adds `4 × Z/2` and `4 × Z/3`, the four spurious `Z/6`. The rational rank (56) and the free
rank (4) cannot see them (R7).

**(E) What changes as the order grows (FINITE-EXACT).**

| order | three-run observer (k = 3) | two-arc star law (k = n−2) | torsion vs Hamiltonian cycles |
|---|---|---|---|
| 4 | strong `S4` unimodular, non-strong lose one free orbit, no torsion, no glue loss | universal | the unique-HC class is torsion-free |
| 5 | first torsion; glue loss in 4 of 12 classes, `(spur_2, spur_3) = (4,4)` T5, `(2,0)` 12223/H13, `(0,5)` 12223/H15, `(3,3)` source→C3→sink; T5 the unique boundary killer | universal | unique-HC classes (3) torsion-free; HC ≥ 2 classes have torsion |
| 6 | glue loss in 30 of 56 (k = 3) and 18 of 56 (k = 2, the first `S_2` losses); no boundary killer | fails for 12 of 56 classes | unique-HC (8) torsion-free; HC ≥ 2 (27) all torsion for k = 2, 3, 4 |
| 7 | glue loss in 360 of 456 (k = 3); no boundary killer | fails for 168 of 456 | unique-HC: 5 of 21 carry torsion (k = 2, 3); HC ≥ 2 (332) all torsion |
| 8 | — | fails for 2808 of 6880 | k = 2: unique-HC 19 of 55 carry torsion; HC ≥ 2 (5953) all torsion |

The unique-Hamiltonian-cycle counts `1, 3, 8, 21, 55` (n = 4..8) are Douglas's Fibonacci numbers `F(2n−6)`. A
unimodular observer of a strong tournament has occurred only with a unique Hamiltonian cycle. The one-arc observer
`k = n − 1` is too coarse for the torsion law: at `n = 5` (one class) and `n = 6` (two classes), strong classes with
exactly two Hamiltonian cycles have singular but torsion-free `M_{n−1}`.

## 2. Proofs

**(A).** A spanning cover by `n − 1` paths has `n − (n−1) = 1` arc. That arc goes from the start `s_j` of its path to
the vertex `w`, which is not a start, and the end tuple is `s` with `s_j` replaced by `w`. In permutation language
this is right multiplication by the transposition `(j n)`. The reverse move needs `w → s_j`, so exactly one of the
two directions occurs. Each edge of `ST_n` is `{π, π(j n)}`, labelled by the pair `{π(j), π(n)}`, and each pair
labels `n!(n−1)/2 / C(n,2) = (n−1)!` edges.

The block computations are direct:
- *Trivial block.* A run-symmetric vector supported on the orderings of `V − x` is sent to the counts with start set
  `V − y`, and the cover with arc `x → y` is the only contributor.
- *Sign block.* The end tuple is the sorted start tuple with `x` moved to `y`'s slot, a cyclic shift of the
  `|x − y| − 1` elements between them.
- *Schur.* `M` commutes with right multiplication by `S_{n−1}`. On `Q[S_n]`, viewed as a right `S_{n−1}`-module, each
  irreducible `λ` occurs with multiplicity `n f^λ`. So `M ≅ ⊕_λ (M_λ ⊗ I_{f^λ})` and `det M = Π det(M_λ)^(f^λ)`.

**(B).** With `n − 2` paths there are two arcs. Every non-end vertex has exactly one out-arc in the cover and every
non-start vertex exactly one in-arc. The two arcs therefore form a bijection from non-ends to non-starts along arcs
of `T`. Conversely such a bijection gives a cover: in-degree and out-degree are at most 1, and a cycle would need a
loop or a 2-cycle. The unordered cover is determined by the bijection, so the count is the permanent.

**(C).** By (B), `(P̃ v_w)_α = Σ_y c_w(y) per A[{w,y}, α] = A[w,u1] C_w(u2) + A[w,u2] C_w(u1)` for `α = {u1, u2}`,
where `C_w(u) = Σ_{y ≠ w} c_w(y) A[y,u]` is the signed in-neighbour sum.

The in-neighbourhood of `u` is `{u−1, u−2}`. The 4-element set `Q` in cyclic order has alternating signs, including
the wrap-around pair, because `|Q|` is even. Hence:
- if `w` is not an in-neighbour of `u` (`u ∉ {w+1, w+2}`), the two in-neighbours are cyclically consecutive in `Q`
  and `C_w(u) = 0`;
- `C_w(w+1) = c_w(w−1)` and `C_w(w+2) = c_w(w+1)`, and these are opposite because `w−1` and `w+1` are consecutive
  in `Q`.

So `C_w` is supported on `N⁺(w) = {w+1, w+2}` and is antisymmetric there. For `α = N⁺(w)` the two terms cancel. If
`|α ∩ N⁺(w)| = 1`, the surviving `A[w, ·]` multiplies a `C_w`-value outside `N⁺(w)`, which is 0. If
`α ∩ N⁺(w) = ∅`, both `A[w, ·]` vanish.

The five stars satisfy one relation, the boundary of the 4-simplex, and span a 4-dimensional space. The rational
rank of `P_3(T5)` is 6 (Smith form `1^6 0^4`), so the stars span the kernel.

On `T5 − w = S4`, the cyclic order of `Q` is the Hamiltonian cycle, and alternating signs on a 4-cycle are `±` on its
two antipodal pairs, which are the endpoint pairs of the diagonals.

**(D).** 2-locally the coefficient lattice is `Λ ⊕ M_2(Z_(2))` with `Λ = {(u,v) : u ≡ v (2)}`, and the standard
parts of `M` and `N` coincide. In the basis `u = x`, `v = x + 2y` of `Λ`, `M` acts on the principal part by
`[[P, 0], [(S − P)/2, S]]`, which reduces mod 2 to `[[P̄, 0], [B̄, P̄]]`. For such a block-triangular matrix,
`rank = 2 rank P̄ + rank(ker P̄ → coker P̄ induced by B̄)`. The number of invariant factors divisible by 2 is the
dimension minus the mod-2 rank, and both matrices have the same rational rank, which gives the identity.

For T5, `rank_F2 P = 6`. The star indicators have mod-2 rank 4, so they span `ker P̄`; they are the cut space of
`K5`. The connecting rank is 4 (block_checks.py).

**(E).** Census scripts as listed. The torsion statements are computed Smith forms. The conjectural generalizations
are in HYP-9445.

## 3. Relation to the 4-vertex tournament (the owner's question)

- **The tetrahedron.** Each of the four spurious `Z/6` factors is the mod-2 and mod-3 shadow of one T5 kernel
  direction. The kernel directions are the stars `v_w`, one for each deleted vertex `w`, minus one relation. Each star
  is carried by a tetrahedron `T5 − w`, which is the strong 4-tournament `S4` on 4 vertices with 6 arcs. The arcs
  enter as the Hamiltonian 4-cycle, which fixes the alternating signs, and the two diagonals, whose endpoint pairs
  are the `+` and `−` classes.
- **The 4-vertex tournament on its own.** `S4` is exactly the 4-vertex class whose three-run observer is unimodular:
  it orients the Nauru graph with determinant 1. So the tetrahedron alone hides nothing. The four hidden directions,
  and their ghosts, appear only when five copies of `S4` overlap in `T5`.
- **The 6 in Z/6.** It is `|S_3|`: the primes 2 and 3 are where `Z[S_3]` fails to split, and the index of the chart is
  `54 = 2·3³`. That `S4` also has 6 arcs is the coincidence `3! = C(4,2)` (true only at 4), typed NUMEROLOGY. The
  structural bridge is the 24 ordered triples inside a 4-set: they are the 24 elements of `S_4`, the vertices of the
  Nauru graph, and the 4 triangles × 6 orders over which each boundary is lifted.
