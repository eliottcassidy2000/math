---
id: THM-4624
title: "Finite Beckman-Quarles for the unit-circle Cayley graphs G_q = Cay(F_(q^2), mu_(q+1)): for every prime power q <= 81 except q in {4, 8, 9}, every graph automorphism is an affine semilinear map z -> a sigma(z) + b (a in mu_(q+1), sigma in Gal(F_(q^2)/F_p), b in F_(q^2)), so |Aut G_q| = q^2 (q+1) 2e for q = p^e and the vertex stabiliser is mu_(q+1) x| Gal with 1 + (Frobenius orbits on F_q^x) orbits; the unit circle is a union of lines over a subfield F_(p^f) exactly when p^f - 1 divides q + 1, and it is a basis or a frame (at most dim + 1 lines in general position) exactly for q in {2, 3, 4, 8, 9}: G_2 = K_4, G_3 = K_3 x K_3, G_4 = the Clebsch graph (Aut = 2^4 x| S_5, order 1920), G_8 = the Hamming graph H(3,4) = K_4^(x3) (Aut = S_4 wr S_3, order 82944), G_9 = the ternary five-frame graph Cay(F_3^4, {+-w_0, ..., +-w_4}) with w_0 + ... + w_4 = 0 (Aut = F_3^4 x| (S_5 x C_2), order 19440); the frame structure explains kappa(9) = 3 (a nowhere-zero functional) and kappa(4) = kappa(8) = 4"
status: >
  FINITE-EXACT (nauty/dreadnaut, every prime power q <= 81: |Aut G_q|, vertex-transitivity, the
  order and orbit count of the stabiliser of 0; isomorphism tests and SAT chromatic numbers for
  q = 4, 8, 9; 04-computation/experiments/logic4_20261012/gq_automorphisms.py, gq_exceptional.py).
  PROVED: the affine semilinear group embeds with the stated order (so equality of orders is
  equality of groups), the subfield-line lemma and the classification of the basis/frame cases
  (elementary number theory), the three identifications, the automorphism groups of the three
  exceptional graphs (explicit subgroups of the stated orders, matching nauty's totals; the Clebsch
  and Hamming groups are KNOWN), and the frame colouring of G_9. KNOWN context: Beckman-Quarles
  (1953) for R^n; Kiermaier-Kurz (Discrete Math. 309 (2009)) prove the analogous affine rigidity for
  the INTEGRAL-distance graph on F_q^2 (q = 3 mod 4, or q prime other than 5), a different graph.
  Affine rigidity for all q not in {4, 8, 9} is HYP-9256. Independently audited 2026-10-12 (audit B:
  PASS; every order, ratio, stabiliser order and orbit count re-derived with independent field
  arithmetic and nauty driver for q = 2, 3, 4, 5, 7, 8, 9, 13, 16, 25, 27; the divisibility and frame
  classification brute-forced for p <= 97, e <= 12; two wording nits in the Reading fixed). The audit
  also extended the census: Aut(G_q) = Aff(q) for every prime power 83 <= q <= 256 (38 values;
  audit_B_gq_census_*.out), so the FINITE-EXACT range is now q <= 256.
session: mac-mini-2026-10-12-logic4
source: 05-knowledge/results/logic_lenses_jc2_collatz_hn_20261012.md
scripts:
  - 04-computation/experiments/logic4_20261012/gq_automorphisms.py (+ .out)
  - 04-computation/experiments/logic4_20261012/gq_exceptional.py (+ .out)
related:
  - THM-4558 (kappa(q) = chi(G_q) and the residue-colouring bound for field planes; the table of kappa)
  - THM-4552 (G_7 is the 7 x 7 knight torus)
  - THM-4622 (Beckman-Quarles for R^2 and Q^2)
  - HYP-9256 (affine rigidity for all q not in {4, 8, 9})
  - M. Kiermaier, S. Kurz, Maximal integral point sets in affine planes over finite fields, Discrete Math. 309 (2009) 4564-4575
---

# THM-4624 — finite Beckman–Quarles for the unit-circle graphs

## Setting

`q = p^e`, `F_(q^2) ⊇ F_q`, `μ_(q+1) = {z : z^(q+1) = 1}` = the kernel of the norm `F_(q^2)^× → F_q^×`, and
`G_q = Cay(F_(q^2), μ_(q+1))`: `z ~ z'` iff `z - z'` has norm one. This is decalion89's `G_q`, and
`κ(q) = χ(G_q)` bounds the chromatic number of a field plane at an inert conjugation-fixed prime with
residue field `F_q` (THM-4558, Lemma R). The *affine semilinear group*

    Aff(q) = { z ↦ a σ(z) + b : a ∈ μ_(q+1), σ ∈ Gal(F_(q^2)/F_p), b ∈ F_(q^2) }

acts faithfully on `G_q` by automorphisms; `|Aff(q)| = q^2 (q+1) · 2e`.

## Statement

1. **Census (FINITE-EXACT).** For every prime power `q ≤ 81`:

   | q | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 | 13 | 16 | 17 | 19 | 23 | 25 | 27 | 29 | 31 | 32 |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | Aut/Aff | 1 | 1 | **6** | 1 | 1 | **24** | **6** | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

   | q | 37 | 41 | 43 | 47 | 49 | 53 | 59 | 61 | 64 | 67 | 71 | 73 | 79 | 81 |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | Aut/Aff | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

   Full orders: `|Aut G_q| = q^2(q+1)·2e` except `|Aut G_4| = 1920`, `|Aut G_8| = 82944`,
   `|Aut G_9| = 19440`. Audit B extended the census: ratio `1` for every prime power `83 ≤ q ≤ 256`
   (`83, 89, ..., 241, 243, 251, 256`, including the split-norm cases `q ≡ 1 (mod 4)`), with stabiliser
   orbit counts `1 + #(Frobenius orbits on F_q^×)` throughout. In every case `G_q` is vertex-transitive and the stabiliser of `0` has
   `1 + (number of Frobenius orbits on F_q^×)` orbits in the affine cases (e.g. `q` orbits for prime
   `q`; `15, 6, 11, 28, 14, 24` for `q = 25, 16, 27, 49, 64, 81`) and `3, 4, 5` orbits for `q = 4, 8, 9`.
2. **Affine rigidity.** For every prime power `q ≤ 256` with `q ∉ {4, 8, 9}`, `Aut(G_q) = Aff(q)`: every
   automorphism of `G_q` is an affine semilinear map. (For `q = 2, 3` this holds although the graphs
   are `K_4` and `K_3 □ K_3`: `Aff(2) = S_4` and `Aff(3) = S_3 ≀ S_2`.)
3. **Subfield lines.** Let `F_(p^f) ⊆ F_(q^2)` (`f | 2e`). Then `F_(p^f)^× ⊆ μ_(q+1)` iff `p^f - 1 | q + 1`,
   and in that case `μ_(q+1) ∪ {0}` is a union of `N = (q+1)/(p^f - 1)` distinct `F_(p^f)`-lines through
   `0` in the `F_(p^f)`-space `F_(q^2)` of dimension `D = 2e/f`. This happens only for `p = 2, f ∈ {1, 2}`
   (`f = 2` needs `e` odd) and `p = 3, f = 1`; for `p ≥ 5` only the scalars `±1` lie in `μ_(q+1)`. The
   lines form a basis (`N = D`) or a frame (`N = D + 1`, any `D` independent) exactly for

   | q | 2 | 2 | 3 | 4 | 8 | 9 |
   |---|---|---|---|---|---|---|
   | (p, f) | (2,1) | (2,2) | (3,1) | (2,1) | (2,2) | (3,1) |
   | (N, D) | (3,2) frame | (1,1) basis | (2,2) basis | (5,4) frame | (3,3) basis | (5,4) frame |

   and `N ≤ D + 1` never holds otherwise.
4. **The three exceptional graphs.**
   * `G_4 ≅ Cay(F_2^4, {e_1, e_2, e_3, e_4, e_1+e_2+e_3+e_4})`, the Clebsch graph (`μ_5` sums to zero in
     `F_16`); `Aut = F_2^4 ⋊ S_5`, order `1920`; `χ = 4`.
   * `G_8 ≅ K_4 □ K_4 □ K_4 = H(3, 4)`: `μ_9 = F_4^× · {1, ζ, ζ^2}` with `{1, ζ, ζ^2}` an `F_4`-basis of
     `F_64` (`ζ` a primitive ninth root of unity); `Aut = S_4 ≀ S_3`, order `24^3 · 6 = 82944`; `χ = 4`.
   * `G_9 ≅ Cay(F_3^4, {±w_0, ..., ±w_4})` with `w_j = (-1)^j ζ^j` (`ζ` a primitive tenth root of unity),
     `w_0 + ... + w_4 = 0` and any four of the `w_j` linearly independent: the ternary five-frame
     graph; `Aut = F_3^4 ⋊ (S_5 × C_2)`, order `81 · 240 = 19440`; `χ = 3`, by the nowhere-zero
     functional `λ(w_0, ..., w_3) = (1, 1, 1, 1)`, `λ(w_4) = -1`.
   The graphs `G_2 = K_4` and `G_3 = K_3 □ K_3 = H(2, 3)` are the basis cases with no extra symmetry.
5. **Reading of `κ`.** `κ(4) = κ(8) = 4` are the chromatic numbers of the Clebsch and Hamming graphs
   (a frame over `F_2` admits no nowhere-zero functional, since `λ(w_j) = 1` for all `j` contradicts
   `Σ w_j = 0` with five terms); `κ(9) = κ(3) = 3` come from a frame functional and from `K_3 □ K_3`.
   For the undecided `κ(13), κ(17) ∈ {5, 6}` of THM-4558 the symmetry group is exactly `Aff(13)`
   (order 4732) and `Aff(17)` (order 10404). `κ(13) = 6` is decalion89's external computer proof
   (`notes/g13_chi.md`: largest-class case split with lex-leader symmetry breaking under all 4,732
   automorphisms, cube-and-conquer, kissat 4.0.4, about 470 solver-hours, DRAT-checked, proofs deleted);
   in-repo, `α(G_13) ≥ 36` (FINITE-EXACT) so no counting or spectral bound can decide it, `G_13` is
   triangle-free, and every union of at most ten norm circles is 5-colourable (session note, section 10).

## Proofs

**(1)** is the computation: `F_(q^2)` is built as `F_p[t]/(m)` with `m` irreducible of degree `2e`,
`μ_(q+1)` by exponentiation, the graph is passed to dreadnaut, and the stabiliser of `0` is obtained by
fixing the vertex in the initial partition. The orbit counts of the stabiliser match those of
`μ_(q+1) ⋊ Gal` (`μ_(q+1)` acts transitively on each norm circle `{N(z) = n}`, `n ∈ F_q^×`, and the
Frobenius permutes the circles by `n ↦ n^p`), which confirms the permutation group, not only its order.

**(2).** `Aff(q) ≤ Aut(G_q)`: `z - z' ∈ μ_(q+1)` iff `a σ(z - z') ∈ μ_(q+1)`, since `σ` is a field
automorphism (it preserves `μ_(q+1)`) and `a μ_(q+1) = μ_(q+1)`. Faithfulness: `a σ(z) + b = a' σ'(z) + b'`
for all `z` gives `b = b'` (at `z = 0`), `a = a'` (at `z = 1`), then `σ = σ'`. So `|Aff(q)| = q^2(q+1)·2e`,
and equality of orders in (1) gives equality of groups.

**(3).** `F_(p^f)^×` is the unique subgroup of order `p^f - 1` of the cyclic group `F_(q^2)^×`, and
`μ_(q+1)` the unique subgroup of order `q + 1`; containment is divisibility. Then `μ_(q+1)` is a union of
cosets of `F_(p^f)^×`, each coset `c F_(p^f)^×` being a punctured line `F_(p^f) c ∖ {0}`; distinct cosets
give distinct lines. Divisibility `p^f - 1 | p^e + 1`: write `e = sf + r` with `0 ≤ r < f`; then
`p^e ≡ p^r (mod p^f - 1)`, so `p^f - 1 | p^r + 1`, forcing `p^f - 1 ≤ p^r + 1 ≤ p^(f-1) + 1`, i.e.
`p^(f-1)(p - 1) ≤ 2`: `(p, f) ∈ {(2, 1), (2, 2), (3, 1)}`; for `(2, 2)` one needs `3 | 2^e + 1`, i.e. `e` odd;
for `(2, 1)` and `(3, 1)` the condition is automatic. For `p ≥ 5`, `μ_(q+1) ∩ F_p^× = {±1}` since
`gcd(q + 1, p - 1) = 2`. The inequality `N ≤ D + 1`: for `(2, 1)`, `2^e + 1 ≤ 2e + 1` iff `e ∈ {1, 2}`;
for `(2, 2)`, `(2^e + 1)/3 ≤ e + 1` iff `e ∈ {1, 3}`; for `(3, 1)`, `(3^e + 1)/2 ≤ 2e + 1` iff `e ∈ {1, 2}`.
The table lists `(N, D)` in each case; "frame" means `N = D + 1` with any `D` of the lines spanning,
which holds because the `N` spanning vectors `1, ζ, ..., ζ^(N-1)` satisfy one relation, the minimal
polynomial of `ζ` (irreducible, with all coefficients nonzero in the frame cases: `x^4 + x^3 + x^2 + x + 1`
over `F_2`, `x^4 - x^3 + x^2 - x + 1` over `F_3`, `x^2 + x + 1` over `F_2`).

**(4).** `G_4`: `μ_5 = {1, ζ, ..., ζ^4}` with `ζ^4 + ζ^3 + ζ^2 + ζ + 1 = 0` and `1, ζ, ζ^2, ζ^3` an
`F_2`-basis of `F_16`; in that basis the five generators are `e_1, ..., e_4` and their sum, the Clebsch
graph (folded 5-cube); `|Aut| = 1920` is classical and equals nauty's count. `G_8`: `ζ` a primitive
ninth root of unity has degree `3` over `F_4` (the order of `4` modulo `9` is `3`), `F_4^× = μ_3 ⊂ μ_9`,
and `μ_9 = ∪_(j<3) F_4^× ζ^j` with `{1, ζ, ζ^2}` an `F_4`-basis of `F_64`: the Cayley graph of `F_4^3` on
the nonzero multiples of the basis vectors is `K_4 □ K_4 □ K_4`, whose automorphism group `S_4 ≀ S_3` has
order `82944`, nauty's count. `G_9`: `ζ` a primitive tenth root of unity has degree `4` over `F_3`
(the order of `3` modulo `10` is `4`), `μ_2 = {±1} = F_3^× ⊂ μ_10`, the five lines `F_3 ζ^j` (`j < 5`) are
distinct (`ζ^(j-k) = -1 = ζ^5` is impossible for `0 ≤ j, k ≤ 4`), and `w_j := (-1)^j ζ^j` satisfy
`Σ w_j = Φ_10(ζ) = 0` with any four independent (a relation among four of `1, ζ, ..., ζ^4` would be a
polynomial of degree `≤ 4` with a missing term vanishing at `ζ`, but `Φ_10` has no missing term). The
linear maps preserving the frame `{±w_j}` are the permutations of the `w_j` (each extends uniquely to a
linear map, as any four `w_j` form a basis and the relation is symmetric) times the global sign
(signs `ε_j` must respect `Σ ε_j w_(π j) = 0`, whose only relations are multiples of `Σ w_j = 0`): a group
`S_5 × C_2` of order `240`, which with the translations gives a subgroup of order `19440 = |Aut G_9|`.
The isomorphisms and chromatic numbers were checked by explicit isomorphism tests and SAT (script).

**(5).** For the Clebsch graph a homomorphic `2`-colouring would be a functional `λ : F_2^4 → F_2`
nonzero on all five generators, contradicting `Σ w_j = 0`; `χ(Clebsch) = χ(H(3,4)) = 4` are classical and
re-checked. For `G_9` the functional `λ` is nowhere zero on the ten generators, so `z ↦ λ(z)` is a proper
3-colouring, and `χ ≥ 3` because `G_9` contains triangles: `{0, w_0, -w_0}` is one, since `2w_0 = -w_0 ∈ S`
(every Cayley graph on an elementary abelian 3-group with a symmetric generating set has such triangles).  ∎

## Reading

* **A finite Beckman–Quarles theorem, with exceptions that are geometry.** Over `R^2` every
  unit-distance automorphism is an isometry; over `F_q` the same holds (for `q ≤ 256`) unless the unit
  circle degenerates into a frame or a basis over a subfield (`q ∈ {2, 3, 4, 8, 9}`), in which case the
  frame's symmetric group acts; for `q = 2, 3` that group is already `Aff(q)`, and only for
  `q = 4, 8, 9` is it larger. Subfield lines exist in every characteristic `2` or `3` (`F_2^× = {1}`,
  `F_3^× = {±1}` always lie in `μ_(q+1)`), but a frame or basis needs `e ≤ 3`.
* **Why this matters for the field-plane programme.** THM-4558's bound `χ(G(L)) ≤ κ(q)` is the
  chromatic number of `G_q`; its symmetry group fixes the symmetry-breaking available to SAT searches
  for `κ(13)`, `κ(17)`, and the three exceptional identifications turn the table entries
  `κ(4) = κ(8) = 4`, `κ(9) = 3` into structure.
* **Not claimed.** Nothing is claimed for `q > 256`; a proof of affine rigidity in general (HYP-9256)
  would presumably follow Kiermaier–Kurz's method (directions determined by point sets; Paley
  cliques) adapted to the norm-one circle.
