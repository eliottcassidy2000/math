# The polygon ladder of the finite unit-distance graphs: Kloosterman motives, the four-bar surface, and the real plane

2026-10-17, opus-2026-10-17-S39. **PROVED + FINITE-EXACT + CITED + NUMERICAL** (typed per item below; canon:
[THM-4663](../../01-canon/theorems/THM-4663-polygon-ladder-of-the-finite-unit-distance-graphs-closed-walks-are-kloosterman-moments-pentagons-level-15-hexagons-level-6.md);
the cross-problem reading is [HYP-9286](../hypotheses/HYP-9286-local-system-dictionary-for-tournaments-unit-distances-collatz-and-keller-maps.md)
and [local_systems_dictionary_20261017.md](local_systems_dictionary_20261017.md)). Nothing here decides `χ(R^2)`, a
`κ(q)` value, or any field-plane chromatic number.

## 0. The question and what was inherited

The owner asked (2026-10-17) to keep working Hadwiger–Nelson together with the repository's machinery for computing
tournament and sequence data efficiently, and to look for more fundamental relations. The closest inherited results:

* [THM-4628](../../01-canon/theorems/THM-4628-residue-colourings-of-the-unit-distance-graph-of-q-zeta55-are-kloosterman-forced-at-every-conjugation-fixed-prime.md):
  the finite Euclidean graph `G_q = Cay(F_(q^2), μ_(q+1))` (residue colourings of field planes land in it,
  [THM-4558](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md))
  has eigenvalues `-K_q(a)`, Kloosterman sums; [THM-4552](../../01-canon/theorems/THM-4552-exceptional-knight-tori-six-is-a-two-place-product-seven-is-a-paley-torus.md):
  `G_7` is the 7 x 7 knight torus.
* The Paley pattern-sum calculus of S34–S35 ([THM-4652](../../01-canon/theorems/THM-4652-paley-pattern-sums-vanishing-and-series-parallel-polynomiality.md),
  [HYP-9278](../hypotheses/HYP-9278-rigidity-circuits-calabi-yau-cores-wheel-jacobi-ladder.md),
  [the run ladder note](paley_least_resistance_runs_one_modular_form_20261014.md), [the H(P_31) note](paley_h31_exact_rate_and_apery_wheel_20261014.md)): pattern counts in Paley tournaments
  reduce to wheels, and wheels are power sums of Jacobi sums (cost `O(q)` per prime instead of `O(q^n)`), modular at
  level 8; the first non-polynomial rung is a "core".
* The exact-integer device of the `H(P_31)` computation (S34): word-size primes and the Chinese remainder theorem.

The transfer: the unit circle plays the role of the quadratic character, Kloosterman sums the role of Jacobi sums.

## 1. Exact Kloosterman moments (FINITE-EXACT)

`M_k(p) = Σ_(a ∈ F_p^×) K(a)^k` is a rational integer (Galois permutes the `K(a)`). Choosing primes `l ≡ 1 (mod p)`
below `2^31` with a primitive `p`-th root of unity in `F_l`, the Kloosterman vector is a cyclic convolution of
length `p - 1` over `F_l`; four such primes and the Chinese remainder theorem give `M_k` exactly, the range certified
by Weil's bound. All primes `5 ≤ p ≤ 2999`, `k ≤ 8`, in 62 s
([script](../../04-computation/experiments/kloosterman_moments_exact_20261017.py),
[data](kloosterman_moments_exact_p3000_20261017.json)).

## 2. The ladder (THM-4663 statements 1–3, 6)

`W_k(q)`, the closed `k`-walks from a vertex (equilateral `k`-gons of unit steps, degenerate ones included), is
`q^(-2)[(q+1)^k + (q+1) Σ_a (-K_q(a))^k]`. Closed forms, exact on every prime `5 ≤ p ≤ 2999`:

| k | `W_k(p)` | arithmetic |
|---|---|---|
| 2 | `p + 1` | degree |
| 3 | `(p + 1)(1 - (p/3))` | triangles iff `3 | p + 1` (cube roots of unity in `μ`) |
| 4 | `3p(p + 1)` | rhombi; codegree lemma `1 - η(N(d)(N(d) - 4))` |
| 5 | `(p + 1)(p^2 + 4(1 - (p/3))p + 1 - a_p)` | weight-3 CM newform of level 15, `(η3η5)^3 + (η1η15)^3` (Peters–Top–van der Vlugt) |
| 6 | `p^4 + 11p^3 + 10p^2 + p + 1 - (p + 1) b_p` | `(η1η2η3η6)^2`, weight 4 level 6 (Hulek–Spandaw–van Geemen–van Straten) |
| 7 | via `c_7(p) = (p/105)(|λ_p|^2 - p^2)` | weight 3, level 525, quartic nebentypus (Evans; Yun) — range test only |
| 8 | `(p + 1)(p^5 + 7p^4 + 35p^3 + 7p + 1 - 7p b_p - d_p)` | `b` and the weight-6 level-6 newform (Evans; Yun) |

The level-15 coefficients were produced from the Hecke character `ψ((α)) = α^2` of `Q(sqrt(-15))` (class number 2,
`ψ(𝔭_2) = ω`), found by fitting `M_5` on the 218 inert primes and reading the residual `p^2 a_p` on the 209 split
ones; they then matched Broadhurst's eta-quotient at every prime below 3001 and satisfy the functional equation with
root number `+1`. The weight-6 newform was built inside `S_6(Γ_0(6)) = f_(4,6) · M_2(Γ_0(6))` as the
multiplicity-one `T_5`-eigenvector (`a_5 = -66`, the old eigenvalue is 6). In symmetric-power form:
`S_3 = -1 - (p/3)p^2`, `S_4 = -1 - p^2`, `S_5 = -1 - p^2 a_p`, `S_6 = -1 - p^2 b_p`, `S_8 = -1 - p^4 - p^2 d_p`.
For `q = 49, 121, 169` the same formulas hold with Frobenius-squared traces.

## 3. Pentagons and hexagons (THM-4663 statements 4, 5; PROVED given §2)

* `C_5(G_p) = p^2 (p + 1)(p^2 + 1 - (1 - (p/3))p - a_p)/10` (`p ≥ 7`): 1960 pentagons in the 7 x 7 knight torus,
  14520 in `G_11`, 40222 in `G_13`, all confirmed by search.
* `C_6(G_p) = p^2 (p^4 + p^2 + 2p - (p + 1) b_p - 2(1 - (p/3))(p + 1)(3p - 4))/12` (`p ≥ 5`). The derivation classifies
  closed 6-walks by free reduction (tree-like walks; 4-cycles with an inserted backtrack; lollipops; for
  `p ≡ 2 mod 3` doubled triangles and figure-eights, with three explicit double counts) and was confirmed by
  depth-first search at all 21 primes `5 ≤ p ≤ 83`.

## 4. The first core: the four-bar linkage (THM-4663 statements 8, 9)

Pattern counts that are series-parallel reduce to radial walk functions `R_k(N(d))` (walks of length `k` from 0 to
`d`). Testing the theta graphs at 106 primes below 600:

* `θ(2,2,2) = p^3 + 7p^2 + 4p - 2`, `θ(1,2,3) = 3(1 - (p/3))p(p + 1)`, `θ(1,3,3) = 9p^2(p + 1)`: polynomial.
* `θ(2,2,3)` (cycles 4, 5, 5) is **not** a polynomial plus `a_p, b_p` terms. Exactly,
  `θ(2,2,3) = 2W_5 + (p^2 - 1)W_3 - (p + 1)R_3(4)`, and `R_3(4)`, the configurations of the planar four-bar linkage
  with sides `1, 1, 1, 2`, is `p + 1 - a_p(30a1)`: the configuration curve has `j = 71^3/2160 = j(30a1)` and the
  traces agree.

More generally (PROVED by a two-line character sum): for `c ∉ {0, 1, 9}`, `R_3(c) = #E'_c(F_p)` with
`E'_c : y^2 = (t^2 - 4c)(t - c - 1)(t - c + 3)`. The family is semistable with fibres `I_1, I_6, I_2, I_3` at
`c = 0, 1, 9, ∞`, i.e. Beauville's elliptic modular surface of `Γ_1(6)`; the linkage's own configuration curve is
the same surface with the cusps permuted (fibres `I_2, I_3, I_1, I_6`; equal point counts at every prime
`7 ≤ p ≤ 101` and every `c`), and `R_3(c) ≡ 0 (mod 6)`. The singular lengths `0, 1, 3` are exactly the singular
points of the planar three-step walk density `p_3(x)` (Borwein–Straub–Wan–Zudilin; CITED from memory). So:

* **closed chains** (polygons) are governed by the Kloosterman motives `Sym^k Kl_2`;
* **open three-chains** with prescribed endpoints are governed by the four-bar surface, the `Γ_1(6)` family;
* a pattern count mixes them, as `θ(2,2,3)` does (level 15 through `W_5`, and the `c = 4` fibre 30a1);
* the two are tied together: since `Σ_d R_3(d)^2 = W_6`, the Frobenius traces `a_c` of the four-bar fibres satisfy
  `Σ_c a_c = -2 - 2(p/3)` and `Σ_c a_c^2 = (p + 1)(p - 4) - b_p` (THM-4663 (10); checked at every prime below 1000):
  the second moment of the four-bar surface is the hexagon form `(η1η2η3η6)^2`, the form of the Barth–Nieto
  quintic and of the Kloosterman sixth moment.

This is the HN analogue of the Paley run ladder's cores: there too, the first patterns beyond the wheels needed a new
modular input.

## 5. The real plane (THM-4663 statement 7; NUMERICAL to 25 digits, CITED proofs)

With `L_15 = L(f_(3,15), s)` and `L_6 = L(f_(4,6), s)`:

* `∫_0^∞ x J_0(x)^5 dx = 6 L_15(1)/π^2 = (30/π^4) ∫_0^∞ t I_0 K_0^4 dt`: the five-step planar unit random walk's
  density slope at the origin is the central-edge `L`-value of the same CM form that counts the finite pentagons.
* `∫ t I_0 K_0^4 = π^2 L_15(1)/5` and `∫ t I_0^2 K_0^3 = 3 L_15(2)/4` (Broadhurst eqs. 104–105, proved by Bloch–Kerr–Vanhove
  with Samart, and Bailey–Borwein–Glasser with Rogers–Wan–Zucker), `∫ t I_0 K_0^5 = π^2 L_6(2)/2`,
  `∫ t I_0^2 K_0^4 = π^2 L_6(1)/2 = 3 L_6(3)/2`.

The Fourier transform of the unit circle is `J_0` over `R` and `-K` over `F_q`; the Bessel moments are the de Rham
periods, and the Kloosterman moments the Frobenius traces, of one family of motives (Fresán–Sabbah–Yu; CITED from
memory). Hoffman's bound in `G_q` (THM-4628) and the linear-programming bounds for unit-distance-avoiding sets in
`R^2` (which use the minimum of `J_0`) are the same computation at different places; neither the moments nor the
four-bar data improve either bound.

## 6. What this changes for Hadwiger–Nelson, honestly

* Structural: the finite planes in which residue colourings live have a complete, modular census of short polygons,
  and their first non-moment pattern is a linkage curve. The six-colour programme's gadgets (the fish curve `C3`
  of THM-4662, EI21's one-parameter family) are configuration curves of larger flexible unit-distance frameworks;
  the four-bar curve is the smallest member of the same class, and the `Γ_1(6)` family is its universal form.
* Not changed: `κ(13)` and `κ(17)` (still `κ(13) ≥ 5` by the repo's SAT, `= 6` CITED from decalion89;
  `5 ≤ κ(17) ≤ 6`); the Hoffman floors; `χ(R^2)`.

## 7. Open

* The seventh rung independently (coefficients of the level-525 form), heptagon and octagon cycle counts.
* `R_4(c)`: open four-chains (five-bar linkages) are surfaces; which motive?
* The theta moments `θ(2,3,3) = Σ R_2 R_3^2`, `θ(3,3,3) = Σ R_3^3`: second and third moments of the `Γ_1(6)` family.
  The prediction (HYP-9286(b)) is that they are modular in the way the Paley wheels are moments of the Legendre family.
* Does any four-bar fibre with rational `c` carry an HN-relevant forced angle (as EI19's do)?

## Reproduction

```
python -B 04-computation/experiments/kloosterman_moments_exact_20261017.py 3000 8 05-knowledge/results/kloosterman_moments_exact_p3000_20261017.json
python -B 04-computation/experiments/unit_distance_polygon_ladder_20261017.py 05-knowledge/results/kloosterman_moments_exact_p3000_20261017.json
python -B 04-computation/experiments/unit_distance_hexagons_20261017.py 05-knowledge/results/kloosterman_moments_exact_p3000_20261017.json 83
```

Outputs: [ladder](unit_distance_polygon_ladder_20261017.out), [hexagons](unit_distance_hexagons_20261017.out).
Citations read on 2026-10-17: Broadhurst, *Feynman integrals, L-series and Kloosterman moments*, arXiv:1604.03057
(sections 5 and 7: the attributions for `k = 5, 6, 7, 8` and the Bessel identities); Yun, arXiv:1308.3920 (Compositio
151 (2015)). From memory, not re-fetched: Beauville (1982) on semistable elliptic surfaces with four singular fibres;
Borwein–Straub–Wan–Zudilin (2012) on short planar walks; Fresán–Sabbah–Yu on Kloosterman connections; Harary–Manvel
(1971) on 5-cycle counts.
