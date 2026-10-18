---
id: THM-4668
title: "The degree-192 golden direction e = Zs/h is Cartesian on the golden linear contact space: in K + uK + vK + eK (K = Q(zeta15, sqrt(-11)), clock_duality_oct15's golden contact u, v and escape direction e), every unit vector with a nonzero e-coefficient is a pure e-axis unit q e, |q| = 1; hence every unit-distance graph on A + uC + vC + eL (A, C, L in K) is the old source graph Cartesian-multiplied along e, every boundary colouring of the old source extends, and no unit edge joins different port and e-coefficient fibres"
status: >
  PROVED. The proof is a case analysis on the norm equation in E' = F(s, h), together with a degree argument (N is not
  in K) and one exact computation: the minimal polynomial u^3 = m2 u^2 + m1 u + 1 of u over Q(j), and the nonvanishing
  of an explicit element of Q(zeta5, j), done in sympy with exact reduction to the basis zeta5^i j^k.
  FINITE-EXACT control: an independent multi-prime sieve plus exact norm census of all 9.1 x 10^8 differences with
  coefficients in D(A), D(C), D(C) u D(L) of clock_duality's golden source finds exactly the 26 pure e-axis units.
  The direction e and the field E' are clock_duality_oct15's (golden_square_contexts_20261017h.md, G2); this theorem
  closes their stated next target inside the linear space.
session: opus-2026-10-17-S41
source: 05-knowledge/results/golden_additive_boundary_20261017.md
scripts:
  - 04-computation/experiments/golden_e_direction_proof_check_20261017.py (+ .out): minimal polynomial and exact nonvanishing
  - 04-computation/experiments/golden_e_direction_census_20261017.py (+ .out, .json): finite census control
related: [clock_duality_oct15 golden_square_contexts_20261017h.md (G1, G2, attachment P + eL), joint_contact_field_20261017f.md (G5, section 5), shared_golden_contact_20261017g.md, HYP-9289]
---

# THM-4668 — the golden direction e is Cartesian on the golden linear contact space

## Setting (clock_duality_oct15)

- **Base fields.** `K = Q(zeta15, j)` with `j^2 = -11`, and `F = K(w)` with `w^3 + 2w^2 - 4w - 4 = 0`, `1 < w < 2`;
  `[F:Q] = 48`, so `[F:K] = 3`.
- **Contact.** `d = 1 + zeta5`, `u = (w + 3 - j)/(w + 3 + j)` (so `conj(u) = 1/u`), `Z = d + u`,
  `N = Z conj(Z) ∈ (37/15, 59/20)`, `D = 1 - 4/N < 0`, `R = 4 - N > 0`.
- **Escape field.** `s^2 = D`, `h^2 = R`, `E' = F(s, h)` of degree 192. Complex conjugation fixes `h` and sends `s` to
  `-s`.
- **Directions.** `v = Z(-1 + s)/2` and `e = Z s/h = (Z/R) sh`; `|v| = |e| = |d + u + v| = 1`.

## Statement

Let `x = a + ub + vc + eg` with `a, b, c, g ∈ K` and `g ≠ 0`. If `x conj(x) = 1`, then `a = b = c = 0` and
`g conj(g) = 1`.

**Corollary.** Let `A, C, L ⊂ K` be finite, with distinct addresses. The unit graph on `A + uC + vC + eL` has exactly:
- the edges of `A + uC + vC` in each `e`-layer (classified by G5: axes and the contact family `q(d, 1, 1)`);
- the unit edges of `L` along `e`, with the other coordinates fixed.

Consequently:
- every proper `k`-colouring `c` of the old source with `k >= χ(L)` extends by `c + c_L (mod k)` (additive), so the
  boundary relation on the old ports is unchanged;
- no unit edge joins different port fibres and different `e`-coefficients.

## Proof

**The norm equation.** Write elements of `E'` as `c0 + c1 s + c2 h + c3 sh` over `F`. Then `x` has
`c0 = a + ub - Zc/2`, `c1 = Zc/2`, `c2 = 0`, `c3 = (Z/R) g`. Since conjugation sends `(c0, c1, c2, c3)` to
`(c̄0, -c̄1, c̄2, -c̄3)`, the equation `x conj(x) = 1` is equivalent to:

    (i)   c0 c̄0 - D c1 c̄1 - DR c3 c̄3 = 1,     (ii)  c1 c̄0 = c0 c̄1,
    (iii) c1 c̄3 + c3 c̄1 = 0,                    (iv)  c3 c̄0 = c0 c̄3.

**N is not in K.** `N = |d|^2 + 1 + d/u + d̄u`. If `N ∈ K`, then `d̄ u^2 + (|d|^2 + 1 - N) u + d = 0` would be a
quadratic over `K` with leading coefficient `d̄ ≠ 0`, contradicting `[K(u):K] = [F:K] = 3`. So `R = 4 - N ∉ K`.

**Case `c ≠ 0`.**
- (iii) reads `(N/(2R))(c ḡ + g c̄) = 0`, so `c ḡ` is purely imaginary.
- If `c0 ≠ 0`, then (ii) and (iv) make `c1 c̄0` and `c3 c̄0` real. Hence `c1 c̄3 = (c1 c̄0)(c3 c̄0)/|c0|^2` is real, so by
  (iii) it is 0. That contradicts `c1, c3 ≠ 0`. So `c0 = 0`.
- Using `-D = R/N`, (i) becomes `R|c|^2/4 + |g|^2 = 1`. With `|c|^2, |g|^2 ∈ K` and `c ≠ 0`, this would put `R` in
  `K`, which is impossible.

**Case `c = 0`.**
- (iv) says `c0 Z̄ ḡ` is real. Hence `c0 = a + ub = t Z g` with `t ∈ F` real.
- (i) becomes `|g|^2 (1 + t^2 N) = 1`. So `t^2 N = κ := 1/|g|^2 - 1 ∈ K`, i.e. `c0^2 N = κ Z^2 g^2`.
- Cancelling `Z` and using `conj(u) = 1/u` gives `(a + ub)^2 (1 + d̄ u) = κ' u (d + u)` with `κ' = κ g^2 ∈ K`.
- Reduce modulo the exact minimal polynomial `u^3 = m2 u^2 + m1 u + 1`, where `m2 = -(20 + 11j)/19` and
  `m1 = (20 - 11j)/19` (computed, and checked numerically at the intended embedding to `2·10^(-16)`). Comparing
  coefficients of `1, u, u^2`:

      a^2 + d̄ b^2 = 0,     2ab + d̄ a^2 + d̄ m1 b^2 = κ' d,     b^2 + 2 d̄ ab + d̄ m2 b^2 = κ'.

- If `b = 0`, then `a = 0`, `κ' = 0`, `t = 0` and `|g| = 1`.
- If `b ≠ 0`, put `a = λ b`. The last two equations give `λ = (d̄^2 - d̄ m1 + d + |d|^2 m2)/(2 - 2|d|^2)`, where
  `|d|^2 = φ^2 ≠ 1`. The first requires `λ^2 + d̄ = 0`. But `λ^2 + d̄` is an explicit nonzero element of
  `Q(zeta5, j)`: its numerator in the basis `zeta5^i j^k` is
  `-44 j z^3 + 396 j z^2 + 880 j z + 902 j - 3770 z^3 - 5788 z^2 - 5195 z - 893` with `z = zeta5`
  (golden_e_direction_proof_check_20261017.out). So there is no solution.

**Corollary.** A difference of two addresses lies in `K + uK + vK + eK` with `e`-coefficient `l - l'`. By the
statement it is a unit only if `l ≠ l'` forces the other coordinates to agree and `|l - l'| = 1`, or if `l = l'` and
the difference is a unit of `K + uK + vK` (G5). The additive extension is the usual Cartesian one. ∎

## Reading

- clock_duality's degree-192 escape made `σ`-mates unit-distance pairs (their G1/G2). Their first attachment
  `G ∪ (P + eL)` was Cartesian (their formula (3)). They named "a unit edge joining different port and `e`-coefficient
  fibres" as the next target.
- This theorem shows that such an edge does not exist anywhere in `K + uK + vK + eK`. A forcing construction through
  `e` must use coefficients outside `K`, or products such as `ue` and `ve`. Any such coordinate should first be tested
  with the census method (golden_vsquare_unit_census_20261017.py, golden_e_direction_census_20261017.py).
- Together with HYP-9289 (additive freedom of `K + uK + vK`, OPEN), this confines six-colour forcing in the golden
  contact field to genuinely nonlinear coordinates.
