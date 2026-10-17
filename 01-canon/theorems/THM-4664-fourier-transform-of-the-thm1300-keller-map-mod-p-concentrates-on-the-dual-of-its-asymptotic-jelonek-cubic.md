---
id: THM-4664
title: "The exponential sums of THM-1300's Keller map modulo p: S(alpha, beta, gamma) = sum over F_p^3 of e_p(alpha F_1 + beta F_2 + gamma F_3) equals p e_p(-4 beta^2/alpha)(p [alpha^2 gamma = 4 beta^3] - 1) minus two explicit two-variable sums; the cone alpha^2 gamma = 4 beta^3 is the projective dual of the cuspidal cubic 27 a^2 c + b^3 = 0 at infinity of the Jelonek surface, and it carries 92-98% of the Fourier energy for 11 <= p <= 41"
status: >
  PROVED: the exact decomposition (statement 1), by THM-4621/THM-4625's parametrization of the (F_2, F_3)-fibres (a
  bijection for every b, checked exhaustively for p <= 13 by the audit) and a one-line orthogonality in b; the
  projective duality (statement 3, elementary). FINITE (float64 FFT, margins at least 0.25 p^2): statement 2 (the
  decomposition at every frequency for p = 5, 7, 11, 13, both remainders at most 2.5p; verified exactly in Z[zeta_p]
  at every alpha != 0 frequency for p = 5, 7, 11 by the audit; the large-coefficient set equal to the cone for every
  prime 11 <= p <= 41; Parseval and the fibre frequencies). Independently audited 2026-10-17 (s39_audit_20261017,
  MISTAKE-624: the Gao counts below were corrected). PROVED successor (clock_duality_oct15, scoped, self-audited):
  K1-K3 in keller_marked_cubic_fourier_20261017d.md evaluate both remainders exactly, give |S| = q^2 + O(q) on the
  cone and |S| <= 3q off it, and prove the exact cone-energy fraction for every finite field of characteristic at
  least 5. This closes the former conditional remainder target; the original numerical evidence below is retained.
  No planar Jacobian statement is claimed; HYP-9286 records the planar prediction.
session: opus-2026-10-17-S39
source: 05-knowledge/results/local_systems_dictionary_20261017.md (section 5)
scripts:
  - 04-computation/experiments/keller_thm1300_fourier_modp_20261017.py (+ 05-knowledge/results/keller_thm1300_fourier_modp_20261017.out)
related:
  - THM-1300 (the map), THM-1310 (fibre sizes 0, 1, 3 at odd p), THM-4562 (fibre counts over F_q)
  - THM-4621 (the coordinate-pair fibrations; the Jelonek set 27a^2c^2 - 18abc + 16a + b^3c - b^2 = 0), THM-4625 (the tangent-sweep family; statement (1) used here)
  - HYP-9286 (the local-system dictionary; row J), THM-4663 (the unit-distance analogue: Fourier transforms of the circle)
---

# THM-4664 — where the exponential sums of THM-1300's Keller map live

**2026-10-17 successor.** [K1-K3, the marked cubic Fourier proof](../../05-knowledge/results/keller_marked_cubic_fourier_20261017d.md)
close this theorem's former uniform remainder obligation by summing the source variable `z` first and retaining
the coefficient `(alpha t + 2 beta)/x`. For every odd finite field the transform reduces to a sum over the distinct
roots of `alpha t^3 + 3 beta t^2 = gamma`; for characteristic at least 5 the cone carries exactly the fraction
`1 - (q^2 - 2q + 3)/(q(q^2 + 2))` of the nonzero Fourier energy. Its large-coefficient set is exactly the nonzero
cone for `q >= 7`. The original proof and finite receipts below preserve the route by which this target arose.

## Setting

`F = (F_1, F_2, F_3)` is THM-1300's map, `u = 1 + xy`,
`F_1 = u^3 z + y^2 u (4 + 3xy)`, `F_2 = y + 3x u^2 z + 3x y^2 (4 + 3xy)`, `F_3 = 2x - 3x^2 y - x^3 z = xΓ`,
`Γ = 2 - 3xy - x^2 z`, `det JF = -2`. For an odd prime `p` and `(α, β, γ) ∈ F_p^3` put
`S(α, β, γ) = Σ_(v ∈ F_p^3) e_p(α F_1(v) + β F_2(v) + γ F_3(v))`, `e_p(t) = exp(2πi t/p)`. An automorphism would have
`S = 0` off the origin; here `Σ_(a ≠ 0) |S(a)|^2 = p^3 #{(v, v') : v ≠ v', F(v) = F(v')}` (Parseval), which is
`p^6 (1 + O(1/p))`: the defect `F_* 1 - 1` is the standard representation of `S_3` on the generic fibre (THM-1310).

## Statements

1. **(Exact decomposition, PROVED.)** With `p(w) = -3w^2 + 4w`, `q(w) = -w^3 + w^2` (THM-4625's data for THM-1300):
   for `α ≠ 0`,

       S(α, β, γ) = p e_p(-4β^2/α) (p [α^2 γ = 4β^3] - 1) - E + S_Γ,
       E   = Σ_(c ≠ 0) Σ_(w ∈ F_p) e_p(α q(w)/c^2 + β p(w)/c + γ c),
       S_Γ = Σ_(x ≠ 0) Σ_(y ∈ F_p) e_p(α F_1 + β F_2)  at  z = (2 - 3xy)/x^2   (the surface Γ = 0, where F_3 = 0);

   for `α = 0`: `S(0, 0, γ) = p^2` (`γ ≠ 0`) and `S(0, β, γ) = -E + S_Γ` (`β ≠ 0`).
2. **(Sizes, FINITE-EXACT.)** For `p = 5, 7, 11, 13` the decomposition holds at every frequency (floating-point
   residual below `2·10^(-13)`) with `max |E| ≤ 2.47p` and `max |S_Γ| ≤ 1.45p`. For every prime `11 ≤ p ≤ 41`, the set
   `{a ≠ 0 : |S(a)| > p^2/2}` is exactly the `p^2 - 1` nonzero points of the cone `α^2 γ = 4β^3`; on it
   `|S| ≥ 0.83 p^2` (rising to `0.95 p^2` at `p = 41`), off it `|S| ≤ 0.81 p^(3/2)` (falling to `0.45 p^(3/2)`, i.e.
   about `3p`), and the cone carries the fraction 0.925 (`p = 11`) to 0.977 (`p = 41`) of the Fourier energy.
3. **(Duality, PROVED.)** The top-degree part of the Jelonek surface `27a^2c^2 - 18abc + 16a + b^3c - b^2 = 0`
   (THM-4621) is `c (27 a^2 c + b^3)`, so its curve at infinity is the line `c = 0` and the cuspidal cubic
   `27 a^2 c + b^3 = 0`. Parametrizing the cubic by `(1 : -3t : t^3)`, its tangent line at `t` is `(2t^3 : t^2 : 1)`,
   so its projective dual is `α^2 γ = 4β^3`, a cuspidal cubic whose cusp `(0 : 0 : 1)` is the dual of the line `c = 0`.
   The exponential sums of the Keller pencil `αF_1 + βF_2 + γF_3` are of size `p^2` exactly in the directions dual to the
   asymptotic non-properness locus.

## Proof

**1.** `F_3 = xΓ`. Points with `F_3 ≠ 0` are those with `x ≠ 0, Γ ≠ 0`; by THM-4625 (1) (THM-4621 for THM-1300) they
are in bijection with the triples `(b, c, w)`, `c ≠ 0`, `p(w) ≠ cb`, via `x = 2c/(cb - p(w))`, `y = w/c - 1/x`,
`z = (2 - 3xy - c/x)/x^2` (i.e. `Γ = c/x`), with `F_2 = b`, `F_3 = c`, `F_1 = H_(b,c)(w)/c^2`,
`H_(b,c)(w) = q(w) + w(cb - p(w))/2 = w^3/2 - w^2 + (cb/2) w`. Summing first over all `b ∈ F_p`,
`Σ_b e_p(b(αw/(2c) + β)) = p [αw = -2βc]`, so the unrestricted sum is `p Σ_(c ≠ 0) e_p(α(w_c^3/2 - w_c^2)/c^2 + γc)` with
`w_c = -2βc/α`; since `α(w_c^3/2 - w_c^2)/c^2 = -4β^3 c/α^2 - 4β^2/α`, this is
`p e_p(-4β^2/α) Σ_(c ≠ 0) e_p(c(γ - 4β^3/α^2)) = p e_p(-4β^2/α)(p[α^2γ = 4β^3] - 1)`. The excluded term `b = p(w)/c`
has `H = q(w)`, which gives `E`. Points with `F_3 = 0` are the plane `x = 0`, on which `F = (z + 4y^2, y, 0)` and the sum
over `z` vanishes for `α ≠ 0` (for `α = 0` it is `p^2 [β = 0]`), and the surface `Γ = 0, x ≠ 0`, which gives `S_Γ`. For
`α = 0, β ≠ 0` the `b`-sum vanishes identically. For `(0, 0, γ)`: on `x ≠ 0` the phase `γ x Γ` is linear in `z` with
coefficient `-γx^3 ≠ 0`, so only `x = 0` contributes, `p^2`.

**2.** Direct computation (script): the fibre-count function `N` by evaluating `F` on `F_p^3`, `S` by a
three-dimensional FFT of `N`, `E` and `S_Γ` by direct summation.

**3.** As stated; `(1, -3t, t^3) × (0, -3, 3t^2) = (-6t^3, -3t^2, -3)`.

## Reading

The exponential sum of a polynomial without critical points is controlled by its behaviour at infinity (Suzuki's
formula `χ(generic fibre) = 1 - μ - λ` with `μ = 0` for every member of a Keller pencil; CITED): it can only be large
in the directions where the level sets meet the non-properness locus at infinity non-transversally. THM-4664 makes this
exact for the dimension-3 counterexample: the defect's spectrum is the dual of the Jelonek surface's curve at infinity.
For a planar Keller pair the Jelonek set is a curve with finitely many points at infinity, so the analogous prediction
(HYP-9286, row J) is that a planar counterexample would have its exponential sums `Σ_(x ∈ F_p^2) e_p(aP + bQ)` large only
on finitely many lines of frequencies; by the finite-field criterion of local_systems_dictionary_20261017.md (section 5),
the planar conjecture is equivalent to these sums vanishing identically.

## Closed remainder target and remaining open work

The uniform bounds for `E` and `S_Γ` are closed by the successor above (clock_duality_oct15): for `α != 0`,
`|E| <= 3q + sqrt(q)` and `|S_Γ| <= q + sqrt(q)`, with their Gauss terms cancelling in `-E + S_Γ`. Open: the phases
`e_p(-4β^2/α)` on the cone as a function on the dual cubic. For the general sweep family the same orthogonality gives
the main term `p Σ_(c ≠ 0) e_p(γc - (α/2) P(-2βc/α)/c^2)`, `P = ∫_0^w p`, which is linear in `c` only for `d = 2`; for
`d ≥ 3` it is `p` times a Weil sum of a degree-`(d - 1)` polynomial in `c`. For Gao's degree-four map (`det = -2`) the
audit (a8b, exact substitution) finds exactly `2(p - 1)` coefficients with `|S| > p^2/2` (20, 24, 32 at
`p = 11, 13, 17`): the lines `(0, 0, γ)` (dual of `c = 0`) and `(α, 0, 0)`, as the main term predicts (`β = 0` forces
`w_c = 0`); a session scratch run had reported 340, 600, 384 (wrong, MISTAKE-624). Whether this is the dual of Gao's
asymptotic Jelonek locus, the generic coefficient distribution for `d ≥ 3`, and the planar prediction of HYP-9286
remain OPEN.
