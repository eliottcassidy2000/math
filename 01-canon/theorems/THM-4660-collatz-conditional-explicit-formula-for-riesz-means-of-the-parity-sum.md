---
id: THM-4660
title: "Conditional explicit formula for the Collatz parity sum: if the mod-3 discrepancy H has dyadic energy E_X = O(X^(2+epsilon)) for every epsilon > 0, then F(s) = sum t(n) n^(-s) is meromorphic in Re s > 1/2 with at most simple poles, at zeros rho of D(s) = 1 + 2^(-s) + (3/2)^s/3, residues r_rho = O(|rho|) (this continuation and the identity F D = G are Codex's P4 and renewal_shapes_oct11's Lemma B1 / Theorem B3, there under stronger budgets), and the second Riesz mean satisfies sum_(n<=x) t(n)(1 - n/x)^2 = sum_(Re rho > 1/2 + delta) 2 r_rho x^rho/(rho(rho+1)(rho+2)) + O_delta(x^(1/2+delta)) with an absolutely convergent sum; unconditionally, D has uniformly simple zeros, bounded zero density and a uniform lower bound off its zeros"
status: >
  CONDITIONAL (complete proof under the stated energy hypothesis, which is weaker than every budget of HYP-9283).
  Lemmas A-C on D are PROVED unconditionally (compactness of the torus of phases; Hurwitz; Jensen). The hypothesis
  itself is OPEN (it is a square-root-cancellation statement for H; HYP-9283 calibrates even the stronger budgets as
  beyond Krasikov-Lagarias). NUMERICAL context: the S38 fits of the unsmoothed sum (HYP-9285).
  Independently audited 2026-10-16 (04-computation/experiments/s38_audit_20261016/REPORT.md): the identity F D = G
  checked exactly on cut-off sums (to 1.4e-12 at 24 points), Lemmas A, B and part 2 CONFIRMED, Lemma C and part 1
  CONFIRMED WITH FIXES (enlarged rectangle; bounded strips), credit added (E2), the residue remark corrected (E5);
  MISTAKE-623.
source: opus-2026-10-16-S38
depends_on: [THM-4658 (class-sum recursion), collatz_amplitude_boundary_20261013c.md (renewal_shapes_oct11: Lemma B1, Theorem B3)]
prior_art: [collatz_signed_resonance_pair_mass_20261013b.md (Codex P4, (9a)-(9d): the same identity F D = G with the same E and R, and the conditional meromorphic continuation with residue G(rho)/D'(rho)), collatz_amplitude_boundary_20261013c.md (renewal_shapes_oct11, Lemma B1 / Theorem B3 under a constant budget)]
note: 05-knowledge/results/collatz_explicit_formula_dense_census_20261016.md
related: [THM-4659, HYP-9283, HYP-9285, THM-4632]
---

# THM-4660: a conditional explicit formula for the Riesz means of the parity sum

## Setting

`t(n) = (−1)^τ(n)` on the basin of 1 and 0 elsewhere, `S(x) = Σ_(n<=x) t(n)`, `F(s) = Σ t(n) n^(−s)` (absolutely convergent
for `Re s > 1`), `w_m = 3[m ≡ 2 (3)] − 1`, `H(x) = Σ_(m<=x) w_m t(m)`, `F_H(s) = Σ w_m t(m) m^(−s)`, and
`D(s) = 1 + 2^(−s) + (3/2)^s/3`.

**Hypothesis (W).** For every `ε > 0`, `Σ_(X<=n<2X) H(n)² = O_ε(X^(2+ε))` over dyadic `X`.

## Statement

**Lemma A (uniformly simple zeros).** There is `c_0 > 0` with `|D'(ρ)| >= c_0` for every zero `ρ` of `D` with
`0 <= Re ρ <= 3`.

**Lemma B (lower bound off the zeros).** For every `η > 0` there is `c(η) > 0` with `|D(s)| >= c(η)` whenever
`1/4 <= Re s <= 3` and `s` has distance at least `η` from every zero of `D`.

**Lemma C (bounded density, good heights).** There is `K` such that every rectangle `[1/4, 3] × [T − 1/2, T + 3/2]`
contains at most `K` zeros of `D`. Hence for every integer `m` there is `T_m ∈ [m, m + 1]` such that the segment
`{σ + iT_m : 1/4 <= σ <= 3}` has distance at least `η_K = 1/(8K + 8)` from every zero.

**Theorem.** Assume (W). Then:
1. `F` continues meromorphically to `Re s > 1/2`, its only possible poles are simple, at zeros `ρ` of `D`, with
   residue `r_ρ = G(ρ)/D'(ρ)` for an explicit analytic `G` (below), and `|r_ρ| <= K_δ |ρ|` when `Re ρ >= 1/2 + δ`.
2. For every `δ ∈ (0, 1/4)` and `x >= 2`,

   `Σ_(n<=x) t(n) (1 − n/x)² = Σ_(D(ρ)=0, Re ρ > 1/2+δ) 2 r_ρ x^ρ / (ρ(ρ+1)(ρ+2)) + O_δ(x^(1/2+δ))`,

   where the sum over zeros converges absolutely.

## Proof

**Lemma A.** For phases `α, β` put `f_(α,β)(s) = 1 + e^(iα) 2^(−s) + e^(iβ)(3/2)^s/3`; then
`D(s + iτ) = f_(α(τ),β(τ))(s)` with `α = −τ ln 2`, `β = τ ln(3/2)` (mod 2π). A common zero of `f` and `∂_s f` would
give, with `A = e^(iα)2^(−s)`, `B = e^(iβ)(3/2)^s/3`, `L = ln 2`, `M = ln(3/2)`: `1 + A + B = 0` and `−LA + MB = 0`, so
`A = −M/(L+M)` and `B = −L/(L+M)`. Then `2^(−σ) = |A| = M/ln 3` forces `σ = log_2(ln 3/ln(3/2)) = 1.438...`, while
`(3/2)^σ/3 = |B| = L/ln 3` forces `σ = 1.573...`: impossible. So no `f_(α,β)` has a multiple zero. If zeros
`ρ_n = σ_n + iτ_n` of `D` in the strip had `D'(ρ_n) → 0`, a subsequence would have `σ_n → σ*` and
`(α(τ_n), β(τ_n)) → (α*, β*)`; uniform convergence of `f_(α,β)` and `∂_s f_(α,β)` in the parameters gives
`f_(α*,β*)(σ*) = ∂_s f_(α*,β*)(σ*) = 0`, a contradiction. ∎

**Lemma B.** Otherwise there are `s_n = σ_n + iτ_n` at distance `>= η` from all zeros with `D(s_n) → 0`. Along a
subsequence `σ_n → σ*` and the phases converge, so `g_n(s) = D(s + iτ_n)` converges uniformly on compact sets to some
`f* = f_(α*,β*)` with `f*(σ*) = 0`. `f*` is not identically zero (its value at `s = 2 + it`... more simply, `f*` has
`|f*(s)| >= 1 − f(1.5) > 0` on `Re s = 1.5`). By Hurwitz's theorem `g_n` has a zero within `η/2` of `σ*` for large `n`,
i.e. `D` has a zero within `η/2 + |σ_n − σ*|` of `s_n`: a contradiction. ∎

**Lemma C.** On `Re s = 3/2`, `|D(s)| >= 1 − (2^(−3/2) + (3/2)^(3/2)/3) = 0.034...`, and `|D| < 5` on the disks used.
Jensen's formula on a disk centred on `Re s = 3/2` bounds the number of zeros in any rectangle `[1/4, 3] × [T, T + 1]` by
`log(5/0.034)/log(2.5/1.6) < 12`, and in the enlarged rectangle `[1/4, 3] × [T − 1/2, T + 3/2]` by a constant `K` (the
audit evaluates the Jensen ratio 11.34 at radius 1.61 for the covering disks). A zero within `η_K` of a segment at
height in `[m, m + 1]` lies in the enlarged rectangle; removing the `<= K` intervals `(Im ρ − η_K, Im ρ + η_K)` from
`[m, m+1]` leaves a point `T_m`, total removed length `2Kη_K < 1`. ∎

**Theorem, part 1.** For `Re s > 1` split `n` into even and odd. Even: `Σ_m t(2m)(2m)^(−s) = −2^(−s) F(s)`. Odd:
`n = 2m' − 1 ↦ m = 3m' − 1 ≡ 2 (3)`, `t(n) = −t(m)`, `n = (2m − 1)/3`, so the odd part is
`−(3/2)^s Σ_(m≡2) t(m) m^(−s) (1 − 1/(2m))^(−s) = −(3/2)^s [F_2(s) + E(s)]` with
`F_2 = Σ_(m≡2) t(m)m^(−s) = (F + F_H)/3` and `E(s) = Σ_(m≡2) t(m) m^(−s)((1 − 1/(2m))^(−s) − 1)`. Hence

`F(s) D(s) = G(s) := −(3/2)^s (F_H(s)/3 + E(s))`.

`E` converges absolutely for `Re s > 0` with `|E(s)| <= C_1 |s|` for `Re s >= 1/4` (each term is at most
`C|s| m^(−1−σ)`). `F_H(s) = s ∫_1^∞ H(x) x^(−s−1) dx`; by Cauchy–Schwarz on dyadic blocks, (W) gives
`∫_X^(2X) |H| x^(−σ−1) dx <= (Σ H²)^(1/2) X^(−σ−1/2) = O_ε(X^(1/2+ε−σ))`, summable over dyadic `X` when
`σ >= 1/2 + δ > 1/2 + ε`. So `F_H` continues analytically to `Re s > 1/2` with `|F_H(s)| <= C_2(δ)|s|` on
`Re s >= 1/2 + δ`. Thus `G` is analytic there, with `|G(s)| <= K'_δ |s|` on every bounded strip `1/2 + δ <= Re s <= 2`
(the factor `(3/2)^σ` grows with `σ`; only `σ <= 2` is used below), `F = G/D` is meromorphic with poles only at
zeros of `D`, all simple by Lemma A, and `r_ρ = G(ρ)/D'(ρ)`, `|r_ρ| <= K'_δ |ρ|/c_0`. ∎

**Theorem, part 2.** For `c = 2` the Riesz-mean Perron formula
`Σ_(n<=x) t(n)(1 − n/x)² = (1/2πi) ∫_((c)) F(s) x^s 2 ds/(s(s+1)(s+2))` holds (absolutely convergent). Shift the
segment `[c − iT_m, c + iT_m]` to the contour consisting of the horizontal segments at heights `±T_m` (Lemma C) and the
vertical line `Re s = σ' = 1/2 + δ`, deformed to the left around every zero within `η` of that line by a circular arc
of radius `η` (with `η <= δ/2` small: zeros are separated by at least `c_0/sup|D''|`, so the arcs are disjoint), so that
those zeros lie inside. Away from the arcs the distance to the zeros is at least `min(η, η_K)`, so
`|D| >= c(min(η, η_K))` (Lemma B); on an arc around `ρ`, `|D(s)| >= c_0 η − sup|D''| η² >= c_0 η/2` (Lemma A). The
integrand is `O(|G(s)| x^σ/(|D||s|³)) = O(x^σ/|s|²)`. Hence the horizontal segments contribute `O(x^c/T_m²) → 0` as
`m → ∞`; the vertical line and the arcs (where `Re s <= σ'`) contribute `O(x^(1/2+δ)(1 + Σ_ρ η/(c_0 η |ρ|²))) =
O(x^(1/2+δ))`, using at most `K` zeros per unit height. The residues inside are `2 r_ρ x^ρ/(ρ(ρ+1)(ρ+2))`, absolutely
summable because `|r_ρ|/|ρ|³ = O(1/|ρ|²)` and the zeros have bounded density; those with `1/2 + δ − η < Re ρ <= 1/2 + δ`
contribute `O(x^(1/2+δ))` in total. Letting `m → ∞` proves the displayed formula. ∎

## Remarks

- **What is new.** The identity `F D = G` and the conditional continuation with residues `G(ρ)/D'(ρ)` are Codex's P4
  and renewal_shapes_oct11's Lemma B1 / Theorem B3 (there under a constant budget, with one residue certified). New
  here: the much weaker hypothesis (W), Lemmas A-C, and part 2, the full modal expansion of the second Riesz mean.
  The unsmoothed sum needs considerably more than `|G| = O(|s|)` (audit E4): the vertical Perron segment then
  contributes `x^(1/2+δ) T`; one also needs `∫ |G(σ' + it)| dt/|t| = O(T^ε)` on the vertical line, horizontal bounds of
  order `T^(1/2+δ)` at `T ≍ x`, and the sum `Σ r_ρ x^ρ/ρ` is not absolutely convergent if `|r_ρ| ≍ 1`. The explicit formula
  for `S(x)` itself (HYP-9285) remains a conjecture.
- **Consistency with the data.** The S38 dense fits, the S37 pre-registered out-of-sample test (and the S38 one, see
  the S38 note) and the boundary residues are the numerical shadow of this statement for `S(x)` itself (HYP-9285).
- **Lemma A** also gives `|r_ρ| = O(|ρ|)` uniformly. The boundary estimates of `|r_ρ|` over the 295 zeros below height
  20000 have maximum 11.05 (at `0.501924 + 16556.222i`) and range 0.0052-0.184 on the Mercator family (s_0 0.184,
  720.74 0.070, 5307.3 0.057, 8561.6 0.047, 2533.5 0.037, ...). These finite-N estimates do not test the `O(|ρ|)`
  bound: at `N = 2^40` under the budget `(2, 0)` only 11 of the 295 zeros have an estimate larger than the certified
  radius (median radius/estimate about 47,000 for `Re ρ <= 0.6`).
- **Constants.** `|D'(ρ)|` depends only on `Re ρ`; its minimum over `1/2 < Re ρ < 1` is `κ/2 = 0.143841` (approached as
  `Re ρ → 1`), and over `0 <= Re ρ <= 3` the infimum is `0.1308` (as `Re ρ → 2⁺`). In `1/2 < Re < 1` the zeros are at
  least 44.5 apart in height, but a rectangle of Lemma C can hold two.
