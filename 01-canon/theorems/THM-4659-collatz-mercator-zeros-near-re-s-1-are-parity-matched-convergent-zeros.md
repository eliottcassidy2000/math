---
id: THM-4659
title: "The zeros of the Collatz parity characteristic function D(s) = 1 + 2^(-s) + (3/2)^s/3 near Re s = 1 are exactly the parity-matched near-returns of log2 3: for integers p even, q odd with eps = q log2 3 - p small, D has a unique zero rho = 1 + i pi q/ln 2 + w(eps) with w analytic and independent of q, 1 - Re rho = (pi^2 ln^2 2/ln^3(4/3)) eps^2 + O(eps^4), Im rho = pi q/ln 2 + pi eps/ln(4/3) + O(eps^3), and every zero with Re rho close to 1 arises this way; the members have frequencies within O(eps) of the half-integers q/2 cycles per octave, so they alias pairwise on every log-uniform grid with an integer number of points per octave"
status: >
  PROVED (holomorphic implicit function theorem; the key point is that D(1 + i pi q/ln 2 + w) depends on (p, q) only
  through eps). NUMERICAL controls: all 295 zeros of D with 1/2 < Re < 1 and 0 < Im < 20000 (argument principle in 40
  rectangles); (1 - Re rho)/eps^2 = 197.3 at 84/53 and 199.17 at 1054/665 against the limit 199.16; the odd multiple
  (3162, 1995) gives 1 - Re rho = 7.1e-6 = 9 x 7.9e-7.
source: opus-2026-10-16-S38
depends_on: [THM-4658 (the transfer operator and D), collatz_amplitude_boundary_20261013c.md (D as the Mercator characteristic function)]
verification:
  - 04-computation/experiments/collatz_D_zeros_and_boundary_residues_20261016.py (zeros to height 20000, argument principle)
note: 05-knowledge/results/collatz_explicit_formula_dense_census_20261016.md
related: [THM-4629, THM-4632, HYP-9260, HYP-9285, THM-4658]
---

# THM-4659: the zeros of D near Re s = 1 are the parity-matched near-returns of log2 3

## Setting

`D(s) = 1 + 2^(−s) + (3/2)^s/3` is the characteristic function of the class-sum recursion for the parity
`t(n) = (−1)^τ(n)` (THM-4658 statement 4(a): Haar measure is an eigen-distribution with eigenvalue 1 exactly at its
zeros). Put `L = ln 2`, `M = ln(3/2)`, `κ = L − M = ln(4/3)`, `θ = log_2 3` and

`C = π² L²/κ³ = π² ln² 2/ln³(4/3) = 199.16...`

## Statement

1. **No zeros on or right of the line up to Re s = 2.** `D(s) ≠ 0` for `1 <= Re s <= 2`.
2. **The family.** There are `ε_0 > 0`, `r_0 > 0` and an analytic function `w` on `|ε| < ε_0` with
   `Re w(ε) = −Cε² + O(ε⁴)` and `Im w(ε) = πε/κ + O(ε³)`, such that for all integers `p` even, `q` odd with
   `ε = qθ − p`, `|ε| < ε_0`, the point `ρ_(p,q) = 1 + iπq/L + w(ε)` is a zero of `D`, and it is the only zero in the
   disk `|s − (1 + iπq/L)| < r_0`. In particular `1 − Re ρ_(p,q) = Cε² + O(ε⁴)`.
3. **Converse.** There is `x_0 > 0` such that every zero `ρ` of `D` with `1 − x_0 < Re ρ < 1` equals `ρ_(p,q)` for a
   unique pair `p` even, `q` odd with `|qθ − p| = O(sqrt(1 − Re ρ))`.
4. **Zero-free region.** If `|θ − p/q| >= c q^(−μ)` for all `q >= 1` (an irrationality exponent `μ` of `θ` with
   constant `c`), then every zero with `Re ρ < 1` and large `|Im ρ|` has `1 − Re ρ >= c' |Im ρ|^(−2(μ−1))`.
5. **Aliasing on integer grids.** For a grid `N_j = 2^(j/P)` with `P` a positive integer, the modes of `ρ_(p,q)` and
   `ρ_(p',q')` with `q + q' = 2P k` (k integer) satisfy `N_j^(i Im ρ_(p',q')) = conj(N_j^(i Im ρ_(p,q))) · e^(i j O(ε + ε')/P)`: up to
   the detunings, the two modes are indistinguishable in least-squares fits of `Re Σ c_ρ N^ρ` on that grid.

## Proof

**1.** `f(σ) = 2^(−σ) + (3/2)^σ/3` is strictly convex with `f(1) = f(2) = 1`, so `f(σ) < 1` and
`|D(s)| >= 1 − f(σ) > 0` for `1 < σ < 2`. At `σ = 1` and `σ = 2` a zero needs equality in the triangle inequality, i.e.
both terms equal to negative reals. At `σ = 1`: `D(1 + it) = 1 + (a/2)(1 + b)` with `a = 2^(−it)`, `b = 3^(it)`
vanishes iff `a = −1` and `b = 1`. At `σ = 2`: it needs `2^(−it) = −1` and `(3/2)^(it) = −1`, hence again `b = 1`. In both
cases `t ln 2 ∈ π(2Z + 1)` and `t ln 3 ∈ 2πZ`, which would make `θ` rational. ∎

**2.** Let `t_0 = πq/L`, so `2^(−it_0) = −1`, and `3^(it_0) = e^(iπ(p+ε)) = e^(iπε)` because `p` is even. Then for every
complex `w`,

`D(1 + it_0 + w) = 1 − (1/2) 2^(−w) − (1/2) e^(iπε) (3/2)^w =: Φ(w, ε)`,

which depends on `(p, q)` only through `ε`. `Φ` is entire in `(w, ε)`, `Φ(0, 0) = 0` and `∂_w Φ(0, 0) = (L − M)/2 = κ/2 ≠ 0`.
By the holomorphic implicit function theorem there are `ε_0, r_0 > 0` and an analytic `w(ε)`, `w(0) = 0`, such that
for `|ε| < ε_0` the zeros of `Φ(·, ε)` in `|w| < r_0` are exactly `w(ε)`. Expanding
`Φ = (κw − iπε)/2 − (L²w² + (iπε + Mw)²)/4 + O(3)` gives `w_1 = iπε/κ`; then `iπε + M w_1 = iπεL/κ`, and the second
order equation `(κ/2) w_2 = (L² w_1² + (iπε + M w_1)²)/4` gives `w_2 = −π² L² ε²/κ³ = −Cε²`, which is real. Since `conj(Φ(conj w, ε)) = Φ(w, −ε)`, uniqueness gives
`w(−ε) = conj(w(ε))`: `Re w` is even and `Im w` is odd in `ε`, so `Re w = −Cε² + O(ε⁴)` and `Im w = πε/κ + O(ε³)`. The neighbourhood sizes do not depend on `(p, q)`
because `Φ` does not. ∎

**3.** Let `D(ρ) = 0`, `ρ = σ + it`, `x = 1 − σ` small. With `A = 2^(−ρ)`, `B = (3/2)^ρ/3` we have `A + B = −1`, and
`|A| + |B| = 2^(−σ) + (3/2)^σ/3 = 1 + (κ/2)x + O(x²)` (the derivative at `σ = 1` is `−κ/2`). Hence
`2|A||B|(1 − cos ∠(A, B)) = (|A| + |B|)² − 1 = O(x)`, so `A` and `B` are within angle `O(sqrt x)` of the negative real
axis: `2^(−it) = −e^(iα)` and `3^(it) = e^(iβ)` with `α, β = O(sqrt x)`. Choose `q` odd with `tL = πq + O(sqrt x)` and
`p` even with `t ln 3 = πp + O(sqrt x)`. Then `qθ − p = (tL/π)θ − t ln 3/π + O(sqrt x) = O(sqrt x)`, and
`|ρ − (1 + iπq/L)| = O(sqrt x)`. For `x` small these are inside the neighbourhoods of statement 2, whose uniqueness
gives `ρ = ρ_(p,q)`. ∎

**4.** By statements 2 and 3, `1 − Re ρ >= (C/2) ε²` for zeros close to the line, with `|ε| = |qθ − p| >= c q^(1−μ)`
and `q = |Im ρ| L/π + O(1)`. Zeros with `1 − Re ρ >= x_0` satisfy the bound trivially for large `|Im ρ|`. ∎

**5.** `N_j^(iπq'/L) = 2^(iπq'j/(LP)) = e^(iπq'j/P)`; with `q' = 2Pk − q` this is `e^(2πikj) e^(−iπqj/P)`, the conjugate
of `N_j^(iπq/L)`. The zeros differ from `1 + iπq/L` by `w(ε)`, giving the detuning factor. ∎

## Remarks

- **Which pairs occur.** The parity constraint (`p` even, `q` odd) is THM-4629's selection rule for the parity
  character `(−1)^τ`: that character hears only the near-returns `2^p ≈ 3^q` with this parity. Among the convergents of
  `log_2 3` the matched ones are `2/1`, `8/5`, `84/53`, `1054/665`, `50508/31867`, ... The pair `84/53` gives the
  Mercator zero `s_0 = 0.99820922 + 240.248i` of the amplitude note; `1054/665` (`ε = 6.30e-5`) gives
  `ρ* = 0.99999921 + 3014.0202i`; `50508/31867` (`ε = 1.05e-5`) predicts a zero with `1 − Re ρ ≈ 2.2e-8` near height
  144,432.7 (numerically `0.999999978121 + 144432.720754i`, by Newton on `Φ`). The pair `8/5` (`ε = −0.075`) lies
  outside the perturbative range, but Newton on `Φ` from the perturbative guess still lands on the zero
  `0.693918 + 22.2532i` (the analytic branch continues there; `1 − Re = 0.306` against `Cε² = 1.13`): the source of the
  2.45-cycle line of `S(N)/N` is the `8/5` member of the family. Odd multiples `(mp, mq)` give `ε' = mε` and
  the zeros `ρ_(3p,3q)` (`9042.06`, `1 − Re = 7.1e-6`) and `ρ_(5p,5q)` (`15070.10`).
- **Relation to THM-4632.** THM-4632 measures the same phenomenon on the line: the transfer function peaks at the
  convergents with height proportional to `θ_f^(−2)`, and its boundary growth exponent is `2(μ − 1)`. Statement 4 is
  the zero-side version; the constant `C` makes the peak law explicit for the parity character.
- **Consequence for the amplitude programme.** Under any energy budget of HYP-9283, `F(s) = Σ t(n) n^(−s)` is
  meromorphic in `Re s > 1/2` with poles only at zeros of `D`. A certified nonzero residue at `ρ_(p,q)` would give
  `limsup log B(X)/log X >= Re ρ_(p,q)`. The boundary estimate of the residue at `ρ*` from the exact block moments
  through 2^40 is `0.00993 + 0.00113i`, but the conditional Arb radius under `(2, 0)` at `N = 2^40` is 0.026, so it is
  not yet certified (see the S38 note).
- **Classical context (CITED, no novelty claim).** Locating zeros of exponential polynomials near a boundary line by
  Diophantine approximation of the frequency ratios is classical in the theory of almost periodic functions and
  exponential sums (Bohr; Langer 1931). The uniformity in `q` here comes from the exact identity for `Φ`.
