---
id: THM-4554
title: "The backward (3-adic) sieve of a minimal Collatz counterexample keeps a positive proportion: the Haar fraction of 3-adic units with no smaller predecessor has a limit in [0.28821, 0.29912] (first moment alone: >= 0.22425), answering D62; its generating function is (3/2) rho(theta)^d with rho(theta) = 3^(theta-1)/(2^theta - 1), whose zero set is that of the Moran function g, and min rho = 3^-(1-h), h = H(log_3 2): the 2 <-> 3 mirror of the forward glide law"
status: "PROVED (i), (ii), (iv), (v); FINITE-EXACT (iii), (vi); independent audit: results note mod18_mod19_seven_sixtythree_fractal_20261006.md, section 9"
session: mac-mini-2026-10-06-mod1819
source: 05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md
scripts:
  - 04-computation/experiments/mod1819_20261006_backward_sieve.py (+ .out, ALL CHECKS PASSED)
  - 04-computation/experiments/mod1819_20261006_clock_tower.py (+ .out; Mersenne and residual-seed sections)
related:
  - 05-knowledge/results/collatz_connectivity_from_rigidity_20261001.md (Proposition 3, direction D62)
  - 01-canon/theorems/THM-4504-families-of-27-moran-function-of-the-inverse-tree.md (Moran function g, glide law)
  - 05-knowledge/results/collatz_last_dip_cst_alphabet_exotic_20261005.md (the set D of targets with a smaller ancestor)
  - 05-knowledge/results/checked_switch_phase19_20261004.md (the rewrite compiler's residual seeds)
---

# THM-4554 — the backward sieve keeps a positive proportion

**Setting.** `U(x) = oddpart(3x+1)` on odd integers. From a 3-adic unit `n`, an inverse word is a composition
`w = (k_1, ..., k_d)` (`k_i >= 1`) with `x_0 = n`, `x_i = (2^(k_i) x_(i-1) - 1)/3`. It is **legal** iff each
`x_(i-1)` (`i <= d`) is prime to 3 and `k_i` is odd exactly when `x_(i-1) ≡ 2 (mod 3)`. Its predecessor `x_d` is
**smaller** (multiplier `2^K/3^d < 1`, `K = Σ k_i`) iff `2^K < 3^d`. A **first-passage** word is a smaller one all of
whose proper prefixes have `2^(K_i) > 3^i`. Let `s(r)` be the Haar fraction of units with no legal smaller word of
length `<= r` (equivalently: no legal first-passage word of length `<= r`); `s(r)` is non-increasing, and
`s_∞ = lim s(r)` is the Haar measure of the units with no smaller predecessor at all. Note 18 (Proposition 3) shows a
minimal counterexample must lie in this set; its direction D62 asked whether `s_∞ > 0`.

**(i) Legality weight (PROVED).** Every word of length `d` is legal on exactly one unit class mod `3^d`: Haar weight
`(1/2) 3^(1-d)`.

**(ii) Positivity (PROVED).** With `F_d` the number of first-passage compositions of length `d`,

    1 - s_∞ = P(n has a smaller predecessor) <= Σ_d (1/2) 3^(1-d) F_d = 0.7757408148...,

so `s_∞ >= 0.2242591852 > 0`. The sum is exact in rationals to `d = 1500`; the dropped window (excess `>= 200`)
and the depth tail are bounded by `R(θ) 2^(-θh)` with `θ = 3/2`, `R = ρ/(1-ρ)`, `ρ = ρ(3/2) = 0.9474`
(Chernoff over all continuations), contributing `< 10^-29`.

**(iii) Bracket (FINITE-EXACT + PROVED tail).** `s(r)` for `r = 1..16` is
`0.5, 1/3, 1/3, 17/54, 17/54, 0.31070, 0.30521, 0.30521, 0.30308, 0.30308, 0.30203, 0.30033, 0.30033, 0.29955,
0.29955, 0.29912`. These come from a vectorised exact recursion over `x mod 3^r`. An independent brute force agrees
for `r <= 8`, and so does the repo's own T-depth table. Hence

    0.28821 <= s(16) - Σ_(d > 16) (1/2) 3^(1-d) F_d <= s_∞ <= s(16) = 0.29912.

**(iv) Moran duality (PROVED).** The Chernoff generating function of the legal words is

    Σ_(|w| = d) (1/2) 3^(1-d) 2^(-θ (K - d log_2 3)) = (3/2) ρ(θ)^d,     ρ(θ) = 3^(θ-1)/(2^θ - 1),

and `ρ(θ) - 1 = 2^θ (g(θ) - 1)/(2^θ - 1)` with `g(s) = 2^(-s) + (1/3)(3/2)^s`, the Moran function of THM-4504
(roots `θ = 1, 2`). Its minimum is attained at `θ* = log_2(log_2 3/(log_2 3 - 1)) = 1.43803` and equals

    min ρ = 3^-(1-h),        h = H(log_3 2) = 0.949956

(Cramér at slope `log_2 3`). The forward glide census decays like `2^-(1-h)` per step (THM-4504: `W_k = Θ(2^(hk) k^(-3/2))`).
The backward first-passage mass decays like `3^-(1-h)`. The two places count the same lattice paths (binary words, odd steps as ones) with
weights `2^-K` and `3^-d`, and these weights are **equal on the critical line `2^K = 3^d`**.

**(v) Mersenne sampling (PROVED).** For odd `a`, `a mod 2·3^(k-1) -> 2^a mod 3^k` is a bijection onto the units
`≡ 2 (mod 3)` (2 is a primitive root mod `3^k`). Consequences:

* For every `k`, the set of odd `a` such that `2^a - 1` has a legal smaller inverse word of length `<= k` is
  periodic with period `2·3^(k-1)` and density exactly `1 - 2 s(k)`.
* `a ≡ 5 (mod 6)` always gives `(8n - 5)/9 < n` (`2^a - 1 ≡ 4 mod 9`).
* Even `a` gives multiples of 3, which have no ancestors at all.

**(vi) Census facts (FINITE-EXACT).**

* The first-passage masses fit `C(d) d^(-3/2) 3^-(1-h)d`, with bounded oscillating `C(d)`: slope within `10^-5` of `-(1-h) ln 3` over the nonzero terms to `d = 1500` (and `3000`).
* Among odd `a <= 120`, `2^a - 1` has no smaller ancestor (exact integers, depth 60) for 35 of 59 exponents (`0.593`; compare `2 s_∞ ∈ [0.5764, 0.5982]`).
  * The exceptions with `a ≡ 1, 3 (mod 6)` are `a = 13, 67, 69, 103`.
  * `13` and `67` lie in the class `13 mod 54 = ord_81(2)` and share the word `(2,2,1,1)`.
* All 239 residual seeds of the rewrite compiler (`checked_switch_phase19_20261004`, `n <= 10000`) have no smaller ancestor: 153 multiples of 3 plus 86 checked to depth 40. The base rate among comparable sources is `0.294`.

## Proofs

(i) Residues of `x_0, x_1, ...` mod 3 are determined successively. Given `x_(i-1) mod 3^(m+1)`, the next
`x_i mod 3^m` is an affine bijection, so each prescribed residue (fixed by the parity of `k_(i+1)`) has conditional
probability `1/3`, and `1/2` at the root among units. ∎

(ii) The union bound over first-passage words (every unit with a smaller predecessor has a shortest smaller prefix,
which is first-passage) and (i). The exact rational partial sums and the explicit Chernoff bounds for the dropped
mass are in the script. ∎

(iii) `s(16)` is computed exactly (comparisons `K - d log_2 3 >= 0` have margin `>= 0.0196` for `d <= 16`). The tail is
the union bound of (ii) restricted to `d > 16`. ∎

(iv) `Σ_K C(K-1, d-1) 2^(-θK) = (2^(-θ)/(1 - 2^(-θ)))^d`, times `3^(θd)` and `(1/2)3^(1-d)`. Then
`d/dθ log ρ = ln 3 - 2^θ ln 2/(2^θ - 1) = 0` gives `2^θ = log_2 3/(log_2 3 - 1)`, i.e. the tilted mean of a part equals
`log_2 3`. At that tilt, the exponential rate of the compositions of `K ≈ d log_2 3` into `d` parts is
`2^(log_2 3 · H(log_3 2))` per part. ∎

(v) `(Z/3^k)^*` is cyclic, generated by 2, and its index-2 subgroup is the classes `≡ 1 (mod 3)`, i.e. the even
powers. Legality of a word of length `<= k` depends only on `n mod 3^k`. For large `n` the integer condition
"smaller ancestor via a word of length `<= k`" coincides with `2^K < 3^d`, because supercritical words give smaller
ancestors only below a finite carry threshold. ∎

## Remarks

* **What it settles.** D62: the backward sieve is not a contraction to measure zero. A minimal counterexample
  cannot be excluded 3-adically from below. All exclusion power sits on the forward side (2-adic survivors: Haar
  measure 0, dimension `h`) and in integrality/size, where the two places couple only through `n < 6^L`
  determining `(n mod 2^L, n mod 3^L)`.
* **Sheet-blind.** The minus sheet (`3n - 1`) has the same backward sieve under `x -> -x`. The theorem is a local
  statement and cannot see the sheet (note 18, Theorem B).
* **Relation to the last-dip set.** `1 - s_∞ ∈ [0.70088, 0.71179]` is the 3-adic measure, among units, of the
  last-dip note's `D` (targets with a smaller ancestor).
