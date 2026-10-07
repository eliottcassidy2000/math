---
id: THM-4554
title: "The backward (3-adic) sieve of a minimal Collatz counterexample keeps a positive proportion: the Haar fraction of 3-adic units with no smaller predecessor has a limit in [0.28820, 0.29912] (first moment alone: >= 0.22425), answering D62; its generating function is (3/2) rho(theta)^d with rho(theta) = 3^(theta-1)/(2^theta - 1), whose zero set is that of the Moran function g, and min rho = 3^-(1-h), h = H(log_3 2): the 2 <-> 3 mirror of the forward glide law"
status: "PROVED (i), (ii), (iv), (v); FINITE-EXACT (iii), (vi) (the census law in (vi) is a numerical fit); INDEPENDENTLY AUDITED (2026-10-06; corrections in MISTAKE-569; results note mod18_mod19_seven_sixtythree_fractal_20261006.md, section 9)"
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
minimal counterexample either lies in this set or is divisible by 3; its direction D62 asked whether `s_∞ > 0`.

**(i) Legality weight (PROVED).** Every word of length `d` is legal on exactly one unit class mod `3^d`: Haar weight
`(1/2) 3^(1-d)`.

**(ii) Positivity (PROVED).** With `F_d` the number of first-passage compositions of length `d`,

    1 - s_∞ = P(n has a smaller predecessor) <= Σ_d (1/2) 3^(1-d) F_d = 0.7757408148...,

so `s_∞ >= 0.2242591852 > 0`. The sum is exact in rationals to `d = 1500` (an untruncated recomputation gives
`0.77574081479760749602...`).

* The dropped window is bounded by `Σ_d N_(d-1) (1/2) 3^(1-d) R(θ) 2^(-θ e_d)/(1 - 2^(-θ)) = 4.1·10^-30`. Here `N_(d-1)` counts the in-window prefixes, `e_d` is the excess at the window edge, `θ = 3/2`, `R = ρ/(1-ρ)` and `ρ = ρ(3/2) = 0.94729` (Chernoff over all continuations).
* The depth tail is at most `(3/2) ρ^1501/(1-ρ) = 1.4·10^-34`.

**(iii) Bracket (FINITE-EXACT + PROVED tail).** `s(r)` for `r = 1..16` is
`1/2, 1/3, 1/3, 17/54, 17/54, 151/486, 445/1458, 445/1458, 3977/13122, 3977/13122, 35669/118098, 106405/354294,
106405/354294, 955147/3188646, 955147/3188646, 4292002/14348907` (the last `= 0.299117`).

These come from a vectorised exact recursion over `x mod 3^r`; the exact rationals are from the audit's independent class marking. Brute force agrees for `r <= 8`, and a per-class search agrees for `r <= 10`.

The repo's own table (`collatz_connectivity_from_rigidity_20261001.out`, section D) is the same sequence on a different clock. It uses standard-map depth `D = K + d`, so its value at `D` is `s(max{d : d + floor(d log_2 3) <= D})`, because every first-passage word of length `d >= 2` has `K = floor(d log_2 3)`. Hence

    0.28820 <= s(16) - Σ_(d > 16) (1/2) 3^(1-d) F_d <= s_∞ <= s(16) = 0.29912

(the lower value is `0.2882095614`).

**(iv) Moran duality (PROVED).** The Chernoff generating function of the legal words is

    Σ_(|w| = d) (1/2) 3^(1-d) 2^(-θ (K - d log_2 3)) = (3/2) ρ(θ)^d,     ρ(θ) = 3^(θ-1)/(2^θ - 1),

and `ρ(θ) - 1 = 2^θ (g(θ) - 1)/(2^θ - 1)` with `g(s) = 2^(-s) + (1/3)(3/2)^s`, the Moran function of THM-4504
(roots `θ = 1, 2`). Its minimum is attained at `θ* = log_2(log_2 3/(log_2 3 - 1)) = 1.43803` and equals

    min ρ = 3^-(1-h),        h = H(log_3 2) = 0.949956

(Cramér at slope `log_2 3`).

* The forward glide census decays like `2^-(1-h)` per step (THM-4495: `W_k = Θ(2^(hk) k^(-3/2))`; THM-4504 gives `log_2 W_m / m -> h`).
* The backward first-passage mass satisfies the proved bound `m_d <= (3/2) 3^-(1-h)d`.
* The two are **mirror families on opposite sides of one critical line.** Glide words keep their density of ones above `log_3 2` on every prefix. Backward first-passage words keep it below until the last step. The shared exponent is the Cramér rate at the line, where the 2-adic weight `2^-K` equals the 3-adic weight `3^-d`.

**(v) Mersenne sampling (PROVED).** For odd `a`, `a mod 2·3^(k-1) -> 2^a mod 3^k` is a bijection onto the units
`≡ 2 (mod 3)` (2 is a primitive root mod `3^k`). Consequences:

* For every `k`, the set of odd `a` such that `2^a - 1` has a legal smaller inverse word of length `<= k` is
  periodic with period `2·3^(k-1)` and density exactly `1 - 2 s(k)`, which tends to `1 - 2 s_∞`. A density statement at infinite depth would need an exchange of limits and is not claimed.
* `a ≡ 5 (mod 6)` always gives `(8n - 5)/9 < n` (`2^a - 1 ≡ 4 mod 9`).
* Even `a` gives multiples of 3, which have no ancestors at all.

**(vi) Census facts (FINITE-EXACT; the first bullet is a numerical fit).**

* The first-passage masses fit `C(d) d^(-3/2) 3^-(1-h)d`, with `C(d)` staying in `[0.22, 0.62]` (so far as computed): slope within `10^-5` of `-(1-h) ln 3` over the nonzero terms to `d = 1500` (and `3000`). A free fit gives the exponent `-1.52`. Only the upper bound in (iv) is proved.
* Among odd `3 <= a <= 119`, `2^a - 1` has no smaller ancestor (exact integers, depth 60) for 35 of 59 exponents (`0.593`; compare `2 s_∞ ∈ [0.5764, 0.5983]`).
  * The exceptions with `a ≡ 1, 3 (mod 6)` are `a = 13, 67, 69, 103`.
  * `13` and `67` lie in the class `13 mod 54` (`54 = ord_81(2)`) and share the word `(2,2,1,1)`.
* All 239 residual seeds of the rewrite compiler (`checked_switch_phase19_20261004`, `n <= 10000`) have no smaller ancestor: 153 multiples of 3 plus 86 checked to depth 40 (and 60). The base rate among comparable sources is `245/834 = 0.294`. The residuals are a proper subset of the backward-minimal sources: 153 of the 416 comparable multiples of 3, and 86 of the 245 comparable others.

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
"smaller ancestor via a word of length `<= k`" coincides with `2^K < 3^d`, because a supercritical word `w` (`2^K > 3^d`) gives
`x_d < n` only for `n < B_w/(2^K - 3^d)`. Here `x_d = (2^K n - B_w)/3^d`, and the threshold is the rational fixed point of the
word's affine map. Over words of length `<= 12` it is at most `381.3`, attained at `(d, K) = (10, 16)`. A direct search finds
no odd unit `n < 3000` with a smaller ancestor through a supercritical word of length `<= 12`. ∎

## Remarks

* **What it settles.** D62: the backward sieve is not a contraction to measure zero. No argument using only the
  3-adic backward sieve can exclude a minimal counterexample. Any exclusion must use more: the forward (2-adic)
  side, whose survivors have Haar measure 0 and dimension `h`, or integrality and size, which couple the two
  places (for `n < 6^L`, `n` is determined by `(n mod 2^L, n mod 3^L)`).
* **Sheet-blind.** The minus sheet (`3n - 1`) has the same backward sieve under `x -> -x`. The theorem is a local
  statement and cannot see the sheet (note 18, Theorem B).
* **Relation to the last-dip set.**
  * The last-dip note's `D` (targets with a smaller ancestor) is a set of integers.
  * `1 - s_∞ ∈ [0.70088, 0.71180]` is the Haar measure, among units, of its 3-adic analogue.
  * The observed proportion among units below `2^28` is `0.703004`, which is consistent, but equality of densities is not proved.
