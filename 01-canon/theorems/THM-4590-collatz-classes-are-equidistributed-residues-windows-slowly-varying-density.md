---
id: THM-4590
title: "Collatz classes are equidistributed: every union B of grand-orbit classes of the Terras map on the positive (or negative) integers has o(x) boundary points n (n in B, n+1 not, or conversely) in [x, 2x); B is asymptotically equidistributed among residue classes mod every M and uniform across multiplicative windows [cx, c(1+eta)x); its dyadic density varies slowly (beta(cx) - beta(x) -> 0 for each fixed c). On the negative integers the basins of the cycles -1, -5, -17 have densities 0.3268, 0.3248, 0.3484 at 2^28, residue deviations <= 0.004 up to mod 1024, and a slowly falling cut fraction 0.41 -> 0.34"
status: >
  PROVED from THM-4581 (3), almost-sure coalescence: statements 1-4 (elementary, given density-one coalescence of n and n+1).
  FINITE-EXACT / NUMERICAL: negative-integer basins and positive-integer entry classes to 2^28 (scripts below). The rate is not proved:
  the cut density is at most q(c log x), i.e. O((log x)^(-1/2) (log log x)^2) at THM-4581 (4)'s sketch level, and the approach
  to multiplicative uniformity has a Diophantine resonance spectrum (HYP-9230). The theorem says nothing about whether a non-1 class exists;
  if the Collatz conjecture holds, it is vacuous on Z_{>0} except for "end-saturated" partitions such as entry classes (statement 5).
session: mac-mini-2026-10-07-golden
source: 05-knowledge/results/golden_collatz_resonance_20261007.md
scripts:
  - 04-computation/experiments/golden_20261007_lead/ (negbasins.c, negwindows.c, negfourier.c, negfourier70.c, posentry85.c, posentry85_70.c, odd-step variance; + .out)
related:
  - THM-4581 (Haar coalescence; (6c) density-one coalescence of n and n+1 is the only input)
  - THM-4569 (7) (index [R_A : R_C] = 1; statement 3 below is its integer shadow)
  - HYP-9230 (the resonance spectrum of the approach to uniformity: 5/12/41/53-TET)
  - Kontorovich-Miller (2005) and Lagarias-Soundararajan (2006), Benford's law for 3x+1 (the same irrational-rotation mechanism for leading digits of iterates)
---

# THM-4590 — Collatz classes are equidistributed

## Setting

* `T` is the Terras map on `Z`. A set `B` of positive integers (or of negative integers) is **saturated** if `n ∈ B ⟺ T(n) ∈ B`. Equivalently, `B` is a union of grand-orbit classes. Examples:
  * the basin of 1;
  * the basin of any hypothetical other cycle or divergent class;
  * on `Z_{<0}`, the basins of the cycles −1, −5 and −17.
* The **cut count** is `c_B(x) = #{n ∈ [x, 2x) : 1_B(n) ≠ 1_B(n+1)}`.
* For a window `W`, the density is `b(W) = |B ∩ W|/|W|`.
* The dyadic density is `β(x) = b([x, 2x))`.

## Statements

1. **Few cuts (PROVED).** `c_B(x) = o(x)`. More precisely, `c_B(x) ≤ q(K)·x + O(2^K)` for every `K`, where `q(K) → 0` is the Haar probability that the pair chain of `(n+1, n)` is not absorbed within `K` steps (THM-4581).
2. **Residues (PROVED).** For every `M ≥ 1` and residue `r`, `|B ∩ [x,2x) ∩ (r + MZ)| = |B ∩ [x,2x)|/M + O(M·c_B(x) + M)`. The same holds inside every window of length `≥ ηx`.
3. **Multiplicative windows (PROVED).** For every fixed `η > 0`,

       lim_(x→∞) sup_(1 ≤ c ≤ 2) | b([cx, c(1+η)x)) − b([x, (1+η)x)) | = 0 ,

   and `b([2^j x, 2^j(1+η)x)) − b([x, (1+η)x)) → 0` for each fixed `j ∈ Z`. In particular `β(cx) − β(x) → 0` for each fixed `c > 0`: the density is slowly varying.
4. **Affine invariance (PROVED; the integer shadow of index 1).** For every affine `γ(n) = 2^a 3^b n + c` that maps an arithmetic progression `P` into the integers, `B ∩ P` and `γ^(−1)(B) ∩ P` differ by a set of density 0. For example, `n ∈ B ⟺ 3n ∈ B` and `n ∈ B ⟺ n + 7 ∈ B` hold for density-one `n`.
5. **End-saturated partitions (PROVED).** Statements 1–4 hold for any partition of `Z_{>0}` that is `T`-saturated away from a density-0 set of "endpoints". An example is the **entry class** `e(n)` = the last odd number before the orbit reaches a power of 2 (5, 85, 341, …).

## Proofs

**1.**
* Fix `K` and let `E_K ⊆ Z/2^K` be the residues `r` for which the pair chain from `(0, 1)`, driven by the first `K` parities of `v = r`, is not absorbed within `K` steps. Those parities depend only on `r mod 2^K` (Terras), and `|E_K|/2^K = q(K) → 0` by THM-4581 (3).
* If `n mod 2^K ∉ E_K` and `|n| > 2^(K+1)`, then `T^t(n) = T^t(n+1)` for some `t ≤ K`. So `n` and `n + 1` lie in one class, and every saturated set contains both or neither.
* Counting gives `c_B(x) ≤ |E_K|(x/2^K + 1) + 2^(K+1)`. Hence `limsup c_B(x)/x ≤ q(K)` for every `K`. ∎

**2.**
* Put `f = 1_B`. Then `Σ_(n∈I, n≡r+1) f(n) − Σ_(n∈I, n≡r) f(n) = Σ_(m∈I, m≡r) (f(m+1) − f(m)) + O(1)`, which is at most the number of cuts in `I`, plus `O(1)`.
* So adjacent residue classes differ by `≤ c + 1`, and every class lies within `M(c + 1)` of the average. ∎

**3.** Let `W = [y, y')` with `y' − y ≥ ηx` at scale `x`. Two exact identities hold, because `T(2n) = n` and odd `n = 2j + 1` maps to `T(n) = 3j + 2`:

    |B ∩ 2W ∩ 2Z| = |B ∩ W| ,
    |B ∩ W'' ∩ (2 + 3Z)| = |B ∩ W ∩ (1 + 2Z)| ,   with W'' = [(3y+1)/2, (3y'+1)/2).

* By 2 (inside windows, cuts `o(x)`): `b(2W) = b(W) + o(1)` and `b(W'') = b(W) + o(1)`, with `|W''| = (3/2)|W|`. Reading the identities backwards gives the inverse maps.
* `log_2(3/2)` is irrational, so `{a + b log_2(3/2) : a, b ∈ Z}` is dense in `R`.
* Given `δ > 0`, compactness gives finitely many pairs `(a, b)` such that every `c ∈ [1, 2]` is within a factor `2^(±δ)` of some `2^a (3/2)^b`.
* Composing the bounded number of exact maps changes `b` by `o(1)`. Moving from `2^a (3/2)^b W` to `[cx, c(1+η)x)` changes it by `O(δ/η)`.
* Let `x → ∞`, then `δ → 0`. ∎

**4.**
* For `γ = +1` this is 1.
* For `γ(n) = 3n`, use THM-4581 (6b) at `(1, 0)`: the pair chain of `(3n, n)` is absorbed within `K` steps outside a residue set of density `q_(1,0)(K) → 0`, and the counting is as in 1.
* A general `γ` reduces to these through THM-4581 (6b)'s routing.
* A composition of `m` maps costs `m` density-0 sets. ∎

**5.**
* The entry class is constant along an orbit until the orbit reaches a power of 2.
* An equal-time merge of `n` and `n+1` before that point preserves it. The exceptions (orbits that reach a power of 2 within `K` steps of `n ~ x`) have density 0. ∎

## Numbers

**Negative integers: the three cycle basins** (`m = −n ∈ [2^(k−1), 2^k)`, i.e. the `3m − 1` map; `negbasins.c`, all `m ≤ 2^28`).

| `k` | basin of −1 | basin of −5 | basin of −17 | cut fraction | max residue deviation mod 3, 4, 9 |
|---|---|---|---|---|---|
| 10 | 0.3281 | 0.3164 | 0.3555 | 0.4141 | 0.028, 0.078, 0.093 |
| 16 | 0.3330 | 0.3215 | 0.3455 | 0.3998 | 0.006, 0.007, 0.014 |
| 22 | 0.3277 | 0.3243 | 0.3480 | 0.3611 | 0.0004, 0.0002, 0.0020 |
| 28 | 0.3268 | 0.3248 | 0.3484 | 0.3392 | 0.0000, 0.0001, 0.0001 |

* At `k = 28` (`negwindows.c`), the residue deviations of the −1 basin are about 2–3 sampling standard deviations:

  | modulus | 16 | 64 | 256 | 1024 | 27 | 81 | 243 | 729 | 11 | 121 | 7 | 49 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | deviation | 0.0006 | 0.0008 | 0.0015 | 0.0036 | 0.0004 | 0.0007 | 0.0013 | 0.0032 | 0.0002 | 0.0009 | 0.0002 | 0.0007 |

* The cut fraction falls slowly, as predicted. Unmerged pairs land in basins roughly independently, so cuts ≈ `(1 − Σ d_i²)·q(4.8 log_2 x) ≈ 0.67·q`, which gives 0.41 at `k = 10` and 0.34 at `k = 28`.
* **Multiplicative windows at `k = 28`** (64 windows of width 1.09%). The −1 basin ranges over `[0.257, 0.370]`. This is not yet uniform, because the cut fraction is still 0.34. The non-uniformity is a sharp log-periodic spectrum (HYP-9230), decaying with scale.

**Positive integers: entry classes** (`posentry85.c`, all `n ≤ 2^28`).
* Entry via 5 has density 0.9380; entry via 85 has 0.0236; both are stable from `k = 18` to `k = 28`. (21 is never an entry, since `21 ≡ 0 mod 3`.)
* The cut fraction is 0.0707 at `k = 16` and 0.0593 at `k = 28`, slowly falling.

## Scope

* This is a theorem about how classes are arranged, not about how many there are. It neither proves nor needs the Collatz conjecture.
* It constrains any counterexample. The basin of a hypothetical nontrivial cycle or divergent class must be equidistributed modulo every `M`, uniform across multiplicative windows, slowly varying in density, and almost invariant under `n ↦ 3n` and `n ↦ n + c`.
* It does not give a natural density. A slowly varying `β` need not converge.
* The rates are not uniform: see HYP-9230 for the Diophantine resonance spectrum.
