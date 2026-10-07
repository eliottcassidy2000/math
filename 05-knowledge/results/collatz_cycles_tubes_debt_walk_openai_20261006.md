# Cycles, tubes and debts: the Dumont–Reiter odd-critical-point conjecture via Lygerōs–Rozier's tubes, Chamberland's two cycles, the Mersenne debt walk, and four openai/math preprints (Liénard, Hilbert 16, Ostmann, witnessed choice)

Session opus-2026-10-06-S19. Owner prompt: "work the next steps as they emerge, also mix in and extend ideas from more papers": witnessed symmetric choice vs CPT, the additive indecomposability of the primes, two limit cycles for quintic Liénard systems, and uniform bounds for planar polynomial limit cycles (openai/math, 24 September 2026), "etc at your discretion, let your mathematical hypotheses grow freely into many proofs".

Scripts (each prints ALL CHECKS PASSED; outputs committed next to them):
* `04-computation/experiments/chamberland_tubes_dumont_reiter_20261006.py` (about 30 s; mpmath interval arithmetic), for section 2.
* `04-computation/experiments/chamberland_dumont_reiter_critical_census_20261006.py`: NUMERICAL census, run as `2001 1500`.
* `04-computation/experiments/chamberland_even_critical_fates_20261006.py`: NUMERICAL, even critical points, run as `1200`.
* `04-computation/experiments/mersenne_haar_debt_walk_20261006.py`, for section 1. Default run about 90 s (`.out`); the large run `3000 20000 61 7` about 15 min (`_large.out`).
* `04-computation/experiments/mersenne_certified_density_numba_20261006.py` (the S18 script): run as `9 31`; `.out` committed now, including the K ≤ 27 values the S18 note cites.
* `04-computation/experiments/collatz_papers_misc_20261006.py`, for sections 1.8, 3 and 5.

Status words follow the repo: PROVED, PROVED (computer-assisted, interval arithmetic), FINITE-EXACT, NUMERICAL, HEURISTIC, ANALOGY, DICTIONARY, CITED, KNOWN. The openai/math preprints are AI-written and unrefereed; their results are reported as claims.

**Independently audited** (section 7). The first version missed the decisive prior art, **Lygerōs–Rozier 2014**, for the Chamberland half. Corrections are logged as MISTAKE-580.

## 0. Answers

**Tubes for two real extensions of the Collatz map (THM-4563, PROVED, computer-assisted).**
* `T` is the Terras map. The two extensions are:
  * Chamberland's `C(x) = x + 1/4 − (2x+1)/4 · cos(πx)`;
  * Dumont–Reiter's 3-power extension `D(x) = (3^s x + s)/2`, `s = sin²(πx/2)`.
* For each, every right tube `τ(m) = [m, m + a/(2m+1)]` with `m ≥ 1` is mapped into `τ(T(m))`: `a = 0.8, 0.9` for `C` and `a = 0.6, 0.7` for `D`. The local maximum `c_n` next to an odd `n` lies in `τ(n)`, so the critical orbit shadows the Collatz orbit of `n` forever.
* **For `C` this is KNOWN.** It is Lygerōs–Rozier 2014, Lemma 2.4, with tubes `[n, n + a/(π²n)]`, `a = 7/2`, for all `n > 0`; their small odd `n` were checked in floating point.
* What is ours on `C` is a fully interval-certified proof covering all integers at once. The exact scaled maps `O(ξ, h)`, `E(ξ, h)` with `h = 1/(2m+1)` are verified over the continuum `h ∈ [0, 1/3]`.
* **For `D` it is new**, by the same method.

**Dumont–Reiter's Odd Critical Point Conjecture (2003) is settled as far as it can be.** No prior proof found; Lygerōs–Rozier, footnote 4, mention the conjecture.
* (ii) `c_n` and `n` have the same total stopping time: **PROVED for every odd `n`**.
* (iii) `c_n` lies in `n`'s immediate basin of that stopping time: **PROVED for every odd `n`**, in the real sense.
* (i) "`c_n` is attracted to `(1,2)`" is **PROVED equivalent to the Collatz orbit of `n` reaching 1**. So (i) for all odd `n` is exactly the 3x+1 conjecture.

**Chamberland's map** (KNOWN results re-proved rigorously, plus a few new items).
* KNOWN (Lygerōs–Rozier 2014, Theorem 3.3 and Corollary 3.4), re-proved:
  * every odd `n ≥ 7` whose orbit reaches 1 sends `c_n` to `A1 = {1,2}`;
  * 3x+1 holds iff all those `c_n` go to `A1`.
* New, rigorous:
  * `c_1`, `c_3` and `c_5` go to `A2 = {1.19253…, 2.13866…}` (interval trap; these were numerical);
  * a **flip lemma**: a point within scaled distance `ξ ∈ [−1.3, −0.45]` left of an even integer is thrown into the right tube, and shadows from then on;
  * a **Singer refinement**: an attracting cycle whose immediate basin contains an odd critical point is `A1`, `A2`, or lives in the tubes of a nontrivial positive integer cycle.
* NUMERICAL:
  * On `τ(1)` the scaled return map has a repelling separatrix at `ξ = 0.0710584`, Lygerōs–Rozier's `x1 = 1.023686`, and `A2` at `ξ = 0.5775957`. In these coordinates `A2` is the satellite of `A1` in its own tubes.
  * Even critical points can leave the shadow (`c_54`). Some even ones go to `A2` (`c_382`, `c_496`, `c_502`; as in Lygerōs–Rozier).
* **Fixed points.**
  * The mirror law `C′(−1 − x) = 2 − C′(x)` (Lygerōs–Rozier (5.2)).
  * The attracting fixed points are exactly 0 and `−1.2777338` (interval root isolation; Lygerōs–Rozier Table 1).
  * Correction: the S15 note had `0.2777` and `−1.2777` swapped (MISTAKE-578).
* **Negative Schwarzian on `[0, ∞)`** (Chamberland's claim, re-proved): `2C′C‴ − 3C″² = π²Q(π(x + 1/2))`. It fails at `x = −0.0217160`.

**Next steps 1–2 of S18, the Mersenne switch.**
* **Haar model** (PROVED reduction, via THM-4556 (iv)): switching is a question about Haar-random 2-adic inputs.
* **Certified shares** (FINITE-EXACT): `0.1603, 0.1649, 0.1695, 0.1740` at `K = 28, …, 31`.
* **Monte Carlo** (NUMERICAL, two committed runs).
  * The non-switch share decays like `T^(−α)`: `α = 0.68` (large run) and `0.78` (small run), with bootstrap intervals together spanning about `[0.63, 0.90]`.
  * The lag-1 non-merge share fits `T^(−1/2)`: tail exponent `0.48`–`0.52`, `√T·q₁ ≈ 14–17`.
* **Theorem D** (PROVED; the bookkeeping is mac-mini's HYP-9214 model, made exact).
  * With `x = 2^L y + Δ` and the debt `D = 3Δ + 1 − 2^L`, the next relation is the identity iff `D = 0`, and then `L = 2k ≠ 0` and `x = 4^k y + (4^k − 1)/3` (a sibling pair).
  * A merge with `D ≠ 0` is a value coincidence of Haar measure 0. All 834 sampled merges pass through `D = 0`.
  * For integer debts the debt's odd part runs the Collatz map kicked by a power of two: `D′_odd = U(D_odd) − 2^(L′−w′)`.
* **The two exponent streams** (NUMERICAL).
  * In generic steps with `L ≥ 8`, the debt exponent `w` is `Geom(1/2)`-distributed to four decimals, independent of the partner exponent `a` (`χ² = 7.4` on 9 dof), with autocorrelation `−0.002`.
  * So far from 0, `L` is statistically the symmetric random walk with i.i.d. increments `a − w` (variance 4.00).
* **HYP-9217** (lag-1 debt merges almost surely with a `T^(−1/2)` tail; the switching set has full measure) would give `μ_2(S) = 1` and hence HYP-9213.
* **Against actual exponents** (NUMERICAL).
  * The model predicts the switching share (0.968–0.974 predicted, 0.970 observed) and the `σ`-level counts (`23, 37, 54, 69` observed, `22–24, 40–41, 56, 69–72` predicted).
  * It over-predicts lag-1 merging by 0.02–0.04.

**The four preprints** (claims; ANALOGY / DICTIONARY unless marked).
* **Liénard**, two limit cycles: a half-orbit matching argument with a Schwarzian identity. The rhyme is Chamberland's two attracting cycles and the Singer count.
* **Uniform bounds** (Hilbert 16).
  * The Mersenne run is an exact adelic saddle passage at −1, with Dulac exponent `log₂3` and `g^a(2^a − 1) = 3^a − 1` (PROVED, elementary).
  * Remark H: `3x + d` has unboundedly many cycles (Lagarias 1990, CITED). Our one-run construction gives only non-primitive ones.
* **Ostmann.** The Paley pair `(QR, QR ∪ {0})` is the symmetric local model of the residue partition.
* **Witnessed choice.**
  * The cycle charge `c_w mod (2^A − 3^l)` is a CFI-type global charge (DICTIONARY).
  * Two witnessed choices (a vertex, then an out-neighbour) plus colour refinement canonize every Paley tournament `P_p`, `p ≡ 3 (mod 4)`, `p < 2000` (FINITE-EXACT).
  * The second refinement round reads the Legendre family `Y² = X(X−1)(X−x)` (PROVED identity).

Collatz is OPEN.

---

## 1. The Mersenne switch: a Haar model, a debt walk and a recurrence law

### 1.1 The Haar model (PROVED reduction)

Recall THM-4556 (iv).
* For odd `a`, the post-run start of `M_a = 2^a − 1` is `p = 2·3^(a−1) − 1`, and that of `M_(a−D)` is `q_D = 2·3^(a−D−1) − 1`.
* A **switch at lag `D`** is an equal-time merge `U^i(p) = U^(i+D)(q_D)`.
* Its colliding prefix of template total `K` depends only on `p mod 2^(K+1)`, that is, on `3^(a−1) mod 2^K`, that is, on `a mod 2^(K−2)`.
* As `a` ranges over odd 2-adic integers, `X = 3^(a−1)` ranges over `1 + 8Z_2` with Haar measure: `a ↦ 3^a` is an isomorphism of `Z_2` onto the closed subgroup `{x ≡ 1, 3 (mod 8)}`, and pushes Haar to Haar.

Hence for every `K` the certified share at template total `K` is the Haar measure of an explicit clopen event in `X`. The switching set `S ⊂ Z_2` (exponents that switch at some finite total) has
* `μ_2(S) = lim_K (certified share at K)`;
* `μ_2(S) = 1 ⟹ HYP-9213` (S18, Proposition 6; this gives `o(A)` only).

Sampling `X` uniformly and following the two orbits 2-adically (with precision bookkeeping, 64 guard bits) estimates `P(switch with total ≤ T)` far beyond exhaustive certification.

### 1.2 Certified density to `K = 31` (FINITE-EXACT)

The committed S18 numba search, run as `9 31` (output `mersenne_certified_density_numba_20261006.out`), gives exact shares among odd classes:

| `K` | 27 | 28 | 29 | 30 | 31 |
|---|---|---|---|---|---|
| certified share | 5221376/33554432 = 0.1556 | 10757550/67108864 = 0.1603 | 22135584/134217728 = 0.1649 | 45495372/268435456 = 0.1695 | 93407400/536870912 = 0.1740 |

The increments are about 0.0046 per unit of `K`, consistent with the Haar curve. Exhaustive certification cannot reach the median merge total (several hundred).

### 1.3 Monte Carlo (NUMERICAL; two committed runs)

* **Small run.** Default: seed 2026, `N = 1200`, 12 000-bit precision, odd `D ≤ 41`, totals to 9600.
* **Large run.** `3000 20000 61 7`: seed 7, `N = 3000`, 20 000 bits, odd `D ≤ 61`, totals to 19 924.

| total `T` | 27 | 100 | 400 | 1600 | 3200 | 6400 | 12800 | 19000 |
|---|---|---|---|---|---|---|---|---|
| large `q(T)` | 0.850 | 0.615 | 0.311 | 0.117 | 0.0757 | 0.0443 | 0.0290 | 0.0223 |
| large `q₁(T)` | 0.859 | 0.722 | 0.548 | 0.356 | 0.276 | 0.202 | 0.143 | 0.118 |
| large `√T·q₁` | 4.5 | 7.2 | 11.0 | 14.2 | 15.6 | 16.1 | 16.2 | 16.3 |
| small `q(T)` | 0.848 | 0.604 | 0.323 | 0.120 | 0.0725 | 0.0400 | — | — |
| small `√T·q₁` | 4.5 | 7.1 | 10.9 | 13.9 | 14.7 | 14.5 | — | — |

* At `T = 27` the Haar switch share is `0.150 ± 0.007` (large) and `0.153 ± 0.010` (small), against the exact 0.1556.
* **Any-lag exponent** (least squares on `T ≥ 400`):
  * large: `α = 0.679`, bootstrap 95% `[0.625, 0.745]`;
  * small: `α = 0.782`, bootstrap 95% `[0.689, 0.900]`.
  * We quote `α ≈ 0.7 ± 0.1`. The two runs differ in seed, precision and lag range, and the fit window is pre-asymptotic.
* **Lag 1.**
  * The tail exponent on `T ≥ 3200` is `0.477` (large, `[0.435, 0.520]`) and `0.521` (small, `[0.435, 0.609]`), consistent with 1/2.
  * `√T·q₁` at the largest total is 16.27 (large, `[14.70, 17.87]`) and 14.37 (small, `[12.41, 16.33]`), so `c₁ ≈ 14–17`.
  * Over `T ≥ 400` the exponent is about 0.40–0.42, the value in HYP-9214. Typing it the pre-asymptotic part of a `T^(−1/2)` law is HEURISTIC.
* **Least-total lags** among switching samples: `D = 1` (52% large, 57% small), 3 (17%, 16%), 5 (9%, 10%), 7 (6%, 5%), … .

### 1.4 Theorem D (debt bookkeeping; PROVED)

The relation `x = 2^L y + Δ`, the update `L′ = L + a − b`, `Δ′ = D/2^b`, the walk interpretation and the sibling merge configuration are mac-mini's HYP-9214 model. What we add is (b) in its exact form and (c).

**Setting.** Let `x, y` be odd 2-adic integers with `x = 2^L y + Δ`, `L ∈ Z`, `Δ ∈ Z[1/2]`. Put `D = 3Δ + 1 − 2^L`. Let `a = v_2(3y + 1)`, `y′ = U(y)`, `b = v_2(3x + 1)`, `x′ = U(x)`.

(a) **Transition.** `3x + 1 = 2^(L+a) y′ + D`, so `b = v_2(2^(L+a) y′ + D)`, and `x′ = 2^(L′) y′ + Δ′` with `L′ = L + a − b`, `Δ′ = D/2^b`.

(b) **Merge criterion.**
* `D = 0` iff the next relation is the identity (`L′ = 0`, `Δ′ = 0`), i.e. `x′ = y′` as affine functions of `y`.
  * Then `Δ = (2^L − 1)/3`, so `L` is even and nonzero, and `x = R^k(y) = 4^k y + (4^k − 1)/3` (`k = L/2`), or the mirror for `k < 0`.
* A merge `x′ = y′` with `D ≠ 0` requires the value coincidence `y′ = −Δ′/(2^(L′) − 1)`.
  * It happens: `x = 5`, `y = 1`, `L = 0`, `Δ = 4` gives `D = 12` but `U(5) = U(1)`.
  * When `x` is an affine function of Haar `y`, it has measure 0.
* So in the Haar model, merges occur almost surely through `D = 0`.

(c) **Kicked Collatz law.** Suppose `D ≠ 0` is an integer (as when `L ≥ 0`) with `v_2(D) < L + a`, so that `b = v_2(D)`. Write `D_odd = D/2^(v_2 D)` and `w′ = v_2(3D_odd + 1)`. If `L′ > w′`, then
`v_2(D′) = w′`, `D′/2^(w′) = U(D_odd) − 2^(L′ − w′)`, and `L′ = L + a − v_2(D)`.

**Proof.**
* (a) Substitute `x = 2^L y + Δ` into `3x + 1` and use `3y + 1 = 2^a y′`.
* (b) `x′ − y′ = (2^(L′) − 1) y′ + Δ′`. This vanishes identically in `y′` iff `L′ = 0` and `Δ′ = 0`, i.e. `D = 0`. Then `3Δ = 2^L − 1`, and `Δ ∈ Z[1/2]` forces `3 | 2^|L| − 1`, i.e. `L` even.
* (c) `D′ = 3D/2^b + 1 − 2^(L′) = 3D_odd + 1 − 2^(L′)`. If `L′ > w′`, factor `2^(w′)`. ∎

**Scope.** (c) needs an integer debt. Over the 1000 sampled pairs, steps with `L ≤ 0` (where `D` may be dyadic) are about 45% of all steps.

**Checks** (FINITE-EXACT; `mersenne_haar_debt_walk_20261006.py` B). Over 1000 lag-1 pairs `x_0 = 3·2^v y_0 + 1`, the Mersenne debt states of THM-4556 (ii):
* the transition law holds on all first-60 steps;
* the kicked law holds in all 28 693 generic integer-debt steps tested;
* all 834 merges pass through `D = 0` with `L` even and nonzero.
* `L` one step before a merge is `−2` (489), `+2` (301), `−4` (20), `+4` (11), `±6` and `−8`. About 61% have `k < 0`, matching mac-mini's "about 2/3".

### 1.5 The walk and its two exponent streams (NUMERICAL)

Over all steps of the 1000 pairs (12 000-bit sources, up to 4800 steps each):

| `L` range | `≤ 0` | 1–6 | 7–15 | `≥ 16` | all |
|---|---|---|---|---|---|
| mean increment (s.e.) | +0.0006 (0.0026) | −0.050 (0.010) | +0.012 (0.008) | +0.002 (0.003) | +0.0004 |
| variance | 3.98 | 4.04 | 4.01 | 4.02 | **4.002** |

* Globally the walk has mean 0 and variance 4. In the merge region `L ∈ [1, 6]` it has a small negative drift (about 5 s.e.), which is where the merges and ties act.
* **The two streams.** In the 672 443 generic integer-debt steps with `L ≥ 8`, the increment is `a − w`. Here `a` is the partner's exponent and `w = v_2(D) = b` is the debt's own exponent; `b = w` in every such step, as (c) says.
  * `P(w = k) = 0.5001, 0.2500, 0.1250, 0.0624, 0.0311` for `k = 1, …, 5`, against `2^(−k)`.
  * `a` and `w` are independent: `χ² = 7.4` on 9 dof.
  * `w` has lag-1 autocorrelation `−0.002`.
  * Far from 0, `L` is statistically a symmetric random walk with i.i.d. `Geom − Geom` increments. A HEURISTIC reason: Haar measure is invariant under `U`, and the kicks rewrite only bits at depth `L′ − w′`, so the debt's low bits stay Haar-distributed.
* In steps, `P(no merge by t)` is 0.627, 0.490, 0.324, 0.216 and 0.166 at `t = 100, 300, 1000, 3000, 4800`; `√t·P` is 6.3, 8.5, 10.3, 11.8 and 11.5.

**HEURISTIC mechanism.**
* A symmetric walk with finite variance is recurrent.
* The debt is pinned to the scale `|D| ≈ 2^L` (the kick `2^(L′)` dominates in (c)). So returns of `L` to `{±2, ±4}` bring the debt to small values, where `D = 0` has positive chance.
* That gives `P(no merge by t) ≈ c·t^(−1/2)`.
* Not proved: the asymptotic independence and fairness of the debt's exponent stream, a deterministic Collatz orbit randomized only through the kicks. This is mac-mini's open question (per-visit success bounded below).

### 1.6 HYP-9217 (the debt recurrence law)

`05-knowledge/hypotheses/HYP-9217-mersenne-debt-recurrence-law-t-minus-half.md`. In the Haar model:
1. the lag-1 Mersenne debt merges almost surely, with `P(no merge by template total T) = c₁ T^(−1/2)(1 + o(1))`, `c₁ ≈ 14–17`;
2. the any-lag non-switch probability satisfies `q(T) = T^(−α + o(1))` with `α ≈ 0.7 ± 0.1`.

(2) gives `μ_2(S) = 1` and hence HYP-9213 (`o(A)` new `σ`-levels). The finer count `#levels(A) ≈ A^(1−α)` is HEURISTIC: it also needs the Haar model to govern actual exponents at `T_post(a)`, which 1.7 shows is only approximately true. The reset-2 debt of HYP-9214 is a natural extension; we did not test it.

### 1.7 The model against the actual exponents (NUMERICAL)

For `2 ≤ a ≤ 2001` the script computes `σ(M_a)` and the post-run template total `T_post(a)` (about `7.68a`).

* **Any partner.** For odd `a ∈ [1001, 2001]` the share with a smaller `σ`-partner is 0.970. The Haar prediction `1 − mean q(T_post(a))` is 0.974 with the small run's fit and 0.968 with the large run's. mac-mini's `[10^3, 6·10^3]` census gives 0.980.
* **Level counts.** Distinct `σ` values for `a ≤ 100, 400, 1000, 2001` (counting `a = 1`) are `23, 37, 54, 69`.
  * Predicted `2 + Σ_(odd a) q(T_post(a))`: `24.2, 41.4, 55.8, 68.6` (small fit) and `21.9, 39.5, 56.0, 72.0` (large fit).
  * In this window both fits grow like `A^0.35`, partly through the cap `min(1, ·)` and the offset. So the agreement does not pin `α`.
* **Lag 1.** The share with `σ(M_(a−1))` is 0.824. Predicted: 0.864 with the small run's `c₁ = 14.4`, 0.846 with the large run's 16.3.
  * The model over-predicts lag-1 merging by 0.02–0.04 (1.3–2.3σ). Actual orbits end, so their effective time is shorter than `T_post`.
  * The ratio of the observed lag-1 failure shares in `[1001, 2001]` and in mac-mini's `[10^3, 6·10^3]` (0.176/0.118 = 1.49) matches the `a^(−1/2)` law's 1.43.

### 1.8 The adelic saddle (PROVED, elementary; DICTIONARY with the Hilbert-16 preprint)

Let `g(x) = (3x+1)/2`. Then `g(x) + 1 = (3/2)(x + 1)`, so `g^k(x) + 1 = (3/2)^k (x + 1)` and `g^a(2^a − 1) = 3^a − 1` (THM-4556 (i), binary repunit to ternary repunit).

At the fixed point −1 the multiplier `3/2` has
* `|3/2|_2 = 2`: unstable;
* `|3/2|_3 = 1/3`: stable;
* `|3/2|_∞ = 3/2`;
* and the product formula gives volume neutrality.

So the run is a **saddle passage**.
* It enters at 2-adic distance `r = 2^(−a)` and leaves at 3-adic distance `3^(−a) = r^λ`, after flight time `a`.
* The Dulac exponent is `λ = log 3/log 2 = log₂3`.

The uniform-bounds preprint's local model `ẋ = x, ẏ = −λy` (flight time `L = log(1/r)`, contraction exponent `W = λL`) is this, with the real and `p`-adic places as the two directions. Its "two clocks" are the CRT-independent 2-adic and 3-adic clocks of S18 Proposition 5.

## 2. Chamberland's map and the 3-power map: tubes, shadows, limit cycles

### 2.1 Prior work (CITED)

* **Chamberland 1996** (Dynam. Contin. Discrete Impuls. Systems 2, 495–509; read through Lagarias's annotated bibliography, entry 35, and Chamberland's 2003 survey, §6.5):
  * negative Schwarzian on `R+`;
  * `[0, μ1)` is attracted to 0;
  * `[μ1, μ3]` is invariant, with a.e. point going to `A1` or `A2`;
  * a nontrivial positive cycle of `T` would be attracting;
  * the "Critical Points" conjecture: all `c_n` go to `A1` or `A2`.
  * The survey calls the two-cycle conjecture "equivalent to the 3x+1 problem"; by 2.6 (e) it implies the cycle half.
* **Dumont–Reiter 2003** (Dyn. Contin. Discrete Impuls. Syst. Ser. A 10, 875–893; preprint read):
  * negative Schwarzian of `D` for `x ≥ 0` (Theorem 5);
  * `c_n = n − (−1)^n (2/(π² ln 3))/n + O(n⁻²)` (Theorem 6);
  * every point of `(μ1, μ3)` is attracted to `(1,2)` (Theorem 8);
  * total stopping time of real `x` = least `k` with `μ1 < D^k(x) < μ2`;
  * **Conjecture 1 (Odd Critical Point Conjecture)**: for odd `n`, (i) `c_n` is attracted to `(1,2)`; (ii) `c_n` and `n` have the same total stopping time; (iii) `c_n` is in the immediate basin of total stopping time equal to the value at `n`.
  * Evidence: odd `n < 180 000` and Vyssotsky's `n_V`. The even `c_54` fails to shadow (their Table 6).
* **Lygerōs–Rozier 2014** (Ratio Math. 26, 77–94; arXiv:1402.1979; read), for Chamberland's `C`:
  * Lemma 2.3: localization of `c_n`.
  * **Lemma 2.4: tubes `[n, n + a/(π²n)]`**, `27/8 < a < 6` for large `n`, and **all `n > 0` for `a = 7/2`**. Odd `n = 1, 3, 5, 7, 9` are checked in floating point.
  * **Theorem 3.3**: `c_n → A1` for odd `n ≥ 7` reaching 1 (nodes 13, 16, 40; separatrix `x1 = 1.023686`).
  * **Corollary 3.4**: the 3x+1 reformulation.
  * Divergent orbits and wandering intervals.
  * Census `n ≤ 2000` (1500 digits): `A2` captures `c_n` for `n = 1, 3, 5, 382, 496, 502, …`; the orbits stay close except for `n ≡ −2 (mod 64)`, 54, ….
  * The mirror law (5.2) and the fixed points (5.3), Table 1.
  * Footnote 4 mentions Dumont–Reiter's conjecture; it is not proved there.
* **Letherman–Schleicher–Wood 1999** (Experiment. Math. 8, 241–251; via Lagarias's entry 115): entire interpolations with superattracting integers. No two integers share a Fatou component (except possibly −1, −2); divergence corresponds to a wandering domain.

We found no proof of the Odd Critical Point Conjecture in the sources checked: the three papers, Lagarias's bibliographies, Chamberland's survey, the audit's literature search, and web searches on 2026-10-06/07.

### 2.2 Proposition C1 (mirror law, attracting fixed points; PROVED; KNOWN: Lygerōs–Rozier (5.2)–(5.3))

With `u = x + 1/2`, `C(x) − x = 1/4 − (u/2) sin(πu)`, which is even in `u`. So `C′(−1 − x) = 2 − C′(x)` at mirror fixed points, checked at all 41 fixed points in `[−20, 20]`.

At a fixed point, `cos πx = 1/w` with `w = 2x + 1`, so `C′ = 1 − 1/(2w) ± (π/4)√(w² − 1)`. Then `|C′| < 1` forces `|w| < √(1 + 100/π²) = 3.3365`.

Interval root isolation in that window finds exactly four fixed points. Each candidate cluster has an interval derivative excluding 0 and a sign change. Interval multipliers:
* attracting at `0` (`C′ ∈ [0.4991, 0.5012]`) and `−1.2777338` (`[0.3833, 0.3914]`);
* repelling at `0.2777338` and `−1`.

*Correction (MISTAKE-578).* The S15 note had the two non-integer members swapped. The note, its D51 direction, its INDEX line and the procgen synthesis line are corrected.

### 2.3 Proposition C2 (negative Schwarzian on `[0, ∞)`; PROVED, computer-assisted; Chamberland's claim re-proved)

**Identity.** With `w = π(x + 1/2)`, `s = sin w`, `c = cos w`:

`2C′C‴ − 3C″² = π² Q(w)`, where `Q(w) = −(1/2 + s²/4)w² + c(1 + s)w + (3/2)(s² + 2s − 2)`.

**Proof of `Q(w) < 0` for `x ≥ 0`.**
* For `w > 2√3`: `Q ≤ −w²/2 + (3√3/4)w + 3/2 < 0`.
* On `[π/2, 2√3]`: 4000 interval cells, with maximum `Q(π/2) = −0.3506`.

The sign changes at `x = −0.0217160`.

### 2.4 Theorem T (tube theorem; PROVED, computer-assisted; for `C` KNOWN: Lygerōs–Rozier Lemma 2.4)

**Scaled coordinates.** For `F ∈ {C, D}` and an integer `m ≥ 1`, write `x = m + ξh` with `h = 1/(2m+1)`, and the image coordinate `ξ′ = (2T(m) + 1)(F(x) − T(m))`. Then `ξ′ = O(ξ, h)` for odd `m` and `ξ′ = E(ξ, h)` for even `m`, with:

* `C`, odd `m`: `O = (3+h)/2 · [ξ(1 + cos(πξh)/2) − (π²ξ²/8) S²]`.
* `C`, even `m`: `E = (1+h)/2 · [ξ(1 − cos(πξh)/2) + (π²ξ²/8) S²]`.
* `D`, odd `m`: `O = (3+h)/4 · [−(3(1−h)/2) ln3 (π²ξ²/4) S² φ(σ ln3) + 3^(1−σ) ξ − h(π²ξ²/4) S²]`.
* `D`, even `m`: `E = (1+h)/4 · [((1−h)/2) ln3 (π²ξ²/4) S² ψ(σ ln3) + 3^σ ξ + h(π²ξ²/4) S²]`.

Here `S = sinc(πξh/2)`, `σ = sin²(πξh/2)`, `φ(y) = (1 − e^(−y))/y` and `ψ(y) = (e^y − 1)/y`. These are exact identities at `h = 1/(2m+1)` and entire in `(ξ, h)`. At `h = 0` they give `O → (9/4)ξ − (3π²/16)ξ²` and `E → ξ/4 + (π²/16)ξ²` for `C`, with `ln 3` factors for `D`. Checked against direct evaluation in 48 cases.

**Theorem T.**
* For `C` with `a ∈ {0.8, 0.9}`, and for `D` with `a ∈ {0.6, 0.7}`:
  * `0 < O(ξ, h) ≤ a` on `(0, a] × [0, 1/3]`;
  * `0 ≤ E(ξ, h) ≤ a` on `[0, a] × [0, 1/5]`.
* Hence `F(τ(m)) ⊂ τ(T(m))` for every integer `m ≥ 1`.
* For odd `n`, `F′(n) = 3/2`, `F′(n + a/(2n+1)) < 0` and `F″ < 0` on `τ(n)`, both for `a` and both maps. So `c_n` is the unique critical point in `τ(n)`.

**Proof.**
* Outward-rounded interval arithmetic on box covers, 1200 boxes per constant.
* `sinc`, `φ` and `ψ` are enclosed by monotonicity on the ranges used (`sinc` up to 0.84; `φ`, `ψ` up to 0.23).
* The tube constants are exact decimals: boxes cover up to the interval enclosure's upper end, and images are compared with its lower end.
* Because `h` runs over the continuum `[0, 1/3] ∋ 1/(2m+1)`, all integers are covered at once. ∎

**Numbers.**
* For `C`, `ξ_0 = (2n+1)(c_n − n) → 6/π² = 0.6079`, and the peak image is `6.75/π² = 0.6839`.
* The largest scaled deviation along all orbits of odd `n ≤ 401` is 0.7037 (`C`) and 0.4559 (`D`).
* The audit's independent enclosures give maxima: `O` 0.716 (`C`) and 0.509 (`D`); `E` 0.734 / 0.8965 (`C`, `a = 0.8 / 0.9`) and 0.456 / 0.587 (`D`).
* In the limit the maximal invariant interval for `C` is `[0, 12/π²)`, matching Lygerōs–Rozier's `a < 6`. At `m = 2`, `a = 0.95` fails.
* Lygerōs–Rozier's `[n, n + a_LR/(π²n)]` corresponds to `a ≈ 2a_LR/π²`, so their `7/2` gives `a ≈ 0.709`.

### 2.5 Corollary DR (the Odd Critical Point Conjecture; PROVED, new)

For the 3-power map `D` and every odd `n ≥ 1`:

* **(ii)** The total stopping time of `c_n` equals that of `n`, finite or infinite.
  * `D^k(c_n) ∈ τ(T^k(n))` by Theorem T.
  * `τ(1) = [1, 1.2] ⊂ (μ1, μ2) = (0.3158162, 1.5155526)`, and `τ(m) ⊂ [2, ∞)` for `m ≥ 2`.
  * So `D^k(c_n) ∈ (μ1, μ2)` iff `T^k(n) = 1`.
* **(iii), real sense.** `τ(n)` is connected, contains `n` and `c_n`, and its points all have total stopping time `σ(n)`.
  * So `c_n` lies in the component of the real level set containing `n`.
  * Dumont–Reiter illustrate (iii) in the complex plane, where they do not define the stopping time. We claim only the real statement.
* **(i) ⟺ the Collatz orbit of `n` reaches 1.**
  * On `τ(1)`, `D²` in scaled form is `W_D = E(·, 1/5) ∘ O(·, 1/3)`, and `W_D(ξ)/ξ < 1` on `(0, 0.6]` (interval arithmetic). With `W_D ≥ 0`, every point of `τ(1)` converges to `(1,2)`.
  * Conversely, if `n` never reaches 1, its orbit consists of integers `≥ 3`, and `c_n` stays in their tubes.

So Conjecture 1 (ii) and (iii) (real sense) are theorems, and Conjecture 1 (i) for all odd `n` is equivalent to the 3x+1 conjecture. Dumont–Reiter's remark that the conjecture implies there are no nontrivial cycles is the easy half of this.

### 2.6 Chamberland's map: where the critical points go

* **(a) KNOWN, re-proved: odd `n ≥ 7` go to `A1`** (Lygerōs–Rozier Theorem 3.3).
  * The backward tree: `20 ← {13, 40}`, `32 ← {21, 64}`, `5 ← {10, 3}`, `8 ← {5, 16}`.
  * Odd `n ≥ 5` never meet a multiple of 3 after the start, since `T(x) = 3·2^k` forces `x = 3·2^(k+1)`.
  * So every odd `n ≥ 7` reaching 1 passes through 13, 21, 40 or 64 (checked to 200 001). Lygerōs–Rozier use 13, 16, 40.
  * The interval pushes of those tubes land in `τ(1)` at `ξ ≤ 0.0572, 0.0291, 0.0138, 0.0076`.
  * `W_C = E(·, 1/5) ∘ O(·, 1/3)` satisfies `W_C(ξ) < ξ` on `(0, 0.0705]` with `W_C ≥ 0`, so `[0, 0.0705]` lies in `A1`'s basin.
* **(b) New rigour: `c_1`, `c_3`, `c_5` → `A2`.**
  * `W_C` maps `J = [0.55, 0.60]` into `[0.5702, 0.5820]` with `|W′| ≤ 0.377` (subdivided intervals).
  * Interval enclosures of `c_1`, `c_3`, `c_5` enter `J` after 1, 6 and 5 returns.
* **(c) The separatrix and the satellite.**
  * `W_C(0.3) > 0.3` and `W_C(0.0705) < 0.0705`, so a fixed point lies in `(0.0705, 0.3)` (IVT).
  * NUMERICAL: the fixed points are `0.0710584` (`W′ = 1.2148`; Lygerōs–Rozier's `x1`) and `0.5775957` (`A2`, `W′ = −0.2308`).
  * In the scaled tube coordinates `A2` is a satellite of `A1`. Dumont–Reiter's Table 4 calls it "a cycle shadowing the cycle (1,2)".
* **(d) The equivalence** (KNOWN, Lygerōs–Rozier Corollary 3.4). 3x+1 holds iff every odd critical point `c_n`, `n ≥ 7`, of `C` is attracted to `{1,2}`:
  * if `n` never reaches 1, `c_n` stays in tubes of integers `≥ 3`;
  * if `n` diverges, `c_n → ∞`;
  * if `n` enters a nontrivial cycle `Γ`, `c_n` accumulates in `Γ`'s tubes.
* **(e) Singer refinement (PROVED, using Singer 1978).**
  * Let `Γ*` be an attracting or neutral cycle of `C` in `(0, ∞)` other than `A1`, `A2`, with `[0, μ1)` the basin of 0.
  * Its immediate-basin components are bounded open intervals in `(0, ∞)`. An unbounded one would contain the start of a divergent orbit, and divergent orbits start beyond every bound by the intermediate-value chains through `[2j+2, 2j+3]` (Theorem K(c) of the Kawasaki audit note).
  * So Singer's argument (`C` is `C³`, `SC < 0` by C2, `C([0, ∞)) ⊂ [0, ∞)`) puts a critical point in the immediate basin.
  * If that critical point is odd, `c_n`, then `n` cannot reach 1 (by (a), (b)) and cannot diverge, so `n` enters a nontrivial positive cycle `Γ`, and `Γ*` lies in the closed tubes of `Γ`.
  * Lygerōs–Rozier already note that every positive integer cycle's immediate basin contains a critical point. This is the converse direction for odd critical points.
* **(f) Even critical points** (`c_m` just left of even `m`, `ξ_0 → −2/π²`) are not controlled by tubes.
  * The left side of an odd integer expands: `O(ξ) < (9/4)ξ` for `ξ < 0`.
  * `c_54` leaves the near-integer regime at step 9, at the integer 242.
  * `c_382`, `c_496` and `c_502` go to `A2` (NUMERICAL, 1500 digits, as in Lygerōs–Rozier). All even ones up to 300 go to `A1` (census).
  * **Flip lemma (PROVED, computer-assisted, new).** For every even `k ≥ 2` and `ξ ∈ [−1.3, −0.45]`, `0 ≤ E(ξ, 1/(2k+1)) ≤ 0.8` (4000 boxes). An orbit point at that scaled position left of `k` is thrown into `τ(T(k))` and shadows from then on.
  * Census (NUMERICAL, 400 digits; `chamberland_even_critical_fates_20261006.py`): of the even `m ≤ 1200`, 373 stay left of the integers until 1, 114 flip into a right tube, and 113 escape (`|ξ| > 1.3`).
  * By (e), Chamberland's two-cycle conjecture is equivalent to: no nontrivial positive cycle, and no spurious attracting cycle whose immediate basin holds only even critical points. It does not touch divergence.

### 2.7 Rhymes with the two limit-cycle preprints (ANALOGY / DICTIONARY)

* **Liénard (claims of the preprint).**
  * `ẋ = y − F(x), ẏ = −x`, `deg F ≤ 5`: at most two limit cycles, sharp, e.g. `F = ε(4x − 20x³/3 + 8x⁵/5)`, averaged displacement `−πs²(s² − 1)(s² − 4)`.
  * The proof:
    * keeps the two half-orbits separate (profiles `F(±√(2u))`, even coefficients shared, odd ones flipped);
    * fits quadratic profiles;
    * controls a width coefficient by a Riccati linearization and a Schwarzian identity for the endpoint correspondence;
    * counts zeros of a matching function made monotone by an integrating factor.
  * The parallels here:
    * Chamberland's core interval has exactly two attracting cycles, again under Schwarzian control (Singer).
    * The mirror law is an even/odd split of the displacement.
    * The tube theorem is a one-sided statement: right sides of integers are invariant, left sides of odd integers expand.
  * The coincidence of the two "two"s is NUMEROLOGY. The shared projective mechanism is ANALOGY.
* **Uniform bounds (claims of the preprint).**
  * Return maps are matching systems of passages, with flags of scales; a projection-counting argument gives `B(d)`.
  * Here Singer's critical-point count plays the counting role, and the tubes tie the odd critical points to the integers, so the count reduces to Collatz itself.
  * Remark H (section 3) records the opposite of uniform boundedness in the `3x + d` family. Section 1.8 types the saddle passage.

## 3. Remark H: no uniform bound for `3x + d` (CITED; elementary construction)

Lagarias (1990, Acta Arith. 56, 33–53; via Lagarias's bibliography, entry 108) shows that infinitely many `k` have at least `k^(1−ε)` distinct primitive cycles (`gcd(x, k) = 1`) of period at most `log k` for `T_k(x) = (3x + k)/2` or `x/2`. So the number of cycles of `3x + d` is unbounded in `d`.

That is the arithmetic opposite of the Hilbert-16 uniform bound: the parameter enters through divisibility, not analytically.

The elementary construction we first gave is weaker:
* `d_R = lcm(Δ_0, …, Δ_(R−1))`, with `Δ_r = 2^(r+e_r) − 3^(r+1)`, carries at least `R` one-run cycles, through `x_r = d_R(3^(r+1) − 2^(r+1))/Δ_r`.
* They are **non-primitive**: scaled copies of cycles of `3x + Δ_r`, since `T_(md)(mx) = m·T_d(x)`. Script (H) checks `gcd(x_r, d_R) > 1`.
* Its earlier cycle counts were wrong. `3x + 35` has at least nine positive cycles, with minima `7, 13, 17, 25, 35, 133, 161, 1309, 2429`, two of them primitive (13, 17).

## 4. The additive indecomposability of the primes (claims; DICTIONARY)

**Claims of the preprint.**
* Ostmann's conjecture: no finite modification of the primes is `A + B` with `|A|, |B| ≥ 2`.
* After Laffer–Mann, both summands are infinite. For each prime `p`, the images of `A` and `−B` in `F_p` are disjoint, and completing them gives a partition `S_p ⊔ S_p^c`; a collision estimate makes the two parts about half each.
* Section 3 rules out translated quadratic characters (Poisson summation, the quadratic large sieve). Section 4 rules out higher-order characters (Cauchy–Schwarz transfers). Section 5 handles prime coverage via a tensor `L¹` estimate. Sections 6–9 run a positive statistic on transfer trees.
* No Green–Harper classification is claimed.

**Typed links.**
* **Local model (PROVED, trivial).** A pair with `A_p ∩ (−B_p) = ∅` and `|A_p| + |B_p| = p` has `B_p = −(F_p ∖ A_p)`, and then `A_p + B_p = F_p^×`, since a proper nonempty `A_p` is never translation-invariant.
  * For `p ≡ 3 (mod 4)` the most symmetric choice is the Paley pair `(QR, QR ∪ {0})`.
  * The translated-character configurations the preprint must kill (`S_p = c_p + QR`) are translated Paley out-sets, the local shadows of `A ≈ c + {squares}`.
  * This is the quadratic, Heegner side of S17 (`{1,2,4} = QR_7`) and the Paley thread (THM-4526). ANALOGY.
* **Contrast** (DICTIONARY). The Collatz backward sieve (THM-4554) is deep at one prime: `M_a` by `a mod 2·3^(k−1)`, and `S` by 2-adic classes. The preprint's sieve is shallow at all primes.

## 5. Witnessed symmetric choice versus CPT (claims; DICTIONARY + FINITE-EXACT)

**Claims of the preprint.**
* A `CPT + WSC` sentence (eight binary relations, one WSC occurrence) defines a query no CPT sentence defines.
* The structures:
  * a box `{0,…,L−1}³`;
  * edge states in a Heisenberg group over `F_3`;
  * vertex configurations with prescribed divergence;
  * total charge 0 against 1.
* Witnessed choices use a translation on at most four faces and a circulation around a tree cycle. The lower bound uses short `K`-supports for the central subgroup and a counting transfer.

**Links.**
* **The cycle charge** (DICTIONARY). A cyclic Collatz word `w` (`l` odd steps, total `A`) gives an integer cycle iff `c_w ≡ 0 (mod 2^A − 3^l)`. Every proper window is realized by integers (Terras).
  * Local consistency everywhere, with a global charge: the CFI and WSC shape.
  * Here the base graph is a cycle, so one word's charge is a single walk ("spanning tree plus one non-tree edge"). The difficulty of Collatz is the infinitude of words.
  * This is the cocycle reading of the S15 twelfth note, `H¹(⟨w⟩; Z) = Z/(2^A − 3^l)`.
* **Paley tournaments and witnessed choice** (FINITE-EXACT, script (P)).
  * `Aut(P_p) = {x ↦ ax + b : a ∈ QR}` is arc-transitive. Choosing a vertex and then an out-neighbour are **two** witnessed symmetric choices.
  * For all 155 primes `p ≡ 3 (mod 4)` below 2000, colour refinement after individualizing the arc `0 → 1` is discrete. One individualized vertex leaves 3 classes, the stabilizer's orbits.
* **Round 2 reads the Legendre family** (PROVED identity; script (L) checks it for `p ≤ 47`).
  * The round-1 class of `x` is `(χ(x), χ(x−1))`.
  * Round 2 counts `N_(ε,δ)(x) = #{y : χ(y) = ε, χ(y−1) = δ, χ(y−x) = 1}`, and `8N_(ε,δ)(x) = p − εδ − ε − δ + εδ·S(x) + O(1)`.
  * Here `S(x) = Σ_y χ(y(y−1)(y−x)) = −a_p(E_x)` for `E_x : Y² = X(X−1)(X−x)`, and the `O(1)` term (from `y ∈ {0, 1, x}`) is constant on round-1 classes.
  * So within a class, round 2 separates `x` exactly by the Frobenius trace of the Legendre curve: `O(√p)` values on about `p/4` points. Later rounds carry the discreteness, and a proof for all `p` is an elliptic-curve question (OPEN).
* **Doubly regular tournaments** of equal order are 2-WL-equivalent (classical: the rank-3 skew scheme, intersection numbers fixed by `n`). Like strongly regular graphs, they defeat 2-WL; this is weaker than CFI-type hardness.

## 6. Hostile checks

* **Interval arithmetic.**
  * mpmath `iv` with outward rounding.
  * `sinc`, `φ` and `ψ` are enclosed by monotonicity on the ranges used.
  * The scaled maps agree with direct evaluation in 48 cases.
  * The tube constants are exact decimals (interval-enclosed).
  * The audit re-derived the maps by hand and re-checked the inequalities with independent Taylor enclosures and direct interval evaluation for `m ≤ 150`.
* **The continuum `h ∈ [0, 1/3]`** covers every integer: the identities are exact at `h = 1/(2m+1)`, and the inequalities hold on the whole box.
* **Identification of `c_n`.** `F″ < 0` on odd tubes gives uniqueness. Identification with Dumont–Reiter's `c_n` uses their Theorem 6 (CITED; proof sketched there).
* **Conjecture 1 (i)** is not proved. It is proved equivalent to Collatz.
* **Singer.** The hypotheses are listed in 2.6 (e). Only odd critical points are controlled.
* **The Haar model** is exact for the 2-adic question. For actual exponents it is a model: 1.7 shows where it fits and where it misses.
* **Exponents** come from fits over `[400, 19 000]`. The two runs give `α ≈ 0.68` and `0.78`, quoted as `0.7 ± 0.1`.
* **The debt identities** are algebra, checked exactly. The kicked law applies only to integer debts.

## 7. Audit (independent, adversarial; 2026-10-07)

The audit:
* reran all scripts, which reproduce their outputs exactly, and reran the census with small arguments;
* reproduced the K = 27–31 certified fractions;
* fetched Dumont–Reiter, Chamberland's survey, Lagarias's bibliography, Lygerōs–Rozier and the four preprint sources.

It found no false mathematics in the tube theorem. Its findings and our responses (all applied; MISTAKE-580):

1. **MAJOR, missed prior art.** Lygerōs–Rozier 2014 (Lemma 2.4, Theorem 3.3, Corollary 3.4, (5.2)–(5.3), the census) contain the Chamberland half. *Applied:* credits throughout. For `C` our part is a rigorous re-proof plus the new items listed in §0. No prior proof of the Dumont–Reiter conjecture was found, by us or by the audit.
2. **Critical as stated: Theorem D (b) "merge iff `D = 0`".** It fails for arbitrary relations (`x = 5, y = 1`). *Applied:* restated (`D = 0` iff the next relation is the identity; other merges have Haar measure 0).
3. **MAJOR: Monte Carlo constants not reproduced by committed code.** *Applied:* the script is parameterized with bootstrap and lag shares, both runs are committed, and constants are given as ranges.
4. **MAJOR: the level-count exponent `1 − α` does not follow.** *Applied:* typed HEURISTIC.
5. **MAJOR: Proposition H's cycle counts false, and the construction non-primitive.** *Applied:* now Remark H, citing Lagarias 1990.
6. **Smaller findings, all applied.**
   * The "real sense" qualifier.
   * Interval certification of C1.
   * `W_C` fixed points typed NUMERICAL, with the separatrix's existence by the IVT.
   * The Singer citation (Kawasaki K(c), not LSW).
   * MISTAKE-578 propagated.
   * The off-by-one in 1.1.
   * Credit to mac-mini's HYP-9214 model.
   * The scope of the kicked law, and an unverified tie remark removed.
   * The local drift near the merge region.
   * Consistent lag-1 constants, and an uncomputed "about 0.90" dropped.
   * HYP-9217 restricted to the lag-1 Mersenne debt.
   * The K = 28–31 output committed (with K ≤ 27).
   * Runtimes and dead code.
   * Paley wording ("two witnessed choices") and the DRT "CFI hardness" overstatement.
   * `D″ < 0` checked at `a = 0.7`.
   * The even critical points sent to `A2`.

Added after the audit, not re-audited (all mechanical or FINITE-EXACT):
* the flip lemma (2.6 (f));
* the exponent-stream statistics (1.5);
* the Legendre round-2 identity and the Paley range to 2000 (§5);
* the interval root isolation of C1.

## 8. Next steps

1. **Even critical points.** Classify their fates by the parity sequence of `T(m)`. By the flip lemma, an odd-step run that ends with the left deviation in `[−1.3, −0.45]` captures the point. Escape needs the run to overshoot.
2. **A rigorous HYP-9217 (1).** Prove the `T^(−1/2)` first-passage law for the kicked debt process, or for a model in which the debt's exponents are fresh `Geom(1/2)` draws.
3. **Paley canonization for all `p`.** Joint distribution of the Legendre traces `a_p(E_x)` and their higher-round analogues.
4. **More openai/math.** The CPT noncapture companion over `F_3` and the Weisfeiler–Leman complexity papers (#133), against the repo's tournament automorphism results (THM-4557).

## 9. References

* openai/math preprints (AI-written, unrefereed; read read-only from https://github.com/openai/math):
  * `preprints/Witnessed-symmetric-choice-is-strictly-stronger-than-choiceless-polynomial-time-with-counting-September-24-2026`;
  * `preprints/the-additive-indecomposability-of-the-primes-September-24-2026`;
  * `preprints/two-limit-cycles-for-quintic-lienard-systems-September-24-2026`;
  * `preprints/uniform-bounds-for-planar-polynomial-limit-cycles-September-24-2026`.
* N. Lygerōs, O. Rozier, Dynamique du problème 3x+1 sur la droite réelle, Ratio Mathematica 26 (2014) 77–94, arXiv:1402.1979 (read).
* M. Chamberland, A continuous extension of the 3x+1 problem to the real line, Dynam. Contin. Discrete Impuls. Systems 2 (1996) 495–509 (via J. C. Lagarias, *The 3x+1 problem: an annotated bibliography (1963–1999)*, arXiv:math/0309224, entry 35, and M. Chamberland, *An update on the 3x+1 problem*, Butl. Soc. Catalana Mat. 18 (2003) 19–45, §6.5, read).
* J. P. Dumont, C. A. Reiter, Real dynamics of a 3-power extension of the 3x+1 function, Dyn. Contin. Discrete Impuls. Syst. Ser. A 10 (2003) 875–893 (preprint read: webbox.lafayette.edu/~reiterc/3x+1/w3x+1_pp.pdf).
* S. Letherman, D. Schleicher, R. Wood, The 3n+1 problem and holomorphic dynamics, Experiment. Math. 8 (1999) 241–251 (via Lagarias's entry 115).
* D. Singer, Stable orbits and bifurcation of maps of the interval, SIAM J. Appl. Math. 35 (1978) 260–267 (standard; primary not re-read).
* J. C. Lagarias, The set of rational cycles for the 3x+1 problem, Acta Arith. 56 (1990) 33–53 (via Lagarias's bibliography, entry 108).
* Repo:
  * THM-4554, THM-4555, THM-4556 (mac-mini);
  * HYP-9213, HYP-9214 (mac-mini's model);
  * THM-4526 (Paley);
  * the S18 note `mersenne_switch_parity_f21_compression_20261006.md`;
  * the S15 note `collatz_one_out_edge_systems_20260930.md`;
  * the Kawasaki audit note `collatz_procgen_20260924_fixed_points_kawasaki_audit.md` (Theorem K(b), K(c)).
