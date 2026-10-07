# Cycles, tubes and debts: the Dumont–Reiter odd-critical-point conjecture, Chamberland's two cycles, the Mersenne debt walk, and four openai/math preprints (Liénard, Hilbert 16, Ostmann, witnessed choice)

Session opus-2026-10-06-S19. Owner prompt: "work the next steps as they emerge, also mix in and extend ideas from more papers": witnessed symmetric choice vs CPT, the additive indecomposability of the primes, two limit cycles for quintic Liénard systems, and uniform bounds for planar polynomial limit cycles (openai/math, 24 September 2026), "etc at your discretion, let your mathematical hypotheses grow freely into many proofs".

Scripts (each prints ALL CHECKS PASSED):
* `04-computation/experiments/chamberland_tubes_dumont_reiter_20261006.py` (+ `.out`; about 25 s; mpmath interval arithmetic) for section 2.
* `04-computation/experiments/mersenne_haar_debt_walk_20261006.py` (+ `.out`; about 100 s on 12 cores) for section 1.
* `04-computation/experiments/collatz_papers_misc_20261006.py` (+ `.out`) for sections 1.8, 3 and 5.

The certified densities of 1.2 come from the committed S18 script `04-computation/experiments/mersenne_certified_density_numba_20261006.py`, run with arguments `28 31`.

Status words follow the repo: PROVED, PROVED (computer-assisted, interval arithmetic), FINITE-EXACT, NUMERICAL, HEURISTIC, ANALOGY, DICTIONARY, CITED. The openai/math preprints are AI-written and unrefereed; their results are reported as claims.

## 0. Answers

**Tubes: a theorem about real extensions of the Collatz map (new as far as we found; THM-4563).**
* **Theorem T (PROVED, computer-assisted).** Let `T` be the Terras map. Two real extensions are studied:
  * Chamberland's `C(x) = x + 1/4 − (2x+1)/4 · cos(πx)`;
  * Dumont–Reiter's 3-power extension `D(x) = (3^s x + s)/2`, with `s = sin²(πx/2)`.
* For each, every **right tube** `τ(m) = [m, m + a/(2m+1)]`, for an integer `m ≥ 1`, is mapped into `τ(T(m))`. The constant is `a = 0.8` (and `0.9`) for `C`, and `a = 0.6` (and `0.7`) for `D`.
* The local maximum `c_n` of the map next to an odd integer `n` lies in `τ(n)`.
* So the critical orbit **shadows the Collatz orbit of `n` forever**: `0 ≤ F^k(c_n) − T^k(n) ≤ a/(2T^k(n)+1)` for all `k`.
* The proof reduces all integers at once to inequalities for two explicit analytic functions of `(ξ, h)`, with `h = 1/(2m+1)` ranging over the continuum `[0, 1/3]`. Interval arithmetic then checks them.

**Dumont–Reiter's Odd Critical Point Conjecture (2003) is settled as far as it can be.**
* Parts (ii) and (iii) are **PROVED for every odd `n`**:
  * `c_n` and `n` have the same total stopping time;
  * `c_n` lies in `n`'s immediate basin of that stopping time, in the real sense.
* Part (i), "`c_n` is attracted to `(1,2)`", is **PROVED equivalent to the Collatz orbit of `n` reaching 1**. So (i) for all odd `n` is exactly the 3x+1 conjecture.
* For Chamberland's map the analogue is sharper:
  * for every odd `n ≥ 7` whose orbit reaches 1, `c_n` is attracted to `A1 = {1,2}` (PROVED, computer-assisted);
  * `c_1`, `c_3` and `c_5` go to `A2 = {1.19253…, 2.13866…}`;
  * hence **Collatz ⟺ every odd critical point `c_n`, `n ≥ 7`, of Chamberland's map is attracted to `{1,2}`**.

**Limit cycles in the tubes.**
* On `τ(1)` the scaled return map of `C` has three fixed points (computer-assisted):
  * `0`, which is `A1`;
  * a repelling separatrix at `ξ = 0.0710584`;
  * `ξ = 0.5775957`, which is `A2`.
* So **Chamberland's second attracting cycle `A2` is the satellite of the trivial cycle inside its own tubes**.
* By Singer's theorem and the tubes, an attracting cycle whose immediate basin contains an odd critical point is either `A1`, `A2`, or a cycle living in the tubes of a nontrivial positive integer cycle (PROVED).
* Even critical points can leave the shadow (`c_54`). All even ones up to 300 still end at `A1` (NUMERICAL).

**Other results on Chamberland's map.**
* **Mirror law (PROVED).** The displacement `C(x) − x = 1/4 − (u/2) sin(πu)`, `u = x + 1/2`, is even in `u`. So mirror fixed points have multipliers summing to 2.
* **Attracting fixed points (PROVED).** They are exactly 0 and `−1.2777338` (multiplier `0.3857`). `0.2777338` is repelling (`1.6143`).
* **Correction.** The S15 note `collatz_one_out_edge_systems_20260930.md` had these two swapped; logged as MISTAKE-578.
* **Negative Schwarzian on `[0, ∞)` (PROVED, re-proving Chamberland's claim).** `2C′C‴ − 3C″² = π²Q(π(x + 1/2))` with an explicit `Q`. The sign fails just left of 0, at `x = −0.0217160`.

**Next steps 1–2 of S18, the Mersenne switch.**
* **Haar model.** By THM-4556 (iv), Mersenne switching becomes a question about Haar-random 2-adic inputs.
* **Certified density.** Certified switching shares (FINITE-EXACT) are `0.1603, 0.1649, 0.1695, 0.1740` at template totals `K = 28, …, 31` (S18 reached 0.1556 at `K = 27`).
* **Monte Carlo** (NUMERICAL, `N = 3000`, totals to `19 924`).
  * The non-switch share decays like `T^(−α)` with `α = 0.68`, 95% bootstrap interval `[0.63, 0.75]`.
  * So `1 − α ≈ 0.32` matches HYP-9213's growth exponent 0.36 for the number of `σ`-levels.
  * The single-partner (lag-1) non-merge share satisfies `√T · q₁(T) → 16.3`, i.e. a `T^(−1/2)` law. mac-mini's `B^(−0.4)` in HYP-9214 is its pre-asymptotic phase.
* **Debt bookkeeping (PROVED).** Write the partner relation `x = 2^L y + Δ` and the debt `D = 3Δ + 1 − 2^L`.
  * The two orbits merge at the next step **iff `D = 0`**, which forces `L = 2k ≠ 0` and `x = 4^k y + (4^k − 1)/3` (a sibling pair).
  * In generic steps, **the debt's odd part runs the Collatz map itself, kicked by a power of two**: `D′_odd = U(D_odd) − 2^(L′−w′)`.
  * `L` moves by the difference between the partner orbit's exponent and the debt orbit's exponent.
* **The walk** (NUMERICAL). `L` has increments of mean `0.000` and variance `4.00`. That is exactly the variance of a difference of two independent `Geom(1/2)` exponents: a recurrent symmetric walk.
* **The model against actual exponents** (NUMERICAL).
  * It predicts the switching share of actual exponents `a ∈ [1001, 2001]`: 0.974 predicted, 0.970 observed.
  * It predicts the `σ`-level counts `23, 37, 54, 69` for `a ≤ 100, 400, 1000, 2001`: predicted `24, 41, 56, 69`.
  * It over-predicts lag-1 merging: 0.863 predicted, 0.824 observed.
* **Conjecture HYP-9217** (the debt recurrence law) states all this. It would give `μ_2(S) = 1` and hence HYP-9213.

**The four preprints** (claims; ANALOGY / DICTIONARY unless marked).
* **Two limit cycles for quintic Liénard systems.** It counts cycles by pairing two half-orbits and controlling a matching function, using a Schwarzian identity for the endpoint correspondence. Chamberland's core interval carries exactly two attracting cycles, `A1` and `A2`, also controlled through a negative Schwarzian (Singer). The tube theorem is our half-orbit statement: only the right side of each integer is invariant.
* **Uniform bounds for planar polynomial limit cycles** (Hilbert 16).
  * The Mersenne run is an exact saddle passage at the adelic fixed point −1, with Dulac exponent `log₂3`. It turns 2-adic depth `a` into 3-adic depth `a`: `g^a(2^a − 1) = 3^a − 1`, PROVED, elementary.
  * In contrast to uniform boundedness, `3x + d` has unboundedly many one-run cycles (Proposition H, PROVED, elementary). For example `3x + 35` has three.
* **Additive indecomposability of the primes** (Ostmann). Its residue partition `A mod p` versus `−B mod p` has the Paley pair `(QR, QR ∪ {0})` as its symmetric local model at `p ≡ 3 (mod 4)`. Its Section 3 kills exactly the translated-Paley (quadratic) configurations.
* **Witnessed symmetric choice versus CPT.** Its `F_3` Heisenberg charge is a CFI-type global obstruction. The Collatz cycle condition `c_w ≡ 0 (mod 2^A − 3^l)` is a charge of the same shape.
  * One witnessed choice of an arc plus colour refinement canonizes every Paley tournament `P_p` with `p ≡ 3 (mod 4)`, `p < 400` (FINITE-EXACT).
  * Fixing one vertex leaves three classes.

Collatz is OPEN.

---

## 1. The Mersenne switch: a Haar model, a debt walk and a recurrence law

### 1.1 The Haar model (PROVED reduction)

Recall THM-4556 (iv).
* For odd `a`, the post-run start of `M_a = 2^a − 1` is `p = 2·3^(a−1) − 1`, and that of `M_(a−D)` is `q_D = 2·3^(a−D−1) − 1`.
* A **switch at lag `D`** is an equal-time merge `U^i(p) = U^(i+D)(q_D)`.
* Its colliding prefix of template total `K` depends only on `3^(a−1) mod 2^(K+1)`, i.e. on `a mod 2^(K−2)`.
* As `a` ranges over odd 2-adic integers, `X = 3^(a−1)` ranges over `1 + 8Z_2` with Haar measure, since `a ↦ 3^a` is a measure-preserving isomorphism of profinite groups up to normalization.

Hence for every `K` the certified share at template total `K` is the Haar measure of an explicit clopen event in `X`. The **switching set** `S ⊂ Z_2` (the open set of exponents that switch at some finite total) has
* `μ_2(S) = lim_K (certified share at K)`;
* `μ_2(S) = 1 ⟹ HYP-9213` (S18, Proposition 6).

Sampling `X` uniformly and following the two orbits 2-adically (with precision bookkeeping, 64 guard bits) estimates `P(switch with total ≤ T)` for `T` far beyond exhaustive certification.

### 1.2 Certified density to `K = 31` (FINITE-EXACT)

The committed S18 numba search, run with `28 31`, gives exact shares among odd classes:

| `K` | 27 (S18) | 28 | 29 | 30 | 31 |
|---|---|---|---|---|---|
| certified share | 0.1556 | 10757550/67108864 = 0.1603 | 22135584/134217728 = 0.1649 | 45495372/268435456 = 0.1695 | 93407400/536870912 = 0.1740 |

The increments are about 0.0046 per unit of `K`, consistent with the Haar curve below. Exhaustive certification cannot reach the median merge total (about 250).

### 1.3 Monte Carlo (NUMERICAL)

The large run used seed 7, `N = 3000` samples, 20 000-bit 2-adic precision and odd lags `D ≤ 61`, reaching totals up to 19 924. The committed script reruns a smaller version (seed 2026, `N = 1200`, 12 000 bits, `D ≤ 41`) with the same picture.

| total `T` | 27 | 100 | 400 | 1600 | 3200 | 6400 | 12800 | 19000 |
|---|---|---|---|---|---|---|---|---|
| non-switch share `q(T)` | 0.850 | 0.615 | 0.311 | 0.117 | 0.0757 | 0.0443 | 0.0290 | 0.0223 |
| lag-1 non-merge `q₁(T)` | 0.859 | 0.722 | 0.548 | 0.356 | 0.276 | 0.202 | 0.143 | 0.118 |
| `√T · q₁(T)` | 4.5 | 7.2 | 11.0 | 14.2 | 15.6 | 16.1 | 16.2 | 16.3 |

* At `T = 27` the Haar switch share is `0.150 ± 0.007`, against the exact certified 0.1556.
* Least squares over `T ≥ 400` gives `q(T) ≈ 18.4·T^(−0.683)`, with 95% bootstrap interval `α ∈ [0.628, 0.748]`.
* The lag-1 tail fit over `T ≥ 3200` gives exponent 0.477, and `√T·q₁` levels off at 16.2–16.3.
* Over `T ≥ 400` the lag-1 exponent is 0.40, the value mac-mini reported for HYP-9214. That value is the pre-asymptotic part of a `T^(−1/2)` law.
* Least-total lags among switching samples are `D = 1` (52%), 3 (17%), 5 (9%), 7 (6%), … .

### 1.4 Theorem D (the debt is a kicked Collatz orbit; PROVED)

**Setting.** Let `x, y` be odd 2-adic integers (or odd integers) with `x = 2^L y + Δ`, where `L ∈ Z` and `Δ ∈ Z[1/2]`. Put `D = 3Δ + 1 − 2^L`. Let `a = v_2(3y + 1)`, `y′ = U(y)`, `b = v_2(3x + 1)`, `x′ = U(x)`.

(a) **Transition.**
* `3x + 1 = 2^(L+a) y′ + D`, so `b = v_2(2^(L+a) y′ + D)`.
* `x′ = 2^(L′) y′ + Δ′`, with `L′ = L + a − b` and `Δ′ = D/2^b`.

(b) **Merge criterion.**
* `x′ = y′` iff `D = 0`. Then `Δ = (2^L − 1)/3`, so `L` is even and nonzero, and `x = 4^k y + (4^k − 1)/3 = R^k(y)` (`k = L/2`), or the mirror with `k < 0`. This is the sibling configuration of HYP-9214.
* In the 2-adic model with `x` an affine function of Haar `y`, a merge at a value-equality without `D = 0` has measure 0.

(c) **Kicked Collatz law.** Suppose `D ≠ 0` is an integer with `v_2(D) < L + a` (so `b = v_2(D)`), and write `D_odd = D/2^(v_2 D)` and `w′ = v_2(3D_odd + 1)`. If `L′ > w′`, then
`v_2(D′) = w′`, `D′/2^(w′) = U(D_odd) − 2^(L′ − w′)`, and `L′ = L + a − v_2(D)`.

So the debt's odd part is a Collatz orbit (of a mostly negative integer), kicked at 2-adic depth `L′ − w′`. The walk `L` moves by `a − w`: the partner's exponent minus the debt's own exponent.

**Proof.**
* (a) Substitute `x = 2^L y + Δ` into `3x + 1` and use `3y + 1 = 2^a y′`.
* (b) `x′ = y′` iff `(2^(L′) − 1) y′ = −Δ′`. In the identity case this needs `L′ = 0` and `Δ′ = 0`, i.e. `D = 0`. Then `3Δ = 2^L − 1`, and `Δ ∈ Z[1/2]` forces `3 | 2^|L| − 1`, i.e. `L` even.
* (c) With `b = v_2(D)`, `D′ = 3D/2^b + 1 − 2^(L′) = 3D_odd + 1 − 2^(L′)`. If `L′ > w′`, then `v_2(D′) = w′` and `D′/2^(w′) = (3D_odd + 1)/2^(w′) − 2^(L′−w′)`. ∎

**Checks** (FINITE-EXACT, `mersenne_haar_debt_walk_20261006.py` B): over 1000 lag-1 pairs `x_0 = 3·2^v y_0 + 1` (the Mersenne debt states of THM-4556 (ii)):
* the transition law holds on all first-60 steps;
* the kicked law holds in all 28 693 generic steps tested;
* all 834 merges pass through `D = 0` with `L` even and nonzero and `Δ = (2^L − 1)/3`.
* `L` one step before a merge is `−2` (489), `+2` (301), `−4` (20), `+4` (11), `±6`, `−8`. About 61% have `k < 0`, matching mac-mini's "about 2/3".

### 1.5 The walk (NUMERICAL)

Over all steps of the 1000 pairs (12 000-bit sources, up to 4800 steps each):

| `L` range | `≤ 0` | 1–6 | 7–15 | `≥ 16` | all |
|---|---|---|---|---|---|
| mean increment | +0.0006 | −0.050 | +0.012 | +0.002 | +0.0004 |
| variance | 3.98 | 4.04 | 4.01 | 4.02 | **4.002** |

* `Var(G − G′) = 2 + 2 = 4` for independent `Geom(1/2)` exponents. The partner's exponent stream and the debt's exponent stream behave like two independent fair streams.
* In steps, `P(no merge by t)` is 0.627, 0.490, 0.324, 0.216 and 0.166 at `t = 100, 300, 1000, 3000, 4800`; `√t · P` is 6.3, 8.5, 10.3, 11.8 and 11.5.

**HEURISTIC mechanism.**
* A mean-zero walk with finite variance is recurrent.
* The debt is pinned to the scale `|D| ≈ 2^L` (by (c), the kick `2^(L′)` dominates). So a return of `L` to `{±2, ±4}` brings the debt back to small integers, where `D = 0` has a fixed positive chance. (`D_odd = −1` is the negative fixed point of `U`; the `k = −1` merges are tie steps from `D_odd = −1`.)
* Then `P(no merge by t) ≈ (expected returns)⁻¹ ≈ c·t^(−1/2)`.
* What is not proved is the independence of the debt's exponent stream: it is a deterministic Collatz orbit randomized only by the kicks. This is precisely mac-mini's open question (per-visit success bounded below).

### 1.6 HYP-9217 (the debt recurrence law)

`05-knowledge/hypotheses/HYP-9217-mersenne-debt-recurrence-law-t-minus-half.md`. In the Haar model:
1. the lag-1 Mersenne debt (and the reset-2 debt of HYP-9214) merges almost surely, with `P(no merge by template total T) ~ c₁ T^(−1/2)`, `c₁ ≈ 16.3`;
2. the any-lag non-switch probability satisfies `q(T) = T^(−α + o(1))` with `α ≈ 0.68` (bootstrap `[0.63, 0.75]`).

Hence `μ_2(S) = 1` and HYP-9213 holds. The number of `σ`-levels among `a ≤ A` is about `A^(1−α)`, consistent with the observed 0.36.

### 1.7 The model against the actual exponents (NUMERICAL)

For `2 ≤ a ≤ 2001` the script computes `σ(M_a)` and the post-run template total `T_post(a)` (about `7.68a`).

* **Any partner.** For odd `a ∈ [1001, 2001]`, the share with a smaller `σ`-partner is 0.970. The Haar prediction `1 − mean q(T_post(a))` is 0.974 (mac-mini's `[10^3, 6·10^3]` census: 0.980).
* **Level counts.** Distinct `σ` values for `a ≤ 100, 400, 1000, 2001` (counting `a = 1`): `23, 37, 54, 69`. Predicted `2 + Σ_(odd a) q(T_post(a))`: `24.2, 41.4, 55.8, 68.6` with the committed script's fit, and `22.0, 39.5, 55.9, 71.7` with the large run's fit (share 0.968).
* **Lag 1.** The share with `σ(M_(a−1))` is 0.824, against a predicted 0.863. The model over-predicts lag-1 merging by about 0.04 (about 2σ). mac-mini's `[10^3, 6·10^3]` figure is 0.882, predicted about 0.90.
  * The likely cause is that actual orbits end: the effective time before the values become small is shorter than `T_post`.
  * The ratio of the two observed lag-1 failure shares, 1.49, matches the `a^(−1/2)` law's 1.43.

### 1.8 The adelic saddle (PROVED, elementary; DICTIONARY with the Hilbert-16 preprint)

Let `g(x) = (3x+1)/2`. Then `g(x) + 1 = (3/2)(x + 1)`, so `g^k(x) + 1 = (3/2)^k (x + 1)` and `g^a(2^a − 1) = 3^a − 1`. This is THM-4556 (i), binary repunit to ternary repunit.

At the fixed point −1 the multiplier `3/2` has
* `|3/2|_2 = 2`: unstable;
* `|3/2|_3 = 1/3`: stable;
* `|3/2|_∞ = 3/2`;
* and the product formula gives volume neutrality.

So the Mersenne run is a **saddle passage**.
* It enters at 2-adic distance `r = 2^(−a)` and leaves at 3-adic distance `3^(−a) = r^λ`, after a flight time of `a` steps.
* The Dulac exponent is `λ = log 3/log 2 = log₂ 3`.

The uniform-bounds preprint's local model `ẋ = x, ẏ = −λy`, with flight time `L = log(1/r)` and contraction exponent `W = λL`, is exactly this with the real and `p`-adic places as the two directions. Its "two clocks that may be comparable or widely separated" are the CRT-independent 2-adic and 3-adic clocks of S18 Proposition 5.

## 2. Chamberland's map and the 3-power map: tubes, shadows, limit cycles

### 2.1 Prior work (CITED)

* **Chamberland 1996** (Dynam. Contin. Discrete Impuls. Systems 2, 495–509; read through Lagarias's annotated bibliography, entry 35, and Chamberland's 2003 survey, §6.5):
  * `C` has negative Schwarzian on `R+`;
  * `[0, μ1)` is attracted to 0;
  * `[μ1, μ3]` is invariant, and almost every point of it goes to `A1` or `A2`;
  * a nontrivial positive cycle of `T` would be attracting.
  * The survey reports the conjecture that `A1`, `A2` are the only attracting cycles on `R+`, and calls it equivalent to the 3x+1 problem. As far as we can prove, it implies only the cycle half (2.6).
* **Dumont–Reiter 2003** (Dyn. Contin. Discrete Impuls. Syst. Ser. A 10, 875–893; preprint read):
  * `D` has negative Schwarzian for `x ≥ 0` (their Theorem 5);
  * critical points `c_n = n − (−1)^n (2/(π² ln 3))/n + O(n⁻²)` (Theorem 6);
  * every point of `(μ1, μ3)` is attracted to `(1,2)` (Theorem 8).
  * They define the **total stopping time** of real `x` as the least `k` with `μ1 < D^k(x) < μ2`.
  * Their **Conjecture 1 (Odd Critical Point Conjecture)**: for odd `n`, (i) `c_n` is attracted to `(1,2)`; (ii) `c_n` and `n` have the same total stopping time; (iii) `c_n` is in the immediate basin of total stopping time equal to the value at `n`.
  * Evidence: all odd `n < 180 000`, and Vyssotsky's `n_V` (2565 steps).
  * They also show (their Table 6) that the even critical point `c_54` does not shadow 54.
* **Letherman–Schleicher–Wood 1999** (Experiment. Math. 8, 241–251) study entire interpolations with superattracting integers. There, no two integers share a Fatou component (except possibly −1, −2), and a divergent orbit corresponds to a wandering domain.

We found no proof of the Odd Critical Point Conjecture in the sources we checked:
* the two papers;
* Lagarias's bibliographies;
* Chamberland's survey;
* web searches on 2026-10-06.

### 2.2 Proposition C1 (mirror law, attracting fixed points; PROVED)

With `u = x + 1/2`, `C(x) − x = 1/4 − (u/2) sin(πu)`, which is even in `u`. So the fixed-point set is symmetric under `x ↦ −1 − x`, and `C′(−1 − x) = 2 − C′(x)`. Mirror multipliers sum to 2, checked at all 41 fixed points in `[−20, 20]`. At most one of each mirror pair is attracting.

At a fixed point, `cos πx = 1/w` with `w = 2x + 1`. Hence `C′ = 1 − 1/(2w) ± (π/4)√(w² − 1)`, and `|C′| < 1` forces `|w| < √(1 + 100/π²) = 3.3365`.

An interval sign scan of that window finds exactly the fixed points `−1.2777338, −1, 0, 0.2777338`, with multipliers `0.3857, 1.5, 0.5, 1.6143`. **The attracting fixed points of `C` are exactly `0` and `−1.2777337662`.** This agrees with Dumont–Reiter's Table 4 (`0.385708` at `−1.27773`).

*Correction (MISTAKE-578).* The S15 note said "the only attracting real fixed point is 0.278 (multiplier 0.386); its mirror −1.278 is repelling (1.614)". It is the other way round. The note and its INDEX line are corrected.

### 2.3 Proposition C2 (negative Schwarzian on `[0, ∞)`; PROVED, computer-assisted)

**Identity.** With `w = π(x + 1/2)`, `s = sin w`, `c = cos w`:

`2C′C‴ − 3C″² = π² Q(w)`, where `Q(w) = −(1/2 + s²/4)w² + c(1 + s)w + (3/2)(s² + 2s − 2)`.

**Proof of `Q(w) < 0` for `x ≥ 0`.**
* For `w > 2√3`: `Q ≤ −w²/2 + (3√3/4)w + 3/2 < 0`.
* On `[π/2, 2√3]` (`x ∈ [0, 0.6027]`): 4000 interval cells give `Q < 0`, with maximum `Q(π/2) = −0.3506`.

So `SC < 0` on `[0, ∞)` (Chamberland's claim, re-proved). The sign changes at `x = −0.0217160`, and `Q > 0` just left of it.

### 2.4 Theorem T (the tube theorem; PROVED, computer-assisted)

**Scaled coordinates.** For `F ∈ {C, D}` and an integer `m ≥ 1`, write `x = m + ξh` with `h = 1/(2m+1)`, and the image coordinate `ξ′ = (2T(m) + 1)(F(x) − T(m))`. Then `ξ′ = O(ξ, h)` for odd `m` and `ξ′ = E(ξ, h)` for even `m`, with:

* `C`, odd `m`: `O = (3+h)/2 · [ξ(1 + cos(πξh)/2) − (π²ξ²/8) S²]`.
* `C`, even `m`: `E = (1+h)/2 · [ξ(1 − cos(πξh)/2) + (π²ξ²/8) S²]`.
* `D`, odd `m`: `O = (3+h)/4 · [−(3(1−h)/2) ln3 (π²ξ²/4) S² φ(σ ln3) + 3^(1−σ) ξ − h(π²ξ²/4) S²]`.
* `D`, even `m`: `E = (1+h)/4 · [((1−h)/2) ln3 (π²ξ²/4) S² ψ(σ ln3) + 3^σ ξ + h(π²ξ²/4) S²]`.

Here `S = sinc(πξh/2)`, `σ = sin²(πξh/2)`, `φ(y) = (1 − e^(−y))/y` and `ψ(y) = (e^y − 1)/y`.

These are exact identities at `h = 1/(2m+1)`. They come from `cos(π(m + ε)) = (−1)^m cos πε`, `sin²(π(m+ε)/2) = cos²` or `sin²(πε/2)`, `2T(m) + 1 = 3m + 2` or `m + 1`, and `1 − cos z = 2 sin²(z/2)`. They are analytic in `(ξ, h)`, including at `h = 0`, where `O → (9/4)ξ − (3π²/16)ξ²` and `E → ξ/4 + (π²/16)ξ²` for `C` (with `ln 3` factors for `D`). The script checks them against direct evaluation in 48 cases.

**Theorem T.**
* For `C` with `a ∈ {0.8, 0.9}`, and for `D` with `a ∈ {0.6, 0.7}`: `0 < O(ξ, h) ≤ a` for `(ξ, h) ∈ (0, a] × [0, 1/3]`, and `0 ≤ E(ξ, h) ≤ a` for `(ξ, h) ∈ [0, a] × [0, 1/5]`.
* Hence `F(τ(m)) ⊂ τ(T(m))` for every integer `m ≥ 1`, where `τ(m) = [m, m + a/(2m+1)]`.
* For odd `n`, `F′(n) = 3/2 > 0` and `F′(n + a/(2n+1)) < 0`, and `F″ < 0` on `τ(n)` (analytic for `C`, interval-checked for `D`). So `F` has a unique critical point `c_n` in `τ(n)`, the local maximum. For `D` it is Dumont–Reiter's `c_n`.

**Proof.** The inequalities are verified on boxes covering the stated rectangles (1200 boxes per constant), with outward-rounded interval arithmetic. `sinc`, `φ` and `ψ` are enclosed by monotonicity. Because `h` runs over the continuum `[0, 1/3] ∋ 1/(2m+1)`, all integers are covered at once. Induction on `k` gives `F^k(τ(n)) ⊂ τ(T^k(n))`. ∎

**Numbers.**
* For `C`, odd critical points have `ξ_0 = (2n+1)(c_n − n) → 6/π² = 0.6079`, and the image of the peak is `6.75/π² = 0.6839`.
* The largest scaled deviation along all orbits of odd `n ≤ 401` is 0.7037 for `C` and 0.4559 for `D`.
* In the limit, the maximal invariant scaled interval for `C` is `[0, 12/π²)` (`E` has a repelling fixed point at `12/π²`). At `m = 2`, `a = 0.95` already fails.

### 2.5 Corollary DR (the Odd Critical Point Conjecture; PROVED parts)

For the 3-power map `D` and every odd `n ≥ 1`:

* **(ii)** The total stopping time of `c_n` equals that of `n`, finite or infinite.
  * By Theorem T, `D^k(c_n) ∈ τ(T^k(n))`.
  * `τ(1) = [1, 1.2] ⊂ (μ1, μ2) = (0.3158162, 1.5155526)`, while `τ(m) ⊂ [2, ∞)` for `m ≥ 2`.
  * So `D^k(c_n) ∈ (μ1, μ2)` iff `T^k(n) = 1`.
* **(iii)** `τ(n)` is connected, contains `n` and `c_n`, and every point of it has total stopping time `σ(n)`, by the same argument. So `c_n` lies in `n`'s immediate basin of total stopping time, in the real sense: the component of the real level set containing `n`.
  * Dumont–Reiter illustrate (iii) with complex-plane pictures. If their complex basin is a superset of the real level set near `n`, the real segment `[n, c_n]` puts both points in one component there too. We prove only the real statement.
* **(i) ⟺ the Collatz orbit of `n` reaches 1.**
  * On `τ(1)`, `D²` in scaled form is `W_D = E(·, 1/5) ∘ O(·, 1/3)`. The script shows `W_D(ξ)/ξ < 1` on `(0, 0.6]` (interval arithmetic, 300 cells), so every point of `τ(1)` converges to `(1,2)`. Hence if `n → 1`, then `c_n → (1,2)`.
  * Conversely, if `n` never reaches 1, the orbit of `c_n` stays in tubes of integers `≥ 2`, each of length `< 0.2`. Then `c_n` cannot converge to `{1, 2}`.

So **Conjecture 1 (ii)–(iii) are theorems, and Conjecture 1 (i) for all odd `n` is equivalent to the 3x+1 conjecture** (every odd `n` reaching 1). Dumont–Reiter's remark that the conjecture implies there are no nontrivial cycles is the easy half of this equivalence.

### 2.6 Corollary C (Chamberland's map: where the odd critical points go; PROVED, computer-assisted)

**The return map on `τ(1)`.** It is `W_C = E(·, 1/5) ∘ O(·, 1/3)` on `[0, 0.8]`.
* It is increasing on `[0, ξ_(c1)]`, with `ξ_(c1) = 0.5428161`, and decreasing after.
* Its fixed points are `0` (`A1`, `W′ = 3/4`), `ξ_R = 0.0710584` (repelling, `W′ = 1.2148`) and `0.5775957` (`A2`, `W′ = −0.2308`).
* `W_C(ξ) < ξ` on `(0, 0.0705]` (interval), so `[0, 0.0705]` lies in `A1`'s basin.
* `W_C` maps `J = [0.55, 0.60]` into `[0.5702, 0.5820]` with `|W′| ≤ 0.377`, so `J` lies in `A2`'s basin.

So **`A2` is the satellite of `A1` inside the tubes `τ(1) ∪ τ(2)`**. Dumont–Reiter's Table 4 calls it "a cycle shadowing the cycle (1,2)".

**Odd `n ≥ 7` go to `A1`.**
* The backward tree of 1 under `T` reads `20 ← {13, 40}`, `32 ← {21, 64}`, `5 ← {10, 3}`, `8 ← {5, 16}`.
* Odd `n ≥ 5` never meet a multiple of 3 after the start, since `T(x) = 3·2^k` forces `x = 3·2^(k+1)`.
* So every odd `n ≥ 7` whose orbit reaches 1 passes through `13, 21, 40` or `64` (also checked for `n ≤ 200 001`).
* Pushing the whole tube `τ(m)`, for these four nodes, down the fixed tails with interval arithmetic (400 pieces each) lands in `τ(1)` at `ξ ≤ 0.0572` (from 13), `0.0291` (40), `0.0138` (21) and `0.0076` (64), all `< 0.0705`.

**Therefore.**
* (a) For every odd `n ≥ 7` whose Collatz orbit reaches 1, `c_n → A1`.
* (b) `c_1`, `c_3` and `c_5` → `A2` (interval enclosures of `c_n` iterated into `J`).
* (c) **The 3x+1 conjecture holds iff every odd critical point `c_n`, `n ≥ 7`, of Chamberland's map is attracted to `{1,2}`.** If `n` never reaches 1, `c_n` stays in tubes of integers `≥ 2`. If `n` diverges, `c_n → ∞`. If `n` enters a nontrivial cycle `Γ`, then `c_n` accumulates in `Γ`'s tubes.

**(d) Singer bound (PROVED, using Singer 1978).**
* Let `Γ*` be an attracting or neutral cycle of `C` in `(0, ∞)` other than `A1`, `A2`, with `[0, μ1)` the basin of 0.
* Singer's theorem applies to `C` on `[0, ∞)`: `C` is `C³` with `SC < 0` (C2), `C([0, ∞)) ⊂ [0, ∞)`, and no immediate basin is unbounded, because the monotone divergent orbits of Chamberland and LSW start at arbitrarily large points.
* So `Γ*`'s immediate basin contains a critical point.
* If it is an odd one `c_n`: `n` cannot reach 1 (by (a), (b)) and cannot diverge. So `n` enters a nontrivial positive cycle `Γ`, and `Γ*` lies in the closed tubes of `Γ`.
* Hence **every attracting cycle other than `A1`, `A2` that captures an odd critical point is a nontrivial integer cycle or a satellite in its tubes.**

**Even critical points** (`c_m` just left of even `m`, `ξ_0 → −2/π²`) are not controlled.
* The left side of an odd integer expands: `O(ξ) < (9/4)ξ` for `ξ < 0`.
* `c_54` leaves the near-integer regime at step 9, at the integer 242, which is Dumont–Reiter's Table 6 phenomenon in `C`.
* All even critical points up to 300 end at `A1` for both maps (NUMERICAL, 1500 digits). 7 of 150 (`C`) and 32 of 150 (`D`) leave the regime `|ξ| ≤ 1.2` on the way.
* Chamberland's two-cycle conjecture is therefore equivalent, by (d), to: no nontrivial positive cycle, and no even critical point captured by a spurious attracting cycle. It does not touch divergence.

### 2.7 Rhymes with the two limit-cycle preprints (ANALOGY / DICTIONARY)

* **Liénard (claims of the preprint).**
  * `ẋ = y − F(x), ẏ = −x`, `deg F ≤ 5`: at most two limit cycles, sharp, e.g. `F = ε(4x − 20x³/3 + 8x⁵/5)`, averaged displacement `−πs²(s² − 1)(s² − 4)`.
  * The proof:
    * keeps the two half-orbits separate (`u = x²/2`, profiles `F(±√(2u))`, even coefficients shared, odd ones flipped);
    * fits quadratic profiles;
    * controls a width coefficient by a Riccati linearization and a Schwarzian identity for the endpoint correspondence;
    * bounds the zeros of a matching function whose derivative, after a positive integrating factor, has the sign of a monotone function.
  * The parallels here:
    * Chamberland's core interval has **exactly two** attracting cycles, `A1` and `A2`, also by Schwarzian control (Singer).
    * The mirror law of C1 is the even/odd split of the displacement.
    * The tube theorem is a half-orbit statement: only the right side of each integer is invariant; left sides of odd integers expand.
  * The coincidence that both "two" are two is NUMEROLOGY. The shared mechanism (projective, cross-ratio control of fixed-point counts) is ANALOGY.
* **Uniform bounds (claims of the preprint).**
  * Return maps are matching systems of passages, with flags of scales and exact logarithmic relations. Distinct hyperbolic cycles occupy distinct components of a fibre, and a projection-counting argument gives `B(d)`.
  * Here Singer's critical-point count plays the counting role, and the tube theorem shows the odd critical points are "charged" exactly by the integers: the count reduces back to Collatz.
  * Proposition H (section 3) shows the opposite of uniform boundedness in the `3x + d` family. Section 1.8 types the saddle passage.

## 3. Proposition H: no uniform bound for `3x + d` (PROVED, elementary)

For `r ≥ 0`, let `e_r` be least with `2^(r + e_r) > 3^(r+1)`, and put `Δ_r = 2^(r + e_r) − 3^(r+1)` (odd, prime to 3). Let `d_R = lcm(Δ_0, …, Δ_(R−1))`.

Then `T_(d_R)(x) = x/2` or `(3x + d_R)/2` has the `R` distinct positive cycles through `x_r = d_R(3^(r+1) − 2^(r+1))/Δ_r`, `r < R`. The cycle through `x_r` has one run of `r` odd steps of exponent 1, then a reset of exponent `e_r`, and period `r + e_r`.

**Proof.** The cycle point of the word `(1^r, e_r)` is `d·c_w/(2^(r+e_r) − 3^(r+1))` with `c_w = Σ_(i ≤ r) 3^(r−i) 2^i = 3^(r+1) − 2^(r+1)`. With `Δ_r | d`, all its rotations are odd positive integers (each rotation's `c` is odd and `d/Δ_r` is odd). The odd cycle points have the right 2-adic valuations because each image is odd. ∎

Examples:
* `d = 35` (`Δ = 1, 7, 5`): cycles of periods 2, 4 and 5.
* `d = 21385`: five cycles.
* `d = 2 408 613 935`: seven.
* `d = 1 468 678 841 619 535`: nine.

The unboundedness of the number of `3x + d` cycles is known in that literature (Lagarias 1990; Belaga–Mignotte 1998); we did not check whether this one-run construction appears there. It is the arithmetic opposite of the Hilbert-16 uniform bound: the parameter enters through divisibility, not analytically.

## 4. The additive indecomposability of the primes (claims; DICTIONARY)

**Claims of the preprint.**
* Ostmann's conjecture: no finite modification of the primes is `A + B` with `|A|, |B| ≥ 2`.
* After Laffer–Mann, both summands are infinite. For each prime `p`, the images of `A` and `D = −B` in `F_p` are disjoint, and completing them gives a partition `S_p ⊔ S_p^c`. A collision estimate makes the two parts about half each on prime averages.
* Section 3 rules out persistent correlation with translated quadratic characters (Poisson summation forces a common rational centre, then the quadratic large sieve). Section 4 handles higher-order characters by Cauchy–Schwarz transfers. Section 5 adds prime coverage via a tensor `L¹` estimate. Sections 6–9 build a positive statistic on binary transfer trees and contradict it by averaging over permutations.
* It does not prove the Green–Harper inverse-sieve conjecture.

**Typed links.**
* **Local model (PROVED, trivial).** For `p` prime, a pair with `A_p ∩ (−B_p) = ∅` and `|A_p| + |B_p| = p` is `B_p = −(F_p ∖ A_p)`. Every such pair has `A_p + B_p = F_p^×`, because `A_p + d = A_p` is impossible for a proper nonempty `A_p` and `d ≠ 0`. For `p ≡ 3 (mod 4)` the most symmetric choice is the **Paley pair** `(QR, QR ∪ {0})`. The translated-character configurations the preprint must kill (`S_p = c_p + QR`) are translated Paley out-sets, the local shadows of `A ≈ c + {squares}`.
* This is the quadratic, Heegner side of S17 (`{1,2,4} = QR_7`) and the Paley thread (THM-4526), here as the extremal case of a half-residue large sieve. ANALOGY.
* **Contrast** (DICTIONARY). The Collatz backward sieve (THM-4554) is a deep sieve at one prime: `M_a` is classified by `a mod 2·3^(k−1)`, and the switching set `S` by 2-adic classes. The preprint's sieve is shallow at all primes. Neither transfers to the other as stated.

## 5. Witnessed symmetric choice versus CPT (claims; DICTIONARY + FINITE-EXACT)

**Claims of the preprint.**
* A `CPT + WSC` sentence over eight binary relations with one WSC occurrence defines a Boolean query that no CPT sentence defines.
* The structures:
  * a box `{0,…,L−1}³`;
  * edge states in a Heisenberg group over `F_3`, `(u, z)(u′, z′) = (u + u′, z + z′ + B_e(u, u′))`;
  * vertex configurations with a prescribed divergence `Σ ε_ve z_e = b_v`.
* Zero total charge has a global consistent choice and charge one has none. Witnessed choices use a translation on at most four faces and a circulation around a tree cycle.
* The CPT lower bound uses short `K`-supports for the central subgroup and a counting transfer.

**Links.**
* **The cycle charge** (DICTIONARY). For a cyclic Collatz word `w` (`l` odd steps, total `A`), the cycle point is `x_w = c_w/(2^A − 3^l)`. It is an integer cycle iff the charge `c_w mod (2^A − 3^l)` vanishes. Every proper window is realized by integers (Terras: each parity vector of length `k` is a residue class mod `2^k`).
  * As in CFI and WSC, local consistency holds everywhere and the obstruction is a global charge.
  * Unlike the grids, the base graph here is a cycle, so the charge of a fixed word is computed by one walk around it (the WSC sentence's "spanning tree plus one non-tree edge"). The difficulty of Collatz is the infinitude of words, not the charge of one.
  * This repeats the cocycle reading of the S15 twelfth note, `H¹(⟨w⟩; Z) = Z/(2^A − 3^l)`.
* **Paley tournaments and witnessed choice** (FINITE-EXACT, `collatz_papers_misc_20261006.py` (P)).
  * `Aut(P_p) = {x ↦ ax + b : a ∈ QR}` is arc-transitive, so choosing a vertex and then an out-neighbour are two witnessed symmetric choices.
  * For all 40 primes `p ≡ 3 (mod 4)`, `p < 400`, colour refinement after individualizing the arc `0 → 1` is **discrete**. After individualizing one vertex it leaves 3 classes, the stabilizer's orbits.
  * So CPT + WSC canonizes these tournaments with two choices and counting. A proof for all `p` would need character-sum control of the refinement rounds; it is OPEN.
* **Doubly regular tournaments** of the same order are 2-WL-equivalent, since they all carry the rank-3 skew scheme with intersection numbers fixed by `n` (classical). Non-isomorphic DRTs of equal order exist for many orders, via inequivalent skew Hadamard matrices (Reid–Brown); we did not re-count them. This is the tournament analogue of CFI hardness.

## 6. Hostile checks

* **Interval arithmetic.** Every enclosure uses mpmath `iv` (outward rounding). `sinc`, `φ` and `ψ` are enclosed by monotonicity on `[0, 1]`, where they are monotone. The scaled maps agree with direct evaluation in 48 cases. The tube checks use only these maps.
* **The continuum `h ∈ [0, 1/3]` covers every integer.** It does: the identities are exact at `h = 1/(2m+1)`, and the inequalities are proved on the whole box. No large-`m` asymptotics are needed.
* **Identification of `c_n`.** `F″ < 0` on the odd tube gives uniqueness. For `D`, the tube lies inside `(μ_n, μ_(n+1))`, so this is Dumont–Reiter's `c_n`.
* **Conjecture 1 (i).** We do not claim to prove it. We prove it equivalent to Collatz, so any claimed proof of (i) is a proof of Collatz.
* **Singer on `[0, ∞)`.** The hypotheses are listed in 2.6 (d). Only odd critical points are controlled. The census of even ones is NUMERICAL.
* **The Haar model** is exact for the 2-adic question (THM-4556 (iv)). For actual exponents it is a model: actual orbits end. 1.7 shows where it fits (any-lag share, level counts) and where it misses (lag 1, by 0.04).
* **Exponents.** `α = 0.68` comes from a power-law fit over `T ∈ [400, 19 000]`, with bootstrap interval `[0.63, 0.75]`. No asymptotic claim beyond HYP-9217.
* **The debt identities** are algebra, checked exactly (Fractions) on 28 693 generic steps and all 834 merges.

## 7. Audit

To be filled from the independent adversarial audit (MISTAKE ledger if needed).

## 8. Next steps

1. Prove the tube theorem for the even side in a weaker form. Find the set of even `m` whose critical point lands in a right tube (every even step with `ξ ∈ [−12/π², −4/π²]` throws the point into a right tube). Then estimate the density of escapes.
2. A rigorous version of HYP-9217 (i) for a model in which the debt's exponents are replaced by fresh `Geom(1/2)` draws: the first-passage law for the actual kicked process.
3. Character-sum proof that one witnessed arc plus colour refinement canonizes every Paley tournament.
4. The other two openai/math CPT papers (the noncapture companion over `F_3`; the Weisfeiler–Leman complexity papers 133) against the repo's tournament automorphism results (THM-4557, Aut = F_21 towers).

## 9. References

* openai/math preprints (AI-written, unrefereed; read read-only from https://github.com/openai/math):
  * `preprints/Witnessed-symmetric-choice-is-strictly-stronger-than-choiceless-polynomial-time-with-counting-September-24-2026`;
  * `preprints/the-additive-indecomposability-of-the-primes-September-24-2026`;
  * `preprints/two-limit-cycles-for-quintic-lienard-systems-September-24-2026`;
  * `preprints/uniform-bounds-for-planar-polynomial-limit-cycles-September-24-2026`.
* M. Chamberland, A continuous extension of the 3x+1 problem to the real line, Dynam. Contin. Discrete Impuls. Systems 2 (1996) 495–509 (via J. C. Lagarias, *The 3x+1 problem: an annotated bibliography (1963–1999)*, arXiv:math/0309224, entry 35, and M. Chamberland, *An update on the 3x+1 problem*, Butl. Soc. Catalana Mat. 18 (2003) 19–45, §6.5, read).
* J. P. Dumont, C. A. Reiter, Real dynamics of a 3-power extension of the 3x+1 function, Dyn. Contin. Discrete Impuls. Syst. Ser. A 10 (2003) 875–893 (preprint version read: webbox.lafayette.edu/~reiterc/3x+1/w3x+1_pp.pdf).
* S. Letherman, D. Schleicher, R. Wood, The 3n+1 problem and holomorphic dynamics, Experiment. Math. 8 (1999) 241–251 (via the abstract).
* D. Singer, Stable orbits and bifurcation of maps of the interval, SIAM J. Appl. Math. 35 (1978) 260–267 (standard; primary not re-read).
* J. C. Lagarias, The set of rational cycles for the 3x+1 problem, Acta Arith. 56 (1990) 33–53; E. Belaga, M. Mignotte, Embedding the 3x+1 conjecture in a 3x+d context, Experiment. Math. 7 (1998) (cited for 3x+d cycles; not re-read).
* Repo: THM-4554, THM-4555, THM-4556 (mac-mini), HYP-9213, HYP-9214, THM-4526 (Paley), the S18 note `mersenne_switch_parity_f21_compression_20261006.md`, the S15 note `collatz_one_out_edge_systems_20260930.md`, the Kawasaki audit note `collatz_procgen_20260924_fixed_points_kawasaki_audit.md` (Theorem K(b)).
