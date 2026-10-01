# Drop multiplicities of the Syracuse map: the divisor problem over 2^v − 3, k-step drops, other q

**Lane:** procgen cdiff, 2026-10-01 (resumed lane; owner prompt of 2026-10-01 on `F(2N−1) = 2N−1−K_N`).
**Status:** COMPLETE for this pass. Audit by the orchestrator: OWED.
**Builds on:** opus S15, twentieth note
[`collatz_label_map_paley_snarks_20261001.md`](collatz_label_map_paley_snarks_20261001.md) §1 (canon THM-4527,
Theorems 1–2), and THM-4484 (free and sporadic cycles; Gersonides for every odd q). This note re-verifies
THM-4527 with independent code and then goes beyond it. S15's directions D74 and D75 are the starting point.
**Code:** `04-computation/experiments/procgen_cdiff_20261001_run.py` (95 checks, about 3 minutes, peak RSS
200–300 MB, plus a child process of about 270 MB that runs first, ends `ALL CHECKS PASSED`). **Output:** `05-knowledge/results/procgen_cdiff_20261001.out`.
**Labels used:** PROVED, CITED, FINITE-EXACT, VERIFIED, EMPIRICAL, CONJECTURE (= OPEN, with evidence).
**Nothing here proves anything about Collatz itself.** Section 4.7 explains why drop statistics cannot.

Notation. For odd `A`, `S(A) = oddpart(3A+1)`, `v = v(A) = v_2(3A+1)`, the drop is `K(A) = (A − S(A))/2`, and
`M_v = 2^v − 3` (so `M_1 = −1`, `M_2 = 1`, `M_3 = 5`, `M_4 = 13`, `M_5 = 29`, …). The owner's `K_N` is
`K(4N−3)`. For `d ≥ 0`, `m(d)` is the number of odd `A ≥ 1` with `K(A) = d`. THM-4527 gives
`6K + 1 = (2^v − 3)·S` and `m(d) = #{v ≥ 2 : M_v | 6d+1}`. Write `N(n) = #{v ≥ 3 : M_v | n}`, so that
`m(d) = 1 + N(6d+1)`.

## Status table

| # | Result | Label |
|---|---|---|
| 1 | THM-4527 (S15 Theorems 1–2) re-verified with independent code: `F` rules to `10^6`, the identity for odd `|A| < 2·10^6`, the multiplicity formula for all `0 ≤ d ≤ 10^6`, S15's histogram | FINITE-EXACT |
| 2.1 | The density `δ_j` of `{d : m(d) = j}` exists and equals `P(N(ω) = j−1)` for a Haar-random profinite integer `ω`; truncation at `v ≤ V` costs at most `2^{1−V}` | PROVED |
| 2.2 | `E[u^N] < ∞` for `u < 4`; the inclusion–exclusion series `δ_{j+1} = Σ_T (−1)^{|T|−j} C(|T|,j)/lcm(M_T)` converges absolutely | PROVED |
| 2.3 | The lattice of divisibility events: `gcd(M_v, M_w) | 2^{|w−v|} − 1`; `E_w ⊂ E_v ⟺ ord_{M_v}(2) | w − v`; for `v ≤ 64` exactly ten primes are shared | PROVED + FINITE-EXACT |
| 2.4 | Exact densities `δ_1 … δ_10` to 18 decimals (exact rationals for `v ≤ 64`, error `< 1.1·10^{−19}`): `δ_1 = 0.696174903440132595`, `δ_2 = 0.266238036677738374`, … | PROVED + FINITE-EXACT; VERIFIED on `[0, 10^8]` |
| 2.5 | Maximal order: `√(2 log₂X) − 2.5 ≤ max_{d≤X} m(d) ≤ 1 + ½ log₂(6X+1)`; exact records up to 13 copies (`d ≈ 2.4·10^{27}`) | PROVED + FINITE-EXACT |
| 2.5′ | `max_{d≤X} m(d) = (1+o(1))√(2 log₂X)` | CONJECTURE (records fit `log₂ n_k ≈ k²/2 + 2k`) |
| 2.6 | Over all odd integers (both signs) every `d ∈ Z` is a drop exactly `2 + N(6d+1)` times; the two "unit" preimages are `8d+1` and `−4d−1`. In labels: once among the even `M`, once among `M ≡ 1 (mod 4)`, and the extra copies among `M ≡ 3 (mod 4)`. This is the precise form of the owner's "two copies" | PROVED + FINITE-EXACT |
| 3.1 | The summatory function of `K`: exact recursion, a base-4 digit formula, `Σ_{N≤4^n} K_N = (3·16^n − 4^n − 2)/6`, and `Σ_{N≤X} K_N = X²/2 + X·Φ(X) + O(log X)` with `Φ ∈ [−5/6, 1/6]` (both extremes attained in the limit) | PROVED + FINITE-EXACT |
| 3.2 | `Σ_{A odd < 2^n} K(A) = −J_{n−1}` (Jacobsthal numbers): ascents and descents cancel to `O(X)` | PROVED + FINITE-EXACT (`n ≤ 22`) |
| 4.1 | k-step drop identity `(2^p − 3^k)·Syr^k(A) = 2·3^k·D + c(v)` and a bijection between the points with valuation word `v` and a half-progression of drops with step `|2^p − 3^k|` | PROVED + FINITE-EXACT |
| 4.2 | Multiplicity of a k-step drop, its generating function (for `k = 1` a Lambert series) | PROVED + FINITE-EXACT (`k = 2, 3`, `|D| ≤ 3000`) |
| 4.3 | Mean multiplicities: `ρ_k^+ → 1` and `ρ_k^− → 0`. The error is `O(e^{−0.05498 k} k^8)`, with spikes at the resonances `2^p ≈ 3^k` | PROVED (Rhin's irrationality measure, CITED) |
| 4.4 | "Two copies" at `k = 2`: every `D ≤ −1` is the 2-step drop of `−16D−5` and `−16D−7` (Gersonides `3² − 2³ = 1`) | PROVED |
| 4.5 | Two consecutive drops determine the point: `A ↦ (K(A), K(S(A)))` is injective on odd `A ≥ 1`, for `3x+1` and for `3x−1` | PROVED (EMPIRICAL for other `q`) |
| 4.6 | Joint law along orbits: normalized drops are asymptotically i.i.d., with law `(2^V − 3)/2^{V+1}`, `V ~ Geom(1/2)` | PROVED (Terras/Everett, CITED) + FINITE-EXACT |
| 4.7 | `D = 0` is the cycle equation. The `3x+1` and `3x−1` sheets have the same drop statistics but different cycles, so drop statistics do not constrain cycles. The unit ("two copies") words are exactly the free cycles of THM-4484 | PROVED; answer to the cycle question: **NO** |
| 4.8 | Exact k-step multiplicity laws for `k = 2, 3`: on `D ≥ 0`, `P(μ_2 = 0, 1, 2, …) = 0.4066556, 0.4111170, 0.1527614, …` and `P(μ_3 = 0, 1, 2, …) = 0.1000254, 0.3699371, 0.2527803, …`; on `D < 0`, `μ_2 = 2 + [5 | D]` and `μ_3 = [19 | D] + [D mod 11 ∈ {1, 7, 8}]` | PROVED (existence, error `< 3·10^{−12}`) + VERIFIED (float64 elimination, four cross-checks) |
| 5.1 | Other q: `2qK + s = (2^v − q)·T` for `T = oddpart(qA+s)`; multiplicity `#{v : (2^v − q) | 2qd + s, quotient > 0}`; k-step `(2^p − q^k)T^k = 2q^kD + s·c_q(v)` | PROVED + FINITE-EXACT |
| 5.2 | The drop map covers all `d ≥ 0` iff `q = 2^a − 1`, all `d < 0` iff `q = 2^a + 1`, all of `Z` iff `q = 3`; otherwise a side misses a set of density `> 0.38` | PROVED + FINITE-EXACT (`q ≤ 1025`) |
| 5.3 | Exact multiplicity laws for `q = 1, 3, 5, 7, 9, 11, 13`; for example `5x+1` misses 59.27% of all `d ≥ 0` | FINITE-EXACT + VERIFIED |
| 5.4 | Drift dichotomy: `ρ_k^+ + ρ_k^− → 1` for `q ∈ {1, 3}` and `→ 0` for `q ≥ 5` | PROVED (Baker-type bounds, CITED) + FINITE-EXACT |
| 6 | Mod `3^k`: `S ≡ (6K+1)(2^v−3)^{−1}`, a function of `(K mod 3^{k−1}, v mod 2·3^{k−1})`; `Syr^k(A) ≡ 2^{−p} c(v) (mod 3^k)` depends only on the word; `d mod 3^j` is independent of `m(d)` | PROVED + FINITE-EXACT |

## 0. Headlines

1. **The "two copies" are exact over Z.** For the Syracuse map on all odd integers, every difference
   `d ∈ Z` has exactly two unit preimages, `8d+1` (branch `v = 2`, `M_2 = 1`) and `−4d−1` (branch `v = 1`,
   `M_1 = −1`). They have opposite signs. On top of these, `d` has one more preimage for each
   `v ≥ 3` with `(2^v − 3) | 6d+1` (Prop. 2.6). On the positive integers this gives: every `d < 0` exactly
   once, and every `d ≥ 0` exactly `1 + N(6d+1)` times. Among odd `q`, only `q = 3` has a unit branch on both
   sides (`3 = 2 + 1 = 4 − 1`), so it is the only `q` whose drop map hits every integer (Thm 5.2).
2. **The multiplicity law is a divisor problem with sparse dependence.** `m(d) − 1` is the number of moduli
   `2^v − 3` (`v ≥ 3`) that divide a Haar-random integer. The moduli are almost pairwise coprime: among
   `v ≤ 64` only ten primes are shared. This makes the densities exactly computable:

   | `j` | 1 | 2 | 3 | 4 | 5 | 6 |
   |---|---|---|---|---|---|---|
   | `δ_j` | 0.6961749034 | 0.2662380367 | 0.0353895588 | 0.0021345852 | 6.2061·10⁻⁵ | 8.4915·10⁻⁷ |

   Full precision is in §2.4. The dependence is positive (`5 | 125`, and `5` divides every fourth modulus). It
   lifts `P(m = 1)` from 0.69023 (independent model) to 0.69617.
3. **Two consecutive drops determine the point** (Thm 4.5). One drop determines `A` only up to `m(K(A))`
   choices, but the pair `(K(A), K(S(A)))` determines `A`.
4. **k-step drops.** `(2^p − 3^k)·Syr^k(A) = 2·3^k·D + c(v)`. The mean number of points with a given k-step
   drop tends to 1 for `D → +∞` and to 0 for `D → −∞`. The spikes sit exactly at the resonances
   `2^p ≈ 3^k` and are damped by `e^{−0.05498 k}` (Prop. 4.3). The full k-step laws are determined for
   `k = 2, 3` (Thm 4.8, to 12 digits): 40.7% (resp. 10.0%) of all `D ≥ 0` are never a 2-step (resp. 3-step) drop, although
   every `D ≥ 0` is a one-step drop.
5. **Drop statistics cannot see cycles** (Thm 4.7). The `3x−1` sheet has the same drop statistics as `3x+1`,
   at every `k`. It also has the same pair-injectivity. Yet it has the cycles `{1}`, `{5, 7}` and `{17, …}`. The
   only arithmetic special to `q = 3` is the units `|2^p − 3^k| = 1`, and at `D = 0` these produce exactly the
   free cycles `{1}`, `{−1}`, `{−5, −7}` of THM-4484.

## 1. Independent re-verification of THM-4527 (FINITE-EXACT)

Our code shares no code with S15's script. It confirms the following (checks A1–A6):
- `F(2N) = 3N`, `F(4j+1) = 3j+1`, `F(4n−1) = F(n)` (to `10^6`);
- `F(2N−1) = 2N−1−K_N`, the owner's `K = 0,2,1,4,2,10,3,9,4,15,5`, the examples `S(1), S(5), S(9), S(13) = 1, 1, 7, 5`,
  and the 2-regular recursion `K_{2j+1} = j`, `K_{4m} = 5m−1`, `K_{4m−2} = K_m + 6m − 4`;
- `6K + 1 = (2^v − 3)S` for every odd `A` with `|A| < 2·10^6`, of both signs;
- `m(d) = #{v ≥ 2 : M_v | 6d+1}` for **every** `0 ≤ d ≤ 10^6`, by brute force over all `A ≡ 1 (mod 4)` up to
  `8·10^6 + 1` (enough, because `K ≥ (A−1)/8`);
- S15's histogram `{1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6}` on `[0, 10^5]`;
- every negative `d ≥ −10^5` occurs exactly once;
- the mean `Σ_{v≥2} 1/(2^v − 3) = 1.34367343318176901854448…`.

This number is irrational. For integer `q ≥ 2` and rational `r ≠ 0, −q^n`, the sum `Σ_{n≥1} 1/(q^n + r)` is
irrational and not a Liouville number (P. Borwein, J. Number Theory 37 (1991) 253–259; CITED). Apply it with
`q = 2`, `r = −3`.

Verdict on THM-4527 §§1–2: **confirmed**. No discrepancy was found.

## 2. The divisor problem over 2^v − 3 (D74)

### 2.1 The profinite model

**Theorem 2.1 (PROVED).**
- (a) For a finite `T ⊂ {3, 4, …}`, the set `{d ≥ 0 : M_v | 6d+1 for all v ∈ T}` is one residue class modulo
  `lcm(M_T)`, so its density is `1/lcm(M_T)`.
- (b) Let `ω` be Haar-random in `Ẑ` and `N(ω) = #{v ≥ 3 : M_v | ω}`. For every `j ≥ 1`, the natural density
  `δ_j` of `{d ≥ 0 : m(d) = j}` exists and equals `P(N = j−1)`.
- (c) If `N_V` counts only `3 ≤ v ≤ V`, then `|δ_{j+1} − P(N_V = j)| ≤ Σ_{v>V} 1/M_v < 2^{1−V}`.

*Proof.*
- (a) Every `M_v` is prime to 6, since `2^v − 3 ≡ (−1)^v (mod 3)` and it is odd. So 6 is invertible modulo
  `L = lcm(M_T)`, and the conditions all say `6d ≡ −1 (mod L)`.
- (b), (c) The event `{m_V(d) = j+1}` depends only on `d mod L_V`. The map `d ↦ 6d+1` is a bijection of
  `Z/L_V`, and `ω mod L_V` is uniform, so its density is `P(N_V = j)`. Moreover
  `#{d ≤ X : m(d) ≠ m_V(d)} ≤ Σ_{v>V, M_v ≤ 6X+1} (X/M_v + 1) ≤ X·ε_V + log₂(6X+4)` with
  `ε_V = Σ_{v>V} 1/M_v`. Hence the upper and lower densities of `{m = j+1}` lie within `ε_V` of `P(N_V = j)`.
  The same estimate holds in the probability space. Let `V → ∞`. Finally `M_v ≥ 2^{v−1}` for `v ≥ 3`. ∎

So **the law of the multiplicities is the law of the number of moduli `2^v − 3` dividing a random integer.**
Coprimality with 6 means the progression `6d+1` plays no role. In particular `6d−1` (the `3x−1` sheet,
§4.7) gives the same law.

### 2.2 Absolute convergence

**Lemma 2.2 (PROVED).**
- `E[u^N] = Σ_T (u−1)^{|T|}/lcm(M_T) < ∞` for `1 ≤ u < 4`.
- Hence `δ_{j+1} = Σ_{T finite} (−1)^{|T|−j} C(|T|, j)/lcm(M_T)`, and this series converges absolutely.
- `P(N ≥ k) ≤ E[u^N]·u^{−k}` for every `u < 4`.

*Proof.* Let `s = u − 1 < 3`. Split off `T = ∅` and the singletons, whose sum `Σ s/M_v` is finite. For
`|T| ≥ 2`, let `w = max T` and `v` the second largest element. Then:
- `1/lcm(M_T) ≤ 1/lcm(M_v, M_w) = gcd(M_v, M_w)/(M_v M_w)`;
- the rest of `T` is any subset of `[3, v−1]`, with total weight `(1+s)^{v−3}`;
- by Prop. 2.3(a), `gcd ≤ min(M_v, 2^{w−v})`, so `gcd/(M_v M_w) ≤ min(2^{1−w}, 2^{2−2v})`.

The double sum `Σ_{v<w} (1+s)^v min(2^{−w}, 4^{−v})` converges when `1 + s < 4`. Both regimes `v ≤ w/2` and
`v > w/2` give geometric series with ratio `√(1+s)/2` or `(1+s)/4`.

Write `V(ω) = {v ≥ 3 : M_v | ω}`. The identity `1[N = j] = Σ_{T ⊂ V(ω)} (−1)^{|T|−j} C(|T|,j)` holds
pointwise. Since `Σ_{T ⊂ V(ω)} C(|T|,j) = C(N,j)2^{N−j} ≤ 3^N`, which is integrable, we may integrate term by
term. ∎

### 2.3 The lattice of divisibility events

Let `E_v = {ω : M_v | ω}`.

**Proposition 2.3.**
- (a) `gcd(M_v, M_w) | 2^{w−v} − 1` for `v < w`, because `M_w − 2^{w−v} M_v = 3(2^{w−v} − 1)` and the gcd is
  prime to 3. (PROVED)
- (b) For a prime `p ≥ 5`, `p | M_v ⟺ 2^v ≡ 3 (mod p)`. This happens for some `v` iff `3 ∈ ⟨2⟩ ⊂ F_p^×`, and
  then exactly for `v ≡ ind_2 3 (mod ord_p 2)`. (PROVED)
- (c) `E_w ⊊ E_v ⟺ M_v | M_w ⟺ w > v` and `ord_{M_v}(2) | w − v`. So each `E_v` contains infinitely many
  `E_w`. (PROVED; FINITE-EXACT for `3 ≤ v ≤ 40`, `v < w ≤ 400`) For example:

  | `M_v` | divides `M_w` exactly when |
  |---|---|
  | `M_3 = 5` | `w ≡ 3 (mod 4)` |
  | `M_4 = 13` | `w ≡ 4 (mod 12)` |
  | `M_5 = 29` | `w ≡ 5 (mod 28)` |
  | `M_6 = 61` | `w ≡ 6 (mod 60)` |
  | `M_7 = 125` | `w ≡ 7 (mod 100)` |
  | `M_8 = 253` | `w ≡ 8 (mod 110)` |

- (d) **For `3 ≤ v ≤ 64` exactly ten primes divide two of the `M_v`.** (FINITE-EXACT; found from the pairwise
  gcds alone, which are all `≤ 431`; cross-checked with full factorizations)

  | Prime | Divides `M_v` exactly when | Remarks |
  |---|---|---|
  | 5 | `v ≡ 3 (mod 4)` | `5² | M_v` iff `v ≡ 7 (mod 20)`; `5³ | M_v` iff `v ≡ 7 (mod 100)` |
  | 11 | `v ≡ 8 (mod 10)` | |
  | 13 | `v ≡ 4 (mod 12)` | |
  | 19 | `v ≡ 13 (mod 18)` | |
  | 23 | `v ≡ 8 (mod 11)` | |
  | 29 | `v ≡ 5 (mod 28)` | |
  | 37 | `v ≡ 26 (mod 36)` | |
  | 47 | `v ≡ 19 (mod 23)` | |
  | 71 | `v ≡ 16 (mod 35)` | `M_16 = 13·71²` |
  | 431 | `v ≡ 13 (mod 43)` | |

  These primes touch 38 of the 62 moduli. The private parts of all 62 are pairwise coprime and prime to the
  ten primes.

So the events are functions of the independent `p`-adic valuations of `ω`. The ten shared primes have 3072
joint (truncated) states, and given the state, the private events are independent Bernoulli variables with
parameter `1/(private part)`.
This is the "structure theorem": **an exact finite formula for `P(N_64 = j)`**. Two exact identities confirm
it:
- `E[N_64] = Σ_{v≤64} 1/M_v`;
- `E[C(N_64, 2)] = Σ_{v<w≤64} 1/lcm(M_v, M_w)`.

### 2.4 Exact densities (PROVED error bound, FINITE-EXACT rationals)

Absolute error `< 2^{−63} ≈ 1.08·10^{−19}` (Thm 2.1(c) with `V = 64`):

| `j` | `δ_j = dens{d : m(d) = j}` | independent model | count on `[0, 10^8]` |
|---|---|---|---|
| 1 | 0.696174903440132595 | 0.69023469 | 69617492 |
| 2 | 0.266238036677738374 | 0.27724918 | 26623796 |
| 3 | 0.0353895588393158824 | 0.031149657 | 3538961 |
| 4 | 0.00213458515130198837 | 0.001341173 | 213471 |
| 5 | 6.20614225702492·10⁻⁵ | 2.5082·10⁻⁵ | 6192 |
| 6 | 8.491492917053·10⁻⁷ | 2.169·10⁻⁷ | 89 |
| 7 | 5.302819767·10⁻⁹ | | |
| 8 | 1.6800745·10⁻¹¹ | | |
| 9 | 2.8666·10⁻¹⁴ | | |
| 10 | 2.70·10⁻¹⁷ | | |

- The mean is `E[m] = 1.34367343318176901854448…` and the variance is `Var(m) = 0.3099105129…`.
- The counts on `[0, 10^8]` agree with `δ_j` to within `1.5·10^{−7}` (check B3).
- The shared primes matter: compared with the independent model they raise `P(m=1)` by 0.0059 and `P(m=4)` by
  a factor of 1.59. This is positive dependence: once one event holds, nested ones are likelier.
- `E[2^{N}] ≈ 1.3883` and `E[3^{N}] ≈ 1.8763`, both finite as Lemma 2.2 says.

### 2.5 Maximal order and tails

**Proposition 2.5.**
- (a) (PROVED) `m(d) ≤ 1 + ½ log₂(6d+1)` for every `d ≥ 0`.
  *Proof.* Let `n = 6d+1` and `T = {v ≥ 3 : M_v | n}`, with largest element `w`. For `v ∈ T \ {w}`,
  `n ≥ lcm(M_v, M_w) ≥ M_v M_w/2^{w−v} ≥ 2^{v−1} 2^{w−1}/2^{w−v} = 2^{2v−2}`. So every element of `T` except `w`
  is `< 1 + ½ log₂ n`, and `|T| ≤ ½ log₂ n`. ∎
- (b) (PROVED) `max_{0≤d≤X} m(d) ≥ √(2 log₂X + 6.25) − 2.5`.
  *Proof.* Take `T = {3, …, k+2}`, with `lcm(M_T) < 2^{(k+2)(k+3)/2 − 3}`. The least `n ≡ 1 (mod 6)` divisible
  by it is at most `5·lcm(M_T)`, so some `d < lcm(M_T)` has `m(d) ≥ k+1`. ∎
- (c) (FINITE-EXACT; exhaustive search over all 2262581 sets `T` with `lcm(M_T) ≤ 2^{100}`) Let `n_k` be the least
  `n ≡ 1 (mod 6)` with `N(n) ≥ k`, and `d_k = (n_k − 1)/6` the first `d` with `m(d) ≥ k+1`:

  | `k` | `d_k` | `n_k` (its moduli) | `log₂ n_k` |
  |---|---|---|---|
  | 1 | 2 | 13 | 3.70 |
  | 2 | 24 | 145 = 5·29 | 7.18 |
  | 3 | 314 | 1885 = 5·13·29 | 10.88 |
  | 4 | 7854 | 47125 = 125·13·29 | 15.52 |
  | 5 | 479104 | 2874625 = 125·13·29·61 | 21.45 |
  | 6 | 121213354 | 727280125 | 29.44 |
  | 7 | 49576261854 | 297457571125 | 38.11 |
  | 8 | 50617363353104 | 303704180118625 | 48.11 |
  | 9 | 115043252496711229 | 690259514980267375 | 59.26 |
  | 10 | 117459160799142164979 | 704754964794852989875 | 69.26 |
  | 11 | 480760345150888881259729 | 2884562070905333287558375 | 81.25 |
  | 12 | 2423512899905630850430294729 | 14541077399433785102581768375 | 93.55 |

  The records use the nesting `M_3 | M_7` (`lcm(5, 125) = 125`). From `k = 10` on they also take `v = 19`
  and `v = 16`, whose moduli `M_19 = 5·23·47·97` and `M_16 = 13·71²` reuse the primes 5, 23 and 13. So `v = 19`
  costs 12.2 bits instead of 19. For example, `n_12` is divisible exactly by `M_v` for `v = 3, …, 12, 16, 19`.
  Otherwise the records stay close to the consecutive products `Π_{v=3}^{k+2} M_v ≈ 2^{k²/2 + 2.5k}`.
- (d) (CONJECTURE, OPEN) `max_{d≤X} m(d) = (1 + o(1))·√(2 log₂X)`. Equivalently, any `k` distinct moduli have
  `log₂ lcm ≥ (½ − o(1))k²`.
  - Evidence: `log₂ n_k − k²/2 = 3.2, 5.2, 6.4, 7.5, 8.95, 11.4, 13.6, 16.1, 18.8, 19.3, 20.8, 21.6` for
    `k = 1..12`, i.e. about `1.8k` to `2k`.
  - A proof would need lower bounds for the "primitive part" of `2^n − 3`, a Zsigmondy-type statement that
    is not available for `a^n − b`.
- (e) Tails. `P(N ≥ k) ≥ 1/lcm(M_3, …, M_{k+2}) ≥ 2^{−(k²+5k)/2}` (PROVED), and `P(N ≥ k) = O((4−ε)^{−k})`
  (PROVED, Lemma 2.2). Numerically `log₂ δ_{k+1} = −1.91, −4.82, −8.87, −13.98, −20.17, −27.49, −35.79, −44.99`
  for `k = 1..8`, which is `≈ −(k²/2 + 1.6k)` (EMPIRICAL): the tail is Gaussian in `k`.

### 2.6 Exactly two copies over Z

**Proposition 2.6 (PROVED; FINITE-EXACT for `|d| ≤ 2·10^4`).** Extend `S` and `K` to all odd integers
`A ≠ 0`. Every `d ∈ Z` is `K(A)` for exactly `2 + N(|6d+1|)` odd integers `A`:
- `A = 8d+1`, branch `v = 2`, with `S = 6d+1`;
- `A = −4d−1`, branch `v = 1`, with `S = −6d−1`;
- one `A = (2^v s − 1)/3` with `s = (6d+1)/M_v` for each `v ≥ 3` with `M_v | 6d+1`. This `A` has the sign of
  `6d+1`.

*Proof.* This is THM-4527's bijection `(d, v) ↔ (A, v)` without the positivity condition. `s = (6d+1)/M_v` is
an odd integer exactly when `M_v | 6d+1`, and then `A` is an odd integer automatically. ∎

**In the owner's labels** `M = (A+1)/2`, now allowed to be any integer, this reads as follows (PROVED; checked for
`|d| ≤ 5000`). Every integer `d = M − F(M)` occurs:
- exactly once among the even labels (`M = −2d`, since `F(2N) = 3N`);
- exactly once among the labels `M ≡ 1 (mod 4)` (`M = 4d + 1`, since `F(4j+1) = 3j+1`);
- `N(|6d+1|)` more times among the labels `M ≡ 3 (mod 4)`, the microcosm `F(4n−1) = F(n)`.

So **the owner's "exactly 2 copies of each difference" is exactly true for the two unit branches, counted over
both signs.** One copy is a positive odd number and one is negative; which is which depends on the sign of
`d`. The `3x+1` map on the positive integers sees one of the two copies. The `3x−1` map, which is `3x+1` on
the negatives, sees the other. The extra copies are the divisor events of §2.1–2.5.

## 3. The summatory function of K (D74)

**Theorem 3.1 (PROVED; FINITE-EXACT for `X ≤ 2^20`).** Let `𝒮(X) = Σ_{N=1}^{X} K_N`.
- (a) `𝒮(X) = J(J+1)/2 + 5U(U+1)/2 − U + 3W(W+1) − 4W + 𝒮(W)`, where `J = ⌊(X−1)/2⌋`, `U = ⌊X/4⌋`,
  `W = ⌊(X+2)/4⌋` and `𝒮(0) = 0`. This comes from the three rules of THM-4527: the odd positions, the
  positions `≡ 0 (mod 4)`, and the positions `≡ 2 (mod 4)` that replay the sequence.
- (b) **Digit formula.** Let `X_0 = X` and `X_{i+1} = ⌊(X_i + 2)/4⌋`. Then
  `𝒮(X) = X²/2 + Σ_{i≥0} ℓ(X_i)`, where

  | `x` | `4t` | `4t+1` | `4t+2` | `4t+3` |
  |---|---|---|---|---|
  | `ℓ(x)` | `−t/2` | `−(5t+1)/2` | `(t+1)/2` | `−(3t+2)/2` |

  Equivalently, write `X = Σ_i e_i 4^i` with digits `e_i ∈ {−2, −1, 0, 1}`. Then
  `𝒮(X) = X²/2 + X·Φ(X) + O(log X)` with `Φ(X) = Σ_i c(e_i) 4^{−i}`, where
  `c(0) = −1/8`, `c(1) = −5/8`, `c(−2) = 1/8`, `c(−1) = −3/8`.
- (c) `𝒮(4^n) = (3·16^n − 4^n − 2)/6`.
- (d) `−5/6 ≤ Φ ≤ 1/6`. Both bounds are attained in the limit: `Φ → 1/6` along `X = (4^n+2)/3` (all digits
  −2) and `Φ → −5/6` along `X = (4^n−1)/3` (all digits 1). The first family indexes the trunk
  `A = 4X − 3 = (4^{n+1}−1)/3`. The Cesàro mean of `Φ` is `−1/3` (VERIFIED on `[10^3, 2^20]`; the digits of a random `X` are equidistributed).

*Proof of (b).*
- Write (a) as `𝒮(X) = Q(X) + 𝒮(X_1)`. Then `ℓ(X) = Q(X) − (X² − X_1²)/2`, computed in each residue class
  mod 4.
- Then `X_i = X/4^i + O(1)` and `ℓ(X_i) = c(e_i)X_i + O(1)`.
- The range in (d) follows from `max c = 1/8`, `min c = −5/8` and `Σ 4^{−i} = 4/3`. ∎

The leading term `X²/2` says that on average `K_N ≈ N`. In the owner's indexing `K_N` is linear in `N` on each
valuation class, with slope `(2^v − 3)/2^{v−1}`. So the slopes are `1/2, 5/4, 13/8, 29/16, … → 2`, and the
class with `v` halvings has density `2^{1−v}`.

The counting function of the multiplicities is exact as well:
`Σ_{d≤Y} m(d) = Σ_{v≥2} (⌊(Y − d_v)/M_v⌋ + 1) = 1.34367…·Y + O(log Y)`, where `d_v` is the least solution
(check C5). Its generating function is the Lambert-type series (check C6)

    Σ_{d≥0} m(d) q^{6d+1} = Σ_{v≥2 even} q^{M_v}/(1 − q^{6M_v}) + Σ_{v≥3 odd} q^{5M_v}/(1 − q^{6M_v}).

**Proposition 3.2 (PROVED; FINITE-EXACT `n ≤ 22`).** Over all odd `A < 2^n`, ascents and descents together,
`Σ K(A) = −J_{n−1} = −(2^n + 2(−1)^n)/6`, where `J` is the Jacobsthal sequence `0, 1, 1, 3, 5, 11, 21, 43, …`.

*Proof.* The descents contribute `𝒮(2^{n−2})` and the ascents `A = 4N−1` contribute `−Σ_{N≤2^{n−2}} N`. Use
3.1(c) for `n` even. For `n` odd, use `𝒮(2·4^m) = 2·16^m − (4^m−1)/3`. This follows from `E(4t) = E(t) − t/2`,
where `E(X) = 𝒮(X) − X²/2` (the case `X ≡ 0 (mod 4)` of (b)), and `E(2) = 0`. ∎

So `Σ_{A odd < 2^n} S(A) = Σ_{A odd < 2^n} A + 2J_{n−1}`, and in general `Σ_{A odd ≤ X} K(A) = O(X)`
(observed range `[−X/3, X/6]`, with the extremes next to trunk elements). This is the drop-language form of
`E[S(A)/A] = E[3·2^{−v}] = 1`. The Syracuse map is a **fair game in the arithmetic mean**: up-steps and
down-steps cancel exactly to first order. It contracts only in the geometric mean (by `3/4`).

## 4. k-step drops and cycles (goal 3)

### 4.1 The k-step identity

Let `v = (v_1, …, v_k)` be the valuation word of `A`, `p = Σ v_i`, `p_i = v_1 + … + v_i`, and
`c(v) = Σ_{i=0}^{k−1} 3^{k−1−i} 2^{p_i}`.

**Theorem 4.1 (PROVED; FINITE-EXACT).**
- (a) With `D = (A − Syr^k A)/2`: `(2^p − 3^k)·Syr^k(A) = 2·3^k·D + c(v)`. For `k = 1` this is THM-4527's
  `6K + 1 = (2^v − 3)S`.
- (b) For each word `v`, the map `A ↦ s = Syr^k(A)` is a bijection from the odd `A ≥ 1` with word `v` onto the
  odd `s ≥ 1` with `2^p s ≡ c(v) (mod 3^k)`, a single class mod `2·3^k`. The inverse is
  `A = (2^p s − c(v))/3^k`.
- (c) Hence the drops of the points with word `v` form a half-progression `D_0(v) + (2^p − 3^k)t`, `t ≥ 0`,
  with step `|2^p − 3^k|`. It runs to `+∞` if `2^p > 3^k` and to `−∞` if `2^p < 3^k`.

*Proof.*
- (a) Induction gives `2^p Syr^k(A) = 3^k A + c(v)`; substitute `A = Syr^k A + 2D`.
- (b) Induction on `k`. Write `v = (v_1, v')` and `c(v) = 3^{k−1} + 2^{v_1}c(v')`. The congruence modulo
  `3^k` forces `2^{p−v_1}s ≡ c(v') (mod 3^{k−1})`, so `A' = (2^{p−v_1}s − c(v'))/3^{k−1}` has word `v'` and
  `Syr^{k−1}A' = s`. Then `3A + 1 = 2^{v_1}A'` with `A'` odd, so `A` has first valuation exactly `v_1`.
- (c) Consecutive admissible `s` differ by `2·3^k`, so consecutive drops differ by `2^p − 3^k`. ∎

Checked for `k ≤ 7`, `|A| ≤ 10^5`, and for all words with `k ≤ 4`, `p ≤ 11` (D1–D2).

### 4.2 Multiplicity and generating function

**Proposition 4.2 (PROVED; FINITE-EXACT for `k = 2, 3`, `|D| ≤ 3000`).**
- (a) The number of odd `A ≥ 1` whose k-step drop is `D` is
  `μ_k(D) = #{v ∈ Z_{≥1}^k : (2^p − 3^k) | 2·3^k D + c(v), (2·3^k D + c(v))/(2^p − 3^k) > 0}`.
- (b) Its generating function is
  `Σ_D μ_k(D) x^D = Σ_{v: 2^p>3^k} x^{D_0(v)}/(1 − x^{2^p−3^k}) + Σ_{v: 2^p<3^k} x^{D_0(v)}/(1 − x^{−(3^k−2^p)})`.
  For `k = 1` it is the Lambert series of §3.

For `k ≥ 2` the residues `−c(v)/(2·3^k)` depend on `v`. So, unlike `k = 1`, the k-step law is not a pure
divisor problem: it is a covering-type problem of half-progressions with moduli `2^p − 3^k`.

The exact laws for `k = 2, 3` are in §4.8. **40.7% of the possible 2-step descents and 10.0% of the possible
3-step descents never occur**, whereas every one-step descent occurs.

### 4.3 Mean multiplicities and the resonance

Let `ρ_k^+ = Σ_{p > k log₂3} C(p−1, k−1)/(2^p − 3^k)` and
`ρ_k^− = Σ_{k ≤ p < k log₂3} C(p−1, k−1)/(3^k − 2^p)`. By Thm 4.1(c) these are the densities of the k-step
drops: `Σ_{0≤D≤Y} μ_k(D) = ρ_k^+ Y + o(Y)`, and likewise `ρ_k^−` for `D → −∞`. The words with large `p`
contribute `o(Y)`, because `#{A ≤ X : p_k(A) > P}` is at most `X·P(P_k > P) + O_P(1)`. For fixed `k`, only
finitely many `A` have a descending word (`2^p > 3^k`) but a negative drop: such an `A` has
`A < c(v)/(2^p − 3^k) ≤ (3^k − 2^k)/(2^k(1 − 3^k 2^{−p}))`, which is bounded in terms of `k`. Conversely,
every `A` with an ascending word has `Syr^k A > 3^k A/2^p > A`, so its drop is negative.

**Proposition 4.3 (PROVED, with Rhin's irrationality measure CITED).** Let `λ = log₂ 3` and
`I = log 3 + (λ−1) log(λ−1) − λ log λ = 0.054979…`. This `I` is the large-deviation rate of the negative
binomial `P_k` (the sum of `k` i.i.d. `Geom(1/2)` valuations) at `P_k = kλ`.

Then `|ρ_k^+ − 1| + ρ_k^− = O(e^{−Ik} k^{8})`. The dominant terms are the two resonant ones,
`π_k(p*)/(2^{{kλ}} − 1)` and `π_k(p*+1)/(1 − 2^{{kλ}−1})`, where `π_k(p) = C(p−1,k−1)2^{−p}` and
`p* = ⌊kλ⌋`.

*Proof.*
- Write `1/(2^p − 3^k) = 2^{−p}/(1 − 3^k 2^{−p})`.
- Words with `p ≥ p* + 2`: their total is `P(P_k ≥ p*+2) + O(P''(P''_k > kλ))`, where `P''` is the negative
  binomial with success probability `3/4` (mean `4k/3 < kλ`). Both `P(P_k ≤ kλ)` and `P''(P''_k > kλ)` are
  exponentially small, by Chernoff.
- Words with `p ≤ p* − 1` contribute at most `P(P_k ≤ kλ) ≤ e^{−Ik}`.
- The two resonant terms carry `π_k = O(e^{−Ik})` (Chernoff). Rhin's bound `μ(log 3/log 2) ≤ 8.616` gives
  `{kλ}, 1 − {kλ} ≫ k^{−7.62}` (G. Rhin, *Approximants de Padé et mesures effectives d'irrationalité*,
  Progr. Math. 71 (1987) 155–164; CITED). ∎

Numbers (computed to 40 digits; checks D4):

| `k` | 1 | 2 | 3 | 5 | 7 | 12 | 17 | 29 | 41 | 53 | 100 | 250 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ρ_k^+` | 1.3437 | 0.8077 | 1.8594 | **3.5125** | 0.9229 | 0.9443 | **2.0024** | **1.6187** | **1.5898** | 0.9974 | 1.00027 | 1 + 3·10⁻⁹ |
| `ρ_k^−` | 1 | **2.2** | 0.3254 | 0.1631 | **1.6038** | **4.5097** | 0.0436 | 0.0167 | 0.0070 | **1.4744** | 4.4·10⁻⁴ | 1.6·10⁻⁷ |
| `{kλ}` | .585 | .170 | .755 | .925 | .095 | .020 | .944 | .964 | .984 | .003 | .496 | .241 |

The spikes (bold) sit at the convergents and intermediate fractions of `log₂ 3`
(`k = 2, 5, 7, 12, 17, 29, 41, 53, …`). They are the same resonance that S15's size bound and every cycle
argument meet. The negative-sheet cycles sit on such spikes: `{−5, −7}` has shape `(k, p) = (2, 3)` and
`{−17, …}` has shape `(7, 11)`. **For large `k`, the k-step map covers each large positive drop about once
and almost no large negative drop.**

### 4.4 Two copies at k = 2

**Proposition 4.4 (PROVED; FINITE-EXACT for `|D| ≤ 10^5`).** Every `D ≤ −1` is the 2-step drop of
`A = −16D − 5` (word `(1,2)`) and of `A = −16D − 7` (word `(2,1)`), because `3² − 2³ = 1`. It has a third
preimage iff `5 | D` (word `(1,1)`, modulus `3² − 2² = 5`), and no other. So `μ_2(D) = 2 + [5 | D]` for
`D ≤ −1`, and `ρ_2^− = 11/5`.

*Proof.* The only ascending words for `k = 2` are those with `p ≤ 3`. A descending word with `18D + c > 0` and
`D ≤ −1` needs `v_1 ≥ 4`, and then `s < 1`. ∎

Together with the two `k = 1` units this exhausts the unit moduli `|2^p − 3^k| = 1`: by THM-4484(2) they are
`(k, p) ∈ {(1,1), (1,2), (2,3)}`. So the owner's "two copies" phenomenon has exactly three instances. It
happens only for `q = 3`, and only for `k ≤ 2`.

### 4.5 Consecutive drops determine the point

**Theorem 4.5 (PROVED).** On the odd integers `A ≥ 1`, the map `A ↦ (K(A), K(S(A)))` is injective. The same
holds for the `3x−1` map `A ↦ oddpart(3A−1)`.

*Proof (`3x+1`).* Suppose `A ≠ A'` with the same `D_1 = K(A) = K(A')` and `D_2 = K(S(A)) = K(S(A'))`. Let
their valuations be `v, v'` and those of `S = S(A)`, `S' = S(A')` be `u, u'`.

1. `K < 0` exactly on the branch `v = 1`, and `(D_1, v)` determines `A`. Hence `v ≠ v'`, both are `≥ 2`, and
   `n_1 := 6D_1 + 1 > 0`. Then `S = n_1/M_v ≠ S' = n_1/M_{v'}`. Likewise `u ≠ u'`, both are `≥ 2`, and
   `n_2 = 6D_2 + 1 > 0`. Without loss of generality `u < u'`; put `δ = u' − u`.
2. From `S = 2D_2 + n_2/M_u` we get `3n_1/M_v + 1 = 2^u n_2/M_u`, and the same with primes.
3. Eliminating `n_2` gives `n_1 E = M_v M_{v'}(2^δ − 1)`, where
   `E = 2^δ M_u M_{v'} − M_{u'} M_v = 2^δ M_u (2^{v'} − 2^v) − 3(2^δ − 1)M_v`.
   Hence `E > 0`.
4. Since `lcm(M_v, M_{v'}) | n_1`, we get `E | g(2^δ − 1)` with `g = gcd(M_v, M_{v'})`, so
   `0 < E ≤ g(2^δ − 1)`.
5. `E > 0` forces `v' > v`. Put `e = v' − v`.
6. `g ≤ M_v` gives `M_u(2^e − 1) < 4`, so `u = 2` and `e ∈ {1, 2}`. Then `g | 2^e − 1 ∈ {1, 3}` with `g`
   prime to 3, so `g = 1` and `E ≤ 2^δ − 1`.
7. For `e = 2`, `E = 3·2^v + 9(2^δ − 1)`, which is too big.
8. For `e = 1`, `E = 9(2^δ − 1) − 2^v(2^{δ+1} − 3)`, and `0 < E ≤ 2^δ − 1` puts `2^v` in the interval
   `[8(2^δ−1)/(2^{δ+1}−3), 9(2^δ−1)/(2^{δ+1}−3))`. This interval contains a power of 2 only for `δ = 1`,
   `v = 3`.
9. That case gives `n_1 = 65 ≢ 1 (mod 6)`, a contradiction.

*Proof (`3x−1`).* The same computation gives `n_1 E = −M_v M_{v'}(2^δ − 1)`, so `E < 0`.
- If `v' < v`, then `|E| > 3(2^δ − 1)M_v > g(2^δ − 1)`.
- If `v' > v`, then `e = 1` and `u = 2` again, and `2^v` must lie in
  `(9(2^δ−1)/(2^{δ+1}−3), 10(2^δ−1)/(2^{δ+1}−3)]`. For `δ = 1` this is `(9, 10]`, for `δ = 2` it is
  `(5.4, 6]`, and for `δ ≥ 3` it lies inside `(4.5, 5.4)`. None of these contains a power of 2. ∎

Checks (D6):
- brute force for `A < 10^6` on both sheets;
- an exhaustive search of all quadruples `(v, v', u, u')` with exponents `< 60`. This covers points of
  unbounded size. Its only hit is `(1, 2, 1, 2)`: the fixed points `−1` and `1`, which share the drop pair
  `(0, 0)` across the two sheets;
- the proof's finite case analysis, re-run numerically for all exponents up to 33.

The analogous statement holds empirically for `qx+1` with `q = 1, 5, 7, 9, 11, 13` and for `5x−1`. There is
no repeated pair for `A < 2·10^5` (in the runner) or `A < 4·10^5` (scratch run). EMPIRICAL.

**Meaning.** "Knowing which A values map to each other" is quantitatively this:
- a single drop pins `A` down to `m(K(A))` candidates, with mean 1.3437 and law `δ`;
- two consecutive drops pin it down completely;
- over `k` steps a total drop leaves about `ρ_k^+ ≈ 1` candidates.

### 4.6 The joint law along orbits

**Proposition 4.6.** The valuation words are equidistributed: the odd `A` with a given word of sum `p` form
one class modulo `2^{p+1}`. This is classical (Terras 1976, Everett 1977; CITED), and is checked exactly for
`k = 2` on all odd `A < 2^21` (D9). Hence, under natural density, the normalized drops
`(K(A_i)/A_i)_{i<k} = ((2^{v_i} − 3)/2^{v_i+1} − 1/(2^{v_i+1}A_i))_{i<k}` converge to i.i.d. copies of
`(2^V − 3)/2^{V+1}` with `P(V = j) = 2^{−j}`. The values are `−1/4, 1/8, 5/16, 13/32, … → 1/2`, with mean 0
(Prop. 3.2) and geometric drift `E log(3·2^{−V}) = log(3/4) < 0`. (PROVED)

### 4.7 D = 0 is the cycle equation. Do drop statistics constrain cycles?

**Theorem 4.7.**
- (a) (PROVED) `μ_k(0)` is the number of positive Syracuse periodic points whose period divides `k`. The
  unsigned count `#{v : (2^p − 3^k) | c(v)}` counts all integer periodic points: positive `s` on the `3x+1`
  sheet, negative `s` on the `3x−1` sheet. FINITE-EXACT for `k ≤ 10`, over all words with `k ≤ p ≤ 2k`: the
  fibre is `{1}` and `{−1}`, plus `{−5, −7}` if `2 | k`, plus the seven points of the `−17` cycle if `7 | k`.
- (b) (PROVED) **The `3x+1` and `3x−1` sheets have the same drop-multiplicity densities, for every `k`.**
  - For `k = 1` the laws are `1 + N(6d+1)` and `1 + N(6d−1)`, which are identical by Thm 2.1.
  - For general `k`, the residues of the two sheets are `∓c(v)/(2·3^k)` modulo `2^p − 3^k`. Negation preserves
    Haar measure, so every finite truncation (words with `p ≤ P`) has the same law. The sign conditions remove
    finitely many points per word. The words with `p > P` affect a set of upper density at most
    `2·P(P_k > P)` on either sheet (the argument of Thm 4.8), so the limits agree.
  - Measured `k = 2` laws on `[0, 10^6]`: `(0.40666, 0.41111, 0.15275, 0.02694, 0.00242)` for `3x+1` and
    `(0.40668, 0.41108, 0.15274, 0.02700, 0.00237)` for `3x−1`, against the exact
    `(0.40666, 0.41112, 0.15276, 0.02693, 0.00242)` of Thm 4.8.
  - The sheets also share Theorem 4.5.
- (c) (PROVED, via THM-4484) The unit words, which are what produce the "two copies", are the free shapes of
  THM-4484. At `D = 0` they give exactly the free cycles `{1}`, `{−1}`, `{−5, −7}`. Every other cycle (the
  `−17` cycle, or any hypothetical positive one) needs `(2^p − 3^k) | c(v)` for a non-unit modulus.

**Honest answer: no, drop statistics do not constrain cycles.**
- Every density-level statement about drops (the law `δ`, the k-step laws, the means `ρ_k^±`, the injectivity
  of drop pairs) is identical on the `3x−1` sheet, which has three cycles.
- The drop value `D = 0` is a single point of a sum of periodic functions. Its value there *is* the
  Böhm–Sontacchi cycle equation, not a consequence of the statistics.
- What the drop language does see is the classical size constraint. A cycle balances its drops,
  `Σ_i K(A_i) = 0`, which is `Π(3 + 1/A_i) = 2^p`. This pushes `2^p` into `(3^k, (3 + 1/min A)^k]`, the
  resonance of Prop. 4.3.

Positive sheets of `5x+1` serve as a control. There `ρ_k → 0`, yet the `D = 0` fibres are its three cycles
`{1, 3}`, `{13, 33, 83}`, `{17, 43, 27}` (E7).

### 4.8 Exact k-step multiplicity laws for k = 2 and 3 (D77 for small k)

**Theorem 4.8.** Fix `k ∈ {2, 3}`.
- (a) (PROVED) The natural density `δ_j^{(k)}` of `{D ≥ 0 : μ_k(D) = j}` exists. It differs from the law of
  the truncated system (words with `p ≤ P`) by at most `2·P(P_k > P)`.
- (b) (PROVED) The truncated law is the law of the number of words hit by a Haar-random `x = 2·3^k·D`. The
  word `v` is hit iff `x ≡ −c(v) (mod 2^p − 3^k)`.
- (c) (VERIFIED, float64) For `P = 50` (error `≤ 9·10^{−14}` for `k = 2` and `≤ 2.3·10^{−12}` for `k = 3`):

  | `j` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  |---|---|---|---|---|---|---|---|---|
  | `δ_j^{(2)}` | 0.406655577104 | 0.411117039014 | 0.152761396106 | 0.026927858359 | 0.002419931618 | 1.15124632·10⁻⁴ | 3.027697·10⁻⁶ | 4.5081·10⁻⁸ |
  | `δ_j^{(3)}` | 0.100025419499 | 0.369937051720 | 0.252780315253 | 0.159597308002 | 0.087756944313 | 0.025654295554 | 0.003909015456 | 3.2424097·10⁻⁴ |

- (d) (PROVED) On `D < 0` the laws are finite and explicit:
  - `μ_2(D) = 2 + [5 | D]`, so the law is `4/5` at 2 and `1/5` at 3;
  - `μ_3(D) = [19 | D] + [D mod 11 ∈ {1, 7, 8}]`, so the law is `144/209, 62/209, 3/209` at 0, 1, 2.

*Proof.*
- (a) Words with `p ≤ P` hit `D ≥ 0` in finitely many residue classes, with no sign restriction, since
  `c(v) > 0`. So that part is periodic, and its density is the Haar probability in (b).
- For a word with `p ≥ p* + 2`, `3^k/2^p < 1/2` and `c(v)/2^p < (3/2)^k`. Hence `D > A/4 − (3/2)^k/2`, and the
  pairs `(A, v)` with `p > P` and `D ≤ Y` number at most `#{A < 4Y + C_k : p_k(A) > P}`. That is
  `2Y·P(P_k > P) + O_P(1)` by Prop. 4.6.
- (d) The ascending words are those with `2^p < 3^k`. Their moduli are `9 − 4 = 5` and `9 − 8 = 1` (twice) for
  `k = 2`, and `27 − 8 = 19` and `27 − 16 = 11` (three words, residues `c ≡ 8, 1, 7`) for `k = 3`. ∎

*Computation.* Condition on the `ℓ`-adic digits of `x` at the primes shared by two moduli. There are 10 such
primes for `k = 2` and 11 for `k = 3`; for `k = 2` they are `5, 7, 11, 13, 17, 19, 23, 29, 41, 47`. Given the
digits, the words are independent. The digits are summed out by variable elimination: the largest table is the
5-clique `{5, 7, 23, 29, 47}`, with 1.1·10⁶ cells. Four cross-checks:
1. the same engine reproduces the exact `k = 1` law of §2.4 to `10^{−13}`;
2. its mean equals `ρ_k^+` (`0.807697240794`, `1.85943632458`) to `10^{−12}`;
3. its second factorial moment equals the exact sum over word pairs of `[CRT compatible]/lcm`
   (`0.498524360514`, `3.161071794409`) to `10^{−11}`;
4. sieves over `[0, 10^6]` agree to `5·10^{−5}` (`k = 2`) and `5·10^{−4}` (`k = 3`).

**Reading.** At `k = 1` every `D ≥ 0` is a drop, because of the unit `2² − 3 = 1`. At `k = 2` there is no unit
on the descending side, and 40.7% of all `D ≥ 0` are never a 2-step drop. At `k = 3` the resonance
`2^5 = 32 ≈ 27` (`ρ_3^+ = 1.86`) brings the holes down to 10.0%. The ascending sides are governed by the
small moduli `3^k − 2^p`. They contain the unit `3² − 2³ = 1` at `k = 2` (the "two copies" of Prop. 4.4).

## 5. Other q (D75; goal 4)

### 5.1 The general drop law

Let `T(A) = oddpart(qA + s)` with `q` odd and `s = ±1`, `K = (A − T)/2`, `v = v_2(qA + s)`.

**Theorem 5.1 (PROVED; FINITE-EXACT, checks E1, E2, E6).**
- `2qK + s = (2^v − q)·T`.
- The multiplicity of `d` over odd `A ≥ 1` is `#{v ≥ 1 : (2^v − q) | 2qd + s, (2qd+s)/(2^v−q) > 0}`.
- Over `k` steps, `(2^p − q^k)T^k(A) = 2q^k D + s·c_q(v)`, with `c_q(v) = Σ q^{k−1−i}2^{p_i}`.

The proof is that of THM-4527 and Thm 4.1, with `3 → q` and `+1 → s`. For `q = 3`, `s = −1` this is the
`3x−1` sheet: `(2^v − 3)T = 6K − 1`, so every `d ≤ 0` is a drop exactly once (`v = 1`) and every `d ≥ 1`
exactly `1 + N(6d−1)` times.

### 5.2 Units and surjectivity

**Theorem 5.2 (PROVED; FINITE-EXACT for odd `q ≤ 1025`).** Take `s = +1` and `A ≥ 1`.
- The descending moduli are `2^v − q` with `2^v > q`. They contain a unit iff `q = 2^a − 1`, and then every
  `d ≥ 0` is a drop.
- The ascending moduli are `q − 2^v` with `v ≥ 1`. They contain a unit iff `q = 2^a + 1`, and then every
  `d < 0` is a drop.
- A side without a unit has `Σ 1/|modulus| < 0.62`. So its holes (drop values never taken) have density
  `> 0.38`.

**Hence the drop map hits every integer iff `q = 3`**, the only odd number that is both `2^a − 1` and
`2^b + 1`. For `q = 1` there are no ascents.

*Proof.*
- Descending side: `M_{v_0+1} = 2M_{v_0} + q ≥ 7` whenever `M_{v_0} ≥ 3`, so the sum is at most
  `Σ_{i≥2} 1/(2^i − 1) = 0.6067`.
- Ascending side: the moduli `q − 2^v` are distinct odd numbers. The smallest is `≥ 3`, and the others are
  `≥ 3 + 2^{v_0−2}`, which gives a sum `≤ 1/3 + 2/7`.
- A unit modulus divides everything. ∎

### 5.3 Exact multiplicity laws (FINITE-EXACT, truncation `v ≤ 48`; `q = 1` at `v ≤ 24`)

| `q` | units (desc / asc) | descending `P(m = 0, 1, 2, 3)` | ascending `P(m = 0, 1, 2, 3)` | `ρ^+`, `ρ^−` |
|---|---|---|---|---|
| 1 | 1 / — | 0, 0.548, 0.332, 0.089 | — | 1.6067, 0 |
| 3 | 1 / 1 | 0, 0.69617, 0.26624, 0.03539 | 0, 1, 0, 0 | 1.3437, 1 |
| 5 | — / 1 | **0.59268**, 0.32762, 0.07273, 0.00672 | 0, 2/3, 1/3, 0 | 0.4943, 4/3 |
| 7 | 1 / — | 0, 0.83065, 0.15449, 0.01422 | **8/15**, 2/5, 1/15, 0 | 1.1849, 0.5333 |
| 9 | — / 1 | **0.79945**, 0.18115, 0.01847, 0.00091 | 0, 0.686, 0.286, 0.029 | 0.2209, 1.3429 |
| 11 | — / — | **0.73961**, 0.23985, 0.01850, 0.00194 | **0.5714**, 0.2857, 0.1270, 0.0159 | 0.2831, 0.5873 |
| 13 | — / — | **0.62444**, 0.33094, 0.04250, 0.00208 | **0.6465**, 0.3071, 0.0444, 0.0020 | 0.4224, 0.4020 |

- Bold entries are the hole densities. All of them were confirmed by brute force over all relevant `A`
  (drops up to `3·10^4`), to within `10^{−4}`.
- For `q = 1` the moduli `2^v − 1` form the cyclotomic lattice `gcd(2^a−1, 2^b−1) = 2^{gcd(a,b)} − 1`, so the
  dependence is much stronger.
- For `q ≥ 3` the dependence is sparse, as for `q = 3`: `gcd(2^v − q, 2^w − q) | 2^{w−v} − 1`.
- `ρ^+(1) = 1.60669515…` is the Erdős–Borwein constant.

### 5.4 The drift dichotomy

**Proposition 5.4 (PROVED, with Baker-type lower bounds for `|k log q − p log 2|` CITED; FINITE-EXACT).**
`ρ_k^+(q) + ρ_k^−(q) → 1` for `q ∈ {1, 3}`, and `→ 0` for every odd `q ≥ 5`.

*Proof.* As in Prop. 4.3. The negative binomial `P_k` concentrates at `2k`. For `q ≥ 5`, `2k < k log₂q`, so
both sums are exponentially small: the k-step map expands (drift `log(q/4) > 0`) and its drops thin out. ∎

Numbers for `k = 1, 10, 50, 100, 200`:
- `q = 5`: `1.83, 0.75, 0.246, 0.041, 0.0026`;
- `q = 7`: `1.72, 0.44, 7.8·10⁻⁴, 1.7·10⁻⁶, 8·10⁻¹²`.

### 5.5 What is specific to q = 3

| Feature | Which `q` |
|---|---|
| Units on both sides, so the one-step drop map is onto `Z` (the owner's "two copies") | only `q = 3` |
| Units at `k ≥ 2` (the 2-step "two copies" of Prop. 4.4; Catalan/Gersonides `3² − 2³ = 1`) | only `q = 3` |
| Free cycles from units (`{1}`; `{−1}`; `{−5, −7}`) | `q = 3` (THM-4484) |
| Contracting drift together with `ρ_k → 1` | `q ∈ {1, 3}` |
| The divisor-problem form of the one-step law | every `q` |
| Sparse dependence `gcd | 2^{Δv} − 1` | every `q ≥ 3` |
| Equality of the `qx+1` and `qx−1` statistics | every `q` |
| Injectivity of consecutive drop pairs | every `q` tested (PROVED for `q = 3`) |

The orchestrator's parenthesis "units only for `q = 3` and `q = 1`" should read: **a unit exists iff
`q = 2^a ± 1`, and two units (one on each side) exist only for `q = 3`.**

*Remark on S15's D75 size bound (PROVED, elementary).* The crude bound `x_min ≥ (q^k − 2^k)/(2^L − q^k)`
exceeds 1 whenever `2^L < 2q^k − 2^k`. With `L = ⌈k log₂ q⌉` this holds for every `k` outside a set of density
zero, for every odd `q ≥ 3` (because `k log₂ q` is equidistributed mod 1). So the bound fails infinitely often
for every such `q`. The known cycles are classified by
THM-4484 (free versus sporadic).

## 6. The drop identity mod 3^k (goal 5; short, since the tcpc lane owns the clock)

**Proposition 6 (PROVED; FINITE-EXACT, checks F1–F5).**
- (a) 2 is a primitive root mod `3^k`, so `2^v mod 3^k` has period `2·3^{k−1}`. From
  `6K + 1 = (2^v − 3)S`: `S ≡ (6K+1)·(2^v − 3)^{−1} (mod 3^k)`. This is a function of
  `(K mod 3^{k−1}, v mod 2·3^{k−1})`. Also `K ≡ 2(A − (−1)^v) (mod 3)`, and `S ≡ (−1)^v (mod 3)`.
- (b) Mod 9, `S` is the following function of `(K mod 3, v mod 6)`. The rows are multiplications by
  `6K + 1 ∈ {1, 7, 4}`, i.e. by `2^0, 2^4, 2^2`. In row 0 the entries alternate between the classes
  `2 (mod 3)` and `1 (mod 3)`, as they must.

  | `K mod 3` | `v ≡ 1` | `v ≡ 2` | `v ≡ 3` | `v ≡ 4` | `v ≡ 5` | `v ≡ 0 (mod 6)` |
  |---|---|---|---|---|---|---|
  | 0 | 8 | 1 | 2 | 7 | 5 | 4 |
  | 1 | 2 | 7 | 5 | 4 | 8 | 1 |
  | 2 | 5 | 4 | 8 | 1 | 2 | 7 |

- (c) Over `k` steps, `Syr^k(A) ≡ 2^{−p} c(v) (mod 3^k)` **depends only on the valuation word**. The drop
  enters only at higher 3-adic order: `D = ((2^p − 3^k)Syr^k A − c(v))/(2·3^k)`, so `D mod 3^j` is fixed by
  the word and `Syr^k A mod 3^{k+j}`. (This is the k-step form of the tcpc clock "`Syr(A) mod 9` depends on
  `(A mod 3, v mod 6)`".)
- (b′) Cross-reference. The tcpc lane's note
  [`procgen_tcpc_20261001_tournament_clock_prime_collatz.md`](procgen_tcpc_20261001_tournament_clock_prime_collatz.md)
  (its G1) writes residues prime to 3 as `x ≡ (−1)^c 4^s (mod 9)`. In those coordinates one Syracuse step is
  `(c, s) ↦ (v mod 2, 1 + c + v mod 3)`. Check F5 confirms this for all odd `A < 2·10^5` prime to 3. It is
  consistent with the table above, which is the same map read through the drop, since `A = S + 2K`.
- (d) **Drop residues carry no multiplicity information.** Every `M_v` is prime to 3, so by CRT the density
  of `{d ≡ r (mod 3^j) : m(d) = i}` is `δ_i/3^j`. In particular every residue `6d+1 mod 9` admits every `v`.
  Measured to `10^{−3}` on `[0, 10^7]`.

## 7. For the owner (plain language)

Your corrected map is the Syracuse map (S15). Here is what the differences `K` do.

1. **Every whole number is a difference.** Over the positive odd numbers, every negative difference occurs
   exactly once (the steps that go up). Every difference `d ≥ 0` occurs:
   - exactly once for 69.6% of the numbers `d`;
   - twice for 26.6%;
   - three times for 3.5%;
   - four times for 0.2%.

   The average is 1.3437, not 2. Each extra copy comes from one of the numbers `5, 13, 29, 61, 125, 253, …`
   (that is, `2^v − 3`) dividing `6d+1`. These numbers are almost unrelated to each other. The only overlaps
   come from ten small primes; for instance 5 divides 125 and every fourth one. That is why the percentages can
   be given exactly, to 18 digits.
2. **Your "exactly two copies" is exactly right once negative numbers are allowed.** Every difference `d`,
   positive or negative, is produced by the two numbers `8d+1` and `−4d−1`, one positive and one negative.
   In your labels `M`, these are exactly one even label and exactly one label `≡ 1 (mod 4)`. All the extra
   copies sit in the labels `≡ 3 (mod 4)`.
   The extra copies come on top of these. The two basic copies come from `4 − 3 = 1` and `2 − 3 = −1`.
   Among all multipliers `q` in `qx+1`, only `q = 3` has such a pair, so it is the only one for which every
   whole number is a difference. For `5x+1`, 59% of the possible differences never occur.
3. **Records.**
   - The first `d` with 6 copies is 479104.
   - The first with 10 copies is `115 043 252 496 711 229`.
   - The first with 12 copies is about `4.8·10^{23}`, and the first with 13 copies is about `2.4·10^{27}`.

   The record count grows at least like the square root of the number of binary digits of `d`. That much is
   proved; by the numbers it grows no faster.
4. **Two steps identify the start.** One difference leaves a few candidates, but the two consecutive
   differences `K(A)` and `K(S(A))` determine `A` uniquely.
5. **Ups and downs cancel on average.** Summing all the differences for the odd numbers below `2^n` gives
   exactly `−1, −1, −3, −5, −11, −21, −43, …` (the Jacobsthal numbers). The total is tiny compared with the
   numbers involved. Collatz shrinks numbers only in the typical (multiplicative) step, `×3/4`, not on
   average.
6. **The differences cannot rule out cycles.** Over `k` steps a total drop of `D` comes from about one
   starting number on average. A cycle is exactly `D = 0`. But the `3x−1` map, which has the cycles
   `1`, `5 → 7 → 5` and `17 → … → 17`, has exactly the same difference statistics as `3x+1`. So no
   statistic of the differences can decide whether `3x+1` has other cycles.

## 8. Directions (new)

- **D76 (OPEN).** Prove `log₂ lcm(2^v − 3 : v ∈ T) ≥ (½ − o(1))|T|²`, or even `≥ c|T|²`. This gives
  Conjecture 2.5(d). It is a Zsigmondy-type question for `2^n − 3`: is the part of `2^n − 3` shared with
  earlier terms `2^{o(n)}`? Data: for `3 ≤ v ≤ 64` the overlap `log₂ Π M_v − log₂ lcm M_v` is 141.8 of
  2075.7 bits (6.8%).
- **D77 (PARTLY ANSWERED).** `k = 2, 3` are done (Thm 4.8, numerically to `10^{−12}`). Open: closed forms or
  an exact rational description, and general `k`.
  - The same engine handles `k = 4` in seconds, but it needs about 600 MB, so it is not in the runner. Scratch
    runs give `P(μ_4 = 0) ≈ 0.32797` from the engine and `0.3281` from a sieve on `[0, 10^5]`.
  - The moduli `2^p − 9 = (2^{p/2} − 3)(2^{p/2} + 3)` for even `p` contain the `k = 1` moduli, so the laws for
    different `k` are coupled.
  - **Poisson question (OPEN).** As `k → ∞` along non-resonant `k`, does the law of `μ_k` converge to
    `Poisson(1)`? Scratch sieves on `[0, 10^5]` give laws close to, but measurably different from,
    `Poisson(ρ_k^+)`: `P(μ_4 = 0) = 0.328` against `0.361`, `P(μ_6 = 0) = 0.300` against `0.314`, and
    `P(μ_7 = 0) = 0.385` against `0.397` (EMPIRICAL).
  - Caveat: for `k ≳ 8`, such small windows are contaminated by orbits that have already reached 1 (they make
    `μ = 1` too likely), so the densities must come from the engine and not from sieves.
- **D78 (OPEN).** Prove Theorem 4.5 for every odd `q`; it is empirically true for `q ≤ 13`. Characterise
  `Z`-collisions: across the two sheets, `±1` share the pair `(0,0)`.

## 9. Reproduction and files

```bash
python3 04-computation/experiments/procgen_cdiff_20261001_run.py   # about 3 min, prints ALL CHECKS PASSED
```

- Note: `05-knowledge/results/procgen_cdiff_20261001_drop_multiplicities.md` (this file).
- Runner: `04-computation/experiments/procgen_cdiff_20261001_run.py` (95 checks; sections A–F as above). The
  elimination-engine computations of §4.8 run first, in a child process (`--engine`), to keep the memory low.
- Output: `05-knowledge/results/procgen_cdiff_20261001.out`.

Independence: the section-B densities use only trial factorization of the pairwise gcds (all `≤ 431`). sympy
is used for a cross-check of the shared primes, for multiplicative orders, and for factoring gcds in the
other-`q` laws.
