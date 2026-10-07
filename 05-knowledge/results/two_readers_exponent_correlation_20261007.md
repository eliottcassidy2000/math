# Two readers of one tape: how the exponents of the coupled Collatz orbits correlate, at every gap

Session opus-2026-10-07-S20.
* **Owner prompt:** "figure out how the two orbits' exponents correlate over long gaps. consider all the openai math results to be sufficiently correct and verified, no need to double check, just explore heavily for possible connections to leverage, more papers are attached, but also pick out a few of your own from that repo to also explore and merge together with our ideas".
* **The two orbits** are S19's coupled pair `y` and `x = 3·2^v·y + 1` (HYP-9217).
* **Scripts:**
  * `04-computation/experiments/two_readers_correlation_20261007.py` (+ `.out`; ALL CHECKS PASSED, about 1.5 min on 12 cores);
  * `two_readers_merge_side_20261007.py` (+ `.out`).
* **Canon:** THM-4565 and MISTAKE-582 (this note's audit corrections).

**Concurrent work on the same prompt.**
* mac-mini: THM-4564, HYP-9218, THM-4569, THM-4581, and the note `oai3_two_orbits_twos_and_threes_20261007.md`.
* codex-tiling: the integration audit `collatz_overlap_kernel_integration_20261007.md`. It independently confirmed the kernel and covariance identity of THM-4565, refuted HYP-9218's full-past formulation, and supplied the hostile family quoted in §3.
* Section 4 states how the pieces fit.

**Status.**
* PROVED (elementary, 2-adic): Theorems 1–5; `M = ∞` only at `d = −1`; the merge remark (almost surely).
* Conditionally PROVED: Corollary 6, given mac-mini's THM-4581 (audit pending at the time of writing).
* FINITE-EXACT: exact rational kernel algebra, and the pointwise identities at 8.1M overlaps and 13.0M non-overlapping neighbours.
* VERIFIED: the kernel by truncated enumeration.
* NUMERICAL: the depth law by lag, the covariance tables, the near-synchronous coupling, the merge side.
* CITED: the openai/math preprints, taken as correct per the owner.
* DICTIONARY / ANALOGY: the connections in section 5.
* Collatz OPEN.

## 0. The answer

1. **One tape, two readers (Theorems 1–2).** Both exponent streams read the same 2-adic digits.
   * At every pair of steps `(s, k)`, the reader whose window starts later on the tape is Haar given the joint past.
   * The other reader's exponent, minus the offset `|λ|` between the two window starts, is `min(F, M)` away from ties. Here `F` is the fresh exponent, and `M` is a **saturation depth** that the joint past has already written.
2. **Exact covariance at every pair (Theorem 4).**
   * `Cov(A_s, B_k) = E[(2 − 6·2^−M); windows overlap]`, with `|Cov| ≤ 2·P(overlap)`.
   * Non-overlapping windows never contribute.
   * Overlapping windows contribute only through the one moment `E[2^−M | overlap]`. That moment is `1/3` for `Geom(1/2)`, but other depth laws have it too.
3. **Independence needs the whole depth law (Theorem 3(d)).**
   * Within an offset class, the fresh and re-read exponents are independent iff `M` is `Geom(1/2)` there.
   * Covariance, agreement `P(R = F)` and the average conditional information all read the same moment `E[2^−M]`. They cannot certify independence (codex-tiling's two-depth law, §3).
   * At the level of whole processes the streams are **totally dependent**: `y`'s exponent sequence determines `y`, and `v` then determines `x`. So the `B`-stream is a function of `(v, A`-stream`)`.
4. **Long gaps (NUMERICAL, ensemble average over the natural law of `(y, v)`, `s ≥ 10`).**
   * Every single time lag with `d ≥ 6` or `d ≤ −9` has a depth law on overlaps consistent with `Geom(1/2)`: Pearson `χ² < 28` on 6 dof, out to `|d| = 40`.
   * Pooled over `13.5 ≤ |d + ½| ≤ 40.5`, `P(M = 1..4) = .4993/.2502/.1252/.0625` (2.03M overlaps, `χ² = 8.2`).
   * The single-lag defect roughly halves per lag: `√(χ²/n) = .17, .085, .043, .030, .020` at `d = −4…−8`.
   * Covariances at bit-gap `|L| ≥ 9` are below detection at every tested lag.
   * This is an *ensemble* statement. On explicit cylinders (`v = 4k`), codex-tiling exhibits non-Geom depths at arbitrarily large gaps.
5. **Long times (Corollary 6, conditional on THM-4581).** Merging happens almost surely at rate `T^(−1/2)`. Then the cross-covariance function converges to `2·1{d = −1}` at rate `O(s^(−1/2) log² s)`: in the long run the two streams coincide up to an index shift, and nothing else survives.
6. **Near synchrony** (`|L| ≲ 16`, i.e. `|d| ≲ 8`). This is where all detectable coupling sits.
   * It is parity-signed and peaks at the re-reading lag, at sizes up to 0.18.
   * Merges pass, almost surely, through sibling relations at even nonzero `L`.
   * Early merges come almost all from the `L = −2` side (1128 : 24 for `τ ≤ 5`). Later ones are nearly balanced (1.14–1.27).

## 1. Setting

* `e(z) = v_2(3z+1)` and `U(z) = (3z+1)/2^e(z)` on odd 2-adic integers.
* `y` is Haar on the odd 2-adic integers; `v ≥ 2` with `P(v = k) = 2^−(k−1)`; `x = 3·2^v y + 1`.
  * Because `v` is geometric, `x` is Haar on `1 + 4Z_2`. So `B_0 = 1 + Geom(1/2)`, and `B_1, B_2, …` are i.i.d. `Geom(1/2)`, independent of `B_0`.
  * `y_s = U^s(y)` and `x_k = U^k(x)`, with exponents `A_s = e(y_s)` and `B_k = e(x_k)`.
  * Partial sums: `S_y(s) = Σ_(i<s) A_i` and `S_x(k) = Σ_(i<k) B_i`.
* S19's debt relation: `x_t = 2^(L_t) y_(t+1) + Δ_t`, with `L_t = v + S_y(t+1) − S_x(t)` and `D_t = 3Δ_t + 1 − 2^(L_t)`.

**The tape.** In `x`'s digit frame:
* `y_s` reads `W_y(s) = (S_y(s) + v, S_y(s+1) + v]`;
* `x_k` reads `W_x(k) = (S_x(k), S_x(k+1)]`.

**For a pair `(s, k)`:**
* offset `λ = v + S_y(s) − S_x(k)`, `j = k + 1 − s`, and lag `d = k − s`;
* comparison constant `E = (3x_k+1) − 2^λ 3^j (3y_s+1)`, computed 2-adically. Equivalently, `E = (3^(k+1) + c'_(k+1) − 2^v 3^j c_(s+1))/2^(S_x(k))`, where `c_s = Σ_(i<s) 3^(s−1−i) 2^(S_y(i))` and `c'_k` is the same sum built from the `B`'s.
* `μ = v_2(E)`;
* saturation depth `M = μ − max(λ, 0)`;
* joint past `G = σ(v, A_(<s), B_(<k))`.

**Fresh and re-read exponents.**
* `F` is the exponent of the later-starting reader. `R` is the other exponent minus `|λ|`.
* If `λ ≥ 0`: `F = A_s`, `R = B_k − λ`, `w = 3^j(3y_s + 1)/2`, `c = −E/2^(λ+1)`.
* If `λ < 0`: `F = B_k`, `R = A_s + λ`, `w = (3x_k + 1)/2`, `c = E/2`.
* In both cases `F = 1 + v_2(w)`, `R = 1 + v_2(w − c)` and `v_2(c) = M − 1`.

## 2. Exact theory (PROVED)

**Theorem 1 (the later reader is fresh).** Given `G`:
* `y_s` is Haar on the odd 2-adic integers when `λ ≥ 0`, and `x_k` is when `λ < 0`.
* The number `w` above is Haar on `Z_2`.
* `E` is `G`-measurable.

*Proof.*
* `{B_(<k) = b}` is one residue class of `x` mod `2^(S_x(k)+1)`, which is one class of `y` mod `2^(S_x(k)+1−v)`. Likewise `{A_(<s) = a}` is one class of `y` mod `2^(S_y(s)+1)`.
* **Case `λ ≥ 0`.** `A_(<s)` decides `B_(<k)`, and `y_s` is an odd constant plus `2·`Haar on the `A`-class.
* **Case `λ < 0`.** `x` is Haar on the coset of `v`. The `B`-class lies inside that coset and decides `A_(<s)`, because `S_y(s)+1+v ≤ S_x(k)`.
* The formula for `E` involves only `v`, `A_(<s)` and `B_(<k)`. ∎

*Credit.* For one orbit this is the classical parity-vector bijection:
* Terras 1976 and Everett 1977;
* Lagarias 1985 (measure-preserving on `Z_2`);
* the Bernstein–Lagarias conjugacy map (1996);
* Tao 2022 (Syracuse valuations approximately i.i.d. geometric).

Theorem 1 is its two-orbit extension. The proof has the same shape as the Thorp-shuffle path-flow lemma (§5.1).

**Theorem 2 (reading identity).**
* (a) `R = min(F, M)` if `F ≠ M`, and `R ≥ M + 1` if `F = M`. In the tie case `R − M ~ Geom(1/2)` given `G` and `F = M`.
* (b) The windows overlap iff `M ≥ 1`. If `M ≤ 0`, the earlier reader's exponent is the `G`-measurable number `μ + |λ|·1{λ<0}`.
* (c) `M = ∞` (that is, `E = 0`) is possible only at `j = 0`, i.e. `d = −1`.

*Proof.*
* **(a, b)** `F = 1 + v_2(w)` and `R = 1 + v_2(w − c)`, with `w` Haar and `v_2(c) = M − 1`. The overlap condition is `B_k > λ` (case `λ ≥ 0`) or `A_s > |λ|` (case `λ < 0`).
* **(c)** `E = 0` means `3^(k+1) + c'_(k+1) = 2^v 3^j c_(s+1)`. Mod 3, `c'_(k+1) ≡ 2^(S_x(k))` and `c_(s+1) ≡ 2^(S_y(s))`.
  * For `j ≥ 1` the left side is `≢ 0` while the right side is `≡ 0`.
  * For `j ≤ −1`, multiply by `3^(−j)`: the left side becomes `≡ 0` and the right side does not. ∎
* (c) is due to the audit. FINITE-EXACT: no `M = ∞` at `d ≠ −1` among 8.1M overlaps.

**Theorem 3 (memoryless competition).** On `{M ≥ 1}`, given `G`:
* (a) `Cov(F, R | G) = κ(M) = 2 − 6·2^−M` (`κ(1) = −1`, `κ(2) = 1/2`, `κ(∞) = 2`), and `E[R | G] = 2`.
* (b) `1{R = 1} = 1{F = 1} XOR 1{M = 1}`.
* (c) `I(F; R | G) = 2(1 − 2^−M)` bits and `P(F = R | G) = 1 − 2^(1−M)`.
* (d) Given the overlap event and `λ`, `(F, R)` is `Geom(1/2) ⊗ Geom(1/2)` iff `M` is `Geom(1/2)` given the overlap event and `λ`.

*Proof.*
* **(a)** `E[F min(F,m)] = 6 − 2^(1−m)(m+3)`, and the tie term is `2m·2^−m`. FINITE-EXACT for `m ≤ 30`.
* **(d) `(⇐)`** The joint law depends on `c` only through `M`, and a fair independent `c` makes `w − c` fair and independent of `w`. FINITE-EXACT: the `Geom(1/2)` mixture equals `2^−(a+r)` for `a, r ≤ 25`.
* **(d) `(⇒)`** `P(F=a, R=b) = 2^−a P(M=b)` for `b < a`, and `P(F=a, R=a) = 2^−a P(M>a)`. Independence then forces `P(M=a) = P(M>a)` for all `a`. ∎
* The three summaries in (a) and (c) are affine in `E[2^−M]`. Only the full law of `M` decides (d).

**Theorem 4 (covariance at every pair).** For all `s, k ≥ 0`:

    Cov(A_s, B_k) = E[κ(M_(s,k)); M_(s,k) ≥ 1] = 6·E[(1/3 − 2^−M); overlap],

with `−P(overlap) ≤ Cov(A_s, B_k) ≤ 2·P(overlap)`.

*Proof.*
* Given `G`, the fresh reader's conditional mean is exactly 2. So `E[(A_s−2)(B_k−2) | G]` is `κ(M)` on overlaps and 0 otherwise.
* `E[A_s] = 2` exactly. ∎

Codex independently derived the same identity and bound. The script checks it at all small `(s, k) ∈ [0,6] × [0,8]`, including `S_x(k) < v` and `λ < 0`: max `|z| = 2.41`, and `Cov(A_0, B_0) = −0.503` against the exact `P(v = 2)·κ(1) = −1/2`.

**Corollary (variance).** `Var(A_s) = Var(B_k) = 2` for all `s, k`, and both streams are independent sequences. So, exactly,

    Var(L_(t+K) − L_t) = 4K − 2 Σ_(t<s≤t+K, t≤k<t+K) Cov(A_s, B_k).

**Theorem 5 (two-point function of `v_2`; elementary, likely folklore).**
* For `n` Haar on `Z_2` and fixed `h`: `Cov(v_2(n), v_2(n+h)) = 2 − 3·2^−v_2(h)`.
* Averaging over a Haar shift gives 0.
* Theorems 2–4 are this two-point function, read along two affine forms of one Haar variable whose shift `E` is written by the dynamics.

**Corollary 6 (long times; conditional on THM-4581).**
* Let `τ` be the merge time, the first `t` with `x_t = y_(t+1)`.
* After the merge the windows coincide (`x_k = y_(k+1)`). So for `d ≠ −1` an overlap needs `τ > min(k, s−1)`.
* By Theorem 4:
  * `|Cov(A_s, B_k)| ≤ 2·P(τ > min(k, s−1))` for `k − s ≠ −1`;
  * `|Cov(A_s, B_(s−1)) − 2| ≤ 4·P(τ > s−1)`.
* mac-mini's THM-4581 (checkpoint; audit pending at the time of writing) states `P(τ > T) ≤ C·T^(−1/2)(log T)²` for this pair. Given that, the cross-covariance function converges to `2·1{d = −1}` at rate `O(s^(−1/2) log² s)`.
* NUMERICAL: the pooled `d = −1` covariance is 0.9697. The merged pairs alone contribute `2 × 3,215,968/6,600,000 = 0.975` over the same window (12,000 pairs × 550 values of `s`).

**Remark (merges; S19 Theorem D, both signs; almost surely).**
* In the Haar model, almost surely `D_(τ−1) = 0`. One step before the merge the relation is then a sibling relation, the classical `4n+1` predecessor structure:
  * `x = 4^k y' + (4^k − 1)/3` with `L = 2k > 0`, or
  * `y' = 4^|k| x + (4^|k| − 1)/3` with `L = −2|k| < 0`.
* Integer value coincidences without `D = 0` exist. Example: `y = 1, v = 2` gives `τ = 2`, `L_1 = 3`, `D_1 = −16`. So the statement is almost sure, not pointwise. This is MISTAKE-580's lesson, re-learned (MISTAKE-582).
* In the sample, all 7444 merges have `D = 0` (FINITE-EXACT).

## 3. Numbers

All numbers come from `two_readers_correlation_20261007.py` (`.out`): 12,000 exact integer pairs, 3000-bit sources, 600 steps; 62.0% merged.

**(X), (K), (Z), (R).**
* Exact rational kernel algebra (FINITE-EXACT).
* Kernel enumeration at 20 maps, including `m = ∞` (VERIFIED, error < 3e-3).
* Small-index Theorem 4 (above).
* The leader's exponent is `Geom(1/2)` in both regimes.
* The lagger's prediction from the pinned debt holds at every non-tie step, for both signs of `L`.

**(G) Overlaps, `|k − s| ≤ 40`, `10 ≤ s < 560`.**

The pointwise identities hold at all 8,121,017 overlaps and 12,966,272 non-overlapping neighbours. The law of `M` on overlaps by lag class:

| class | n | P(M=1..4) | `E[κ]` | max dep. (F,R) |
|---|---|---|---|---|
| `d = −1` merged | 3,215,968 | `M = ∞` | +2 | — |
| `d = −1` pre-merge (incl. merging steps) | 109,590 | .524 .235 .142 .024 | −0.044 | .012 |
| `d = 0` | 126,803 | .472 .270 .180 .040 | +0.022 | .014 |
| `d ∈ {−2, 1}` | 250,807 | .485 .261 .136 .058 | +0.021 | .0076 |
| `d ∈ [−6,−3] ∪ [2,5]` | 980,023 | .509 .249 .118 .060 | −0.018 | .0043 |
| `d ∈ [−13,−7] ∪ [6,12]` | 1,404,628 | .5002 .2498 .1249 .0628 | −0.0002 | .0003 |
| `d ∈ [−40,−14] ∪ [13,40]` | 2,033,198 | .4993 .2502 .1252 .0625 | +0.0015 | .0003 |
| Geom(1/2) | | .5 .25 .125 .0625 | 0 | 0 |

**Single lags** (Pearson `χ²` against `Geom(1/2)` over cells `M = 1..6, ≥ 7`, 6 dof; effect `√(χ²/n)`):

| lag | χ² | effect |
|---|---|---|
| `d = −3, −4, −5, −6, −7, −8` | 6023, 3202, 703, 180, 80, 35 | .23, .17, .085, .043, .030, .020 |
| `d = −9 … −16` | ≤ 13.3 | noise |
| `d = +2, +3, +4, +5` | 3677, 733, 224, 71 | — |
| `d = +6, +7, +8` | 22.5, 15.3, 11.4 | below threshold |

* Lags with `χ² ≥ 28` are exactly `d ∈ [−8, 5] \ {−1}`. For `s ≥ 100` the set is `[−7, 5]`, with smaller `χ²`, so part of the excess is early-time.
* Pooling hides single-lag structure (the audit's finding 9). The table above is the correct long-gap statement.

**Offsets.**
* Among pre-merge overlaps, exact alignments (`λ = 0`) are a third: 0.331 near, 0.333 far. A renewal heuristic gives 1/3: each `y`-window meets 3/2 `x`-windows on average, and one is aligned with probability 1/2.
* Counting merged pairs (all `λ = 0`), offset overlaps are about 40% of all overlaps.
* Offsets carry the same conditional coupling: `E|κ(M)| ≈ 1.00`. At long gaps their depth law matches the aligned one: `P(M=1..3) = .4994/.2504/.1251` against `.5003/.2492/.1252`.

**Covariance at lag `d`, pooled over `s`.**
* Full covariance, overlap part and kernel prediction agree at all 23 listed lags.
* `d = −1` gives +0.9697. That is the merged pairs plus the merging steps.
* Off `d = −1`, `|Cov| ≤ 0.0024` out to `|d| = 40`.

**(S) Same step, strictly before the merging step.**

| `L_t` | predicted `E[κ]` | measured Cov |
|---|---|---|
| 0 | −0.2800 | −0.2759 |
| 1–4 | −0.0190 | −0.0182 |

* At `L = 0`, `P(m = 1 | m ≥ 1) = 0.586`.
* `P(m ≥ 1 | L = l) ≈ 2^−l` with a parity modulation.
* Merging steps by `L`: `−8: 1, −6: 48, −4: 164, −2: 4548, +2: 2553, +4: 97, +6: 32, +8: 1`. All are even and nonzero.

**Merge side by merge time** (`two_readers_merge_side_20261007.out`).

| `τ` | `−` side : `+` side |
|---|---|
| ≤ 5 | 1128 : 24 |
| 6–20 | 887 : 355 |
| 21–60 | 751 : 593 |
| 61–200 | 1034 : 904 |
| 201–600 | 961 : 807 |

* The early excess is a start effect: `L_0 = v + A_0 ≥ 3`, and the Mersenne start produces `y`-sibling relations quickly.
* The late ratio, 1.14–1.19, is a residual about 3–4 s.e. from 1. It is unexplained.

**(P) Near-synchronous coupling by single `L_t`.** Cov(`A_(t+1), B_(t+d')`), not merged by `t`. Here `k − s = d' − 1`, and the re-reading lag is `d' ≈ L/2`.

| `L` | peak entries | max \|Cov\|/s.e. |
|---|---|---|
| −4 | `d'=−1`: −0.068 | 6.6 |
| −3 | `d'=0`: −0.071 | 6.9 |
| −2 | `d'=−1`: +0.181 | 17.8 |
| −1 | `d'=−2`: −0.079 | 7.8 |
| +1 | `d'=1`: −0.086 | 8.8 |
| +2 | `d'=1`: +0.149 | 15.7 |
| +3 | `d'=0`: −0.071 | 7.6 |
| +4 | `d'=2`: +0.070 | 7.5 |
| +5 | `d'=3`: −0.046 | 5.1 |
| +6 | `d'=4`: −0.087 | 9.6 |
| +7 … +12 | | 3.5, 3.2, 2.2, 1.6, 1.2, 1.6 |

* The pooled `L ∈ [1,8]` row of (C) at `d' = 4` is −0.030. Up to the between-group term, that is the `n`-weighted average of these entries.

**codex-tiling's hostile family** (cited from `collatz_overlap_kernel_integration_20261007.md`; not recomputed here).
* Take raw `y ≡ 3 mod 4` and `v = 4k`. The initial alignment sits at gap `2k` with depth law `2^d/(2^ν − 2)`, `1 ≤ d < ν`, `ν = 3 + v_2(k)`. For odd `k` this is the two-depth law `1/3, 2/3`.
* Those pairs are uncorrelated (`E[2^−M] = 1/3`) but dependent, at arbitrarily large gaps.
* Under the natural law of `v` the family has weight `2^−4k`, and my sample starts at `s ≥ 10`. So it does not contradict the ensemble tables. It does rule out any statement that the depth is Geom uniformly, or conditionally on coarse data such as `v`.

## 4. How the pieces fit, and what remains

**THM-4564 (mac-mini).**
* Index translation: their `y_t` is our `y_(t+1)`, their pair `(s, t)` is our `(s+1, t)`, and their `k` is our `j`.
* Their statements 2–3 are THM-4565 (2)–(3) at `λ = 0` with `d ≥ 0`, the case their code detects. So `λ = 0` pairs at `d ≤ −1` are also outside their stated scope. These include the same-step `L = 0` pairs with `Cov = −0.28`.
* Their lockstep minimum rule holds away from an equality/cancellation branch (codex-tiling's repair).
* Their title's "coupled only at tape alignments" should read "correlated only at window overlaps". Pairs that do not overlap are conditionally uncorrelated given `G`; they are not shown to be independent.

**HYP-9218.**
* As first stated, it conditioned the depth on the full joint past, which determines it. It is refuted (codex-tiling; this note reached the same conclusion).
* A replacement must name a coarser σ-field, or an averaged offset/lag class, and must cover offset overlaps too.
* The per-lag table above is evidence for an averaged version only. Whether such a version suffices for anything is HEURISTIC.

**HYP-9217 and THM-4581.**
* mac-mini's THM-4581 attacks the merge directly: a Lyapunov weight on the Terras-clock pair chain, almost-sure merging, and `c T^(−1/2) ≤ P(no merge by T) ≤ C T^(−1/2)(log T)²`. That bypasses the depth-law route. If it survives its audit, HYP-9217's almost-sure part and exponent are settled there.
* What this note adds is the exact correlation structure the owner asked about:
  * the covariance at every pair;
  * the depth profile by lag;
  * the near-synchronous band;
  * with THM-4581, Corollary 6's long-time limit.
* What stays open:
  * the multi-point (`k`-reader) structure;
  * a selection-free test of it (naive pairing of consecutive overlaps selects on the exponents);
  * the late merge-side asymmetry;
  * any coarse-conditional depth law.

## 5. Connections

The openai/math preprints are taken as correct, per the owner. Each connection is labelled with its status.

### 5.1 Thorp sweeps ("Conditional information under deterministic coordinate sweeps"; "Optimal-order mixing of the Thorp shuffle"; Sep 26) — DICTIONARY

* **Their path-flow lemma.** Observed paths fix a cylinder of the coin array, so unvisited coins stay fair. The conditional law of a further card is *averaged* on pairs with two available positions and *transported* on pairs with one.
* **Theorem 1 is the 2-adic counterpart.** The proof shape is the same, and classical. The fresh reader averages. The re-reader is transported while it reads revealed digits; where it reaches the fresh region, its exponent is `min(F, M)`.
* **A Thorp pair is depth `M = 1`.** One shared coin forces opposite outputs; at `M = 1`, the "exponent = 1" indicators are exactly complementary (`κ(1) = −1`).
* **Where the analogy stops.** Their block-variance contraction (Prop. 7.8) needs a sign: all pairwise covariances are `≤ 0`. Ours have none: `κ(1) < 0 < κ(m ≥ 2)`, and they cancel only on average. Their variance route does not transfer as is.
* **Template (HEURISTIC).** Their "one sweep of memory" corresponds to one re-read of the pinned region.

### 5.2 Two-point Chowla ("Ordinary two-point correlations of multiplicative functions", Sep 24) — DICTIONARY

* The preprint proves `Σ λ(a_1 n + b_1) λ(a_2 n + b_2) = o(N)` for nonproportional affine forms, with a log-power saving, and the binary corrected Elliott conjecture under non-pretentiousness.
* The multiplicative function `n ↦ z^(v_2(n))` is pretentious. Its two-point structure along `n, n + h` is the local factor of Theorem 5, which never vanishes at a fixed shift.
* The Collatz pair is `v_2` along two nonproportional affine forms with a *dynamically written* shift.
* Long-gap decorrelation is "Chowla for `v_2` on average over the shift". The averaging must come from the dynamics, and §3's hostile family shows it can fail on thin cylinders.

### 5.3 Courtade–Kumar ("Sharp binary-information contraction on the discrete cube", Sep 24) — elementary + ANALOGY

* By Theorem 3(b), the first re-read bit is the fresh dictator bit through a binary symmetric channel with crossover `α = P(M = 1 | class)`. Its information `1 − h_2(α)` is elementary (one bit).
* Measured values:

  | class | bits |
  |---|---|
  | same step, `L = 0` | 0.0214 |
  | `d = −1` pre-merge | 0.0017 |
  | `d = 0` | 0.0022 |
  | long lags | 0.0000 |

* Courtade–Kumar (cited) enters only counterfactually. If the comparison digits were i.i.d. `Bernoulli(α)`, no Boolean summary of the fresh window could carry more than `1 − h_2(α)`, and the exponent-1 indicator would be optimal.

### 5.4 Rokhlin's multiple mixing; pointwise multiple averages (Sep 23, Oct 4) — DIRECTION

* `(v−1, A_0, A_1, …) ↦ (B_0 − 1, B_1, …)` is a measure-preserving bijection of `Geom(1/2)^N`. In exponent coordinates it is the conjugate of `y ↦ 3·2^v y + 1` (Bernstein–Lagarias).
* The long-gap question asks for asymptotic independence of `f(σ^s a)` and `g(σ^k Ψ a)` as `|k − s| → ∞`.
* The preprint's Theorem 1.1, which answers Rokhlin's question, upgrades 2-mixing to every order for a *single* invertible mixing transformation. Our coupling is not one.
* The relation chain's only stationary probability is the point mass at the merged state.
* The needed multi-point statement (the `k`-reader Theorem 3) is OPEN.

### 5.5 Thompson's group F is nonamenable (Sep 23) — ANALOGY

* Relations live in `BS(1,2)`, the affine germs of F.
* Amenability alone gives no polynomial return rate: symmetric random walks on `BS(1,2)` return with probability `≈ exp(−c n^(1/3))` (Varopoulos; Pittet–Saloff-Coste). The relation chain is not a random walk, and its polynomial merge rate is a 2-adic normalisation effect (THM-4581).
* `Ψ` belongs to the full group of the lagged shift-orbit relation (`Ψ(a)_k = a_(k+1)` after `τ`), which is hyperfinite. It needs infinitely many pieces because merge times are unbounded. That alone explains it; F's nonamenability (cited) and Kesten 1959 play no role.

### 5.6 Bernoulli convolutions (Oct 3) — DICTIONARY

* For fixed `(s, k)` the law of `E` is atomic, and the atom at 0 occurs only at `d = −1` (Theorem 2(c)).
* The meaningful analogue of absolute continuity is weak convergence of the digit laws of `c` (beyond its first 1) toward Haar as the gap grows, along averaged classes.

### 5.7 The five attached papers

* **Uniform matrix hitting points (Oct 4) — ANALOGY.** Relations are 2×2 matrices `[[2^L, Δ], [0, 1]]`, updated by `g ↦ F_b g F_a^(−1)` with `F_c = [[3, 1], [0, 2^c]]`. A merge is a word identity. No further transfer.
* **Universal tensor squares (Sep 24) — DICTIONARY.** Under Geom depth, the self-coupling `(F, R)` charges every pair. In general the support is full iff `M` has full support: `(a, b)` needs `M = b`, `M > a`, or `M = a`.
* **Beauville splitting (Sep 23) — ANALOGY.** Theorem 1 gives a local product (fresh × transported) at each overlap. Theorem 3(d) is the condition for the global product.
* **Einstein zero-plane rigidity (Oct 4) — DICTIONARY, with a reversal.** One zero mixed plane forces a product. One zero comparison constant fuses the readers: `κ(∞) = 2`, the merge.
* **Yau uniformization (Sep 23) — ANALOGY (weak).**

## 6. Status and next steps

| item | status |
|---|---|
| Theorems 1–5; `M = ∞` only at `d = −1`; merge remark (almost surely) | PROVED |
| Corollary 6 (long-time limit `2·1{d=−1}`, rate) | PROVED conditional on THM-4581 |
| rational kernel algebra; pointwise identities at 8.1M + 13.0M pairs; 7444 merges with `D = 0` | FINITE-EXACT |
| kernel enumeration | VERIFIED |
| depth law by lag (Geom for `d ≥ 6` or `d ≤ −9`, defect halving per lag), offsets, covariances, near-synchronous band, merge side | NUMERICAL |
| openai/math preprints | CITED |
| §5 | as labelled |
| HYP-9218 (as first stated: refuted; replacement OPEN), HYP-9217 (see THM-4581) | — |

**Next steps.**
1. A `k`-reader Theorem 3, with a selection-free multi-point test whose pairs are chosen by events fixed before the first reading.
2. Explain the late merge-side ratio of 1.14–1.19.
3. A coarse-conditional depth law that names its σ-field.
4. Port the Thorp "one sweep of memory" to the pinned-region re-read.

## References

**openai/math preprints** (`https://github.com/openai/math/tree/main/preprints/`):
* `Conditional-information-under-deterministic-coordinate-sweeps-September-26-2026`
* `Optimal-order-mixing-of-the-Thorp-shuffle-September-26-2026`
* `Ordinary-two-point-correlations-of-multiplicative-functions-September-24-2026`
* `Sharp-binary-information-contraction-on-the-discrete-cube-September-24-2026`
* `Rokhlins-multiple-mixing-problem-for-one-transformation-September-23-2026`
* `Pointwise-Multiple-Ergodic-Averages-for-Mixing-Transformations-October-4-2026`
* `Thompsons-group-F-is-nonamenable-September-23-2026`
* `Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026`

**Owner-attached OpenAI preprints:**
* *Uniform Matrix Hitting Points in Every Positive Characteristic* (Oct 4)
* *Universal Tensor Squares for Symmetric Groups* (Sep 24)
* *Universal-cover splitting for compact Kähler manifolds* (Sep 23)
* *Uniformization of complete Kähler manifolds with positive bisectional curvature* (Sep 23)
* *Zero-Plane Rigidity for Einstein Four-Manifolds* (Oct 4)

**Classical:**
* R. Terras, *Acta Arith.* 30 (1976); C. J. Everett, *Adv. Math.* 25 (1977); J. C. Lagarias, *Amer. Math. Monthly* 92 (1985).
* D. J. Bernstein and J. C. Lagarias, "The 3x+1 conjugacy map", *Canad. J. Math.* 48 (1996) 1154–1169.
* T. Tao, *Forum Math. Pi* 10 (2022).
* L. E. Garner, "On heights in the Collatz 3n+1 problem", *Discrete Math.* 55 (1985) 57–64; M. Elia and A. Tucker, "Consecutive integers and the Collatz conjecture", *Integers* (2015), arXiv:1511.09141. Coalescence of related orbits.
* A. Kontorovich and J. C. Lagarias, arXiv:0910.1944.
* H. Kesten, *Trans. AMS* 92 (1959); N. Varopoulos; C. Pittet and L. Saloff-Coste (return probabilities on `BS(1,2)`).

**Repo:**
* THM-4564, THM-4569, THM-4581, HYP-9218 and `oai3_two_orbits_twos_and_threes_20261007.md` (mac-mini)
* `collatz_overlap_kernel_integration_20261007.md` (codex-tiling)
* S19's `collatz_cycles_tubes_debt_walk_openai_20261006.md` (Theorem D)
