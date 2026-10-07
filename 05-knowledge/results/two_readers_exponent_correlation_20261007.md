# Two readers of one tape: how the exponents of the coupled Collatz orbits correlate, at every gap

Session opus-2026-10-07-S20.
* **Owner prompt:** "figure out how the two orbits' exponents correlate over long gaps. consider all the openai math results to be sufficiently correct and verified, no need to double check, just explore heavily for possible connections to leverage, more papers are attached, but also pick out a few of your own from that repo to also explore and merge together with our ideas".
* **The two orbits** are S19's coupled pair `y` and `x = 3·2^v·y + 1` (HYP-9217).
* **Script:** `04-computation/experiments/two_readers_correlation_20261007.py` (+ `.out`; ALL CHECKS PASSED, 15 checks, about 85 s on 12 cores).
* **Canon:** THM-4565.
* **Concurrent work.** mac-mini worked the same prompt (checkpoint `9dc16269c`: THM-4564, HYP-9218). Their statements 2–3 are the `λ = 0` case of THM-4565 (2)–(3). The two derivations are independent. Section 4 states what each adds.

**Status.**
* PROVED (elementary, 2-adic): Theorems 1–5 and the merge remark.
* FINITE-EXACT: kernel enumeration; pointwise identities at 8.1M overlaps and 13.0M non-overlapping neighbours.
* NUMERICAL: the depth law by gap, covariance tables, near-synchronous coupling.
* CITED: the openai/math preprints, taken as correct per the owner.
* DICTIONARY / ANALOGY: the paper connections in section 5.
* Collatz OPEN; HYP-9217 and HYP-9218 OPEN.

## 0. The answer in five lines

1. **One tape, two readers.** Both exponent streams read the same 2-adic digits. At every pair of steps `(s, k)`, the reader whose window starts later on the tape is fresh (Haar) given the joint past. The other re-reads, and its exponent is `min(fresh exponent, M)` with a past-written **saturation depth** `M` (Theorems 1–2).
2. **Exact covariance at every gap.** `Cov(A_s, B_k) = E[κ(M); windows overlap]` with `κ(M) = 2 − 6·2^−M` (Theorem 4). Non-overlapping windows never correlate. Overlapping windows correlate exactly as much as the depth `M` deviates from `Geom(1/2)`: since `E_Geom[κ] = 0`, a Haar depth means no correlation.
3. **Haar depth means independence, not just zero correlation (Theorem 3(d)).** The re-read exponent is the first disagreement of the fresh digits with a past-written string. It is independent of the fresh exponent iff that string's first 1 is `Geom(1/2)`. Even then, *given the past* the two exponents share `2(1 − 2^−M)` bits, `4/3` bits on average. The streams are pairwise independent but conditionally dependent.
4. **Long gaps.** At time lags `|d + 1/2| ≥ 6` the depth on overlaps is `Geom(1/2)` to four decimals: 3.7M overlaps, never infinite, joint law within 0.0006 of the product. Every covariance at bit-gap `|L| ≥ 9` is zero within noise, at all time lags tested. Two thirds of these coupling events are offset overlaps (`λ ≠ 0`) that THM-4564's alignment law does not cover; they behave the same way.
5. **All correlation is near synchrony** (`|L| ≤ 8`). It sits at the re-reading lag `d ≈ L/2`, with parity-dependent sign and size up to 0.18. Plus the merge itself: one step before a merge the pair is a sibling pair, `L` even and nonzero, 4548 merges at `L = −2` against 2553 at `L = +2`.

## 1. Setting

* `e(z) = v_2(3z+1)` and `U(z) = (3z+1)/2^e(z)` on odd 2-adic integers.
* `y` is Haar on the odd 2-adic integers; `v ≥ 2` with `P(v = k) = 2^−(k−1)`; `x = 3·2^v y + 1`.
* `y_s = U^s(y)`, `x_k = U^k(x)`, with exponents `A_s = e(y_s)` and `B_k = e(x_k)`.
* Partial sums: `S_y(s) = Σ_(i<s) A_i` and `S_x(k) = Σ_(i<k) B_i`.
* **S19's debt relation:** `x_t = 2^(L_t) y_(t+1) + Δ_t`, with `L_t = v + S_y(t+1) − S_x(t)` and `D_t = 3Δ_t + 1 − 2^(L_t)`.

**The tape.** `x` is a function of the digits of `y`, so both orbits read one digit string. In `x`'s frame:
* `y_s` reads `W_y(s) = (S_y(s) + v, S_y(s+1) + v]`;
* `x_k` reads `W_x(k) = (S_x(k), S_x(k+1)]`.

**For a pair `(s, k)`:**
* offset `λ = v + S_y(s) − S_x(k)`, and `j = k + 1 − s`;
* comparison constant `E = (3x_k+1) − 2^λ 3^j (3y_s+1)`, computed 2-adically;
* `μ = v_2(E)`;
* saturation depth `M = μ − max(λ, 0)`;
* joint past `G = σ(v, A_(<s), B_(<k))`.

**Fresh and re-read exponents.**
* `F` is the exponent of the later-starting reader. `R` is the other exponent minus `|λ|`.
* If `λ ≥ 0`: `F = A_s`, `R = B_k − λ`. If `λ < 0`: `F = B_k`, `R = A_s + λ`.

## 2. Exact theory (PROVED)

**Theorem 1 (the later reader is fresh).**
* Given `G`, `y_s` is Haar on the odd 2-adic integers when `λ ≥ 0`, and `x_k` is when `λ < 0`.
* `E` is `G`-measurable.

*Proof.*
* The event `{B_(<k) = b}` is one residue class of `x` mod `2^(S_x(k)+1)`: each `{e(z) = c}` is one class of `z` mod `2^(c+1)`. Since `x = 3·2^v y + 1`, that is one class of `y` mod `2^(S_x(k)+1−v)`. Likewise `{A_(<s) = a}` is one class of `y` mod `2^(S_y(s)+1)`.
* **Case `λ ≥ 0`**, i.e. `S_x(k) ≤ S_y(s) + v`. Then `A_(<s)` already decides `B_(<k)`, so conditioning on `G` is conditioning on one `y`-class mod `2^(S_y(s)+1)`. On it `y_s = (3^s y + c_s)/2^(S_y(s))` is an odd constant plus `2·`Haar, with `c_s = Σ_(i<s) 3^(s−1−i) 2^(S_y(i))`.
* **Case `λ < 0`.** Here `x` is Haar on the coset `1 + 2^v + 2^(v+1) Z_2`. The `B`-cylinder is one class mod `2^(S_x(k)+1)` inside that coset. The `A`-cylinder is one class of `x` mod `2^(S_y(s)+1+v)`, and since `S_y(s)+1+v ≤ S_x(k)` it is decided by the `B`-cylinder. So `x_k` is Haar.
* **Measurability.** `Δ = x_k − 2^λ 3^j y_s = (3^k + c'_k − 2^v 3^(k+1−s) c_s)/2^(S_x(k))` involves only `v`, `A_(<s)` and `B_(<k)`. ∎

This is the 2-adic twin of the conditional path-flow lemma for the Thorp shuffle (section 5.1). Revealed digits are a cylinder in an independent digit array, so the unrevealed digits stay fair.

**Theorem 2 (reading identity; overlap iff `M ≥ 1`).**
* `R = min(F, M)` if `F ≠ M`, and `R ≥ M + 1` if `F = M`. In the tie case `R − M ~ Geom(1/2)` given `G` and `F = M`.
* The windows overlap iff `M ≥ 1`. If `M ≤ 0`, the earlier reader's exponent is the `G`-measurable number `μ + |λ|·1{λ<0}`.

*Proof.*
* **Case `λ ≥ 0`.** `3x_k + 1 = 2^λ 3^j (3y_s + 1) + E`, and the first term has valuation `λ + A_s`. So `B_k = min(λ + A_s, μ)` unless `λ + A_s = μ`. In the tie, `3x_k + 1 = 2^μ (3^j y_(s+1) + E/2^μ)` with both terms odd, and `y_(s+1)` is Haar given `G` and `A_s`, so the extra valuation is `Geom(1/2)` on `{1, 2, …}`.
* **Overlap.** `W_x(k)` meets `W_y(s)` iff `B_k > λ`, i.e. iff `R ≥ 1`. If `M ≤ 0`, then `B_k = μ ≤ λ`. If `M ≥ 1`, then `R ≥ min(F, M) ≥ 1`.
* **Case `λ < 0`.** Symmetric: `3y_s + 1 = 2^|λ| 3^−j ((3x_k + 1) − E)`. ∎

**Theorem 3 (memoryless competition).** On `{M ≥ 1}`, given `G`:
* (a) `Cov(F, R | G) = κ(M) = 2 − 6·2^−M`, with `κ(1) = −1`, `κ(2) = 1/2`, `κ(3) = 5/4` and `κ(∞) = 2`.
* (b) `1{R = 1} = 1{F = 1} XOR 1{M = 1}`.
* (c) `I(F; R | G) = 2(1 − 2^−M)` bits. `R` determines `min(F, M+1)` and nothing else about `F`.
* (d) Given the overlap event and `λ`, `(F, R)` is `Geom(1/2) ⊗ Geom(1/2)` iff `M` is `Geom(1/2)` given the overlap event and `λ`.

*Proof.*
* **(a)** Use `E[F min(F, m)] = 6 − 2^(1−m)(m+3)` and the tie term `2m·2^−m`. Then `E[FR] = 6 − 6·2^−m` and `E[R] = 2`.
* **(b)** Check the four cases of `(F = 1?, M = 1?)`. A tie at 1 sends `R ≥ 2`.
* **(c)** `H(min(F, m+1)) = 2 − 2^(1−m)` bits. The tie noise is independent.
* **(d), the clean picture.** Write the fresh reader's relevant digits as a fair string `u`, so `F = ` first 1 of `u`. Then `R = ` first 1 of `u ⊕ c`, where `c` is the past-written digit string of `−E/2^(λ⁺+1)`; `R` is the first disagreement of `u` with `c`.
  * `(⇐)` If `c` is replaced by an independent fair string, `u ⊕ c` is fair and independent of `u`. The joint law of `(F, R)` depends on `c` only through `M` = (first 1 of `c`), so a `Geom(1/2)` depth gives the same joint law: product.
  * `(⇒)` For `b < a`: `P(F = a, R = b) = 2^−a P(M = b)`, so `P(R = b) = P(M = b)`. For `a = b`: `P(F = a, R = a) = 2^−a P(M > a)`, so `P(M = a) = P(M > a)` for all `a`. Hence `M ~ Geom(1/2)`. ∎

**Theorem 4 (covariance at every pair).** For all `s, k ≥ 0`:

    Cov(A_s, B_k) = E[κ(M_(s,k)); M_(s,k) ≥ 1] = 6 E[(1/3 − 2^−M); overlap],

with `−P(overlap) ≤ Cov(A_s, B_k) ≤ 2 P(overlap)`.

*Proof.*
* By Theorem 1, the fresh reader's conditional mean is exactly 2. So `E[(A_s − 2)(B_k − 2) | G] = Cov(F, R | G)`, which is `κ(M)` on `{M ≥ 1}` (Theorems 2–3) and 0 on `{M ≤ 0}`, where the other exponent is `G`-measurable.
* `E[A_s] = 2` exactly. Also `κ ∈ [−1, 2]`. ∎

**Corollaries.**
* **Window-overlap identity.** `Cov(A_s, B_k) = E[(A_s − 2)(B_k − 2); overlap]`.
* **Variance of the debt walk.** Exactly, `Var(L_(t+K) − L_t) = Σ Var(A) + Σ Var(B) − 2 Σ_(t<s≤t+K, t≤k<t+K) Cov(A_s, B_k)`, with `Var(A_s) = 2` and `Var(B_k) = 2` once `S_x(k) ≥ v`. The walk is diffusive with variance 4 per step exactly when the saturation-depth defects `E[2^−M; overlap] − P(overlap)/3` sum to `o(K)`.
* **Not claimed.** Zero covariance does not by itself make the *raw* pair `(A_s, B_k)` independent. The offset class is past-measurable but correlates with both exponents. Theorem 3(d) is the exact independence statement, for `(F, R)` within an offset class.

**Theorem 5 (the static two-point function of `v_2`).**
* For `n` Haar on `Z_2` (or uniform on `[1, N]` with `N → ∞`) and fixed `h`: `Cov(v_2(n), v_2(n+h)) = 2 − 3·2^−v_2(h)`.
* *Proof.* Apply Theorem 3(a) to `F = v_2(n) + 1` and `R = v_2(n+h) + 1` with depth `v_2(h) + 1`.
* So `v_2` is a "pretentious" function whose two-point correlation at a fixed shift never vanishes. It averages to 0 exactly over a Haar shift. Theorems 2–4 say the coupled Collatz pair is `v_2` read along two affine forms of one Haar variable. The shift `E` is not fixed: the past writes it.

**Remark (merges; S19 Theorem D, both signs).**
* Let `τ` be the first `t` with `x_t = y_(t+1)`. Then `g_(τ−1) = F_b^−1 F_a` with `F_c(z) = (3z+1)/2^c`, i.e. `g_(τ−1)(z) = 2^(b−a) z + (2^(b−a) − 1)/3`.
* That lies in `BS(1,2) = Z[1/2] ⋊ 2^Z` iff `b − a` is even, and it is not the identity because `τ` is minimal.
* So `L_(τ−1) = 2k ≠ 0`, and the relation is a sibling pair: `x = 4^k y' + (4^k − 1)/3` for `k > 0`, or `y' = 4^k x + (4^k − 1)/3` for `k < 0`.

## 3. Numbers

All numbers come from `two_readers_correlation_20261007.py` (`.out`): 12,000 exact integer pairs, 3000-bit sources, 600 steps; 62.0% merged.

**(K) Kernel, exact enumeration over `z mod 2^22` (FINITE-EXACT).**
* `Cov(e(z), e(2^l 3^j z + Δ)) = κ(m)` at 20 maps, including `m = ∞` for `z ↦ z`, `16z + 5` and `64z + 21` (sibling maps).
* Truncation error is at most `2.2·10^−3`.

**(R) The leader.**
* The leader's exponent is `Geom(1/2)` in both regimes: `0.5002, 0.2500, 0.1248, 0.0623` and `0.5006, 0.2500, 0.1244, 0.0629`.
* The lagger's exponent equals its prediction from the pinned debt at every non-tie step.

**(G) Every overlapping pair `(s, k)`, `|k − s| ≤ 40`, `10 ≤ s < 560`.**

The statements "overlap iff `M ≥ 1`", the reading identity and the XOR identity hold at all 8,121,017 overlaps and 12,966,272 non-overlapping neighbours (FINITE-EXACT). The law of `M` on overlaps (NUMERICAL):

| class (`d = k − s`) | n | P(M=1..4) | P(7≤M<∞) | P(M=∞) | `E[κ]` | max dep. of (F,R) |
|---|---|---|---|---|---|---|
| `d = −1`, merged | 3,215,968 | — | 0 | 1 | +2 | — |
| `d = −1`, pre-merge | 109,590 | .524 .235 .142 .024 | .0013 | .0515 (merging steps) | −0.044 | .012 |
| `d = 0` | 126,803 | .472 .270 .180 .040 | .0026 | 0 | +0.022 | .014 |
| `d ∈ {−2, 1}` | 250,807 | .485 .261 .136 .058 | .0043 | 0 | +0.021 | .0076 |
| `2 ≤ \|d+½\| ≤ 5` | 741,631 | .510 .250 .118 .059 | .0188 | 0 | −0.022 | .0048 |
| `6 ≤ \|d+½\| ≤ 12` | 1,472,757 | **.501 .249 .124 .063** | .0157 | 0 | −0.002 | **.0006** |
| `13 ≤ \|d+½\| ≤ 40` | 2,203,461 | **.499 .250 .125 .063** | .0158 | 0 | +0.002 | **.0003** |
| Geom(1/2) | | .5 .25 .125 .0625 | 1/64 = .0156 | 0 | 0 | 0 |

**Offsets.**
* Among pre-merge overlaps, the exact alignments (`λ = 0`, the case of THM-4564) are a third: 0.331 near and 0.333 far.
* The offset overlaps are two thirds.
* Both carry the same conditional coupling: `E|κ(M)|` = 1.007 / 0.984 near and 1.001 / 1.000 far; `E[I(F;R | past)]` = 1.34 / 1.32 / 1.333 / 1.334 bits.
* Both average to zero at long gaps: `E[κ] = +0.0003` (aligned) and `+0.0011` (offset).

**Covariance at lag `d`, pooled over `s`.**
* Full covariance, overlap part and kernel prediction agree at every listed lag, e.g. `d = −1`: +0.9697 / +0.9699 / +0.9738 (pair-level s.e. 0.0016). The `d = −1` value is the merged pairs.
* For `d ∉ {−1}`, `|Cov| ≤ 0.0024` at all 22 listed lags out to `|d| = 40`, while `P(overlap)` is 0.001–0.022.

**(S) Same step `(A_(t+1), B_t)`, strictly before the merging step.**

| `L_t` | n | `P(m ≥ 1)` | `P(m=1..4 \| m≥1)` | predicted `E[κ]` | measured Cov |
|---|---|---|---|---|---|
| 0 | 31,943 | 1 | .586 .307 .062 .039 | −0.2800 | −0.2759 |
| 1–4 | 196,201 | 0.2125 | .538 .209 .201 .019 | −0.0190 | −0.0182 |
| 5–12 | 427,796 | 0.0074 | | −0.0014 | +0.0021 |
| ≥ 13 | 1,650,475 | 0.0000 | | 0.0000 | −0.0002 |

* `P(m ≥ 1 | L = l)` is about `2^−l` with a parity modulation, odd `l` higher. For `l = 1..12`: `.479 .196 .164 .045 .033 .0093 .0081 .0031 .0030 .0008 .0007 .0002`.
* **Merging steps (`D_t = 0`) by `L_t`:** `−8: 1, −6: 48, −4: 164, −2: 4548, +2: 2553, +4: 97, +6: 32, +8: 1`. All are even and nonzero (FINITE-EXACT). The `L = −2 : +2` ratio is 1.78 (NUMERICAL; unexplained).

**(P) Near-synchronous coupling by single `L_t`** (NUMERICAL; `Cov(A_(t+1), B_(t+d))`, `d = −2..8`, not merged by `t`).

| `L` | largest entries | max `\|Cov\|`/s.e. |
|---|---|---|
| −2 | `d=−1`: +0.181, `d=0`: +0.107 | 17.8 |
| −1 | `d=−2`: −0.079 | 7.8 |
| +1 | `d=1`: −0.086 | 8.8 |
| +2 | `d=1`: +0.149, `d=0`: +0.083 | 15.7 |
| +3 | `d=0`: −0.071, `d=3`: −0.054 | 7.6 |
| +4 | `d=2`: +0.070, `d=3`: +0.056, `d=4`: −0.043 | 7.5 |
| +6 | `d=4`: −0.087 | 9.6 |
| +7, +8 | | 3.5, 3.2 |
| +9 … +12 | | 2.2, 1.6, 1.2, 1.6 |

* The lag-4 entry of the pooled `L ∈ [1, 8]` row in (C), −0.030, is (up to the between-group term) the `n`-weighted average of the `L = 1..8` entries at `d = 4`. It is dominated by `L = 4..7`, and `L = 6` gives −0.087.
* For `|L| ≥ 9` the pooled rows are within ±0.006 at every lag `−8..8`. The binned covariance at `L ≥ 13` (1.65M steps) is −0.0002.

## 4. What this settles, and what it does not

**Relation to THM-4564 (mac-mini, same prompt, concurrent).**
* Their statement 2 (ultrametric law at exact alignments) and statement 3 (Haar depth iff independence) are THM-4565 (2)–(3) at `λ = 0`.
* Their lockstep cap `δ' = min(δ − a, v_2(3^k − 1))`, their one-sided causality, and their Kesten perpetuity (tail index 1, from the mean-one martingale `3^j/2^(A_j)`) are not in THM-4565. They are complementary.
* THM-4565 adds:
  * the general offset (two thirds of all coupling events);
  * the exact identity `Cov = E[κ(M); overlap]` at every pair;
  * the XOR proof;
  * the conditional-information and BSC statements;
  * the static form.
* **Scope correction for THM-4564.** Its title "coupled only at tape alignments" should read "coupled only at window overlaps". With exact alignments meaning `B_t − B_s = L_s`, the offset overlaps carry the same conditional coupling (`E|κ| = 1`, 4/3 bits).

**HYP-9218 should be stated for all overlaps.**
* HYP-9218's reduction "Haar depth at alignments ⇒ asymptotically uncorrelated increments" needs the depth to be Haar at offset overlaps too. Theorem 4 sums `κ(M)` over *all* overlaps.
* The data support the widened statement: offset overlaps at long gaps have `P(M = 1..3) = .4994, .2504, .1251` and `E[κ] = +0.0011`.
* With that change, HYP-9218 gives exactly what Theorem 4 needs: the covariance defects vanish beyond the near-synchronous band. That leaves the variance-4 law of S19's block variances (4.01–4.12) and of mac-mini's (4.00), with the near-synchronous negative same-step kernel as the only visible correction.

**What the exact theory buys for HYP-9217.**
* The increments `A_(t+1) − B_t` of `L` have exactly computable correlations: `κ`-averages of past-written depths.
* The pairwise structure is solved. The remaining content is a 2-adic equidistribution statement for the past-written comparison strings: their first 1 should be `Geom(1/2)` conditionally on the past, at bit-gap `≥ 9`.
* Conditional dependence is real (4/3 bits per overlap). A proof must therefore run in the filtration where the later reader is fresh (Theorem 1), not by independence of the streams.

## 5. Connections

The openai/math preprints are taken as correct, per the owner. Each connection is labelled with its status.

### 5.1 Thorp sweeps (openai/math, Sep 26: "Conditional information under deterministic coordinate sweeps"; "Optimal-order mixing of the Thorp shuffle") — DICTIONARY, with one exact transfer

* **Their path-flow lemma.** Conditional on observed card paths, the coins on unvisited edges stay independent and fair. The conditional law of a tagged card is *averaged* on pairs with two available positions and *transported* on pairs with one.
* **Our Theorem 1 is the same mechanism on the 2-adic tape**, with the same proof (revealed exponents form a cylinder in an independent digit array):
  * the fresh reader averages: its exponent is `Geom(1/2)`;
  * the re-reader is transported: its exponent is a function of the revealed digits, until it saturates.
  * "Observing more paths changes an averaging process into a random mixture of averaging and transport." Here that is the lagger's alternation between pinned reads and ties.
* **Their exact noise/covariance identities** carry geometric weights `2^−j`: the present energy is noise from the last sweep. Ours carry the geometric moment: `Cov = 2 − 6 E[2^−M]` on overlaps.
* **A Thorp pair is depth `M = 1`.** One fair coin orders two paired cards oppositely. At `M = 1`, Theorem 3(b) makes the readers' "exponent = 1" indicators exactly complementary (`κ(1) = −1`). Depths `M ≥ 2` are partial synchronizations with positive `κ`. Under `Geom(1/2)` they cancel the Thorp anti-correlation exactly.
* **A template for HYP-9218 (HEURISTIC).** Their full-deck proof controls *conditional* laws averaged over observed paths, with one sweep of memory. Here a "sweep" is one re-read of the pinned region (about `L/2` steps). HYP-9218 asks that the comparison string written during that re-read be fresh in its first 1. That is the 2-adic analogue of "the present energy consists entirely of noise generated during the preceding sweep".

### 5.2 Two-point Chowla (openai/math, Sep 24: "Ordinary two-point correlations of multiplicative functions") — DICTIONARY

* The preprint proves `Σ λ(a_1 n + b_1) λ(a_2 n + b_2) = o(N)`, with a log-power saving, for nonproportional affine forms. It also proves binary corrected Elliott under non-pretentiousness.
* Theorem 5 is the opposite extreme. For the maximally pretentious additive function `v_2`, the two-point function along `n, n + h` is exactly the local factor `2 − 3·2^−v_2(h)`, which never vanishes at fixed shift.
* The two Collatz orbits are `v_2` along two nonproportional affine forms of one Haar variable: `3y_s + 1` and `3x_k + 1 = 2^λ 3^j (3y_s + 1) + E`. Their correlation is that local factor at the dynamically written shift `E` (Theorems 2–4).
* **Long-gap decorrelation = "Chowla for `v_2` on average over the dynamical shift".** Pretentiousness is defeated not by the function but by the randomness of the shift. That randomness is what HYP-9218 must supply.

### 5.3 Courtade–Kumar (openai/math, Sep 24: "Sharp binary-information contraction on the discrete cube") — CITED + exact first bit

* The theorem (cited): for any Boolean `f` of uniform bits passed through i.i.d. `BSC(α)` noise, `I(f(X); Y) ≤ 1 − h_2(α)`, with equality for a dictator.
* Theorem 3(b) makes the first re-read bit exactly the fresh dictator bit `1{F = 1}` through a `BSC(α)`, with `α = P(M = 1 | overlap class)`. The fresh string is `u`, and the re-reader sees `u ⊕ c`.
* If the comparison digits `c` were i.i.d. `Bernoulli(α)`, Courtade–Kumar would make the exponent-1 indicator the most informative Boolean summary of the fresh window, with sharp budget `1 − h_2(α)`.
* Measured first-bit budgets (NUMERICAL):

  | class | bits |
  |---|---|
  | `d = −1` pre-merge | 0.0017 |
  | `d = 0` | 0.0022 |
  | `d ∈ {−2, 1}` | 0.0007 |
  | `2 ≤ \|d+½\| ≤ 5` | 0.0003 |
  | same step at `L = 0` (`α = 0.586`) | 0.0214 |
  | long gaps | 0.00000 |

* The comparison digits are not i.i.d. beyond the first 1. Only `M` matters for `(F, R)`, and Theorem 3(d) is the exact replacement. So the Courtade–Kumar statement is an ANALOGY beyond the first bit.
* The contrast is the point: zero first-bit information at long gaps, yet 4/3 bits of conditional information per overlap.

### 5.4 Rokhlin multiple mixing and pointwise multiple averages (openai/math, Sep 23 and Oct 4) — DIRECTION

* The coupling map `Ψ: (A_s)_s ↦ (B_k)_k` is a measure-preserving map of `Geom^N`.
* Long-gap pairwise independence (Theorem 3(d) + the depth law) is 2-mixing of the pair (shift, `Ψ`).
* Rokhlin's theorem (cited) upgrades 2-mixing to every order for a *single* invertible mixing transformation. The two-reader system is not one: the relation chain on `BS(1,2)` has no stationary probability, since `L` is null-recurrent and absorbed by merges.
* The needed multi-point statement is the `k`-reader version of Theorem 3 (several overlaps sharing digits). Its proof shape is the same XOR picture with several comparison strings. OPEN.

### 5.5 Thompson's group F is nonamenable (openai/math, Sep 23) — ANALOGY

* Relations `g_t` live in `BS(1,2)`, the affine germs of F. That group is amenable, and the merge is a return to the identity at polynomial rate (S19: `T^−α`).
* The coupling map `Ψ` is a countable-piece prefix replacement on merged cylinders: an element of the full group of the tail relation, which is hyperfinite.
* Finite-piece truncations live in Thompson-type groups containing F, which is nonamenable (cited). By Kesten, symmetric walks there return at exponential rates.
* Merges are therefore germ events with unbounded merge times. No finite-piece description of `Ψ` exists, consistent with the heavy merge-time tail.

### 5.6 Bernoulli convolutions (openai/math, Oct 3: "Arithmetic classification and non-Pisot singularity") — DICTIONARY

* The comparison constant is a 2-adic random sum: `Δ = (3^k + c'_k − 2^v 3^(k+1−s) c_s)/2^(S_x(k))`, with `c_s = Σ 3^(s−1−i) 2^(S_y(i))`.
* Its law is an atom at `E = 0` plus a diffuse part. The atom is the merge, an arithmetic coincidence (the sibling relation), like the Pisot singularity of a Bernoulli convolution.
* Geom depth means the diffuse part is Haar in its first 1. The preprint's one-sided approximation by finite unit sets has the shape a 2-adic singularity criterion for these sums would need. No transfer is attempted.

### 5.7 The five attached papers

* **Uniform matrix hitting points (Oct 4) — ANALOGY.**
  * Relations are 2×2 matrices `[[2^L, Δ], [0, 1]]`, and Collatz steps act by conjugation with `F_c = [[3, 1], [0, 2^c]]`.
  * A merge is a word identity. One fixed substitution detecting every nonzero formula of bounded size parallels one 2-adic valuation `μ = v_2(E)` that detects every non-identity relation at finite depth (`E ≠ 0 ⇒ μ < ∞`).
* **Universal tensor squares (Sep 24) — DICTIONARY.**
  * Under Haar depth, the self-coupling `(F, R)` is the full product law, with every pair `(a, b)` charged. In general the support is full iff `M` has full support: `(a, b)` needs `M = b` (`b < a`), `M > a` (`b = a`) or `M = a` (`b > a`).
  * The anti-diagonal Thorp pair `M = 1` is the "sign" constituent.
* **Beauville splitting (Sep 23) — ANALOGY.** Theorem 1 gives a local product (fresh × transported) at every overlap. Theorem 3(d) says the local product continues to a global product (independence) iff the "holonomy" (the depth law) is Haar.
* **Einstein zero-plane rigidity (Oct 4) — DICTIONARY, with a reversal.**
  * There, one zero mixed plane forces a product, `S^2 × S^2`. Here, one zero comparison constant (`E = 0`, the sibling relation) forces the two readers to *fuse* forever: `κ(∞) = 2`, the merge.
  * Both are rigidity from a single zero. The zero splits a manifold but fuses the readers.
* **Yau uniformization (Sep 23) — ANALOGY (weak).** Positivity everywhere forces the flat model. Haar depth at every overlap would force the "flat" model of two independent streams. Nothing transfers.

## 6. Status and next steps

| item | status |
|---|---|
| Theorems 1–5 (later reader fresh; reading identity; memoryless competition; covariance at every pair; `v_2` two-point function) | PROVED (elementary) |
| merges are sibling relations at even nonzero `L` | PROVED (S19 Theorem D, both signs) |
| kernel enumeration; pointwise identities at 8.1M + 13.0M pairs | FINITE-EXACT |
| depth law Geom at `\|d+½\| ≥ 6`; no covariance at `\|L\| ≥ 9`; offset overlaps 2/3 of coupling events | NUMERICAL |
| near-synchronous coupling, parity signs, `L = −2 : +2` merge ratio 1.78 | NUMERICAL |
| openai/math preprints | CITED (taken as correct per owner) |
| sections 5.1–5.7 | DICTIONARY / ANALOGY / DIRECTION as labelled |
| HYP-9217, HYP-9218 (to be read for all overlaps) | OPEN |

**Next steps.**
1. A `k`-reader version of Theorem 3 (several comparison strings), for higher cumulants of the `L`-walk.
2. A conditional (filtration) version of the depth law at offset overlaps: the widened HYP-9218, tested with consecutive offset depths as mac-mini did for alignments.
3. Explain the `L = −2 : +2` merge asymmetry (1.78) from the sibling-relation debt law.
4. Port the Thorp sweep's "one sweep of memory" argument to the re-read of the pinned region.

## References

**openai/math preprints** (`https://github.com/openai/math/tree/main/preprints/`, taken as correct per the owner):
* `Conditional-information-under-deterministic-coordinate-sweeps-September-26-2026`
* `Optimal-order-mixing-of-the-Thorp-shuffle-September-26-2026`
* `Compatibility-entropy-and-the-spectrum-of-a-Thorp-sweep-September-26-2026`
* `Ordinary-two-point-correlations-of-multiplicative-functions-September-24-2026`
* `Sharp-binary-information-contraction-on-the-discrete-cube-September-24-2026`
* `Rokhlins-multiple-mixing-problem-for-one-transformation-September-23-2026`
* `Pointwise-Multiple-Ergodic-Averages-for-Mixing-Transformations-October-4-2026`
* `Triple-ergodic-averages-with-distinct-integer-slopes-October-4-2026`
* `Thompsons-group-F-is-nonamenable-September-23-2026`
* `Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026`

**Owner-attached PDFs (OpenAI preprints):**
* *Uniform Matrix Hitting Points in Every Positive Characteristic* (Oct 4, 2026)
* *Universal Tensor Squares for Symmetric Groups* (Sep 24, 2026)
* *Universal-cover splitting for compact Kähler manifolds* (Sep 23, 2026)
* *Uniformization of complete Kähler manifolds with positive bisectional curvature* (Sep 23, 2026)
* *Zero-Plane Rigidity for Einstein Four-Manifolds* (Oct 4, 2026)

**Repo:**
* THM-4564 and HYP-9218 (mac-mini, checkpoint `9dc16269c`)
* HYP-9217 and S19's note `collatz_cycles_tubes_debt_walk_openai_20261006.md` (Theorem D)
* H. Kesten, *Acta Math.* 131 (1973); C. M. Goldie, *Ann. Appl. Probab.* 1 (1991) — via THM-4564
