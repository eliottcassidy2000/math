# Cycle gates: exponential sums and the equidistribution of Collatz carries modulo 2^p − 3^a

**Status: FINITE-EXACT + PROVED (elementary) + CITED. As a proof angle the exponential-sum/equidistribution route ranks LOW: every counting statement it reaches is already classical (Belaga's perigee bound), and a square-root barrier stops it short of excluding any cycle. Its value is diagnostic. Naive equidistribution of the gate residues is false, and the heuristic built on it fails for 3x+1: it predicts `ln P + 0.11` positive cycles of period `≤ P`. The repair is side-aware — only rational cycles whose least point is `≥ 1` may count. With that single input the model matches every control, including the published 3x+d cycle counts of Belaga–Mignotte (9 of 11 exactly, two off by one; **RESOLVED 2026-09-26**: 11 of 11). The two missing cycles are long primitive cycles found by the [sporadic lane](procgen_sporadic_20260926_free_and_sporadic_cycles.md): `d = 14303` with least element 101 and clock `(2155,1092)`, and `d = 17021` with least element 5 and clock `(2140,1088)`. Both sit near `a/p ≈ 1/2`, beyond this lane's scan (`p ≤ 250`, main family `p ≤ 600`). No HYP or THM file was created.**

* **FINITE-EXACT (census; SHEET and DRIFT controls).** Every clock `(p,a)` with `p ≤ 40` was checked for `q = 3` (860 clocks) and `q = 5` (860), by two independent methods: a sorted join of residues and orbit iteration over the size range. They agree. The integral periodic points are exactly the known cycles.
  * `3x+1`: `0`, `{1,2}`, `{-1}`, `{-5,-7,-10}`, `{-17,…,-136}`. The last three are the `3n−1` cycles on the side `3^a > 2^p` (SHEET).
  * `5x+1`: `0`, `{1,3,8,4,2}`, `{13,…}`, `{17,…}`, `{-1,-2}` (DRIFT).
* **PROVED / certified machinery.**
  * `S(h) = Σ_w e(h c_w/|G|)` by an `O(pa)` DP, certified against brute force (`6·10^{-13}`) and a 40-digit DP (`2·10^{-17} C`).
  * The counting identity `N = |G|^{-1} Σ_h S(h)` holds *exactly* in `F_ℓ` on 164 clocks.
  * The rotation law `c_{Rw} = c_w/2` or `(q c_w + G)/2`. It gives `c_{R^j w} ≡ 2^{-j} q^{k_j(w)} c_w` and a transport identity for `S`.
* **The gate residues are not equidistributed** (Proposition A; PROVED for `a/p ∉ [1/2, 2/3]`, numerical at every ratio tested).
  * At fixed `θ = a/p ≠ log_3 2`, `S(h)/C → E e(hY_θ)`, where `Y_θ` is the real law of the rational periodic points.
  * The limits are non-zero: `|E e(Y_θ)| = 0.24, 0.07, 0.015, 0.09, 0.50` at `θ = 0.30, 0.45, 0.55, 0.75, 0.90`.
  * The excess mass sits next to integers the perigee window forbids (0 on the dyadic side, 1 on the 3-adic side).
* **Structure of the large spectrum** (FINITE-EXACT, `p ≤ 30`). The top frequency is `h ≡ u 2^{-j} q^k (mod |G|)` with `|u| ≤ 3` on 99 of 100 (`q = 3`) and 66 of 67 (`q = 5`) clocks with a computed full spectrum. For a random `h` the minimal `|u|` is typically `10`–`10^4`, depending on `M`.
  * The dominant peak is usually `h ≡ ±3^{-k}`, because the carry's ternary digit `c_w ≡ 2^{s_{a-1}} ≢ 0 (mod 3)`.
  * `h = 1` is the archimedean peak and `h ≡ ±2^{-j}` the binary one. They are one family under rotation, since `2^p ≡ 3^a`.
  * Peak height falls from `0.45 C` (far from the critical line) to `0.05 C` (`amp ≥ 10`).
  * Off the structured set `|S(h)|` has RMS `0.97 √C` (median), which is the uniform value, but a heavy tail.
  * Coincidences run below the uniform model: 0.43 of it far from the line, rising to 0.96 at the critical line.
* **The naive heuristic fails; the side-aware repair works** (Propositions N and P).
  * PROVED: `Σ_{p≤P} Σ_{2^p>3^a} L(p,a)/|G| ≥ ln P − O(1)`; numerically it equals `ln P + 0.111`. This predicts 26 sporadic positive `3x+1` cycles below the verified period `2.18·10^11`, where there are none (Poisson probability `4·10^{-12}`).
  * All of the divergence sits on clocks whose rational cycles dip below 1.
  * Counting only necklaces with least point `≥ 1` gives finite totals: 1.81 (`3x+1`), 1.21 (`3x−1`), 2.77 (`5x+1`), 0.65 (`5x−1`). Observed: 0, 1, 3, 0.
* **Dense divisor regime against published data (FINITE-EXACT).** For the eleven `d` of Belaga–Mignotte 2006, table (20), we enumerated exactly the primitive cycles of `T_d(y) = y/2, (3y+d)/2` on all clocks `p ≤ 250`, plus the main family to `p ≤ 600`.
  * The totals reproduce their `ω(d)` for 9 of the 11 `d`. For `d = 14303` and `17021` we find one fewer (943, 257). **RESOLVED 2026-09-26** by the sporadic lane: each has one more primitive cycle with period `p > 2000` (least elements 101 and 5), outside this scan. With them, all 11 match.
  * Near the critical line, equidistribution modulo `G/d` predicts every clock within noise.
  * On far clocks the multiples are depleted by factors up to 10.
* **Proof angle.** The following are all dominated or classical:
  * The perigee bound `m ≤ 1/(2^{p/a} − 3)` (Belaga 2003, via Baker/Rhin) already gives polynomially many cycles per clock.
  * A Parseval/second-moment bound cannot certify `N` below about `√C = 2^{0.475p}`.
  * The switching (Weyl-type) certificate is rigorous but decays only like `e^{-0.1p}`.
  * What is provable: equidistribution modulo fixed `m` coprime to `2q` (PROVED here, rate `E cos(π/m)^K`), and plausibly Tao-type decay off the structured set. None of it reaches the gate.

Session `collatz-procgen-20260922`, lane "gates" (mac-mini), 2026-09-25.
* Scripts (`04-computation/experiments/`): `procgen_gates_20260925_{core,expsum,stats,heuristic,switching,run}.py`.
* Output: [procgen_gates_20260925.out](procgen_gates_20260925.out), parts A–D.
* Inputs read:
  * THM-4471 `§4`;
  * [atlas](procgen_atlas_20260924_collatz_implication_atlas.md) `§§3.4–3.5`;
  * [mod-192 note](collatz_procgen_20260924_inverse_tree_mod192.md) `§§2.1–2.4`;
  * [gates note](collatz_mod6_20260921_pillai_convergents_cycle_gates.md).

---

## 0. Setting and conventions

* A **word** `w` of length `p` has `a` ones, at positions `s_0 < … < s_{a-1}`.
* Its **carry** is `c_w = Σ_i q^{a-1-i} 2^{s_i}`, and the **gate** is `G = 2^p − q^a`, `M = |G|`.
* For `T(x) = x/2, (qx+1)/2` one has `2^p T^p(x) = q^a x + c_w`.
  * Proof in one line: `y_t = 2^t T^t(x)` satisfies `y_{t+1} = q y_t + 2^t` at an odd step.
  * So the unique rational periodic point with word `w` is `x_w = c_w/G`.
  * This is THM-4471 `§4`, and the gates note's `B(w)` for words that start with 1.
* For `G < 0` the point `−x_w > 0` is the periodic point of `x/2, (qx−1)/2` with the same word (mod-192 note, Theorem 6).
* **The divisibility `G | c_w` is the same statement on both sheets.** The sign of `G` only decides whose positive cycle an integral point is:
  * `qx+1`'s when `2^p > q^a` (*dyadic* clocks);
  * `qx−1`'s when `q^a > 2^p` (*`q`-adic* clocks).
* Further notation:
  * `C = C(p,a)` is the number of words, `L(p,a)` the number of primitive necklaces (Lyndon words);
  * `N(p,a) = #{w : M | c_w}`;
  * `amp = max(2^p, q^a)/|G|` measures closeness to the critical line `a/p = log_q 2`.

## 1. The exponential sum (A1, A2, A5)

**The DP.**
* `e(h c_w/M) = Π_i e(h q^{a-1-i} 2^{s_i}/M)`, so `S(h)` is the last entry of a product of `p` bidiagonal `(a+1)×(a+1)` matrices, one per position; the state is the number of ones placed so far.
* This costs `O(pa)` per `h` and is vectorised over `h`.
* Certification:
  * against brute force, 896 (clock, `h`) pairs with `p ≤ 12`, `q ∈ {3,5}`: max error `5.7·10^{-13}`;
  * against a 40-digit mpmath DP, 63 pairs up to `(84,53)`: max error `1.8·10^{-17} C`.
* In the normalised version (`§4`) the double-precision values of `|S(1)|/C` down to `5·10^{-17}` at `(1054,665)` agree with a 40-digit DP to 6 digits (C4).

**Counting identity.** `N(p,a) = M^{-1} Σ_{h mod M} S(h)`.
* It was checked *exactly* on 164 clocks with `M ≤ 20000`. The same DP was run in `F_ℓ` for a prime `ℓ ≡ 1 (mod M)`; the reduction `Z[ζ_M] → F_ℓ` is a ring map, so the identity holds there exactly.
* In floats it was checked on 236 clocks with `M ≤ 2·10^5`: error `4·10^{-15}`.

**Rotation law (PROVED).** Let `R` move the first letter to the end.
* If `w_0 = 0` then `c_{Rw} = c_w/2`. If `w_0 = 1` then `c_{Rw} = (q c_w + G)/2`.
  * *Proof:* if `s_0 = 0`, `Rw` has its ones at `s_1 − 1, …, s_{a-1} − 1, p − 1`, so `c_{Rw} = (q(c_w − q^{a-1}) + 2^p)/2`. ∎
* Hence `c_{R^j w} ≡ 2^{-j} q^{k_j(w)} c_w (mod M)`, where `k_j(w)` is the number of ones among the first `j` letters.
* Since `R^j` permutes `W(p,a)`,

  `S(2^{-j} q^k h) = Σ_{v : k ones among the last j letters} e(h c_v/M) + Σ_{w : k_j(w) ≠ k} e(2^{-j} q^k h c_w/M).`

  A large coefficient at `h` is therefore transported to `2^{-j} q^k h`.
* Checked on all 65476 words with `p ≤ 14`, `q ∈ {3,5}`, and on 2160 instances of the identity (error `2·10^{-12}`).
* A coincidence `c_w ≡ c_{w'} ≢ 0` survives rotation exactly as long as the words agree letter by letter. At the first disagreement the residues become `2^{-1}r` and `2^{-1}qr`, which differ because `gcd(q−1, M) = 1`.

## 2. Census and controls (A3, A4)

`N(p,a)` was computed for every clock with `p ≤ 40`, `q ∈ {3,5}`, in two ways:
* a sorted meet-in-the-middle join of prefix and suffix residues, on every clock with `M < 2^62`: 859 clocks for `q = 3`, 755 for `q = 5`;
* orbit iteration of every integer in `[c_min/M, c_max/M]`, where `c_min = (q^a − 2^a)/(q−2)` and `c_max = 2^{p-a} c_min`, on all 860 clocks.
  * An orbit that leaves this interval is discarded: every point of an integral cycle on the clock is itself a periodic point of the clock.

The two methods give identical counts and points. Every `N(p,a)` is the sum, over the cycles below, of their periods at the multiples of their clocks:

| map | clock | `G` | cycle |
|---|---|---|---|
| `3x+1` | (2,1) | 1 | `{1,2}` (forced, `|G| = 1`) |
| | (1,1) | −1 | `{-1}` = the `3n−1` cycle `{1}` |
| | (3,2) | −1 | `{-5,-7,-10}` = `{5,7,10}` of `3n−1` |
| | (11,7) | −139 | `{-17,…,-136}` = `{17,…}` of `3n−1`, the only *sporadic* cycle (`|G| > 1`) |
| `5x+1` | (2,1) | −1 | `{-1,-2}` (forced) |
| | (5,2) | 7 | `{1,3,8,4,2}` |
| | (7,3) | 3 | `{13,33,83,208,104,52,26}` and `{17,43,108,54,27,68,34}` |

Plus the fixed point `0` on every `(p,0)`.

**Main term against counts.**
* The clocks with the largest word-level main term `C/|G|`:
  * `(19,12)`: 7.04, `N = 0`;
  * `(8,5)`: 4.31, `N = 0`;
  * `(11,7)`: 2.37, `N = 11`;
  * `(5,3)`: 2.00, `N = 0`;
  * `(27,17)`: 1.66, `N = 0`.
* Integral words come in multiples of the period, so the natural unit is the necklace level `L/|G|` (`(11,7)`: 0.216).
* Summed over `p ≤ 40`, the word-level main term is 51.7 on the dyadic side and 21.2 on the 3-adic side. The dyadic sum gains about 1 per period; the reason is in `§5`.

## 3. Equidistribution statistics, `p ≤ 30` (part B)

**What was computed.** For every clock with `p ≤ 30` (`q = 3`) or `p ≤ 24` (`q = 5`), the whole multiset `{c_w mod M}` was generated, streamed and sorted in residue buckets (up to `1.55·10^8` words per clock). From it:
* the zeros `N`;
* the coincidences `Coll − C = #{(w ≠ w') : c_w ≡ c_{w'}}`, against the uniform mean `C(C−1)/M`;
* the star discrepancy `D*` of `{c_w/M}`;
* the near-integer ratios `λ_δ = #{x_w ≥ 1/2, ||x_w|| ≤ δ}/(2δC)`;
* the share of words with `x_w < 1`;
* the perigee statistics of `§5`;
* where feasible, the whole spectrum: an FFT of the residue histogram for `M ≤ 2^20`, or the blocked DP for near-critical clocks (`amp ≥ 1.5`) with `M ≤ 2^23`;
* otherwise the DP on the structured set `Σ_3 = {±u 2^{-j} q^k : 1 ≤ u ≤ 3, j ≤ p, k ≤ a}` and on 1000 random `h`.

**Aggregates for `q = 3`, clocks with `C ≥ 200` (B3).**

| amp | clocks | median `max|S|/C` | median `|S(1)|/C` | median `D*` | top frequency with `|u| ≤ 3` | `Σ(Coll−C) / Σ C(C−1)/M` |
|---|---|---|---|---|---|---|
| [1, 1.05) | 221 | 0.447 | 0.251 | 0.220 | 44/44 | 0.43 |
| [1.05, 1.5) | 69 | 0.238 | 0.015 | 0.0090 | 27/27 | 0.70 |
| [1.5, 3) | 25 | 0.203 | 0.0018 | 0.0018 | 17/17 | 0.83 |
| [3, 10) | 11 | 0.158 | 0.0012 | 0.0023 | 8/8 | 0.86 |
| ≥ 10 | 5 | 0.047 | 0.0021 | 0.0049 | 3/4 | 0.96 |

For `q = 5` the same columns are: `max|S|/C` 0.57, 0.23, 0.15, 0.07, 0.05; top frequency 66 of 67; coincidence ratios 0.38, 0.68, 0.87, 0.95, 0.98.

What the numbers show:

1. **The top of the spectrum is structured.**
   * Write each top frequency as `h ≡ u 2^{-j} 3^k (mod M)` with `|u|` minimal over the transport box `0 ≤ j ≤ p`, `0 ≤ k ≤ a` (`§1`). Then `|u| ≤ 3` on all but one full-spectrum clock.
   * The most frequent top frequency has `u = ±1` and `j = p`, i.e. `h ≡ ±3^{k-a}`, a power of `3^{-1}`.
   * Example: at `(16,5)` the top is `h = 21764 = (M−1)/3 ≡ −3^{-1}` with `|S|/C = 0.41`; `h = 1` is second (0.23).
   * *Mechanism:* for `h = (M−1)/3`, `e(hc_w/M) = e(c_w/3) e(−x_w/3)`. Now `c_w ≡ 2^{s_{a-1}} ≢ 0 (mod 3)` (all words, D1), so `|E e(c_w/3)| ≈ |e(1/3) + e(2/3)|/2 = 1/2`, and far from the critical line `x_w` is small. This is the coarse 3-adic irregularity that Tao observes for `Syrac(Z/3^n)`: it never takes values divisible by 3 (`§8`).
   * `h = 1` is the archimedean peak (the sizes `x_w`).
   * `h ≡ ±2^{-j}` is the binary peak: `c_w ≡ w_0 (mod 2)`, and `c_w mod 2^j` is fixed by the first `j` letters.
   * `2^p ≡ 3^a (mod M)` and the transport identity make the three one family.
   * Among the top four frequencies, a minority have large `|u|` (e.g. `u = ±18431` at `(24,15)`). We did not identify their structure; a sum-of-two test is uninformative because `|Σ_3|² > M`.
2. **Energy.**
   * `Σ_3` carries a median 9% of the Fourier energy `Σ_{h≠0}|S(h)|² = M·Coll − C²` (far and near alike). `Σ_30` carries 16–17%, on about 1% of the frequencies.
   * Large coefficients remain outside `Σ_30`: the maximum there is a median `4.9` (far) and `7.9` (near) times `√(C ln M)`, up to 21.8 at `(27,17)`, where the uniform model gives about 1 (B4).
   * Example: at `(27,17)` the second peak has `u = −37`.
   * Reason: the archimedean law has slowly decaying Fourier coefficients (`§4`), so structured frequencies reach fairly large `|u|`.
3. **Off the structured set the spectrum is random-sized.**
   * On the 231 clocks without a full spectrum, the RMS of `|S(h)|` over 1000 random `h` has median `0.97 √C` (range 0.71–1.94). On the 12 near-critical ones it is 0.81–1.25.
   * The maximum over those `h` has median `7.6 √C` (up to 52), against about `2.6 √C` for iid residues. In several cases it sits at a partially structured `h` (`u = −95` at `(25,15)`, `u = −175` at `(29,18)`, against a typical `|u|` of `2·10^4`–`10^5`).
4. **Coincidences are below uniform.**
   * `Coll − C` is 0.43 of the uniform mean far from the line and 0.96 at `amp ≥ 10` (table).
   * So the average of `|S(h)|²` over `h ≠ 0` is 0.8–0.96 `C` near the critical line: the residues are slightly *more* regular than random.
   * Why: words with `c_w < M` cannot coincide, and agreements propagate along common letters (`§1`).
   * Uniform-model controls with the same `(C, M)`: coincidences equal to the mean within noise, `D* ≈ 1/√C`, and `max|S| ≈ √(C ln M)` (B4).
5. **Discrepancy.** `D* ≈ 0.22` far from the line (macroscopic non-uniformity), `2·10^{-3}` near it.

**Near-integer ratios.**
* `λ_δ ≈ 1` on near-critical clocks.
* On far dyadic clocks `λ_δ < 1`: 0.90 at `(30,15)`, where 27% of the words have `x_w < 1`.
* On far 3-adic clocks `λ_δ ≫ 1` for small `δ` (36 at `(28,26)`, `δ = 0.003`): the pile-up just above the integer 1 (`§4`).

## 4. The archimedean law (C4)

**Proposition A (limit law).** Fix `θ ∈ (0,1) \ {log_3 2}` and let `a/p → θ`.
* On the dyadic side, `x_w ⇒ Y_θ := Σ_{n≥0} 3^n 2^{-D_n}`, where `D_0 < D_1 < …` are the positions, counted from the end, of the ones of an iid Bernoulli(`θ`) sequence.
* On the 3-adic side, `|x_w| ⇒ Z_θ := Σ_{n≥0} 2^{σ_n}/3^{n+1}`, with `σ_n` the positions of the ones counted from the start.
* Hence `S(h)/C → E e(hY_θ)` (resp. `E e(hZ_θ)`) for every fixed `h`.

*Proof (sketch).*
1. `x_w = (c_w/2^p)/(1 − 3^a/2^p)`, and `c_w/2^p = Σ_n 3^n 2^{-(p − s_{a-1-n})}`.
2. The last `K = O(log p)` letters of a uniform word in `W(p,a)` are within total variation `O(K/p)` of iid Bernoulli(`a/p`).
3. The tail beyond them is small in probability. A fractional moment gives this: `E(3·2^{-g})^s < 1` for small `s > 0`, since `E log(3·2^{-g}) = log 3 − θ^{-1} log 2 < 0` (`g` a geometric gap). Conditioning on the total `a` costs a factor `O(√p)`, absorbed by taking `K = C log p`.
4. The 3-adic side is the same argument read from the start of the word. ∎

**Non-uniformity (PROVED for `θ < 1/2` and `θ > 2/3`).**
* *Dyadic, `θ < 1/2`.*
  * Write `Y = 2^{-D_0} Y'` with `Y'` independent of `D_0` and finite almost surely.
  * Then `P(Y ≤ ε) ≥ c ε^γ` with `γ = log_2(1/(1−θ)) < 1`.
  * So `Y mod 1` has unbounded density at `0^+`, some `E e(hY) ≠ 0`, and `D*` does not tend to 0.
* *3-adic, `θ > 2/3`.* A word that starts with `L` ones has `Z − 1 = (2/3)^L (Z'' − 1)` with `Z''` a copy of `Z`, so `P(Z − 1 ≤ ε) ≥ c ε^{γ'}` with `γ' = log(1/θ)/log(3/2) < 1`.
* In both regimes the extra mass sits next to an integer that is forbidden:
  * dyadic `θ < 1/2` means `p > 2a`, where every rational cycle dips below 1, so 0 is the only integer near the pile-up and it is never attained (`x_w > 0`); see `§5`;
  * on the 3-adic side `|x_w| > 1` strictly for `p > a`, because `c_min > M`.

**Numerics.**
* The normalised DP at `p = 50, …, 800` converges to Monte Carlo estimates of the limit series (`2·10^5` samples, 160 terms). At `p = 800`, `|S(1)/C − E e(Y)| ≤ 0.0056`.

| θ | `E e(Y)` or `E e(Z)` | law of `Y` or `Z` mod 1 |
|---|---|---|
| 0.30 | `0.148 + 0.192i` | `P(Y < 1) = 0.79`; `P(Y mod 1 < 0.05) = 0.19` (uniform 0.05) |
| 0.45 | `0.006 + 0.068i` | `P(Y mod 1 < 0.05) = 0.070` |
| 0.55 | `−0.008 + 0.013i` | `P(Y mod 1 < 0.05) = 0.052` |
| 0.75 | `0.044 + 0.080i` | `P(Z mod 1 < 0.05) = 0.085` |
| 0.90 | `0.459 + 0.193i` | `P(1 ≤ Z < 1.05) = 0.39`; `P(Z mod 1 < 0.05) = 0.43` |

**Near the critical line.**
* `|S(1)|/C` is about `10^{-3}` at `(19,12)` and at the random level `C^{-1/2}` up to `p ≈ 65`.
* Beyond that it decays only slowly, while the random level falls exponentially:

| clock | `|S(1)|/C` | `C^{-1/2}` |
|---|---|---|
| `(84,53)` | `9.4·10^{-12}` | `2^{-38}` |
| `(149,94)` | `4.0·10^{-10}` | `2^{-69}` |
| `(485,306)` | `1.3·10^{-13}` | `2^{-228}` |
| `(1054,665)` | `4.9·10^{-17}` | `2^{-498}` |

* On convergent clocks the archimedean coefficients therefore dwarf the random scale. **Equidistribution of `c_w mod G` can only be a statement about the arithmetic part, after the size law has been factored out.**

## 5. The random model: naive versus side-aware (C1, C2)

**Proposition P (perigee window; PROVED, elementary; Belaga 2003 CITED).** Take a positive rational cycle of `x/2, (qx+1)/2` with word `w` (`a ≥ 1`) and least point `m`.
* Multiplying `T(x_i)/x_i` around the cycle gives `Π_{odd}(q + 1/x_i) = 2^p`. The least point is odd.
* Hence `q^{a-1}/G ≤ m ≤ 1/(2^{p/a} − q)`.
* On the `q`-adic side, `max(q^{a-1}, c_min)/M ≤ m ≤ 1/(q − 2^{p/a})`.
* So every clock falls into one of three classes:

| class | dyadic side | `q`-adic side |
|---|---|---|
| **none** (no necklace has least point `≥ 1`) | `2^p > (q+1)^a`; for `q = 3`, `p > 2a` | `2^p < (q−1)^a` |
| **all** (every necklace eligible) | `q^{a-1} ≥ G` | `q^{a-1} ≥ M` or `c_min ≥ M`; for `q = 3`, every 3-adic clock |
| **transition** | otherwise: the eligible share is computed | otherwise: the eligible share is computed |

* Verification (B): the perigees were computed on every clock with `p ≤ 26` (`q = 3`) and `p ≤ 22` (`q = 5`), and on the transition clocks up to `p = 30` and `24`.
  * `q = 3`: 156 'none', 147 'all', 60 transition.
  * `q = 5`: 219 'none', 6 'all', 32 transition.
  * The classification held everywhere.
* Integral cycles on a clock are at most the odd integers in the window, which is polynomially many in `p` by Rhin/Baker (Belaga's `d k^{32}`).

**Proposition N (the naive model diverges; PROVED).**
* Statement: `E_naive(P) := Σ_{p≤P} Σ_{3^a<2^p} L(p,a)/(2^p − 3^a) ≥ ln P − O(1)`.
* *Proof:* for `a ≤ 0.6p`, `2^p − 3^a ≤ 2^p` and `L(p,a) ≥ (C(p,a) − p 2^{p/2})/p`. Hoeffding gives `Σ_{a≤0.6p} C(p,a) 2^{-p} ≥ 1 − e^{-0.02p}`. Summing `1/p` over `p` gives the bound. ∎
* Numerically `E_naive(P) − ln P = +0.114, +0.112, +0.111, +0.111` at `P = 100, 300, 1000, 4000`.
* Every nontrivial positive `3x+1` cycle has at least `2.18·10^11` T-steps (atlas `§3.5`: Hercher 2023 Cor. 29 with Barina's `2^71`, CITED). Up to that period the naive model predicts 26.2 sporadic cycles, and there are none. The Poisson probability is `4·10^{-12}`, so **the naive equidistribution model is falsified.**

**Where the divergence lives.**
* The bulk of `W(p,·)` has `a ≈ p/2`. For `q = 3` that lies on the dyadic side (`1/2 < log_3 2`), where `|G| ≈ 2^p`.
* So every period contributes about one word, or `1/p` necklaces.
* By Proposition P those clocks (`p > 2a`) carry no positive integral cycle.
* In the transition window `p/2 ≤ a < 0.63p` the eligible share is non-zero only in a band of about 10 clocks next to the 'all' boundary. It was found by Monte Carlo, validated on 36 transition clocks at `20 ≤ p ≤ 30` (max error 0.019, mean 0.003).
* The refined contribution of the window falls quickly: `1.2·10^{-4}` at `p = 100`, `6.6·10^{-8}` at `p = 300`, `3·10^{-11}` at `p = 500`.

**Refined totals.** Only necklaces with least point `≥ 1` count. Exact eligible counts for `p ≤ 30`, the proven classes elsewhere, and Monte Carlo on the transition band up to `p = 1000`:

| sheet | naive, `P = 40` | refined, `P = 40` | refined, `P = 4000` | observed sporadic (`p ≤ 40`) | Poisson |
|---|---|---|---|---|---|
| `3x+1` (dyadic) | 3.75 | 1.69 | 1.805 | 0 | naive `P(0) = 0.024`; refined 0.185 |
| `3x−1` (3-adic) | 1.15 | 1.15 | 1.210 | 1 (the 17-cycle) | `P(X ≤ 1) = P(X ≥ 1) = 0.68` |
| `5x+1` (dyadic) | 3.66 | 2.56 | 2.772 | 3 | refined `P(X ≥ 3) = 0.47` |
| `5x−1` (5-adic) | 1.12 | 0.39 | 0.653 | 0 | refined `P(0) = 0.68` |

* The refined expectation is concentrated at small `p`. For `3x+1`, the largest contributions are:
  * `(8,5)`: 0.54;
  * `(5,3)`: 0.40;
  * `(7,4)`: 0.11;
  * `(16,10)`: 0.08;
  * `(10,6)`: 0.07;
  * `(27,17)`: 0.07.
* The tail beyond `p = 100` is `0.005`; beyond the verified range it is negligible.
* This is the classical expectation that any cycle would be short, with the naive model's divergence removed.
* Transition clocks outside the sampled band carry naive weight 2.04 (dyadic, `q = 3`). They are counted as 0: the Monte Carlo share is 0 at the band edge and decreases away from the critical line.
* **DRIFT.** For `q = 5` the bulk `a ≈ p/2` lies on the expanding side (`1/2 > log_5 2`), where `G ≈ 5^a`, so the naive dyadic sum converges (4.03). The naive model's failure is thus tied to the sign of the drift `½ log(q/4)`: the typical word is contracting exactly when `q < 4`.

## 6. The dense divisor regime (D1, D3)

**Proposition SW (switching bound; PROVED).**
* Pair the positions as `(2k, 2k+1)`. In a *mixed* pair (`01` or `10`), swapping the two letters keeps the index `i_k` of its one and changes `c_w` by `δ_k = q^{a-1-i_k} 2^{2k}`.
* Group the words by skeleton; within a group the choices at the mixed pairs are independent. For every modulus `m` and every `h`:

  `|Σ_w e(h c_w/m)| ≤ Σ_w Π_{k mixed in w} |cos(π h δ_k/m)| =: C·B_m(h).`
* If `gcd(m, 2q) = 1` and `h ≢ 0 (mod m)`, every factor is at most `cos(π/m)`. So `max_{h≠0}|S_m(h)|/C ≤ E cos(π/m)^K`, where `K` is the number of mixed pairs (exact formula in D1), and this is `≈ exp(−c θ(1−θ) p/m²)`.
* **So carries are equidistributed modulo every fixed `m` coprime to `2q`, at an exponential rate.** For `m | G` this counts the cycles of `qx + G/m`.
* Modulo 2 and modulo `q` the carries are never uniform: `c_w ≡ w_0 (mod 2)` and `c_w ≡ 2^{s_{a-1}} ≢ 0 (mod q)` (checked on all words with `p ≤ 14`, `q ∈ {3,5}`).
* **D1** (exact residue counts modulo small primes `ℓ`, dividing the gate or not; DP over residues):
  * At `(60,30)` the deviations `max_r |n_r ℓ/C − 1|` are `10^{-7}`–`10^{-10}`, and at `(149,94)` modulo 7 they are `2·10^{-23}`. All are far inside the proven bound (which is weak: 0.98 for `ℓ = 61`), and the bound holds on every row.
  * The slowest case is `q = 5` modulo 7, where 2 has order 3: deviation `6·10^{-4}` at `p = 60`.

**D3 (Belaga–Mignotte 2006, table (20); primary text read).** The table lists the eleven `d ≤ 19999` for which `T_d(y) = y/2, (3y+d)/2` has more than 160 primitive cycles.
* A primitive `T_d`-cycle with clock `(p,a)` is a word with `(G/d) | c_w` and `gcd(c_w/(G/d), d) = 1`. So it lives on a clock with `d | G`.
* For every such clock with `p ≤ 250` (and `2^p > 3^a`), and along the main family up to `p ≤ 600`, we enumerated the odd least points in the perigee window `d 3^{a-1}/G ≤ y ≤ d/(2^{p/a} − 3)` exactly.
  * The window comes from the product identity `Π(3 + d/y_i) = 2^p`.
  * The map is the same as theirs, eq. (6).
* Per clock, the table compares the exact count with the equidistribution prediction `L(p,a)·(d/G)·Π_{ℓ|d}(1 − 1/ℓ)`:

| `d` | `ω(d)` | ours | prediction | main-clock `a/p` | found / predicted on the leading clocks |
|---|---|---|---|---|---|
| 14303 | 944 | 943 | 955.6 | 0.630 | `(27,17)` 843/880; `(54,34)` 76/64; `(81,51)` 20/9; `(108,68)` 3/2 |
| 17021 | 258 | 257 | 252.3 | 0.631 | `(65,41)` 254/247; `(130,82)` 3/5 |
| 6487 | 534 | 534 | 528.1 | 0.625 | `(16,10)` 456/458; `(32,20)` 59/52; `(48,30)` 12/12 |
| 13085 | 335 | 335 | 338.4 | 0.600 | `(15,9)` 277/266; `(30,18)` 34/44; `(45,27)` 17/15 |
| 10289 | 214 | 214 | 218.0 | 0.588 | `(17,10)` 166/163; `(34,20)` 26/31; `(51,30)` 10/11 |
| 7727 | 198 | 198 | 192.4 | 0.611 | `(18,11)` 162/161; `(36,22)` 28/22 |
| 9823 | 241 | 241 | 263.2 | 0.571 | `(14,8)` 183/179; `(28,16)` 28/40; `(56,32)` 5/9 |
| 18359 | 164 | 164 | 196.4 | 0.542 | `(24,13)` 112/114; `(48,26)` 23/34; `(96,52)` 5/10 |
| 15655 | 207 | 207 | 273.8 | 0.429 | `(14,6)` 160/163; `(28,12)` 34/49; `(42,18)` 11/23; `(56,24)` 1/13 |
| 14197 | 329 | 329 | 490.5 | 0.500 | `(14,7)` 245/245; `(28,14)` 43/77; `(42,21)` 21/42; `(70,35)` 4/19 |
| 7463 | 162 | 162 | 237.3 | 0.462 | `(13,6)` 126/124; `(26,12)` 23/39; `(39,18)` 6/20; `(65,30)` 1/9 |

* **Near the critical line** (`a/p ≥ 0.59`), equidistribution modulo `G/d` predicts each clock within Poisson noise, and the totals within 1–4%.
* **On far clocks** the main clock is exact: it has cofactor `G/d = 1`, so every necklace is a `T_d`-cycle. The multiples, however, are depleted by factors up to 10.
* This is the dense-regime face of the archimedean non-uniformity of `§4`. It confirms, independently and from published data, that equidistribution holds exactly near the critical line.
* **Resolved (2026-09-26, sporadic lane).** `ω(14303) = 944` and `ω(17021) = 258` exceeded our counts by one each. The missing cycles have clocks `(2155,1092)` (least element 101) and `(2140,1088)` (least element 5), at `a/p ≈ 0.507`. The claim here that the other lattice clocks have `a/p < 0.41` holds only up to `p ≈ 600`; the same lattice families return to `a/p ≈ 1/2` near `p ≈ 2000`. See [procgen_sporadic_20260926_free_and_sporadic_cycles](procgen_sporadic_20260926_free_and_sporadic_cycles.md).
  * No window was skipped and no orbit overflowed.
  * The main families were extended to `p ≈ 600` with no further cycle; the other lattice clocks with `d | G` beyond `p = 250` have `a/p < 0.41` and negligible predictions.
  * We cannot tell whether a cycle sits on a clock we did not scan or the table entry counts something else. Their exhaustiveness argument is in the 62-page report [4] (UNVERIFIED, not read).

## 7. The proof angle (C3, D2)

**The counting question is classical, and much easier than Fourier makes it.**

| bound on `N(p,a)` | source | size at `a/p ≈ log_3 2` |
|---|---|---|
| trivial | `N ≤ C` | `2^{0.950p}` |
| second moment (Parseval + Cauchy–Schwarz) | `|N − C/M| ≤ √((1−1/M)(Coll − C²/M))` | at least `≈ √C = 2^{0.475p}`, since `Coll ≥ C` |
| size range (`w ↦ x_w` injective) | `N ≤ (c_max − c_min)/M + 1` (Böhm–Sontacchi type) | `C^{0.389+o(1)}·amp` |
| perigee window (Belaga 2003) | `N ≤ p (1/|2^{p/a} − 3|/2 + 1)` | `poly(p)` by Rhin/Baker |
| target | `N = 0` | — |

* C3 prints these on convergent and intermediate clocks up to `(50508, 31867)`. There the perigee bound is 45 bits, against `√C` of 23986 bits and a size range of 18658 bits.
* **The second-moment method cannot compete.**
  * Moving `√C` singleton residues onto 0 raises `Coll` only by about `C`. So no bound that depends only on `(C, M, Coll)` certifies `N < √C` in the sparse regime (`C < M`).
  * For `p ≳ 100` the sparse regime is the only regime at the gate modulus: `C/|G|` exceeds 1 on few clocks (`(84,53)`: 2.3; `(485,306)`: `2^{-19}`).
  * Higher moments give `N ≤ (Σ n_r^k)^{1/k}`, which is never below 1.
  * No Fourier or moment method yields `N = 0`. Exclusion needs the exact arithmetic of each individual gate.

**Weyl differencing on the DP (D2).**
* The switching certificate at the gate modulus is rigorous: `|S(h)|/C ≤ B_M(h)` was checked on every frequency tested (7 clocks; 200 random `h` and up to 40 structured `h` each).
* It separates structured from random frequencies by only about 10×. For example at `(84,53)`: structured `|S|/C = 1.6·10^{-6}` with `B = 1.8·10^{-2}`, while random `h` give median `B = 2.8·10^{-4}`, against a true median `|S|/C` of `3·10^{-13}`.
* `B` decays like `e^{-0.1p}`, whereas `C^{-1/2} ≈ e^{-0.33p}`, and exclusion needs far more.
* A Tao-type refinement (pairs of terms, a two-dimensional renewal process, black triangles around `h 2^s 3^t ≡ small`) could plausibly prove "no large coefficient off the `⟨2,3⟩`-orbit of small integers", with polynomial or exponential savings. **Even a complete proof of that leaves the cycle half untouched,** because of the square-root barrier.

**What is provable here, and what it buys:**
1. equidistribution modulo fixed `m` coprime to `2q`, i.e. counts of `qx + G/m` cycles (Proposition SW; proved);
2. the archimedean limit law and its non-uniformity (Proposition A; proved off `[1/2, 2/3]`);
3. the divergence of the naive model (Proposition N; proved);
4. the perigee classification (Proposition P; elementary, known in substance).

None of these bears on whether a large clock carries an integral cycle.

## 8. Literature and priority

**Read in primary text.**
* **Tao, arXiv:1909.03562, `§1.4`.**
  * Syracuse random variables `Syrac(Z/3^n) ≡ Σ_m 3^{m-1} 2^{-a_{[1,m]}}` with Geom(2) gaps.
  * Prop. 1.14 (fine-scale mixing) and Prop. 1.17: `|E e(ξ Syrac/3^n)| ≪_A n^{-A}` for `3 ∤ ξ`.
  * Remark 1.15: the entropy heuristic, with `4^n` tuples mapped into `3^n` residues.
  * Our `S(h)` is the same object with modulus `2^p − 3^a` and the fixed-weight measure. Tao's regime is *supercritical* in entropy, so fine-scale mixing is the right statement there. The cycle gate is *subcritical* (`2^{0.950p}` words into about `2^p` residues), so the right statement is "random-like sparse set", and hitting 0 is a sparse hitting problem.
  * Tao's coarse 3-adic irregularity ("it always avoids the multiples of 3") is our dominant peak `h ≡ ±3^{-k}`.
* **Belaga–Mignotte 2006**, DMTCS proc. AG, 249–260 (open-access text downloaded by the atlas lane).
  * Map (6).
  * Lagarias's conjectures on primitive `3x+d` cycles.
  * 42765 primitive cycles for `d ≤ 19999`.
  * Table (20), reproduced in `§6`.
  * The perigee bound credited to Belaga 2003.

**Read as summaries only**, in Lagarias's annotated bibliographies (arXiv math/0309224 and math/0608208). Primary texts UNVERIFIED:
* **Böhm–Sontacchi 1978:** the rational-cycle formula; every integer in a cycle of length `n` has `|x| < 3^n`.
* **Belaga–Mignotte 1999** (Exp. Math. 7): at most `d k^c` periodic orbits of `3x+d` with at most `k` odd integers, via Baker–Wüstholz.
* **Belaga 2003** (Acta Arith. 106):
  * perigee `≤ d/(2^{l/k} − 3)`;
  * at most `d k^{c_0}` cycles of oddlength `k`, with `c_0 = 32`;
  * largest element `< d k^{c_0} (3/2)^k`.
* **Lagarias 1990** (Acta Arith. 56): rational cycles, i.e. `T_k`-cycles averaged over `k`.
  * Infinitely many `k` have at least `k^{1−ε}` cycles of period at most `log k`.
  * Estimates for the counting function.
* **Eliahou 1993; Halbeisen–Hungerbühler 1997:** bounds on the least element of rational cycles.
* **Belaga–Mignotte 2000** (57-page report on the distribution of `3x+d` cycles): its HAL copy is behind a bot check, which was not bypassed. UNVERIFIED.

**Priority.**
* The counting statement suggested in the task, `N = O(C^{1−δ})` for large `|G|`, is **not new**: the size range gives `δ ≈ 0.61`, and Belaga–Mignotte 1999 and Belaga 2003 give polynomial bounds.
* Not found in the sources read:
  * the divergence of the naive residue model and its sign-aware repair;
  * the archimedean limit law of the gate residues;
  * the `⟨2,3⟩`-orbit structure of the gate spectrum;
  * the per-clock comparison with table (20).
* They may be implicit in the unread Belaga–Mignotte reports, so **no priority is claimed**.
* Seen, abstract only, not used: arXiv:2603.11066 (Chang 2026, an LLM-assisted structural study) and arXiv:1909.00213 (Tremblay 2019).

## 9. Controls: what is sheet-aware

* **SHEET.**
  * Every statistic of `{c_w mod G}` is literally the same for `3n+1` and `3n−1` on a given clock (transport `x ↦ −x`).
  * The census finds the `3n−1` cycles on the 3-adic clocks just as it finds the `3n+1` cycles on the dyadic ones.
  * The residue model is side-blind. The side-aware inputs are exactly two:
    1. the sign law: positive `3n+1` cycles need `G > 0`;
    2. the least-point condition `m ≥ 1` of Proposition P, i.e. the positive cone of the mod-192 note `§2.4`.
  * Without them the model mis-weights the sheets by `ln P`. With them both sheets are Poisson-consistent.
* **DRIFT.**
  * `5x+1` runs through the same code and passes the census.
  * The naive divergence is special to negative drift (`q < 4`): which side of the critical line the typical word lies on is exactly the sign of the drift.
  * The refined model predicts 2.8 sporadic `5x+1` cycles; 3 are observed.
* **DEFECT.** Density-type statements about gate residues tolerate finitely many exceptional clocks, and a cycle is one exceptional clock. So no residue statistic can exclude it.

## 10. Assessment

* **Does equidistribution of `c_w` modulo the gates look true?**
  * *Macroscopically, no.* There is a non-uniform archimedean law plus coarse 2-adic and 3-adic digit biases, all on the `⟨2,3⟩`-orbit of small frequencies. This is proved off `[1/2, 2/3]`.
  * *After removing that part, yes.* Coefficients are random-sized, coincidences approach uniform at the critical line, and published `3x+d` counts are predicted within noise near the critical line.
* **Is it strong enough to exclude cycles heuristically?**
  * The naive version is falsified by its `ln P` divergence.
  * The side-aware version gives finite expectations concentrated at `p ≤ 30`, consistent with all four sheet × drift controls. The expected number of undiscovered sporadic positive `3x+1` cycles beyond `p = 100` is `0.005`.
* **Is it provable in some range?**
  * Only in dense regimes: fixed or slowly growing moduli, divisors of gates, and Tao-type decay off the structured set.
  * At the gate itself, for large `p`, the square-root barrier and the sparse regime make Fourier counting strictly weaker than the classical perigee, Baker and verification route.
* **Rank.** LOW as a proof angle; MEDIUM as a heuristic instrument.
* **What to keep:** the side-aware random model, and the rule never to quote the naive `C/|G|` sum.

## 11. Reproduction

* **Command:** `python3 04-computation/experiments/procgen_gates_20260925_run.py`. It runs four parts, one process at a time, about 18 minutes in all.
* **Output:** [`procgen_gates_20260925.out`](procgen_gates_20260925.out). Every check raises on failure. Each part's peak memory is in the header: A 375 MB, B 433 MB, C 276 MB, D 252 MB.
* **A** (`expsum`): DP certification, exact counting identity, two-method census (SHEET, DRIFT), rotation law.
* **B** (`stats`):
  * full residue statistics for `p ≤ 30` (`q = 3`) and `p ≤ 24` (`q = 5`);
  * spectra;
  * the perigee classification;
  * writes `scratch/procgen_gates/partB_elig.tsv`, which C reads.
* **C** (`heuristic`):
  * naive and refined sums to `P = 4000`;
  * Poisson comparison;
  * the bound hierarchy;
  * the archimedean limit law, with 40-digit certification of the tiny coefficients.
* **D** (`switching`):
  * Proposition SW;
  * residue counts modulo small primes;
  * the certificate at the gate modulus;
  * the Belaga–Mignotte reproduction.
* **Memory (macOS).**
  * numpy's `rfft` needs about 150 bytes per point for lengths with large prime factors (842 MB at `M = 5·10^6`). So FFTs are capped at `M ≤ 2^20`; larger near-critical spectra use the blocked DP.
  * The scripts re-exec with `MallocLargeCache=0`, which returns freed large blocks to the system and roughly halves resident memory.
  * **Disclosure:** an exploratory run before these fixes exceeded the 700 MB budget, peaking at 1.6 GB during an FFT of length about `2^23`. The recorded run stays below 450 MB.
* **Literature texts.** The Lagarias bibliographies and Tao's paper are in `scratch/procgen_gates/` and are not for commit. Belaga–Mignotte 2006 was read from `scratch/procgen_atlas/lit_dioph/bm2006.txt`.
