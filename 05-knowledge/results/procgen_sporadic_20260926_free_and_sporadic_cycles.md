# Free and sporadic cycles of `qx + d`, and the Belaga–Mignotte off-by-one

**Status.**
* **T2 is RESOLVED (FINITE-EXACT): the off-by-one was an omission in the gates lane's scan.** There is no convention difference, and Belaga–Mignotte made no error on these two entries. Each of `d = 14303` and `d = 17021` has one more primitive cycle, on a clock the gates lane never scanned. Both extra cycles are long and lie far from the critical line:
  * `d = 14303`: least element **101**, period `p = 2155` with `a = 1092` odd steps. Its clock is `(2155, 1092) = 75·(27,17) + (130,−183)`.
  * `d = 17021`: least element **5**, period `p = 2140` with `a = 1088` odd steps. Its clock is `(2140, 1088) = 31·(65,41) + (125,−183)`.

  With these two cycles the counts are 944 and 258, as published. A least-element search, combined with the lattice of admissible clocks (§2.3), shows that these are *all* the primitive cycles with period `p < 203 720` (`d = 14303`) and `p < 1 631 245` (`d = 17021`).
* **Whole table (FINITE-EXACT, T2).** One uniform search reproduces all eleven entries of table (20), with no special cases.
  * The search covers every `d ≤ 19999` prime to 6 and every least element `≤ 1200 d`.
  * It also reproduces the published number of systems with `ω = 1` (1481) and `ω = 2` (1507).
  * Every `T_d` has a primitive cycle.
* **UNRESOLVED residual.** We find 42757 primitive cycles in all; Belaga–Mignotte report 42765. We find `ω = 3` for 1004 systems; they report 1005.
  * Extending the search adds nothing. No further cycle exists with period `< 4000` for any `d ≤ 19999`. No further cycle exists with least element `≤ 2·10^7` for any `d < 16667`.
  * A counting argument (§2.5) shows the residual is two-sided: either their list lacks at least one cycle we exhibit, or their summary has a misprint.
  * Settling it needs their Table 1 (HAL `hal-00129727`), which sits behind a bot check. The check was not bypassed.
* **PROVED (T1; hand proofs in §1).** Take `T_{q,d}(y) = y/2` or `(qy+d)/2`, with `q ≥ 3` odd, `d` odd and `gcd(q,d) = 1`.
  * **Proposition 1 (shift criterion).** A mixed shape `(p,a)`, `1 ≤ a ≤ p−1`, is free iff `(2^p − q^a) | d`. The single-word shapes give `0`, and give `−d/(q−2)` when that is an integer.
  * **Proposition 2 (Gersonides for every odd `q`, elementary).** `|2^p − q^a| = 1` with `p, a ≥ 1` iff `a = 1` and `q = 2^p ± 1`, or `(q,p,a) = (3,3,2)`.
  * **Proposition 3 (primitive decomposition).**
    * The cycles of `T_{q,d}` are `⊔_{e|d} e·Prim(T_{q,d/e})`.
    * A primitive `T_{q,f}`-cycle forces `f | D`.
    * On a clock, put `g = gcd(D,d)`. The integral words are exactly those with `(D/g) | c_w`. The clock is free if `g = |D|`, partially free if `1 < g < |D|`, and sporadic if `g = 1 < |D|`.
  * **Corollaries.**
    * The free cycles of `3x+1` are `{0}, {−1}, {1,2}, {−5,−7,−10}`. This confirms wave 13, and no appeal to Mihăilescu is needed.
    * The free cycles of `5x+1` are exactly `{0}` and `{−1,−2}` (wave 13 confirmed).
    * The cycles of `3x−1` are the negatives of those of `3x+1` (wave 13 confirmed).
* **CITED.**
  * Finiteness of the free shapes of each `d`: Pillai 1931 (ineffective), as reported in Ellison 1971 §I; and Ellison 1971 Theorem 1 (effective, via Baker).
  * For `q = 3`, Ellison's explicit bound `|2^x − 3^y| > 2^x e^{−x/10}` for `x ∉ {1,…,11,13,14,16,19,27}`.
  * Belaga 2003's perigee bound (via Lagarias's summary; re-proved in §2.3).
* **FINITE-EXACT (T1, T4).**
  * All `(p,a)` with `|2^p − 3^a| ≤ 20000`.
  * Free-cycle censuses for ten values of `d`.
  * `5x+1` on `Z`, all periods `p ≤ 60`, both signs: exactly `0, {−1,−2}, {1,…}, {13,…}, {17,…}`.
  * `3x+1` on `Z`, all periods `p ≤ 90`: its five known cycles.
  * `3x±1`: least elements `≤ 10^9`.
* **EMPIRICAL (T2, T3).**
  * **Long cycles are common.** 2204 primitive cycles have `p > 600`. Their `a/p` lies in `[0.444, 0.561]`, with mean `0.502`. For 683 values of `d` the *only* primitive cycle is a long one.
  * **Least element.** The least element scales like `m ≈ 5d/p`: over the 18697 cycles with `p ≥ 100`, the quartiles of `mp/d` are 2.03, 4.96 and 10.30, and `p/d ≤ 1.18` throughout. The PROVED ingredient is `m > d/2^{K+1}`, where `K` is the longest run of even steps.
  * **Best approximations.** "Cycles sit at best approximations" holds for `3x±1` and `5x±1`, where every cycle lies on the Stern–Brocot path. It does not hold for `3x+d` in general: only 22.5% of the 42757 primitive cycles lie on the path.
* **ANALOGY (T3).** The Kuratowski-style reading, briefly, in §3.3.
* **Corrections to the gates note** (`procgen_gates_20260925_gate_equidistribution.md`).
  * Its §6 "Unresolved" item is now resolved.
  * Its sentence that the other lattice clocks beyond `p = 250` "have `a/p < 0.41` and negligible predictions" is true only up to `p ≈ 600`. Beyond that, the same lattice families climb back to `a/p ≈ 1/2`. Summed over the unscanned clocks with `p ≤ 2200`, the equidistribution prediction is **1.63** (`d = 14303`) and **0.84** (`d = 17021`). So its own model predicted about one extra cycle each.
  * Its status line "9 of 11 exactly, two off by one" should read 11 of 11.

No HYP or THM file was created, and no novelty is claimed. The Belaga–Mignotte data and the classical Diophantine facts are theirs; the free/partially free/sporadic bookkeeping is elementary.

Session `collatz-procgen-20260922`, lane "sporadic" (mac-mini), 2026-09-26.
* Scripts (`04-computation/experiments/`): `procgen_sporadic_20260926_{run,lib}.py`, `procgen_sporadic_20260926_{traj,sweep}.c`.
* Output: [`procgen_sporadic_20260926.out`](procgen_sporadic_20260926.out). 70 checks; every printed claim is a raising `check`.
* Inputs read:
  * [wave-13 findings](procgen_wave13_20260926_orchestrator_findings.md) §1;
  * [gates note](procgen_gates_20260925_gate_equidistribution.md) §§5–6 and its `procgen_gates_20260925_switching.py` (`td_cycles_on_clock`, `d3_belaga_mignotte`);
  * Belaga–Mignotte 2006, primary text (`scratch/procgen_atlas/lit_dioph/bm2006.txt`);
  * Ellison 1971, primary text (`…/ellison1971.txt`);
  * the atlas lane's `notes_dioph.md`;
  * Lagarias's annotated bibliographies (`scratch/procgen_gates/lagarias_bib{1,2}.txt`).

---

## 1. T1: free, partially free and sporadic cycles of `qx + d`

### 1.1 Setting

**The map.** Fix `q ≥ 3` odd, and `d` odd with `gcd(q,d) = 1`; `d` may be negative. Let `T = T_{q,d}` act on `Z` (and on `Z_2`) by `T(y) = y/2` for even `y` and `(qy+d)/2` for odd `y`.

**Words.** A word `w = w_0…w_{p−1}` has `a` ones, at positions `s_0 < … < s_{a−1}`.
* Its *carry* is `c_w = Σ_i q^{a−1−i} 2^{s_i}`.
* Its *gap* is `D = D(p,a) = 2^p − q^a`, which is odd and non-zero.
* The pair `(p,a)` is its *shape*, or *clock*. The shape is *mixed* if `1 ≤ a ≤ p−1`.

**Lemma 1.1 (affine form).** If `y ∈ Z_2` has first `p` parities `w`, then `2^p T^p(y) = q^a y + d c_w`.

*Proof.* Put `z_t = 2^t T^t(y)`. Then `z_{t+1} = z_t` at an even letter and `z_{t+1} = q z_t + d 2^t` at an odd letter. The odd step at time `s_i` therefore contributes `d 2^{s_i}`, multiplied by `q` once for each of the `a−1−i` later odd steps. ∎

**Lemma 1.2 (periodic point).** The unique `y ∈ Z_2` with itinerary `w^∞` is `y_w = d c_w/D`. This is Böhm–Sontacchi 1978 and Lagarias 1990, known here from Lagarias's summaries; the proof below is self-contained.

*Proof.*
* *Parities determine residues.* The map `Z/2^p → {0,1}^p` sending `y` to its first `p` parities is a bijection. If `y' = y + 2^k u` with `u` odd and `k < p`, and `y, y'` agree in their first `k` parities, then `T^k y' = T^k y + q^{a_k} u`, which has the other parity.
* *The fixed point.* So the cylinder `C_w` is a single residue class mod `2^p`. On it, `T^p` is the map `y ↦ (q^a y + d c_w)/2^p`, a bijection `C_w → Z_2`. Its inverse is a `2^{−p}`-contraction of `Z_2`, so it has a unique fixed point, and that point solves `D y = d c_w`. ∎

**Cycles and freeness.**
* An integer cycle of exact period `p` is the same thing as a necklace of *primitive* words of length `p` with `y_w ∈ Z`. By Lemma 1.2, distinct words give distinct points.
* A shape is **free** (for `T_{q,d}`) if every word of that shape is integral.
* A cycle is **free** if the shape of its primitive word is free.

### 1.2 The shift criterion

**Proposition 1.**
* **(a)** A mixed shape `(p,a)` is free iff `D | d`.
* **(b)** `(p,0)` has the single word `0^p`, which gives `0`. `(p,p)` has the single word `1^p`, which gives `−d/(q−2)`; this is integral iff `D(1,1) = 2−q` divides `d`. So a shape that carries a primitive word is free iff `D | d`.
* **(c)** Sub-shapes of free shapes are free, because `D(p/k, a/k) | D(p,a)`. A free mixed shape carries exactly `L(p,a) = (1/p) Σ_{k | gcd(p,a)} μ(k) C(p/k, a/k)` free cycles of exact period `p`, one per Lyndon word.

*Proof.*
* **(a), "if":** `D | d` implies `D | d c_w` for every `w`.
* **(a), "only if":** For `1 ≤ a ≤ p−1` the words `u = 1^a 0^{p−a}` and `u' = 1^{a−1} 0 1 0^{p−a−1}` both have shape `(p,a)`. They differ only in the last one, which moves from position `a−1` to position `a`; its coefficient is `q^0`. So `c_{u'} − c_u = 2^{a−1}`. If both words are integral, then `D | d 2^{a−1}`, and since `D` is odd, `D | d`.
* **(b):** `c_{1^p} = (q^p − 2^p)/(q−2)`.
* **(c):** `2^{p'k} − q^{a'k}` is divisible by `2^{p'} − q^{a'}`. ∎

Checked (T1.2) exhaustively on every mixed shape with `p ≤ 12`, for `q = 3, 5, 7` and 17 values of `d`: 2706 cases.

**Finiteness.** Free mixed shapes satisfy `|2^p − q^a| ≤ |d|`, which has finitely many solutions for fixed `q`.
* **Pillai 1931.** J. Indian Math. Soc. 19: `|a m^x − b n^y| > m^{(1−δ)x}` for `x > x_0`. The proof uses Thue–Siegel, so it is ineffective. (Statement as reported in Ellison §I.)
* **Ellison 1971.** Sém. Théorie Nombres Bordeaux, exp. 12. Theorem 1 and Corollary 1 make the bound effective via Baker.
* **For `q = 3`, explicitly (Ellison 1971).** The theorem on p. 12-05, numbered 3 in the scan, states `|2^x − 3^y| > 2^x e^{−x/10}` for all positive `x, y` with `x ∉ S = {1,…,11,13,14,16,19,27}`. Both statements were read in the primary text.
* **Stroeker–Tijdeman 1982** (Pillai's conjecture for bases 2 and 3) was not read (UNVERIFIED).

Consequence (T1.5, FINITE-EXACT).
* Since `2^{17}e^{−1.7} = 23946 > 20000`, every `(p,a)` with `|2^p − 3^a| ≤ 20000` has `p ≤ 16` or `p ∈ {19, 27}`.
* The full list (129 pairs with `p, a ≥ 1`; the mixed ones are printed in the `.out`) agrees with brute force up to `p = 400`.
* It contains the free main clocks of six of the eleven table-(20) entries: `D(13,6) = 7463`, `D(14,6) = 15655`, `D(14,7) = 14197`, `D(14,8) = 9823`, `D(15,9) = 13085`, `D(16,10) = 6487`.

### 1.3 Gersonides for every odd multiplier

**Proposition 2.** Let `q ≥ 3` be odd and `p, a ≥ 1`. Then `|2^p − q^a| = 1` iff `a = 1` and `q = 2^p ± 1`, or `(q,p,a) = (3,3,2)`.

*Proof.*
* **Case `q^a − 2^p = 1`.**
  * `p = 1` forces `q^a = 3`.
  * Let `p ≥ 2` and `a` odd. Then `q^a − 1 = (q−1)(q^{a−1}+…+1)`. The second factor is a sum of `a` odd terms, so it is odd; it divides `2^p`, so it is 1, and `a = 1`.
  * Let `a = 2b` be even. Then `q^b − 1` and `q^b + 1` are powers of 2 that differ by 2, so `q^b = 3`.
* **Case `2^p − q^a = 1`.**
  * `p = 1` is impossible.
  * If `a` is even, then `q^a + 1 ≡ 2 (mod 8)`, which is impossible.
  * If `a` is odd, then `q^a + 1 = (q+1)(q^{a−1} − … + 1)` with an odd second factor, which must be 1; so `a = 1`. ∎

Checked for all odd `q < 1000` and `p ≤ 64` (T1.4).

**Corollary (free cycles of `T_{q,±1}`).** Combining Propositions 1 and 2, the free mixed shapes are:
* `(2,1)` and `(3,2)` for `q = 3`;
* `(p,1)` when `q = 2^p ± 1` and `p ≥ 2`;
* none otherwise.

| map | free cycles |
|---|---|
| `3x+1` | `{0}, {−1}, {1,2}, {−5,−7,−10}` (from `2−1, 3−2, 4−3, 9−8`) |
| `3x−1` | `{0}, {1}, {−1,−2}, {5,7,10}` |
| `5x+1` | `{0}, {−1,−2}` (the shape `(1,1)` gives `−1/3`) |
| `5x−1` | `{0}, {1,2}` |

So `3` is the only multiplier with a free cycle of two odd steps. The wave-13 claim "for odd `q ≥ 5` every solution has `a ≤ 1`" is PROVED without Mihăilescu.

### 1.4 Primitive decomposition

**Proposition 3.**
* **(a)** `gcd(T(y), d) = gcd(y, d)`.
* **(b)** Let `C` be an integer cycle and `e = gcd(y,d)` for `y ∈ C` (the same for every `y ∈ C`). Then `C/e` is a *primitive* cycle of `T_{q,d/e}`, with the same necklace. So the integer cycles of `T_{q,d}` are `⊔_{e|d} e·Prim(T_{q,d/e})`, with `e` running over the positive divisors of `d`.
* **(c)** A primitive `T_{q,f}`-cycle with word `w` forces `f | D`, `(D/f) | c_w` and `gcd(c_w f/D, f) = 1`. Conversely, these three conditions give a primitive cycle.
* **(d)** On a clock, let `g = gcd(D,d)`. The integral words are exactly those with `(D/g) | c_w`. They are partitioned by the *primitive level* `|f| = |d|/gcd(y_w,d)`, which divides `g`. The `f`-part consists of the multiples `y_w(d) = (d/f)·y_w(f)` of primitive `T_{q,f}` points.

*Proof.*
* **(a)** For odd `y`, `gcd((qy+d)/2, d) = gcd(qy+d, d) = gcd(qy, d) = gcd(y, d)`.
* **(b)** `T_{q,d}(ez) = e·T_{q,d/e}(z)` because `e` is odd. Also `gcd(y/e, d/e) = 1`, since at most one of `v_ℓ(y) − v_ℓ(e)` and `v_ℓ(d) − v_ℓ(e)` is positive.
* **(c)** From `D y_w = f c_w` and `gcd(y_w, f) = 1` we get `f | D`, and then `y_w = c_w/(D/f)`.
* **(d)** `y_w ∈ Z` iff `(D/g) | (d/g) c_w` iff `(D/g) | c_w`; then apply (b) and (c). ∎

Checked:
* (d) on 455 (composite `d`, clock) cases (T1.7).
* (b) on the whole sweep: `#cycles(T_d, least ≤ 1200d) = Σ_{e|d} #Prim(T_{d/e}, least ≤ 1200d/e)` for all 6667 values of `d` (T2.6).

**The three kinds.** Relative to `d`, a clock is:
* **free** if `g = |D|`;
* **partially free** if `1 < g < |D|` ("`d` pays the factor `g` of the gap");
* **sporadic** if `g = 1 < |D|`. A sporadic cycle of `T_{q,d}` is `d·C_1` for an integer cycle `C_1` of `T_{q,±1}`.

A Belaga–Mignotte (Lagarias) primitive cycle, one with `gcd(y,d) = 1`, lives on the lattice `L_d = {(p,a) : 2^p ≡ 3^a mod d}`, since `d | D` there.
* It is **free** iff `D = d`. Examples are the main clocks of `d = 6487, 13085, 9823, 15655, 14197, 7463`.
* It is **partially free** otherwise. Examples: `d = 14303` on `(27,17)` with `D/d = 355`; `d = 17021` on `(65,41)` with `D/d = 19·29·44835377399`.
* The scaled `3x+1` cycles appear in every `T_{3,d}`: `d·{1,2}`, `{−d}` and `d·{−5,−7,−10}` are free, and `d·{−17,…}` is sporadic unless `139 | d`.

**Example: `T_{3,5}` on `Z`** (census of least elements `|·| ≤ 10^6`; T1.8).

| least | clock | `D` | `g` | kind | level `f` |
|---|---|---|---|---|---|
| 0 | (1,0) | 1 | 1 | free | — |
| 5 (`= 5·{1,2}`) | (2,1) | 1 | 1 | free | 1 |
| −5 | (1,1) | −1 | 1 | free | 1 |
| −25 (`= 5·{−5,…}`) | (3,2) | −1 | 1 | free | 1 |
| 1 (`{1,4,2}`) | (3,1) | 5 | 5 | free | 5 |
| 19, 23 | (5,3) | 5 | 5 | free, `L(5,3) = 2` | 5 |
| 187, 347 | (27,17) | 5·71·14303 | 5 | partially free, `D/g = 1015513` | 5 |
| −85 (`= 5·{−17,…}`) | (11,7) | −139 | 1 | sporadic | 1 |

Two further examples:
* **`d = 5`.** Pillai's two-solution identity `8 − 3 = 32 − 27 = 5` makes both `(3,1)` and `(5,3)` free.
* **`d = 139`.** The whole shape `(11,7)` becomes free: all 30 necklaces are cycles `y = −c_w` of `T_{3,139}` (T1.6). 29 of them are primitive. The thirtieth, the only necklace with `139 | c_w`, is `139·{−17,…}`.

## 2. T2: the off-by-one

### 2.1 What was compared

**Belaga–Mignotte 2006** (DMTCS proc. AG, 249–260), read in the primary text:
* The map is (6): `T_d(n) = n/2` or `(3n+d)/2` on `N`, with `d` prime to 6.
* Definition 3: a trajectory is primitive iff `gcd(n,d) = 1`, and non-primitive cycles are multiples of `T_{d/q}`-cycles (eq. (19)).
* §6 (4) reports:
  * "all 6667 systems" with `1 ≤ d ≤ 19999`;
  * 42765 primitive cycles;
  * `ω = 1` for 1481 systems, `ω = 2` for 1507, `ω = 3` for 1005;
  * table (20): the eleven `d` with `ω(d) > 160`.

**The gates lane** (`td_cycles_on_clock` and `d3_belaga_mignotte` in `procgen_gates_20260925_switching.py`) enumerated the odd least points in the perigee window of:
* every lattice clock with `p ≤ 250`;
* the main family `k·(p_0,a_0)` up to `p ≤ 600`.

This scanned 22 clocks for `d = 14303` and 9 for `d = 17021`, and found 943 and 257 cycles.

### 2.2 The two missing cycles

Exact verification (T2.1): each cycle is primitive, `mD = dc_w`, `m` lies in the perigee window, and the clock has integer coordinates in the Gauss-reduced basis of `L_d`.

| `d` | least `m` | `(p,a)` | `a/p` | max element | `D/d` | window | reduced-basis coordinates | prediction `L(p,a)(d/D)(1−1/d)` |
|---|---|---|---|---|---|---|---|---|
| 14303 | **101** | (2155, 1092) | 0.5067 | 3 863 248 | 2142 bits | `[1, 15428]` | `75·(27,17) + 1·(130,−183)` (index 7151) | 0.094 |
| 17021 | **5** | (2140, 1088) | 0.5084 | 9 122 452 | 2126 bits | `[1, 18718]` | `31·(65,41) + 1·(125,−183)` (index 17020) | 0.101 |

**The lattices.** Both `d` are prime.
* `L_{14303}` has reduced basis `{(27,17), (130,−183)}`; `ord_{14303} 2 = 7151`.
* `L_{17021}` has reduced basis `{(65,41), (125,−183)}`; 2 is a primitive root mod 17021.
* The only lattice clocks with `p ≤ 250` are multiples of the main clock: 9 of them for `d = 14303` and 3 for `d = 17021`.

**Why the scan missed them.**
* **The next family starts far from the line.** The family `i·(main) + (second vector)` begins at `(427,4), (454,21), (481,38), …` for `d = 14303` and at `(450,22), (515,63), (580,104), …` for `d = 17021`. Up to `p = 600` these have `a/p ≤ 0.18`, which is the gates note's "`a/p < 0.41`".
* **It returns to `a/p = 1/2` near `p ≈ 2000`.** As `i` grows, `a/p` increases towards the main ratio and crosses `1/2` at `p ≈ 2000`.
  * At `a/p = 1/2` the window `d/(2^{p/a} − 3)` equals `d`.
  * The per-clock prediction peaks at the maximum-entropy value `≈ √(2/π)·d/p^{3/2}`: `0.125` at `(2020,1007)` for 14303, and `0.151` at `(2010,1006)` for 17021 (T2.3 tables).
* **So one cycle each was about expected.** Summed over the unscanned lattice clocks with `p ≤ 2200`, the prediction is **1.625** (`d = 14303`) and **0.842** (`d = 17021`). Exactly one cycle was found for each.
* **Independent check (T2.4).** The perigee-window method, run on every clock that carries a cycle, reproduces the least-element search clock by clock:
  * `d = 14303`: `(27,17)`: 843, `(54,34)`: 76, `(81,51)`: 20, `(108,68)`: 3, `(135,85)`: 1, `(2155,1092)`: 1.
  * `d = 17021`: `(65,41)`: 254, `(130,82)`: 3, `(2140,1088)`: 1.

**The requested divisibility facts** (T2.2).
* `2^27 − 3^17 = 5077565 = 355·14303`, with `355 = 5·71`.
* `2^65 − 3^41 = 420491770248316829 = 17021·19·29·44835377399`.
* The numbers 14303, 17021 and 44835377399 are prime.

So both values of `d` divide a near-critical gap: `17/27` is a semiconvergent of `log_3 2` and `41/65` a convergent. That explains the *size* of `ω(d)`. It has nothing to do with the missing cycles, which live far from the critical line.

### 2.3 Completeness: the lattice cone

**Lemma 2.1 (perigee window; Belaga 2003, re-proved).** Let `C` be a positive cycle of `T_{q,d}`, `d > 0`, with odd elements `y_1…y_a`, period `p` and least element `m`. Then `Π(q + d/y_i) = 2^p`, and hence `d q^{a−1}/D ≤ m ≤ d/(2^{p/a} − q)`.

*Proof.* Multiply `T(x)/x` around the cycle; the factor is `1/2` at even `x` and `(q + d/x)/2` at odd `x`. The least element is odd, since otherwise `m/2` would be in the cycle. For the lower bound, `m = y_w` with `w_0 = 1`, so `c_w ≥ q^{a−1}`. ∎

**Lemma 2.2 (least-element search).** For odd `y`, follow the odd values `x ↦ (3x+d)/2^{v_2(3x+d)}`. Then `y` is the least element of a cycle iff this orbit returns to `y` before any odd value below `y` appears.

*Proof.* The even values between consecutive odd values are `x'2^j > x'`. ∎

**The C programs** implement Lemma 2.2 with 128-bit arithmetic.
* After 3000 odd steps they switch to Brent's cycle detection. This handles long cycles, and also orbits that enter a cycle with a larger least element.
* In every recorded run there was no overflow and no unresolved start.
* No orbit from `y ≤ X` entered a cycle whose least element exceeds `X`. So every cycle whose basin meets `[1, X]` has been found. The least element of a basin is odd, and its orbit never goes below it.
* Consequently, a Belaga–Mignotte-style search "iterate every `n ≤ N`" with `N ≤ X` cannot find more cycles.

**Lemma 2.3 (cone).** Let `C` be a primitive `T_d`-cycle with least element `m ≥ X+1`. Then its clock satisfies:
* `(p,a) ∈ L_d`;
* `D > 0`;
* `3^a < 2^p ≤ (3 + d/(X+1))^a`. The computation below asserts that equality never occurs.

So if `P_d(X)` is the least `p` over such lattice clocks, every primitive cycle with period `< P_d(X)` has `m ≤ X`.

*Proof.* Proposition 3(c) and Lemma 2.1. ∎

**Computing `P_d(X)`.**
* It is computed exactly: loop over `a`, with `p ≡ π(a) (mod ord_d 2)`.
* The comparisons use 66-digit rationals, asserted to be at least `a·10^{−60}` away from a tie.
* The optimum is re-verified with integers.
* Clocks with `p/a ≥ 65/41` have window `≤ d/(2^{65/41} − 3) = 1192.08 d < 1200 d`. So `1200 d` covers every far clock, and also the whole main family of 17021, whose window `20 290 484` is the same for every multiple.

Final bounds, with `X_d` as in §2.4 (T2.7):

| `d` | `ω(d)` | `X_d` | least uncovered clock | complete for |
|---|---|---|---|---|
| 14303 | 944 | 2·10^7 | (203720, 128533) | `p < 203 720` |
| 17021 | 258 | 20 425 200 | (1631245, 1029201) | `p < 1 631 245` |
| 14197 | 329 | 115 744 533 | (80028, 50492) | `p < 80 028` |
| 6487 | 534 | 2·10^7 | (21004, 13252) | `p < 21 004` |
| 18359 | 164 | 22 030 800 | (15576, 9827) | `p < 15 576` |
| 10289 | 214 | 2·10^7 | (14176, 8944) | `p < 14 176` |
| 13085 | 335 | 61 647 379 | (14100, 8896) | `p < 14 100` |
| 7727 | 198 | 2·10^7 | (10835, 6836) | `p < 10 835` |
| 15655 | 207 | 105 271 617 | (10166, 6414) | `p < 10 166` |
| 9823 | 241 | 2·10^7 | (4974, 3138) | `p < 4974` |
| 7463 | 162 | 2·10^7 | (4245, 2678) | `p < 4245` |

Whether `ω(d)` is finite at all is Lagarias's Conjecture 2(2), which is OPEN for every `d`, `d = 1` included.

### 2.4 The whole table: one uniform sweep

**The sweep (T2.6).**
* Every `d ≤ 19999` prime to 6 (6667 systems) and every odd `y ≤ 1200 d`, run as 2 processes in 287 s wall time.
* No overflow, no unresolved start, and no orbit entering a cycle beyond `1200 d`.
* All 99924 cycles found, primitive or not, were re-verified by exact Python iteration.

**Results.**
* **Table (20):** all eleven `(d, ω)` pairs reproduced exactly.
* **Distribution:** `ω = 1` for 1481 systems and `ω = 2` for 1507, as published. Every `T_d` has a primitive cycle, which is Lagarias's Conjecture 2(1) in this range.
* **Residual:** `ω = 3` for 1004 systems (they report 1005), and 42757 primitive cycles in all (they report 42765).

**The extension (T2.7).**
* For 5663 values of `d`, least elements in `(1200d, X_d]` were searched too: a span of `1.436·10^11` integers, in 519 s.
* The bound is `X_d = max(1200d; 2·10^7 if d < 16667; the perigee windows of the 1252 clocks with p < 4000, log_2 3 < p/a < 65/41 and d | 2^p − 3^a)`.
* **No further cycle.**

**Consequence.** For every `d ≤ 19999`, the 42757 primitive cycles include every primitive cycle of period `p < 4000`. No clock of `L_d` with `p < 4000` has a window reaching `X_d + 1`; this was checked for all 6667 values of `d`. The cycles found also include the long ones beyond that, up to `p = 4686`.

### 2.5 The residual is two-sided

Suppose `ω_BM(d) ≥ ω_ours(d) ≥ 1` for every `d`, i.e. Belaga–Mignotte only ever found extra cycles.
* Then `{ω_BM = 1} ⊆ {ω_ours = 1}`. Both sets have 1481 elements, so they are equal.
* Given that, `{ω_BM = 2} ⊆ {ω_ours = 2}`. Both have 1507 elements, so they are equal.
* Given that, `{ω_BM = 3} ⊆ {ω_ours = 3}`, which would force `1005 ≤ 1004`.

That is impossible (T2.8), so there are two possibilities.
* **At least one `d` has `ω_BM(d) < ω_ours(d)`.** Then their table lacks one of our exactly verified cycles. To still reach the total, their list must contain at least 9 cycles outside ours.
  * By §2.4 such cycles would need period `≥ 4000` and least element `> X_d`, on near-critical clocks.
  * There the equidistribution expectation is about `2^{−0.05p}·poly(p)`, below `2^{−150}` at `p = 4000`. This is a heuristic, not a proof.
* **Their published summary numbers contain a misprint.** We regard this as the more likely explanation (an assessment, not a proof).

Their Table 1 (`hal-00129727`) would settle it. One request returned an Anubis bot-check page, which was not bypassed. So the residual stays UNRESOLVED; it does not affect table (20).

### 2.6 Long cycles (EMPIRICAL, with one proved bound)

**The sweep shows long primitive cycles are not rare.**
* 2204 of them have `p > 600`, 29 have `p ≥ 2000`, and the longest has `p = 4686`.
* They sit at `a/p ∈ [0.444, 0.561]`, with mean `0.502`. That is the maximum-entropy density, not a best approximation.
* For 683 values of `d` the only primitive cycle is such a long one. The four longest:
  * `d = 16819`: least 7, `(4686, 2292)`;
  * `d = 11491`: least 5, `(4531, 2253)`;
  * `d = 13829`: least 13, `(3918, 1954)`;
  * `d = 19427`: least 25, `(3890, 1948)`.

**Scaling (EMPIRICAL).** Over the 18697 primitive cycles with `p ≥ 100`:
* `mp/d` has quartiles 2.03, 4.96 and 10.30, so the least element is `m ≈ 5d/p`;
* the median of `a/p` is `0.500`;
* `p/d ≤ 1.18`.

**The proved half.** Every odd element of a cycle equals `(3y'+d)/2^{j+1} > d/2^{j+1}`, where `j` is the length of the run of even steps before it. So `m > d/2^{K+1}`, with `K` the longest such run; this was checked on 5000 cycles.

**Heuristic cutoff.** A random word has `K ≈ log_2 p`. So a least point `≥ 1` stops being automatic once `p` passes `d` in order of magnitude, which is consistent with `p/d ≤ 1.18`.

**What this means for enumerations.** The gates lane's side-aware model, applied in the dense regime `d ≫ 1`, must scan the clocks at `a/p ≈ 1/2` up to `p ≈ d`. A clock scan capped at `p ≤ 600` misses all 2204 of these cycles.

## 3. T3: Stern–Brocot positions

### 3.1 Paths

Positions are tested exactly: `a/p < log_q 2` iff `2^p > q^a`. `L` marks a best lower approximation, `U` a best upper one.
* **`q = 3`:** `1/1U 1/2L 2/3U 3/5L 5/8L 7/11U 12/19U 17/27L 29/46L 41/65L 53/84U 94/149L …`
* **`q = 5`:** `1/1U 1/2U 1/3L 2/5L 3/7L 4/9U 7/16U …`

### 3.2 Table

The censuses cover `Z`, both signs: `3x+d` by least element `|·| ≤ 10^6`, `5x±1` by all periods `p ≤ 60`. The scaled copies of the `3x+1` cycles inside `3x+d` are omitted; they sit at the same positions as for `3x+1`.

| map | least | `(p,a)` | reduced `a/p` | on path | `D = 2^p − q^a` | kind |
|---|---|---|---|---|---|---|
| `3x+1` | 0 | (1,0) | 0/1 | L | 1 | free |
| | 1 | (2,1) | 1/2 | L | 1 | free |
| | −1 | (1,1) | 1/1 | U | −1 | free |
| | −5 | (3,2) | 2/3 | U | −1 | free |
| | −17 | (11,7) | 7/11 | U | −139 | **sporadic** |
| `3x−1` | 0, −1, 1, 5, 17 | as for `3x+1` | | all on path | | negatives of `3x+1` |
| `5x+1` | 0 | (1,0) | 0/1 | L | 1 | free |
| | −1 | (2,1) | 1/2 | U | −1 | free |
| | 1 | (5,2) | 2/5 | L | 7 | sporadic (1 of 2 necklaces) |
| | 13, 17 | (7,3) | 3/7 | L | 3 | sporadic (2 of 5) |
| `5x−1` | 0 | (1,0) | 0/1 | L | 1 | free |
| | 1 | (2,1) | 1/2 | U | −1 | free |
| | −1, −13, −17 | (5,2), (7,3) | 2/5, 3/7 | L | 7, 3 | sporadic |
| `3x+5` | 1 | (3,1) | 1/3 | **no** | 5 | free |
| | 19, 23 | (5,3) | 3/5 | L | 5 | free |
| | 187, 347 | (27,17) | 17/27 | L | 5·71·14303 | partially free |
| `3x+7` | 5 | (4,2) | 1/2 | L | 7 | free |
| `3x+13` | 1 | (4,1) | 1/4 | **no** | 13 | free |
| | 211, 227, 251, 259, 283, 287, 319 | (8,5) | 5/8 | L | 13 | free (`L(8,5) = 7`) |
| | 131 | (24,15) | 5/8 | L | 13·186793 | partially free |
| `3x+23` | 5, 7 | (5,2) | 2/5 | **no** | 23 | free |
| | −2263, −2359, −2743, −2963, −3091, −3415, −3743, −4819 | (19,12) | 12/19 | U | −23·311 | partially free (8 cycles; `L(19,12)/311 = 8.5` expected) |
| | 41 | (43,26) | 26/43 | **no** | 23·271922921473 | partially free |

**The next upper approximant carries nothing (T3.3).** For `12/19`, `3^12 − 2^19 = 7153`, and none of the 50388 words of `(19,12)` has `7153 | c_w`. So no integer cycle of `3x±1` sits there.

### 3.3 The two patterns, typed

**Pattern 1: "integer cycles sit at best approximations where `|D|` is small relative to `d`."**
* **For `d = ±1` it holds (EMPIRICAL).** All five cycles of `3x±1` and all five of `5x±1` lie on the path. A non-free cycle needs `D | c_w`, which is plausible only where `|D|` is small relative to the number of words, that is, at best approximations.
* **For `|d| > 1` it fails as stated.**
  * Free cycles are forced by `D | d` wherever the shape sits. They land off the path at `1/3` (`d = 5`), `1/4` (`d = 13`) and `2/5` (`d = 23`).
  * Over the 42757 primitive cycles with `d ≤ 19999`, 9614 (22.5%) have their density on the path and 33143 do not (T3.4).
  * The informative statistic is the cofactor `D/d`, the part of the gap that `d` does not pay:

    | cofactor | cycles |
    |---|---|
    | `D = d` (free) | 2566 |
    | `1 < D/d < 2^10` | 4065 |
    | `2^10 ≤ D/d < 2^100` | 18925 |
    | `D/d ≥ 2^100` (the long cycles at `a/p ≈ 1/2`) | 17201 |
* **Corrected pattern (EMPIRICAL).** Cycles of `T_{q,d}` occur roughly in proportion to the equidistribution weight `L(p,a)·|d|/|D|`, restricted to words whose least point can be `≥ 1`.
  * For `d = ±1` this means the best approximations.
  * For large `d` it means `a/p ≈ 1/2` up to `p ≈ d`, together with the near-critical clocks where `d | D`.

**Pattern 2 (ANALOGY), the Kuratowski-style reading: "free obstructions forced by an identity, plus finitely many sporadic ones."**
* **The identity half is PROVED.** The free cycles of `T_{q,d}` are exactly those forced by identities `2^p − q^a = D` with `D | d`. There are finitely many (Pillai, Ellison), and for `d = ±1` there are four (Gersonides: `2−1, 3−2, 4−3, 9−8`).
* **The "finitely many sporadic" half is OPEN.** It is Lagarias's Conjecture 2(2), open even for `d = 1`.
* **The data show the non-free part is not a small, Kuratowski-like list for large `d`.** There are 943 non-free primitive cycles for `d = 14303`, and long cycles up to `p = 4686`.
* **So the reading fits only `d = ±1`.** There the sporadic list is `{−17,…}` alone, as far as is known.

## 4. T4: `5x+1` and `3x−1`

**`5x+1` on `Z`** (FINITE-EXACT: all periods `p ≤ 60`, both signs; T4.1). The cycles are exactly `0`, `{−1,−2}`, `{1,3,8,4,2}`, `{13,33,83,208,104,52,26}` and `{17,43,108,54,27,68,34}`. This agrees with the gates census for `p ≤ 40`.
* **Free:**
  * `{0}` on (1,0);
  * `{−1,−2}` on (2,1), with `D = 4 − 5 = −1` (1 of 1 necklace).
* **Sporadic:**
  * `{1,…}` on (5,2), with `D = 7` (1 of 2 necklaces integral);
  * `{13,…}` and `{17,…}` on (7,3), with `D = 3` (2 of 5).
* **The wave-13 claim "only the free cycles `{0}` and `{−1,−2}`" is PROVED** by Propositions 1 and 2. The only free mixed shape is `(2,1)`, and `(1,1)` gives `−1/3`.

**`3x−1`** (PROVED). `T_{3,−1}(−y) = −T_{3,1}(y)` identically, because `(−3y−1)/2 = −(3y+1)/2`. So the cycles of `3x−1` are the negatives of those of `3x+1`, and the wave-13 claim is confirmed. They are:
* `{1}` on (1,1), `D = −1`: free;
* `{5,7,10}` on (3,2), `D = −1`: free;
* `{17,25,37,55,82,41,61,91,136,68,34}` on (11,7), `D = −139`: sporadic, 1 of 30 necklaces;
* `{−1,−2}` on (2,1): free;
* `{0}`.

**Censuses (FINITE-EXACT).**
* `3x+1` on `Z`, all periods `p ≤ 90`: exactly the five known cycles (T4.5).
* All odd `y ≤ 10^9`: the positive cycles of `3x+1` are `{1,2}` only, and those of `3x−1` are `1, 5, 17` only (T4.6). This is a self-check of the code, far inside known verification bounds.

## 5. Reproduction

* **Command.**

  ```
  python3 -u 04-computation/experiments/procgen_sporadic_20260926_run.py > 05-knowledge/results/procgen_sporadic_20260926.out
  ```

  The flag `--quick` skips T2.6–T2.7 and runs in about 30 s.
* **Compilation and scratch files.** The runner compiles `procgen_sporadic_20260926_traj.c` and `procgen_sporadic_20260926_sweep.c` with `cc -O2` into `scratch/procgen_sporadic/`, and writes its intermediate files there (not for commit).
* **Recorded run.**
  * 70 checks passed, 886 s wall time.
  * Peak RSS: the Python process 183 MiB (191.9 MB, `/usr/bin/time`); the largest C child 45 MiB, which is the compiler. The search programs use about 1 MB each, with at most 2 running at once.
* **sha256** (also printed at the end of the `.out`):
  * `76429df107ff850c120ce3659dc9049616b4eac6df6f56743468601afef7e04f` `procgen_sporadic_20260926_lib.py`
  * `a3217c5d93f6e1733139dc8bb2e4f7d353e468d58426bf61c494a936fe6e426c` `procgen_sporadic_20260926_run.py`
  * `4343b7abbd0c0c42398731ad7ffac25962db4b65726f52befb129046eb7e7a95` `procgen_sporadic_20260926_traj.c`
  * `8a1257a418197c1c0df052047a828f989c6e522a7bfce42dfac8047c9b0042f9` `procgen_sporadic_20260926_sweep.c`
* **Disclosure: the first full run hung.**
  * *Cause:* `ordmod(2, 1)` never terminated for `d = 1`.
  * *Also fixed:* the cone test used `X` where `X + 1` is correct. This made it slightly conservative, and the `p < 4000` check would have failed on window-boundary clocks.
  * *Consequence:* both are fixed in the recorded run. The parts the aborted run completed (T1 to T2.7's extension) gave the same results, up to timings and the wording of the cone statements.
* **Web.** One request to `hal.science/hal-00129727/document`, with a generic User-Agent, returned the Anubis bot-check page. It was not bypassed.
