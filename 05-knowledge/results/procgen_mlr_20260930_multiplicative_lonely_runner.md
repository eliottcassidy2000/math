# The multiplicative lonely runner: ×2×3 runners between the Lonely Runner Conjecture and the Collatz cycle gates

**Status: PROVED (elementary, plus computer-assisted proofs in exact rational arithmetic) + FINITE-EXACT + EMPIRICAL + CITED, with one EMPIRICAL REFUTATION. The runner `procgen_mlr_20260930_run.py` ends with ALL CHECKS PASSED.**

**LRC side.** For multiplicative speed sets the LRC is arithmetic rather than lacunary.
* Every 3-smooth box `{2^j3^k : j<J, k<K}` with `J ≥ 3, K ≥ 2` has gap of loneliness exactly `κ = 1/5`, against the LRC value `1/(n+1)` (PROVED).
* The ×2×3 lonely spectrum `I(t) = inf ||2^j3^kt||` has the discrete top `{1/5, 1/7, 1/10, 1/11, 1/13, 1/14}`, and only 88 explicit rationals reach `≥ 1/14` (PROVED, computer-assisted).
* For discrete times `h/D` and large boxes, a lonely time at level `≥ 1/14` exists iff 5, 7, 11 or 13 divides `D` (PROVED).
* A 6-free counting criterion proves the discrete LRC `L ≥ 1/(n+1)` for all `D ≥ D₀(J,K)` whenever `J ≥ 3, K ≥ 2`, with margin `(3/2)(n+1)/(n+K+J/2+1/2) → 3/2`, and gives complete exception lists for small boxes (PROVED + FINITE-EXACT).

**Crowded side.**
* Near runners form disjoint triangles with exact residues, 6-free apex numerators and disjoint shadows, for thresholds up to `D/4`. This is the coprime-to-6 analogue of Tao's Lemma 7.4.
* It comes with a uniqueness lemma and a coset crowding bound (PROVED).
* The dichotomy "crowded ⇒ numerator `≤ U(θ,δ)` uniformly in `D`" is **REFUTED in boxes** (EMPIRICAL): the numerator of crowded times grows like `D^{0.5}`.

**Collatz side.**
* For the gate sums the analogous dichotomy holds on every tested clock. `|S(h)| ≥ 0.1C` forces an extended numerator `u_ext(h) ≤ 17`. This was tested on the 29 clocks with `p ≤ 24` where `0.1C` is above the random level (Conjecture D).
* The top frequency has `u_ext = 1` on all 67 clocks. This also identifies the gates lane's unexplained peaks.
* A dichotomy cannot exclude cycles (Proposition B, PROVED). Let `λ` be the energy fraction off the structured set. If `|S| ≤ εC` off that set, the minor arcs carry L¹ mass `≥ λ(1−C/M)/ε`, so the triangle inequality cannot certify `N = 0`. Parseval puts the sup off the set at `≥ √(λC(1−C/M))`.
* The data for `p ≤ 40` show the gap: a factor `√C ≈ 10³–3·10⁵`.

Session `collatz-procgen-20260922`, lane "mlr" (multiplicative lonely runner), 2026-09-30. No HYP or THM file was created and no existing file was edited.

---

## 0. Setting, notation, dictionary

* `D ≥ 5` is coprime to 6. A **time** is `h ∈ Z/D` (real time `h/D`). The **runner box** is `B(J,K) = {2^j 3^k : 0 ≤ j < J, 0 ≤ k < K}`, with `n = JK` runners.
* The runner `(j,k)` sits at `2^j 3^k h/D mod 1`. Its **centred residue** `r(j,k) ∈ (−D/2, D/2)` is the representative of `2^j 3^k h mod D`. Negative exponents mean inverses mod `D`, so `r` is defined on all of `Z²`.
* At an integer threshold `R` (`δ = R/D`), the runner is **near** if `|r| < R` and **far** otherwise. The **near set** is `N_R(h) ⊂ Z²`. The time `h` is **θ-crowded** in a box if at least a fraction `θ` of its runners there are near.
* The **best numerator** is `u*(h) = min_{B(J,K)} |r(j,k)|`, so that `h ≡ ±u* 2^{−j}3^{−k}`. This is the "small numerator on the ×2×3 orbit" of the task, measured inside the box.
* The **discrete maximal loneliness** is `L(D,J,K) = max_{h≠0} min_{B(J,K)} |r|/D = max_h u*(h)/D`.
* The **continuous** gap of loneliness is `κ(V) = sup_t min_{v∈V} ||tv||` (notation of the Perarnau–Serra survey). The classical LRC with `n` speeds asserts `κ(V) ≥ 1/(n+1)`.
* The **×2×3 lonely spectrum** is `I(t) = inf_{j,k ≥ 0} ||2^j 3^k t||` for `t ∈ R/Z`. It is the loneliness of `t` against the whole infinite semigroup of 3-smooth speeds.

**Dictionary with the gates lane** ([gate equidistribution](procgen_gates_20260925_gate_equidistribution.md)).
* Gate `M = |2^p − 3^a|`, words `w` with ones at `s_0 < … < s_{a−1}`, carry `c_w = Σ_i 3^{a−1−i} 2^{s_i}`.
* `S(h) = Σ_w e(h c_w/M) = Σ_w Π_i e(runner (s_i, a−1−i))`. Each word picks a **staircase path** through the box `[0,p) × [0,a)` from upper left to lower right, and `S(h)` multiplies the runners on it.
* Their structured set `Σ_3 = {±u 2^{−j}3^k : u ≤ 3, j ≤ p, k ≤ a}` consists of times whose extended runner box contains a residue `≤ 3`: `h ≡ u 2^{−j}3^k ⟺ r(j, −k) = u`.
* So the natural structural parameter for the gates is the **extended numerator** `u_ext(h) = min{|r(s,m)| : 0 ≤ s < p, −a ≤ m < a}`. Via `2^p ≡ 3^a` this window contains both their transport box and the runner box.

---

## 1. LRC side (existence of lonely times)

### 1.1 Continuous time: the boxes are trivially lonely, κ = 1/5

**Proposition K (PROVED).** Let `κ(J,K) = κ(B(J,K))`. Then

| box | κ(J,K) | lonely time | upper bound from |
|---|---|---|---|
| `J = 1` (odd speeds `3^k`) | 1/2 | `t = 1/2` | `{1}` |
| `J ≥ 2, K = 1` (powers of 2) | 1/3 | `t = 1/3` | `{1,2}` |
| `J = 2, K ≥ 2` | 1/4 | `t = 1/4` | `{1,2,3}` |
| `J ≥ 3, K ≥ 2` | **1/5** | `t = 1/5` | `{1,2,3,4}` |

More generally, `κ(V) = 1/5` for every finite set `V` of 3-smooth speeds containing a dilate `c·{1,2,3,4}`.

*Proof.*
* Upper bounds: `κ` is monotone under inclusion, and the AP `{1,…,m}` has `κ = 1/(m+1)` (Dirichlet; the classical tight instance, Perarnau–Serra §4, CITED).
* Lower bounds:
  * `J = 1`: all speeds are odd, and `t = 1/2` gives `1/2`.
  * `K = 1`: all speeds are coprime to 3, and `t = 1/3` gives `≥ 1/3`.
  * `J = 2`: the speeds are `3^k` and `2·3^k`, and `t = 1/4` gives `1/4` and `1/2`.
  * `J ≥ 3, K ≥ 2`: all speeds are coprime to 5, and `t = 1/5` gives `≥ 1/5`. ∎

* Checked exactly (A1) by the classical reduction `κ(V) = max_{N = v+v'} max_l min_v ||lv/N||`, for all `J ≤ 6, K ≤ 4` with speeds `≤ 2000`.
* The set `{t : min_V ||vt|| ≥ 1/5}` is exactly `{1/5, 2/5, 3/5, 4/5}`. This holds for the AP `{1,2,3,4}` and for all 126 boxes `3 ≤ J ≤ 16, 2 ≤ K ≤ 10` (exact rational interval arithmetic).
* The LRC margin `κ(n+1) = (n+1)/5` is unbounded. **For multiplicative speed sets the LRC holds trivially, by a mod-5 obstruction. Lacunarity is irrelevant**: the 3-smooth numbers are not lacunary, since consecutive ratios tend to 1.

### 1.2 The ×2×3 lonely spectrum: discrete at the top

By Furstenberg (1967, CITED), `I(t) > 0` forces `t ∈ Q`, since the ×2×3 orbit of an irrational is dense. So the lonely spectrum `{I(t)}` is a countable set, the ×2×3 analogue of a Lagrange spectrum.

**Theorem S (PROVED, computer-assisted, exact rational arithmetic).**
* (a) `I(t) ≤ 1/5` for all `t`, with equality iff `t ∈ {1,2,3,4}/5 + Z`.
* (b) `{I(t) : t ∈ R} ∩ [1/14, 1/2] = {1/5, 1/7, 1/10, 1/11, 1/13, 1/14}`.
* (c) Exactly 88 points `t ∈ (0,1)` have `I(t) ≥ 1/14`. Their denominators are `5, 7, 10, 11, 13, 14, 26, 28, 33, 52, 56`; the only ones coprime to 6 are 5, 7, 11, 13.
* Below `1/14`, the finite-exact census over all `r/N`, `N ≤ 1200`, continues `5/73, 7/104, 1/15, 1/17, 1/19, 5/97, 5/99, 11/219, 13/259, 1/20, …`.

*Proof.* Let `δ` be one of the levels `1/7, 1/10, 1/11, 1/13, 1/14`.
1. *Box.* `X_δ = {t ∈ [0,1] : ||vt|| ≥ δ ∀ v ∈ B(16,10)}` is computed exactly as a finite union of closed intervals with rational endpoints: 14, 22, 78, 230, 358 components. Any `t ∉ X_δ` has `I(t) < δ`.
2. *Centres.* Every component lies within `10^{−3}` of a rational `t₀ = r/N` with `N ≤ 300`. There are 10, 14, 34, 64, 92 centres.
3. *Local lemma.*
   * Let `N = N₀N'`, with `N₀` the {2,3}-part of `N`, and take a stabiliser `s* ∈ ⟨2,3⟩` with `s* ≡ 1 (mod N')`.
   * For multipliers `u ∈ N₀⟨2,3⟩` one has `(u s*^i) t₀ ≡ u t₀ (mod 1)`.
   * Hence if the open sets `J_u = {x > 0 : ||u t₀ + σux|| < δ}` cover the closed interval `[w/s*, w]`, then every `t = t₀ + σx` with `0 < x ≤ w` has some multiplier `us*^i` with `||us*^i t|| < δ`.
   * The cover is verified in exact arithmetic for every centre and both sides `σ = ±1`, with `w` the half-width of the components assigned to it. The first admissible `s*` that works is used, e.g. `s* = 6` for `t₀ = r/5` at level `1/7`.
4. *Values at the centres.* `I(t₀)` is computed from the finite orbit.

So `I(t) < δ` except at the centres with `I(t₀) ≥ δ`. Their values match the census exactly at every level. (a) follows from `X_{1/5}(B(3,2)) = {r/5}`. ∎ (A2, 17 s in total.)

**Corollary S′ (discrete time; PROVED).** Let `D` be coprime to 6. Suppose the box contains `B(16,10)` and every 3-smooth number `≤ 219·D`; for instance `2^J, 3^K > 219D` with `J ≥ 16, K ≥ 10`. Then:
* `L(D,J,K) ≥ 1/14` iff `D` has a prime factor `q ∈ {5,7,11,13}`;
* in that case `L(D,J,K) = 1/q` for the smallest such `q`.

*Proof.*
* Write `t = h/D`.
* If `t ∉ X_{1/14}(B(16,10))`, some runner of `B(16,10)` is closer than `1/14`.
* Otherwise `t` lies within `w` of a centre `t₀ = r/N`.
  * If `t ≠ t₀`, then `x = |t − t₀| ≥ 1/(ND)`. The local cover supplies a multiplier `us*^i ≤ u_max w/x ≤ 219·D` (A2 prints the worst factor 219), which lies in the box.
  * If `t = t₀`, then `N | D`. The four non-exceptional centres all have denominator 112 (A2 check), and `112` cannot divide `D`. So `N ∈ {5,7,11,13}`.
* For those denominators, `B(16,10)` covers `⟨2,3⟩ = (Z/q)^*`, so the box minimum equals `I(r/q) = 1/q`. ∎

* Checked on 10 instances (A3 ii). Examples: `D = 5·20011` gives `1/5`, and the prime `D = 340201` gives `L = 0.0217`.

**Reading.** For boxes of size `≥ log D + O(1)`, the LRC-type question "is there a time with all runners `≥ δ`?" is decided, for `δ ≥ 1/14`, by divisibility of `D` by 5, 7, 11 or 13.

**Collatz gates** (A3 iii; 117 boxes `(J,K) = (p,a)`, `12 ≤ p ≤ 21`, `D ≤ 2.2·10^6`).
* All gates with `5 | D` have `L = 1/5`; this is rigid, since the box contains `{1,2,3,4}`.
* With `7 | D` and `5 ∤ D`, `L = 1/7` on all 15 such boxes. With `11 | D` (and no smaller lonely prime), `L = 1/11` on all 7.
* With `13 | D`, `L = 1/13` on 3 of 5. The other two are thin boxes `(12,3)` and `(20,2)`, too small for the corollary, where `L` is larger.
* At `(15,8)`, `D = 26207 = 73·359` and `L = 5/73`, the seventh value of the census. The other 73-gate, `(21,4)`, is a thin box.
* The gate boxes are below the size Corollary S′ needs, but the spectrum is visible anyway.
* `L(n+1) ≥ 2.73` on all 117: the discrete LRC holds at every gate box computed.

### 1.3 Discrete time, generic moduli: size of L

**EMPIRICAL (A3 i).** For primes `D = 10007, 100003, 1000003` and boxes `J = c ln D/ln 2`, `K = c ln D/ln 3`:

| `D` | `c` | `n` | `L` | `L(n+1)` | `L·n/ln D` |
|---|---|---|---|---|---|
| 10007 | 1 | 104 | 0.0689 | 7.2 | 0.78 |
| 100003 | 1 | 170 | 0.0472 | 8.1 | 0.70 |
| 1000003 | 1 | 260 | 0.0438 | 11.4 | 0.82 |
| 1000003 | 2 | 1000 | 0.0131 | 13.1 | 0.95 |
| 1000003 | 3 | 2280 | 0.0062 | 14.2 | 1.03 |

* For small boxes (`c = 1/2`), `L → κ = 1/5` (e.g. 0.1876 at `D = 10^6+3`), as Proposition K predicts.
* In matched and larger boxes, `L ≈ (0.7–1.0) ln D/n`. This is 1.4–2 times the independent-runner prediction `ln D/2n`, because near runners cluster in triangles (Theorem T). The LRC margin `L(n+1) ≈ ln D` grows.
* BLMV (CITED via Wang) makes the decay rigorous only at the rate `(log log log D)^{−c}`, for boxes `m, n < 3 log D`.

### 1.4 Discrete multiplicative LRC by 6-free counting

The union bound over runners gives only `L ≳ 1/(2n)`, which is weaker than the LRC value `1/(n+1)`. The multiplicative structure recovers a factor 3/2.

**Theorem L (PROVED).** Let `D` be coprime to 6 and `R ≥ 1`. The number of `h ≠ 0` at which some runner of `B(J,K)` has `0 < |r| < R` is at most
`F_{J,K}(R) = 2 Σ_{1≤u<R, (u,6)=1} [JK + K·A_u + J·B_u + τ⁺_u]`,
where `A_u = #{a ≥ 1 : 2^a u < R}`, `B_u = #{b ≥ 1 : 3^b u < R}`, `τ⁺_u = #{a,b ≥ 1 : 2^a3^b u < R}`. Moreover
`F_{J,K}(R) ≤ (2R/3)(JK + K + J/2 + 1/2) + 4(JK + K log₂R + J log₃R + σ(R))`, with `σ(R) = #{2^a3^b < R} ≤ (log₂R+1)(log₃R+1)`.

Consequences:
* If `F_{J,K}(R) < D − 1` then `L(D,J,K) ≥ R/D`.
* If `(J−2)(K−1) ≥ 1`, i.e. `J ≥ 3` and `K ≥ 2`, then `L(D,J,K) ≥ 1/(JK+1)` for all `D ≥ D₀(J,K)`, and `liminf_D L(D,J,K)·(JK + K + J/2 + 1/2) ≥ 3/2`.

*Proof.*
1. Write a small residue as `r = ±u 2^a 3^b` with `(u,6) = 1`. Then `h ≡ ±u 2^{a−j}3^{b−k}`.
2. The set of exponent differences `{(a−j, b−k)}`, with `(a,b) ∈ T_u = {2^a3^b u < R}` and `(j,k)` in the box, has exactly `JK + K A_u + J B_u + τ⁺_u` elements. `T_u` is a down-set, so `(a,b) = (α⁺, β⁺)` is the cheapest preimage of `(α,β)`.
3. The majorant follows from three counts. Among `u < X` there are at most `X/3 + 2` integers coprime to 6. Next, `Σ_u A_u = Σ_{a≥1} #{u < R/2^a} ≤ R/3 + 2log₂R`, and similarly `Σ B_u ≤ R/6 + 2log₃R` and `Σ τ⁺ ≤ R/6 + 2σ(R)`.
4. Take `R = ⌈D/(n+1)⌉`. The bound is `< (n+1)(R−1) ≤ D − 1` once `cR − (n+1) − H(R) > 0`, where `c = (n+1) − (2/3)(n + K + J/2 + 1/2) > 0` exactly when `(J−2)(K−1) > 0`, and `H` is the explicit smooth lower-order term.
5. `H'` is decreasing, so a single verified point `R₀` with `H'(R₀) < c` gives the inequality for all `R ≥ R₀`. ∎

Checks (A4):
* `#bad ≤ F` exactly on 95 instances, with equality attained, so the bound is sharp as a set count.
* The majorant was checked numerically.
* Explicit `D₀`: `(3,2)`: 12020, `(4,3)`: 4824, `(8,5)`: 3035.
* Asymptotic margin `(3/2)(n+1)/(n+K+J/2+1/2)`: 1.05 at `(3,2)`, 1.24 at `(8,5)`, tending to `3/2`.

**Complete exception lists (FINITE-EXACT below `D₀`; Theorem L above).** The discrete multiplicative LRC `L(D,J,K) ≥ 1/(n+1)` holds for every `D ≥ 5` coprime to 6 except:

| box | `n` | `D₀` | exceptions |
|---|---|---|---|
| `(3,2)` = speeds `{1,2,3,4,6,12}` | 6 | 12020 | 11, 17, 37 |
| `(4,2)` | 8 | 6985 | 11, 13, 17, 19, 29 |
| `(3,3)` | 9 | 8061 | 11, 13 |
| `(4,3)` | 12 | 4824 | 17, 19, 29, 31 |
| `(5,3)` | 15 | 3905 | 17, 19, 23, 29, 31 |
| `(6,4)` | 24 | 3226 | 29, 31, 41 |

For instance, at `D = 11` the speeds `{1,2,3,4,6,12}` reduce to `{1,2,3,4,6}`, and `±{1,2,3,4,6}` covers `(Z/11)^*`, so `L = 1/11 < 1/7`. Continuous-time LRC is not affected (`κ = 1/5`). These are discretisation effects at moduli comparable to `n`.

### 1.5 Comparison with lacunary speed sets (CITED via Perarnau–Serra, arXiv:2409.20160 §7.3)

* Barajas–Serra: `v_{i+1} ≥ 2v_i` gives `κ ≥ 1/(n+1)`.
* Dubickas: `v_{i+1} ≥ (1 + 22 log n/n)v_i`, `n ≥ n₀`, gives `κ ≥ 1/(n+1)`. The variant condition (20), `v_{i+⌈(n+1)/12e⌉} ≥ (n+1)v_i` for `n ≥ 32`, suffices too.
* de Mathan and Pollington: for an ε-lacunary sequence, `sup_t inf_n ||t v_n|| ≥ δ(ε) > 0`. Peres–Schlag: `≥ cε/|log ε|`.

How the 3-smooth speeds compare:
* They are not lacunary. Their ratio tends to 1, and `|p log 2 − a log 3|` is small along the convergents of `log 3/log 2`. Yet `κ = 1/5` for every box, and `I(t) ≤ 1/5` with the discrete top spectrum of Theorem S. The bound comes from arithmetic (coprimality to 5), not from growth.
* Furstenberg rigidity makes the lonely set `{t : I(t) > 0}` countable (rationals). For lacunary sequences it is uncountable: Pollington and de Mathan, Hausdorff dimension 1 (UNCITED-RECOLLECTION).
* Dubickas' condition (20) does hold for very long boxes. In a box the window `[v, (n+1)v)` holds about `K log(n+1)/log 2` speeds, fewer than `(n+1)/12e` once `J ≳ 47 log(JK)`. These results are not needed here.

---

## 2. Crowded side (structure of crowded times)

### 2.1 One generator: runs (the task's lemma, sharp form)

**Lemma R (PROVED).** Let `0 < δ ≤ 1/4` and `y_j ∈ (−1/2, 1/2]` the signed residue of `2^j x`.
* If `|y_j| < δ` for `j₀ ≤ j ≤ j₁`, then `y_j = 2^{j−j₀} y_{j₀}` exactly and `|y_{j₀}| < δ 2^{−(j₁−j₀)}`.
* So `2^{j₀}x` lies within `δ2^{−(j₁−j₀)}` of an **integer**. The "rational with small denominator" of the task is `m/2^{j₀}`, and for `j₀ = 0` it is an integer.
* After a maximal run, the next `λ = max(1, ⌊log₂(1/2δ)⌋)` points are far, with exact residues `2^t y_{j₁}`. Indeed `δ/2 ≤ |y_{j₁}| < δ`, so `δ ≤ 2^t|y_{j₁}| < 1/2` for `1 ≤ t ≤ λ`.
* Hence a row of length `r` with near density `θ` contains a run of length `ℓ ≳ θλ/(1−θ)`. For `x = h/D` this gives `2^{j₀}h ≡ u` with `|u| < δD·2^{1−ℓ}`.

*Proof.* `|2y_j| < 1/2`, so `2y_j` is the signed residue of `2^{j+1}x`. Induct. ∎ (B1: 91639 runs checked.)

### 2.2 Two generators: the Triangle Lemma

**Theorem T (Triangle Lemma; PROVED).** Let `D ≥ 5` be coprime to 6, `h ≢ 0`, and `1 ≤ R ≤ D/4`.
* **(a) Engine.** If `(j,k)` and `(j+1,k)` are both near, then `r(j+1,k) = 2r(j,k)`. If `(j,k)` and `(j,k+1)` are both near, then `r(j,k+1) = 3r(j,k)`. If `(j,k)` is near and `|2r| < R` (resp. `|3r| < R`), the right (resp. upper) neighbour is near with residue `2r` (resp. `3r`).
* **(b) Triangles.** Every 4-connected component of `N_R(h) ⊂ Z²` is finite and equals
  `T(P₀,u) = {P₀ + (x,y) : x, y ≥ 0, 2^x3^y|u| < R}`,
  where `P₀` is its unique componentwise minimum (the **apex**) and `u = r(P₀)`. On `T`, `r(P₀ + (x,y)) = 2^x 3^y u` exactly.
* **(c) Apexes are 6-free.** `gcd(u, 6) = 1`.
* **(d) Shadows.** Let `Sh(P₀,u) = {P₀ + (x,y) : R ≤ 2^x3^y|u| < D/2}`. On it `r = 2^x3^y u` exactly. Shadows avoid the near set and are pairwise disjoint for distinct components.
* **(e)** `N_R(h)` is invariant under the relation lattice `Λ_D = {(j,k) : 2^j3^k ≡ 1 (mod D)}`.

*Proof.*
* **(a)** We have `|2r| < 2R ≤ D/2`. For the ×3 step: if `|3r| > D/2`, the centred representative is `3r ∓ D` with `|3r ∓ D| = D − 3|r| > D − 3R ≥ R`, which is far.
* **(b)** By (a), residues along any path inside a component change by exact factors `2^{±1}, 3^{±1}`. So `r(P') = 2^{j'−j}3^{k'−k} r(P)` for `P, P'` in the same component.
  * *Finite:* integrality of `r` bounds `j, k` from below, and `|r| < R` bounds them from above.
  * *Closed under meets:* take `P = (j,k)`, `P' = (j',k')` with `j ≤ j'`, `k ≥ k'`. Then `3^{k−k'} | r(P)`. The point `Q = (j,k')` satisfies `3^{k−k'} r(Q) ≡ r(P)`, so `r(Q) = r(P)/3^{k−k'}`: it is near, and its vertical segment up to `P` is near.
  * Hence the component has a least element `P₀`, and it lies inside `T(P₀,u)`. Conversely `T(P₀,u)` is a staircase reachable from `P₀` by near steps, using (a).
* **(c)** If `2 | u`, then `r(P₀ − (1,0)) ≡ 2^{−1}u = u/2` is near and adjacent to `P₀`, contradicting minimality. The same argument works for 3.
* **(d)** Exact residues hold by induction along unit steps, since all values stay `< D/2`. If `P` lies in two shadows, then `2^{x₀}3^{y₀}u₀ = 2^{x₁}3^{y₁}u₁`. By (c) and unique factorisation the two apexes coincide.
* **(e)** Clear. ∎

* **Checks (B1).** 161537 interior components in 2568 `(D,h,R)` instances, with `R/D ∈ {1/4, 1/6, 1/16, 1/40}` and `D ∈ {101, 7³, 7·11·13, 1009, 5⁵, 4001, 7553}`. Claims (a)–(d) hold in every instance; (e) is immediate.
* **Sharpness.** At `R = D/4`, none of 19510 adjacent near pairs violates (a). At `R = 0.3D` and `D/3`, the ×3 step wraps onto near points: 5588 of 28595 and 9429 of 34843 adjacent pairs violate the exact relation. So the threshold `1/4` is sharp for the engine.
* **Relation to Tao.** This is the coprime-to-6 analogue of **Tao 2019, Lemma 7.4** (arXiv:1909.03562, "Structure of black set"; CITED, read in primary text).
  * Tao's points with `|θ(j,l)| ≤ ε`, for ×9 and ÷2 acting on `Z/3^n`, form a disjoint union of triangles separated by `≥ (1/10)log(1/ε)`. His "black" = our near.
  * Theorem T adds three things: exact residues, the arithmetic of apexes (c), and disjoint shadows of width `log(1/2δ)` replacing his separation.

### 2.3 Uniqueness (Diophantine window)

**Lemma U (PROVED).** If `|r(P₀)|, |r(P₁)| ≤ U` and `U·2^{|j₁−j₀|}3^{|k₁−k₀|} < D/2`, then `r(P₁) = 2^{j₁−j₀}3^{k₁−k₀} r(P₀)`. The two small numerators describe the same orbit point. In particular, two distinct apexes with numerators `≤ U` are at least `log(D/2U)` apart in the metric `|Δj| log 2 + |Δk| log 3`.

*Proof.* Clear the negative exponents. Both sides are `≡ 2^{…}3^{…}h (mod D)` and have absolute value `< D/2`, so they are equal integers. For two apexes, (c) forces `Δ = 0`. ∎

B2 checks this on 19101 pairs (`D ∈ {1009, 4001, 30011, 100003}`, `U ∈ {3, 12}`, window `|j| ≤ 7, |k| ≤ 5`): no window contains two 6-free numerators.

This is the ×2×3 analogue of Legendre-type uniqueness of good rational approximations.

### 2.4 Coset crowding bound

**Proposition C (PROVED).** Let `H = ⟨2,3⟩ ≤ (Z/D)^*`, `h` a unit, `m = min_{x ∈ hH}|x|` (centred), and `1 ≤ R ≤ D/4`. Then
`ρ(hH) := #{x ∈ hH : |x| < R}/|H| ≤ max_{u ∈ [m,R), (u,6)=1} σ(R/u)/σ(D/2u)`.

*Proof.*
* `(j,k) ↦ 2^j3^kh` induces a bijection `Z²/Λ_D → hH`.
* Each `Λ_D`-class of components contributes `σ(R/|u|)` near elements.
* Its `T ∪ Sh` contributes `σ(D/2|u|)` distinct elements, namely the distinct integers `2^x3^yu` of absolute value `< D/2`.
* Non-equivalent classes contribute disjoint sets, by Theorem T(d).
* Apex numerators are 6-free and lie in `[m, R)`. ∎

* **Checks (B3).** All 3053 `(coset, R)` pairs with `D < 1500` satisfy the bound; the maximum of `ρ/bound` is 0.83.
* **Example (`D = 10^6 + 3`, `R = D/16`).** A coset with near fraction `≥ 1/2` contains an element `|x| ≤ 1157`, and one with near fraction `≥ 1/4` contains an element `≤ 20833`.
* **What the bound is.** It is uniform in `D` and in `|H|`, but only polynomially strong: crowding forces a small element of size `D^{1−c(θ)}`, not a bounded one. For large `H` (`|H| ≥ D^ε`, `D` prime) the truth is much stronger, by Bourgain–Glibichuk–Konyagin-type equidistribution (UNCITED-RECOLLECTION; not used).

### 2.5 The crowded/structured dichotomy in boxes (EMPIRICAL; the uniform version is REFUTED)

**Single-triangle law (PROVED lower bound, EMPIRICAL equality).** The time `h = u` has the near runners `{(j,k) : 2^j3^k|u| < R}`, so `ρ_R(u) ≥ #{(j,k) ∈ B : 2^j3^k|u| < R}/n`.
* With equality: in boxes smaller than the triangle the equality is exact. At `D = 10^6+3`, box `(10,6)`, `R = D/16`: `u = 1, 5, 25` give `ρ = 0.983, 0.883, 0.767` against `0.983, 0.883, 0.750`.
* In matched boxes the excess is the background `≈ 2δ` plus wrap-around triangles.
* For fixed `u`, the fraction tends to the `u = 1` value as `D → ∞`.

**Uniform-U test (B4).** For each `θ` above background, take the largest best numerator `u*` among θ-crowded times in the matched box (`c = 1`, `R = D/16`):

| `θ` | `D = 10^4` | `D = 10^5` | `D = 10^6` |
|---|---|---|---|
| bg + 0.10 | 41 | 473 | 2353 |
| bg + 0.20 | 13 | 65 | 71 |
| bg + 0.30 | none | none | 1 (only `u* = 1`) |

At `θ = bg + 0.10` the numerator grows like `D^{0.5–0.6}` (`η = log u*/log D ≈ 0.40–0.56`).

**Verdict.**
* In boxes, "θ-crowded ⇒ `h` is on the ×2×3 orbit of `u` with `|u| ≤ U(θ,δ)`, uniformly in `D`" is **REFUTED (EMPIRICAL)**.
* The single-triangle law explains why: crowding decays only like `(1 − log u/log D)²`.
* The correct form is logarithmic (EMPIRICAL): θ-crowded ⇔ `log u*/log D ≤ η(θ) + o(1)`. The PROVED ingredients are the single-triangle lower bound (small `u*` ⇒ crowded) and, for whole cosets, Proposition C (crowded ⇒ a small element).
* Only the most crowded level (`θ` near the `u = 1` value) singles out `u* = 1`.

---

## 3. Collatz side

### 3.1 The gate sum is a runner path product: lonely ⇒ quiet, loud ⇒ crowded along paths

* `S(h) = Σ_w Π_i e(r_i/M)`, where `r_i` is the residue of the runner `(s_i, a−1−i)` at time `h/M`.
* **Switching bound.** Gates lane, Proposition SW (PROVED there): `|S(h)| ≤ C·B(h)` with
  `B(h) = E_w Π_{mixed pairs k of w} |cos(π r_{(2k, a−1−i_k)}/M)|`.
  A pair `(2k, 2k+1)` is mixed when it holds exactly one 1, of index `i_k`.

**Lemma P (PROVED).** Fix `δ ∈ (0, 1/2)`. Let `F_w(h)` be the number of mixed pairs of `w` whose runner is far (`|r| ≥ δM`). If `|S(h)| ≥ εC`, then for a uniformly random word
`P_w(F_w(h) ≤ L) ≥ ε/2`, where `L = ⌈log(2/ε)/log(1/cos πδ)⌉`.

*Proof.* `ε ≤ B(h) ≤ E cos(πδ)^F ≤ P(F ≤ L) + cos(πδ)^{L+1} ≤ P(F ≤ L) + ε/2`. ∎

**Consequences.**
* **Loud ⇒ crowded along paths.** A loud frequency is a time at which a positive proportion of the staircase paths meet only `O_{ε,δ}(1)` far runners in their mixed pairs. By Theorem T, the near runners of such a path lie in triangles with 6-free apex numerators.
* **Lonely ⇒ quiet.** If every runner is far, then `|S(h)| ≤ C·E cos(πδ)^{K_w}`, with `K_w` the number of mixed pairs. This is exponentially small in `p`.
* **The LRC side's lonely times are quiet frequencies.** If `5 | M`, then `h = rM/5` gives `S = Σ_w e(r c_w/5)`, the mod-5 distribution of carries, which is equidistributed at an exponential rate (gates SW).

**Data (C2).** Three clocks, top 10 and 60 random `h`:
* Top frequencies: median `B ≈ 0.53`, path-weighted near fraction `P_near ≈ 0.64–0.77` at `δ = 1/8`.
* Random `h`: median `B ≈ 0.12–0.15` (a weak bound, since the true `|S|/C ≈ 10^{−3}`) and `P_near ≈ 0.25 = 2δ`, the background.
* All 210 checked `h` satisfy `|S|/C ≤ B`.

### 3.2 The large values sit on 6-free orbit numerators (EMPIRICAL, full spectra)

Coverage: 62 clocks with `14 ≤ p ≤ 20`, `M ≤ 2^20`, `C ≥ 1000`, computed for all `h` with the certified DP (C1); plus 5 clocks with `p = 21–24` (C1b).

* The **top frequency has `u_ext = 1` on all 67 clocks**. In the 62 smaller clocks, all of the top 20 have `u_ext ≤ 50`.
* 72280 of the 74214 frequencies (97.4%) with `|S| ≥ 0.05C` and `u_ext < √M` have **`u_ext` coprime to 6**. This is Theorem T(c), apexes are 6-free, seen in the spectrum. The remaining 2.6% are triangle pieces cut by the window.
* The top values run through `u_ext ∈ {1, 5, 7, 11, 13, 17, 31, 35, 37, 47, 59, …}`, all 6-free.
* **The gates lane's unidentified peaks, resolved.** Their decomposition `h = u 2^{−j}3^k` (`j ≤ p`, `k ≤ a`) misses the orbit points `h ≡ ±2^{−s}3^{−m}` with `s, m > 0`.
  * Ten top-4 frequencies with `p ≤ 16` and gates-`|u|` ∈ {56, 125, 265, 1151, 3613} have `u_ext = 1`.
  * The `(24,15)` peaks `u = ±18431` (§3 of their note: "we did not identify their structure") also have `u_ext = 1`: `h = 786952` satisfies `h·2²·3³ ≡ 1` and `h = 1094238` satisfies `h·2·3⁴ ≡ −1 (mod 2428309)` (C1c).
  * Their "partially structured" `u = −37` at `(27,17)`, `−95` at `(25,15)` and `−175` at `(29,18)` are genuine 6-free numerators: `u_ext = 37, 95, 175` (C1c).
* **Two-piece chains (C1c).** At `(20,8)`, `h = 141509` has two path-carrying near triangles: apex `(0,4)` with `u = 64` (36 points, path weight 0.46) and apex `(14,0)` with `u = 81` (34 points, path weight 0.30). They are one orbit point `h ≡ 2^6 3^{−4} ≡ 2^{−14}3^4` (`u_ext = 1`), seen on the two sides of `2^20 ≡ 3^8 (mod M)`.

### 3.3 The S-dichotomy function: bounded on the tested range (EMPIRICAL)

`U_ε = max{u_ext(h) : |S(h)| ≥ εC}`, reported only when `ε ≥ 4√(ln M/C)`, i.e. above the random level:

| clocks | `U_{0.10}` | `U_{0.05}` | `U_{0.03}` |
|---|---|---|---|
| `p = 17–20` (24 clocks) | median 6–7, max 17 | ≤ 175 (7 clocks) | — |
| `(21,13)` | 1 | 13 | — |
| `(21,11)` | 5 | 47 | 781 |
| `(22,13)` | 1 | 29 | 295 |
| `(23,14)` | 1 | 13 | 295 |
| `(24,15)` | 1 | 1 | 35 |

* On every tested clock where `0.1` is above four times the random level, `|S(h)| ≥ 0.1C ⇒ u_ext(h) ≤ 17`, and `|S| ≥ 0.05C ⇒ u_ext ≤ 175`.
* There is no growth with `p`.
* This contrasts with box crowding (§2.5), where the bound grows like `D^{0.5}`.

**Why the two dichotomies differ (heuristic, consistent with all data).**
* By Lemma P, a loud frequency has typical staircase paths almost entirely near.
* A staircase crosses the whole box, which has log-length `≈ log M`. Its near cover must therefore be one triangle of log-size `≈ log M`, or a short chain of pieces of one orbit point. Either way the numerator is bounded.
* Box crowding at a fixed fraction `θ` is produced by a triangle of log-size `(1−η) log M`, i.e. by numerators `M^η`.

**Conjecture D (OPEN; tested for `p ≤ 24`).** For every `ε > 0` there is `U(ε)` such that for all clocks `(p,a)` and all `h`, `|S(h)| ≥ εC ⇒ u_ext(h) ≤ U(ε)`.
* A proof would adapt Tao's renewal argument (arXiv:1909.03562, Prop. 7.3 and Lemmas 7.4–7.10) to fixed-weight words at the gate modulus.
* `U(ε)` should be set by the limit profile of the structured values `S(u2^{−s}3^{−m})/C`: the archimedean factor `E e(uY_θ)` of the gates lane's Proposition A, times 2-adic and 3-adic digit factors.

### 3.4 What a dichotomy buys for cycles: nothing beyond the square-root barrier

**Proposition B (PROVED, elementary).** Let `Σ ⊂ Z/M` with `0 ∈ Σ`, `E_Σ = (1/M)Σ_{h∈Σ} S(h)`, so that `N = E_Σ + (1/M)Σ_{h∉Σ} S(h)`. Then
`(1/M) Σ_{h∉Σ} |S(h)| ≥ LB_Σ := (M·Coll − Σ_{h∈Σ}|S(h)|²) / (M·max_{h∉Σ}|S(h)|)`,
where `Coll = Σ_r n_r² ≥ C`. Two consequences follow.
* Suppose `Σ∖0` carries at most a fraction `1 − λ` of the non-zero energy, and a dichotomy gives `max_{h∉Σ}|S(h)| ≤ εC`. Then `LB_Σ ≥ λ(1 − C/M)/ε`. So the triangle-inequality certificate `N ≤ |E_Σ| + (1/M)Σ_{h∉Σ}|S(h)| < 1` is **impossible as soon as the dichotomy has content** (`ε < λ(1 − C/M)`). The better the dichotomy, the larger the minor-arc L¹ mass.
* Parseval forces `max_{h∉Σ}|S(h)| ≥ √(λC(1 − C/M))`, which is the square-root barrier. A certificate needs `< 1`. The gap is a factor `√C ≈ 2^{0.47p}` near the critical line.

*Proof.* `Σ|S| ≥ Σ|S|²/max|S|`, together with `Σ_{h∉Σ}|S|² = M·Coll − Σ_{h∈Σ}|S|²`. ∎

**Checks (C3a; 62 full-spectrum clocks, `U ∈ {1, 7, 50}`, `Σ_U = {u_ext ≤ U}`).**
* `L1_off ≥ LB_Σ` holds on every row, as an exact identity using `Coll` from the histogram.
* In all 183 cases where `Σ_U` carries `≤ 90%` of the energy, `LB_Σ ≥ 1` (range 2.64–21.48), so no certificate of `N = 0` is possible.
* The minor-arc L¹ mass is `0.64·√C·M` (median at `U = 7`), between `12.5M` and `261M`.

**The structured set is not a good major arc.**
* `E_U` fluctuates in sign and size.
* At `(16,10)`: `C/M = 1.23`, `N = 0`, but `E_1 = −3.72`.
* At `(19,12)`: `C/M = 7.04`, `N = 0`, but `E_1 = +4.01` and `E_7 = −1.04`.
* The reason: the kernel `(1/M)Σ_{h∈Σ} e(hx/M)` of a ×2×3-orbit set is neither positive nor localised. Even the "main term" is not a singular-series prediction.

**Sparse clocks `p = 24–40` (C3b).** `Σ` = numerators `u ∈ {1, 5}` on the extended window, `|Σ| = 2304–8320`, with `E_Σ` computed exactly by DP. All these clocks carry `N = 0` (gates-lane census).

| `(p,a)` | `log₂ M` | `log₂ C` | `E_Σ` | minor-arc RMS / √C | RMS of `|S|` |
|---|---|---|---|---|---|
| (24,12) | 24.0 | 21.4 | +0.288 | 1.06 | 1.8·10³ |
| (24,15) | 21.2 | 20.3 | −0.836 | 0.71 | 8.1·10² |
| (28,17) | 27.1 | 24.4 | −0.053 | 0.96 | 4.4·10³ |
| (32,20) | 29.6 | 27.8 | −0.345 | 1.30 | 2.0·10⁴ |
| (36,22) | 35.1 | 31.8 | −0.113 | 0.84 | 5.2·10⁴ |
| (40,25) | 37.9 | 35.2 | −0.156 | 0.80 | 1.6·10⁵ |
| (40,26) | 40.4 | 34.4 | +0.098 | 0.72 | 1.1·10⁵ |

* `E_Σ` is an exactly computable `O(1)` number, and `N = E_Σ + (signed minor-arc sum)`.
* A cycle would be excluded by `|minor| < 1 − |E_Σ|`, but the minor-arc terms have RMS `0.6–1.3·√C`, i.e. `|S| ≈ 10³–3·10⁵`.

**Answer to "what crowded-time dichotomy would exclude cycles at the gate (p,a)?".**
* None of sup-norm type. By Proposition B, a dichotomy, however explicit its `U` and `θ`, only bounds `|S|` off the structured set, and the triangle-inequality route then fails by `√C`.
* What would exclude a cycle is `|Σ_{h∉Σ} S(h)| < M(1 − |E_Σ|)`. Up to the computable `E_Σ`, this cancellation statement is `N = 0` itself.
* The transfer-matrix structure computes any single `S(h)` in `O(pa)`, but the minor arcs comprise `M − |Σ| ≈ 2^p` frequencies.
* Remark (heuristic, not proved or computed): take the Fourier-majorant variant, `N ≤ Σ_w φ(c_w)` with `φ ≥ 0`, `φ(0) ≥ 1`, `supp φ̂ ⊂ Σ'`. By LP duality its optimum is the largest atom at 0 of a non-negative measure with the same Fourier data on `Σ'`. A dimension count suggests this atom is `≳ C/|Σ'|` unless `|Σ'| ≳ C`, i.e. exponentially many frequencies.

---

## 4. Consequences each way, and assessment

**Runner picture → Collatz gates.**
* Theorem T and Lemma P make "crowded" precise. A loud frequency is a time whose typical staircase paths run through near triangles, and their apexes carry 6-free numerators (PROVED).
* The data say these numerators are bounded (Conjecture D, `p ≤ 24`). The top of the spectrum is now completely accounted for, including the gates lane's unidentified peaks.
* Proposition B shows that no statement about `|S|` reaches the cycle count. The obstruction is cancellation over about `2^p` minor-arc frequencies.
* This is the gates lane's verdict, sharpened: the dichotomy route is self-defeating, not merely insufficient.

**Gate picture → LRC.**
* The quiet frequencies of the gate spectrum are lonely times.
* Lonely times at a fixed level exist only for arithmetic reasons (Corollary S′): `5 | M` gives `L = 1/5` (29 of 117 computed gate boxes), and similarly 7, 11, 13 and, at `(15,8)`, 73.
* For large boxes, a gate with no prime factor `≤ 13` has no lonely time at level `1/14`.

**Honest assessment.**
* **Collatz.** Diagnostic only: rank LOW as a proof angle, MEDIUM as an instrument. The hybrid identifies exactly which frequencies are loud and why, but cycles live in the signed minor-arc sum, which no crowding or loneliness statement controls.
* **LRC.** Multiplicative speed sets are the opposite extreme from the tight LRC instances, the arithmetic progressions.
  * The LRC holds continuously with margin `(n+1)/5`, and discretely with a proved margin tending to `3/2` (Theorem L).
  * The tools (mod-`q` obstructions, Furstenberg rigidity, 6-free counting) are specific to multiplicatively closed speed sets.
  * The one transferable idea (a heuristic generalisation, not proved here): for a speed set closed under multiplication by a set of primes `P`, the union bound should improve by `1/Π_{q∈P}(1 − 1/q)`. That is irrelevant for additive tight instances. Rank LOW.
* **What is genuinely new lives inside the hybrid.**
  * The ×2×3 lonely spectrum, whose top is the discrete set `{1/5, 1/7, 1/10, 1/11, 1/13, 1/14}`, and its discrete-time classification.
  * The discrete multiplicative LRC with complete exception lists.
  * The Triangle Lemma with 6-free apexes, disjoint shadows and the coset bound.
  * The refutation of the uniform box dichotomy, contrasted with the bounded S-dichotomy.

---

## 5. Literature (read and cited)

* **Tao, arXiv:1909.03562**, §7 (read in primary text, PDF fetched 2026-09-30).
  * Lemma 7.4, "Structure of black set": `{|θ(j,l)| ≤ ε}` is a disjoint union of triangles `(7.11)` separated by `≥ (1/10) log(1/ε)`. The setting is ×9 and ÷2 acting on `Z/3^n`; his black = our near.
  * Prop. 7.3 and Lemmas 7.9–7.10: the renewal process meets many white points.
  * Theorem T is the coprime-to-6, two-generator analogue, with exact residues, 6-free apexes and disjoint shadows.
* **Perarnau–Serra, "The Lonely Runner Conjecture turns 60", arXiv:2409.20160v3** (read in primary text).
  * `κ` and `κ_N` (§3.1, eqs. (5), (7)); the reduction to `t = ℓ/(v_i+v_j)` (eq. (4)), credited there to Haralambis and to Czerwiński–Grytczuk Thm 6.
  * AP tightness `κ([n]) = 1/(n+1)` (§4).
  * Lacunary results (§7.3): Barajas–Serra Thm 18; de Mathan and Pollington; Peres–Schlag `cε/|log ε|`; Dubickas Thm 19 and condition (20); Czerwiński.
* **Bourgain–Lindenstrauss–Michel–Venkatesh**, "Some effective results for ×a×b", ETDS 29 (2009) 1705–1722. **Not on arXiv.** The statement is quoted from **Z. Wang, arXiv:1004.0035, Thm 1.4** (secondary source):
  * (i) for Diophantine-generic `x`, `{a^m b^n x : 0 < m,n ≤ N}` is `(log log N)^{−c}`-dense for `N ≥ N₀`;
  * (ii) for `x = p/Q` with `Q` coprime to `ab`, `{a^m b^n x : 0 < m,n < 3 log Q}` is `C(log log log Q)^{−c}`-dense.
* **Furstenberg 1967**, Math. Systems Theory 1:1–49: an irrational has a dense ×2×3 orbit. Standard statement; the reference was taken from the bibliographies of arXiv:1607.00670 and arXiv:1004.0035.
* **Gates lane:** switching bound (Prop. SW), census `p ≤ 40`, archimedean law (Prop. A), structured set `Σ_3`.
* **Seen in search results only, not used:** Rosenfeld arXiv:2509.14111 (eight runners); arXiv:2511.22427 (nine and ten runners); arXiv:2511.16636 (Riesz products); Kravitz–Leng arXiv:2510.17744.
* **Priority.**
  * Proposition K is folklore-level: AP tightness plus coprimality.
  * Not found in the sources read: Theorem S (the values of the ×2×3 lonely spectrum); Corollary S′; Theorem L and its exception lists; Theorem T(c),(d) and Proposition C in this form; the S-dichotomy data.
  * The lonely-spectrum question is related to minimal norms of residue classes (Bourgain–Chang; Wang 1004.0035, Appendix). **No priority is claimed.**

## 6. OPEN

1. **Conjecture D** (bounded S-dichotomy, uniform in `p`) and its proof by a Tao-type renewal argument. Also the limit profile of the structured coefficients `S(u2^{−s}3^{−m})/C`.
2. **The lonely spectrum below 1/14.** The census continues `5/73, 7/104, 1/15, 1/17, 1/19, 5/97, …`. Certifying it needs larger boxes. Several of these values come from moduli where `⟨2,3⟩` is a proper subgroup of the units (73, 97, 259, 601, …).
3. **Asymptotics of `L(D,J,K)`** for primes `D` and boxes of size `≍ log D`: empirically `(0.7–1.0) ln D/n`. Proving `L ≫ ln D/n` needs correlation inequalities (Lovász-type, as in Peres–Schlag) beyond Theorem L's counting.
4. **A box version of Proposition C** with explicit boundary terms, and the exact crowding law `η ↦ θ`.
5. **Corollary S′ for matched gate boxes `(p,a)`.** These are below the size the proof needs. The data agree anyway: `L = 1/5` whenever `5 | M`, and `1/7`, `1/11`, `5/73` occur as predicted.

## 7. Reproduction

* **Command:** `python3 04-computation/experiments/procgen_mlr_20260930_run.py > 05-knowledge/results/procgen_mlr_20260930.out`. It runs three parts, one `nice`d process at a time, and ends with `ALL CHECKS PASSED`.
* **Scripts** (`04-computation/experiments/`):
  * `procgen_mlr_20260930_core.py` (library);
  * `_lrc.py` (Part A: §1);
  * `_crowd.py` (Part B: §2);
  * `_gates.py` (Part C: §3; imports the gates lane's `procgen_gates_20260925_core/_stats/_switching` read-only);
  * `_run.py`.
* **Output:** [procgen_mlr_20260930.out](procgen_mlr_20260930.out).
* **Runtime and memory:** about 4.5 minutes in total (268 s). Peak RSS is below 500 MB per part (A 295 MB, B 263 MB, C 415 MB).
* **Disclosure.** An exploratory version of Part C that used an FFT of length about `10^6` peaked at 524 MB. It was replaced by the certified blocked DP of the gates lane.
* **Literature texts** (Tao, Perarnau–Serra, Wang, arXiv:1607.00670) are in `scratch/procgen_mlr/lit/` and are not for commit.
* **Web access.** Two `WebSearch` queries restricted to arxiv.org, to locate the papers. Every download was from arxiv.org or export.arxiv.org with the generic User-Agent `Mozilla/5.0 (research; math-repo)`. No personal data was sent.
