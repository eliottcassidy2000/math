# AMM 12592: golden-zero super-blocks beat the golden constant (HYP-9128 proved, C* ≤ 159/100), the capacity gap extended to C* ≥ 1.377, and why no termwise majorant can extend it further

**Status.**
* **AUDIT and PROMOTION (orchestrator, 2026-09-23).** Lemma S was re-derived. Independent exact checks of the `N = 16` and `64` blocks from the fair-coin definition pass: realizable, deadlines `<= ceil(1.59L)`, `Phi(p) = Phi(1-p)`, and `Phi = w^N S_N(w)` ([`..._hyp9128_orchestrator_check.py`](../../04-computation/experiments/amm12592_procgen_20260923_hyp9128_orchestrator_check.py)). The contour certificate re-runs identically. The `gamma = 0.377` certificate re-verifies. Promoted as [THM-4468](../../01-canon/theorems/THM-4468-golden-zero-superblocks-beat-golden-amm12592.md).
* **PROVED — HYP-9128 (Theorem 1).**
  * **Claim.** There is an exactly fair, deterministic, complement-symmetric extractor with `T(L) ≤ ⌈159L/100⌉` for every `L ≥ 16` and `T(L) ≤ 2L` for `L < 16`. Hence **`C* ≤ 159/100 < C_* = 1 + log_5 φ² = 1.5979874…`**. The golden constant is not the uniform optimum.
  * **Ingredients of the proof:**
    * Lemma R, the generalised lattice rounding. It is re-derived in §3 and not cited.
    * The fold identities, Lemmas F, K and L0 in §4.
    * Interval-arithmetic certificates, uniform in `v = i/N`, for every `N ≥ 4096` (§5; Part 1 of the `.out`).
    * Exact finite certificates for `16 ≤ N ≤ 2048` (§6; Part 2).
    * For the prefix `L < 16`: C. D. Long's three distributed ratio-2 blocks, re-checked exactly (F5).
  * **Classical input CITED:** Robbins' two-sided Stirling bounds, used only for `C(n,k) ≥ e^{nH(k/n)}/√(2n)`.
  * **Rigor caveat.** Part 1 uses `mpmath.iv` interval arithmetic at 60-bit precision with outward rounding, together with the elementary monotonicity-in-`N` argument of §5.3. Everything else is exact integer or rational arithmetic.
* **FINITE-EXACT:**
  * margin certificates for `N = 32, …, 2048`;
  * explicit integral super-blocks for `N = 16, 32, 64, 128, 256`. These come from an independent implementation of Lemma R and are verified exactly: box, parity, boundary `±1`, and class polynomial equal to the palindromic target;
  * exact checks F1–F6 of every identity and structural hypothesis used in the proof, with structure checked for all `N = 16·4^k ≤ 2^22`;
  * the kernel-saturation examples of §9.
* **PROVED, with the FFT caveat of THM-4467's `3/8` certificate — THM-4467 at `γ = 377/1000`.**
  * A new subordination certificate gives `R(W_γ, 0) ≥ b_1 = 1.0001609 > 1`, hence **`C* ≥ 1.377`**. This is `5·10^−4` below the method's reach `γ_1 = 0.3775246`.
  * The proved window is now **`1.377 ≤ C* ≤ 1.59`**.
* **PROVED (negative) — structured majorant (§9).**
  * No termwise majorant `|W_m(p)| ≤ Φ(d_m, p)` valid for all fair extractors can beat THM-4467's `S(p)^{d_m}` by an exponential factor. This holds at any `p ∉ [0,1]` and for any ratio `d_m/m ∈ (0, 0.589)`.
  * The reason: Lemma R's kernel moves preserve fairness, parity and every deadline, and they load `≍ S(p0)^{d}` into a single level (Proposition 9.3).
  * So box + parity + fairness do not enlarge `U_γ`, and no `γ' > γ_1` follows by this route. The loss is in the sum over `m` (cross-level cancellation), not in the Bernstein sum over `k` (which is sharp up to `1/π`).
* **NUMERICAL — HYP-9129.** Asymptotic fold-sufficient thresholds of ratio-4 super-blocks:
  * realizable handoff states: 2-power cyclotomic `1.567`; integer factors of degree `≤ 4`: `1.5685`;
  * the non-realizable real zero-measure relaxation: `≈ 1.50`;
  * the necessary-condition LP (single-block evaluation + mass + Mahler): `≈ 1.435` at `B = 4`, `≈ 1.41` at `B = 8`, and `1.593 ≈ C_*` at `B = 2` (sanity check).
* **OPEN:**
  * `C* < 3/2` (HYP-9129). The realizable constructions stall near `1.567`; the proved obstruction stops at `1.3775`.
  * The exact value of `C*`.
  * Turning H4 into a theorem.

---

## 0. The answer in brief

* **HYP-9128 is a theorem.** Super-blocks `[N, 4N)` carry the handoff state `S_N(w) = (1+w)^{N/16}(1 + w^{N/16} + … + w^{15N/16})`, which has a zero of order `N/16` at the golden point `w = −1`. They give an exactly fair extractor with deadline `⌈1.59 L⌉`. So `C* ≤ 1.59 < C_*`.
  * The constant `159/100` is limited only by a crude level-0 bottom-regime majorant (§5.2, last paragraph). The same family's measure-level threshold is `1.578` (§10).
* **The lower bound moves to `1.377`,** which is essentially the full reach of THM-4467.
* **The owner's structured majorant does not help.**
  * THM-4467 has two triangle inequalities: the Bernstein one inside `W_m`, and the implicit one in the sum over `m`.
  * The first is sharp up to `1/π`.
  * The second loses exactly the cross-level cancellation that the fairness kernel encodes. Kernel moves can bring an individual `W_m` within a constant factor of `S^{d_m}`, at any prescribed point, while leaving `F` unchanged.
  * So termwise information cannot move the domain. The structured object is the block or handoff polynomial, and that is the H2/H4 programme.
* **`3/2` stays open.** The best realizable constructions sit near `1.567`. The relaxation that reaches `3/2` needs zeros spread along an arc of `|w| = 1`, and an elementary discriminant inequality forbids that for integer handoff states (§10.3).

---

## 1. Setting and the super-block reduction

Bits are i.i.d. with `P(0) = p` and `q = 1 − p`.
* `L` is the length of the initial run.
* `T(L) ≥ L+1` is the deadline: the output is a function of the first `T(L)` bits whenever the initial run has length `L`.
* `R_L = T(L) − L − 1` is the number of free bits after the run-ending bit.
* The extractor is **complement-symmetric**: complementing the word complements the output.
* For a 0-run of length `L`, let `e_{L,z}` be the signed count (#outputs 1 − #outputs 0) over the `C(R_L, z)` continuation words with `z` zeros. Put `E_L(p) = Σ_z e_{L,z} p^z q^{R_L−z}` and `Φ(p) = Σ_L p^L q E_L(p)`.

By complement symmetry, `P(1) − P(0) = Φ(p) − Φ(q)`. So the extractor is fair iff `Φ(p) = Φ(q)` on `(0,1)`.

**Realizability.** Integers `e_{L,z}` come from an actual decision rule iff

```text
|e_{L,z}| ≤ C(R_L, z)   and   e_{L,z} ≡ C(R_L, z) (mod 2).
```

Choose `(C + e)/2` of the words to output 1. In particular `e = ±1` at `z ∈ {0, R_L}`.

**Lemma S (super-block reduction).** Let the levels `L = N + i`, `0 ≤ i ≤ n := (B−1)N − 1`, all have `T(L) ≤ BN`. Put `a_i = BN − T(N+i)`, `R_i = R_{N+i}` and `E_i(y) = Σ_z e_{N+i,z} y^z`. Then

```text
Φ_block(p) := Σ_{i} p^{N+i} q E_{N+i}(p) = p^N q^{(B−1)N} P(p/q),     P(y) = Σ_i y^i (1+y)^{a_i} E_i(y).
```

`Φ_block` is `(p ↔ q)`-invariant iff `P` is palindromic of degree `(B−2)N`. Equivalently, `P(y) = (1+y)^{(B−2)N} S(y/(1+y)²)` with `S ∈ Z[w]` and `deg S ≤ (B−2)N/2`. Then `Φ_block = w^N S(w)`, where `w = pq` and `S(0) = P(0) = e_{0,0} = ±1`.

*Proof.* Multiply each term by `(p+q)^{a_i} = 1`. Since `L + 1 + R_L + a_i = BN`, we get `p^L q · p^z q^{R−z} (p+q)^{a} = q^{BN} y^{L+z} (1+y)^{a}` with `y = p/q`.

Invariance under `p ↔ q` reads `y^N P(y) = y^{(B−1)N} P(1/y)`. Palindromic integer polynomials of degree `2D` are exactly `Σ_k s_k y^k (1+y)^{2D−2k}`, with `s_k ∈ Z` by triangularity.

Finally `q^{2D} (1+p/q)^{2D} S(pq/(p+q)²) = S(pq)`. ∎

If the run lengths `[1, ∞)` are partitioned into blocks, each with an invariant `Φ_block`, then the extractor is fair: runs are a.s. finite for `p ∈ (0,1)`.

**Lattice coordinates** (Long's): `f_{i,r} := (−1)^{i+r} e_{i,r}` and `Q(z) := P(−z)`, so that

```text
(A f)(z) := Σ_{i,r} f_{i,r} z^{i+r} (1−z)^{a_i} = Q(z),        c_{i,r} := C(R_i, r).
```

The box and parity conditions read the same for `f` and `e`.

## 2. The construction

Fix `c = 159/100`, `B = 4`, `N = 16·4^k` with `k ≥ 0`, `m = N/16`, `n = 3N − 1`, and

```text
t_i = min(4N, ⌈c(N+i)⌉),   R_i = t_i − N − i − 1,   a_i = 4N − t_i        (0 ≤ i ≤ n),
S_N(w) = (1+w)^m Σ_{j<16} w^{mj},   P_+(y) = (1+y)^{2N} S_N(y/(1+y)²),   Q = P_+(−z).
```

**Theorem 1.** For every `N = 16·4^k` there are integers `f_{i,r}` with `A f = Q`, `f ≡ c (mod 2)` and `|f_{i,r}| ≤ C(R_i, r)`. Consequently the following extractor is exactly fair:

* for `L = 1`, `T(1) = 2` (output 1 on `01`, 0 on `10`);
* for `2 ≤ L < 16`, C. D. Long's distributed separately balanced blocks `[2,4)`, `[4,8)`, `[8,16)`, re-checked exactly in F5 (class polynomial `±1`, box, parity, `T(L) ≤ 2L`);
* the super-blocks `[16·4^k, 64·4^k)` with the deadlines `t_i` above.

It satisfies `T(L) ≤ ⌈159L/100⌉` for `L ≥ 16`, so `T(L) ≤ 1.59 L + 7` for all `L`, and **`C* ≤ 159/100`**.

**Structural facts** (checked exactly in F6 for all `N = 16·4^k ≤ 2^22`; the general argument is one line each):
* **Drops.** `d_i := a_i − a_{i+1} = t_{i+1} − t_i ∈ {0,1,2}`, because `1 < c < 2`.
* **`R` below the cap.** `R_i` is nondecreasing below the cap, since `⌈cL⌉ − L` is nondecreasing.
* **The cap level `h`.** Let `h := min{i : a_i = 0}`. Then `a_i = 0` for `i ≥ h`, and `h < (4/c − 1)N + 1 ≤ n − 2`.
* **Tail levels.** `R_{n−2}, R_{n−1}, R_n = 2, 1, 0`, and `R_i ≥ 2` for `i ≤ n−2`.
* **Parity hypothesis.** `S_N ≡ (1+w^m)(1+w^m+…+w^{15m}) = 1 + w^N (mod 2)`. Hence `Q ≡ (1+z)^{2N} + z^N ≡ 1 + z^N + z^{2N}`. Also `A c ≡ Σ_i z^i (1+z)^{n−i} = (1+z)^{3N} − z^{3N} ≡ 1 + z^N + z^{2N}`, since `N` is a power of 2 and `R_i + a_i = n − i`.

## 3. Lemma R (lattice rounding for palindromic targets), re-derived

This generalises Long's boundary-preserving "five-unit" rounding (his §6) to three situations:
* arbitrary integer targets `Q` subject to a parity condition (Long has `Q = ±1`);
* drops `d_i ∈ {0,1,2}`;
* a capped tail.

The proof below is self-contained. `amm12592_procgen_20260923_hyp9128_finite.py` (`lemma_R_round`) implements it independently of Long's code.

**Lemma R.** Let `R_i + a_i = n − i`, with drops `d_i ∈ {0,1,2}`. Assume:
* **(R1)** `a_{n−2} = a_{n−1} = a_n = 0`, and `R_i ≥ 2` for `i ≤ n−2`.
* **(R2)** `Q ∈ Z[z]` with `deg Q ≤ n` and `Q ≡ A c (mod 2)`.
* **(R3)** `u ∈ R^{sites}` with `A u = Q`, satisfying:
  * `u_{0,0} = ±1`;
  * `|u_{i,0}| ≤ 2` for `i ≥ 1`;
  * `|u_{i,R_i}| ≤ 1` for `i ≤ n−2`;
  * `u = 0` on levels `n−1` and `n`;
  * at every interior site `0 < r < R_i`, `|u_{i,r}| ≤ C(R_i,r) − ε_i`, where `ε_0 = 1` and `ε_i = 1 + 2^{d_{i−1}}`.

Then there is an `f ∈ Z^{sites}` with `A f = Q`, `f ≡ c (mod 2)` and `|f_{i,r}| ≤ C(R_i, r)`. In particular `f = ±1` at all boundary sites, and `|f − u| ≤ ε_i` on level `i`.

*Proof.*
1. **A lattice point in the right coset.** Write `A c − Q = 2H` with `H ∈ Z[z]_{≤n}`. The polynomials `z^i (1−z)^{a_i}` (`0 ≤ i ≤ n`) have lowest term `z^i`, so they are a `Z`-basis of `Z[z]_{≤n}`. Write `H = Σ h_i z^i (1−z)^{a_i}` and put `q := c − 2 Σ_i h_i 1_{(i,0)}`. Then `A q = Q` and `q ≡ c`.
2. **The kernel.** For `0 ≤ i ≤ n−1` and `1 ≤ r ≤ R_i` put

   ```text
   K_{i,r} := 1_{(i,r)} − Σ_{j=0}^{d_i} (−1)^j C(d_i, j) 1_{(i+1, r−1+j)}.
   ```

   The sites are valid because `R_{i+1} = R_i + d_i − 1`. Also `A K_{i,r} = z^{i+r}[(1−z)^{a_i} − (1−z)^{d_i}(1−z)^{a_{i+1}}] = 0`.
   * The `K_{i,r}` are triangular, hence independent.
   * Their number is `Σ_{i<n} R_i = #sites − (n+1) = dim ker A`, because `A` is onto `R[z]_{≤n}` by step 1.
   * So `q − u = Σ γ_{i,r} K_{i,r}` for unique real `γ`.
   * For **any** choice `δ_{i,r} ∈ γ_{i,r} + 2Z`, the vector `f := u + Σ δ_{i,r} K_{i,r} = q + Σ (δ−γ) K` is integral, satisfies `f ≡ q ≡ c`, and satisfies `A f = Q`.
3. **Errors.** The `(i,r)` entry of `Σ δ K` is `δ_{i,r} − Σ_j (−1)^j C(d_{i−1},j) δ_{i−1, r+1−j}`. If all `|δ| ≤ 1`, then `|f − u| ≤ 1 + 2^{d_{i−1}} = ε_i`, and `|f − u| ≤ 1` on level 0.
4. **Boundary sites.**
   * Only `K_{i−1,1}` touches `(i, 0)`. So `f_{i,0} = u_{i,0} − δ_{i−1,1}` and `γ_{i−1,1} = u_{i,0} − q_{i,0}`.
   * Only `K_{i,R_i}` (coefficient `+1`) and `K_{i−1,R_{i−1}}` (coefficient `−(−1)^{d_{i−1}}`, from `j = d_{i−1}`) touch `(i, R_i)`.
   * *Lower chain*, for `1 ≤ i ≤ n−1`. Set `δ_{i−1,1} := u_{i,0} − b_i` with `b_i = sign u_{i,0} ∈ {±1}` (`sign 0 := 1`). Then `|δ| ≤ 1` and `f_{i,0} = b_i`. The coset is right: `δ − γ = q_{i,0} − b_i` is odd minus odd.
   * *Upper chain*, for `0 ≤ i ≤ n−2`. Set `η_i := δ_{i,R_i} := b_i^+ − u_{i,R_i} + (−1)^{d_{i−1}} η_{i−1}` (with `η_{−1} = 0`). Choose `b_i^+ = ±1` so that `|η_i| ≤ 1`; this is possible because `|u_{i,R_i}| + |η_{i−1}| ≤ 2`. Then `f_{i,R_i} = b_i^+`.
     * By induction, `(δ−γ)_{i,R_i} = (b_i^+ − q_{i,R_i}) + (−1)^{d_{i−1}}(δ−γ)_{i−1,R_{i−1}} ∈ 2Z`.
     * The two chains use distinct coordinates, because `R_i ≥ 2` for `i ≤ n−2`.
   * *Last two levels.* `R_{n−1} = 1`, so the single coordinate `δ_{n−1,1}` controls both of the following. Here `u = 0` there, and `d_{n−2} = 0`:
     * `f_{n,0} = −δ_{n−1,1}`;
     * `f_{n−1,1} = δ_{n−1,1} − η_{n−2}`.

     Now `γ_{n−1,1} = −q_{n,0}` is odd, and the site `(n−1,1)` gives `γ_{n−1,1} − γ_{n−2,2} = q_{n−1,1}`, which is odd. So `γ_{n−2,2}` is **even**. Since `η_{n−2} ∈ γ_{n−2,2} + 2Z` and `|η_{n−2}| ≤ 1`, the parity of the lattice forces `η_{n−2} = 0`. Take `δ_{n−1,1} := 1`, which is odd and so lies in `γ_{n−1,1} + 2Z`. Then `f_{n,0} = −1` and `f_{n−1,1} = 1`.
   * *All other coordinates:* the representative of `γ_{i,r} + 2Z` in `[−1, 1]`.
5. **Box.** Boundary sites are `±1 = ±C`. At interior sites, `|f| ≤ |u| + ε_i ≤ C(R_i, r)`. ∎

For the construction, hypothesis (R3) is supplied as follows:
* **Level 0:** `u_{0,r} = (−1)^r e_{0,r}` must satisfy the margin conditions of §4 and §5.
* **Levels `1 ≤ i ≤ h`:** it suffices that `M_i ≤ 1/2` and `R_i ≥ 10`. Indeed `|u_{i,r}| ≤ M_i C(R_i,r)` because `|K_t(r)| ≤ C(R,r)`, and `C/2 + 5 ≤ C` once `C ≥ 10`.
* **Levels `i > h`:** `u = 0` and `d_{i−1} = 0`, so `ε_i = 2 ≤ C(R_i, r)` at interior sites.

## 4. The real center: fold identities

**Fold** (Long's recursion, general target). Put `F_0(u) = ((1+u)/2)^{N−1} S_N((1−u²)/4)`, which has degree `n`, and `y = (1−u)/2`. For `i = 0, 1, …` set

```text
W_i := [F_i]_{<R_i} + (Σ_{s≥R_i} [u^s]F_i) u^{R_i},        F_{i+1} := (F_i − W_i)/y.
```

Since `F_i − W_i = Σ_{s>R_i} [u^s]F_i (u^s − u^{R_i})`, `F_{i+1}` is a polynomial of degree `n − i − 1`. The fold stops at the cap level `h`, where `W_h = F_h` and `F_{h+1} = 0`.

**Lemma F0 (center).** Define `e_{i,r}` by `E_i(x) := Σ_r e_{i,r} x^r := (1+x)^{R_i} W_i((1−x)/(1+x))`, that is, `e_{i,r} = Σ_t w_{i,t} K_t^{(R_i)}(r)`, where

```text
K_t^{(R)}(r) := [x^r](1−x)^t (1+x)^{R−t}.
```

Then `Σ_i x^i (1+x)^{a_i} E_i(x) = P_+(x)`. So `u_{i,r} := (−1)^{i+r} e_{i,r}` solves `A u = Q`.

*Proof.* Substitute `u = (1−x)/(1+x)`. Then `y = x/(1+x)`, `(1+u)/2 = 1/(1+x)` and `(1−u²)/4 = x/(1+x)²`. Telescoping gives `F_0 = Σ_{i≤h} y^i W_i`; multiply by `(1+x)^n` and use `R_i + a_i = n − i`. ∎

**Lemma F1 (closed form; exact check F1).** For `1 ≤ i ≤ h`, `F_i` is supported on `[R_{i−1}, n−i]`. For `s ≥ R_{i−1}`, and for any contour enclosing `0` and `1`, with `c_x = [u^x]F_0`:

```text
[u^s]F_i = (−2)^i Σ_{x>s} c_x C(x−s−1, i−1) = (−2)^i (2πi)^{−1} ∮ F_0(u) u^{−s−1} (u−1)^{−i} du,
T_i := Σ_{s≥R_i} [u^s]F_i = (−2)^i Σ_x c_x C(x−R_i, i) = (−2)^i (2πi)^{−1} ∮ F_0(u) u^{−R_i} (u−1)^{−i−1} du.
```

Since `R_i − R_{i−1} = d_{i−1} − 1 ≤ 1`, `W_i` has at most two nonzero coefficients: one at `R_{i−1}` (only if `d_{i−1} = 2`) and the tail `T_i` at `R_i`. Hence `M_i := ‖W_i‖_1 ≤ |[u^{R_{i−1}}]F_i| + |T_i|`.

*Proof.* Induction on `i`, using the hockey-stick identity `Σ_{x'=s+1}^{x−1} C(x−x'−1, i−1) = C(x−s−1, i)` and the monotonicity of `R` below `h`. For the contour form, expand `(u−1)^{−i} = Σ_k C(k+i−1, i−1) u^{−i−k}` on a large circle and deform. ∎

**The rate function.** For `u ∈ C`, with `a = log|(1+u)/2|`, `b = log|1+w|`, `o = log|w|` and `w = (1−u²)/4`,

```text
log|F_0(u)| ≤ N·G(u) − a + log 16,        G := a + b/16 + (15/16) max(0, o),
```

because `|Σ_{j<16} w^{mj}| ≤ 16 max(1,|w|)^{15m}`. (At the measure level, `G` is the potential of `(1/16)δ_{−1} + (15/16)·Unif(|w|=1)`.)

**Lemma K (pairing; exact check F4, R ≤ 40).** Let `j = min(t, R−t)`. Then

```text
|K_t^{(R)}(r)| ≤ [x^r](1+x²)^j (1+x)^{R−2j} ≤ e^{R H(ρ)} (1 − 2ρ(1−ρ))^j,   ρ = r/R.
```

*Proof.* `(1−x)^t (1+x)^{R−t} = (1−x²)^j (1±x)^{R−2j}`, and absolute values majorise coefficientwise. Then use `[x^r]G ≤ G(s)/s^r` at `s = ρ/(1−ρ)`, noting `(1+s²)/(1+s)² = 1 − 2ρ(1−ρ)`. ∎

Combining this with Robbins' `n! = √(2πn)(n/e)^n e^{r_n}`, `1/(12n+1) < r_n < 1/(12n)` (H. Robbins, Amer. Math. Monthly 62 (1955) 26–29; CITED, not re-read) gives

```text
C(R, r) ≥ e^{R H(ρ)}/√(2R)    whenever min(r, R−r) ≥ 2.
```

Indeed `√(n/(2πk(n−k))) ≥ √(2/(πn))` and `e^{−1/(12k)−1/(12(n−k))} ≥ e^{−1/12} > √π/2`. Hence `|K_t(r)| ≤ √(2R) g^j C(R,r)` with `g = 1 − 2ρ(1−ρ)`.

**Lemma L0 (level 0; exact checks F2, F3).** Let `K := N − 1 − R_0 = 2N − ⌈cN⌉`. Then `E_0 = 𝒜 − 𝒟` as power series, where

```text
𝒜(x) = S_N(x/(1+x)²) (1+x)^{−K},     𝒟(x) = Σ_{t>R_0} c_t [(1−x)^t (1+x)^{R_0−t} − (1−x)^{R_0}].
```

The pieces are:
* (a) `e_{0,0} = 1`.
* (b) *Bottom.* For `r < m`, `|[x^r]𝒜| ≤ C(K+2m+r−1, r)`.

  *Proof.* For `r < m` only the factor `(1+w)^m` of `S_N` matters. Coefficientwise,

  ```text
  |(1+x)^{−K}(1 + x(1+x)^{−2})^m| ≪ (1−x)^{−K}(1 + x(1−x)^{−2})^m = (1−x)^{−K−2m}(1−x+x²)^m.
  ```

  Since `1 − x + x² = (1+x³)/(1+x)`, the right side equals `(1−x)^{−K−m}((1+x³)/(1−x²))^m`. The series `(1+x³)/(1−x²) = 1 + x² + x³ + x⁴ + …` is `≪ 1/(1−x)`. Everything here has nonnegative coefficients, so the right side is `≪ (1−x)^{−K−2m}`, whose `r`-th coefficient is `C(K+2m+r−1, r)`.
* (c) *𝒟 as a contour integral.* The contour encloses `0`, `1` and `u(x)`; take `x` near `0`, and set `λ = (ζ+1)/(ζ−1)`. From

  ```text
  Σ_{t>R_0} c_t (v^t − v^{R_0}) = (2πi)^{−1} ∮ F_0(ζ) ζ^{−R_0} v^{R_0}(v−1)/((ζ−v)(ζ−1)) dζ
  ```

  (residues at `ζ = 0, v, 1`), we get

  ```text
  𝒟(x) = −2 (2πi)^{−1} ∮ F_0(ζ) ζ^{−R_0} (ζ−1)^{−2} · x(1−x)^{R_0}/(1 + λx) dζ.
  ```

  Hence, for `1 ≤ r < r_1`, `|[x^r]𝒟| ≤ ε_D C(R_0, r)`, where `ε_D := ρ_Γ max_Γ 2|F_0| |ζ|^{−R_0} |ζ−1|^{−2} q̄/(1 − |λ| q̄)`. This uses `C(R_0, r−k)/C(R_0, r) ≤ q̄^k` with `q̄ ≥ r/(R_0−r+1)`.
* (d) *Top, by reflection.* For `r' < N−1`,

  ```text
  e_{0,R_0−r'} = −[z^{r'}](1+z)^{R_0} Δ(−u(z)),     Δ(v) = Σ_{t>R_0} c_t (v^t − v^{R_0}),
  ```

  because `(1+z)^{R_0} F_0(−u(z)) = z^{N−1}(1+z)^{R_0−N+1} S_N(z/(1+z)²) = O(z^{N−1})`. With `λ' = (ζ−1)/(ζ+1)` and a contour enclosing `−1`, `0` and `1`,

  ```text
  (1+z)^{R_0}Δ(−u(z)) = 2(−1)^{R_0+1}(1−z)^{R_0} (2πi)^{−1}∮ F_0 ζ^{−R_0} (ζ²−1)^{−1} (1+λ'z)^{−1} dζ.
  ```

  Hence `|e_{0,R_0−r'}| ≤ ε_T C(R_0, r')`, with `ε_T := ρ_Γ max_Γ 2|F_0||ζ|^{−R_0}|ζ²−1|^{−1}/(1 − |λ'| q̄)`.
* (e) *Middle.* For `r_1 ≤ r ≤ R_0 − r_1`, `|e_{0,r}| ≤ √(2R_0) C(R_0, r) Σ_t |w_{0,t}| g_1^{min(t, R_0−t)}`, where:
  * `w_{0,t} = c_t` for `t < R_0`, and `w_{0,R_0} = T_0`;
  * `g_1 = 1 − 2p̄(1−p̄)`, with `p̄ = r_1/((c−1)N) ≤ ρ ≤ 1 − p̄`.

We take `r_1 = 3N/128 < m`, so that (b) applies and `r_1 < N − 1`.

## 5. The analytic certificate (all `N ≥ N_A = 4096`)

### 5.1 What is certified

For `N = 16·4^k ≥ 4096`, with `θ_1 := (2 − c + 1/8 + 3/128)/(c − 1 − 3/128 − 1/N_A) ≥ max_{r<r_1} ((K+2m+r−1)/R_0)`:

* **(I1)** `M_i ≤ 1/2` for `1 ≤ i ≤ h` (and trivially `R_i ≥ R_0 ≥ 10`).
* **(I2)** *Bottom:* `|e_{0,r}| ≤ (θ_1^r + ε_D) C(R_0,r)`, so the margin is `≥ ((c−1)N − 1)(1 − θ_1 − ε_D) ≥ 5`. This uses `C(K+2m+r−1,r)/C(R_0,r) ≤ ((K+2m+r−1)/R_0)^r`.
* **(I3)** *Top:* `ε_T ≤ 1`, so `|e_{0,R_0}| ≤ 1`, and the margin is `≥ ((c−1)N − 1)(1 − ε_T) ≥ 5`.
* **(I4)** *Middle:* `√(2R_0) Σ_t |w_{0,t}| g_1^{min(t,R_0−t)} ≤ 1/2`, so the margin is `≥ C(R_0,r)/2 ≥ 5`.

Together with §3 this is hypothesis (R3).

### 5.2 How

Each quantity is bounded by Cauchy on an explicit circle `Γ = {|u − x_0| = ρ}` through Lemmas F1 and L0 and the rate function:

* **Levels.** `v = i/N` ranges over `[0, 4/c − 1 + 1/N_A]` in 152 cells of width `0.01`. The exponent is

  ```text
  Ψ = G − σ log|u| − v log|(u−1)/2|,    σ ∈ (c−1)(1+v) + [−2/N_A, 1/N_A].
  ```

  This interval covers `R_{i−1}/N` and `R_i/N` at every level, including the capped level `h`, for every `N ≥ N_A`. The `O(1)` factor is `E = max(−a − log|u|, −a − log|u−1|) + log 16`. The bound is `M_i ≤ 2ρ exp max_Γ(NΨ + E)`.
* **Level 0.**
  * The bottom and top circles use the exponent `G − σ_0 log|u|` with `σ_0 = R_0/N ∈ [c − 1 − 1/N_A, c − 1]`. The bottom circle encloses `0` and `1`; the top circle encloses `−1`, `0` and `1`. Each keeps a gap `≥ 1/4` from `±1`, so that `|λ|q̄ < 1`.
  * The middle uses 118 cells of width `0.005` in `τ = t/N`. There, `|c_t| ≤ ρ exp max_Γ(N(G − τ log|u|) − a − log|u| + log 16)`, with the pairing factor `g_1^{N·min(τ, σ_0 − τ)}` and `≤ N·dτ + 1` integers per cell. The tail `T_0` is handled like a level-0 tail.
* **Certification of each circle.** The arc `θ ∈ [0, π]` is certified (the integrand is conjugation-symmetric) by adaptive bisection with `mpmath.iv` intervals. On every final piece, `Ψ.hi ≤ Ψ_float + 0.004 < 0`.
* **Choice of circles.** They were chosen by a float search and then frozen. Only the interval bounds matter.

**Where `159/100` comes from.** The bottom-regime majorant needs `θ_1 < 1`, i.e. `c > 1.586`. Every analytic rate is comfortably negative already at `c = 1.59`, and the family's measure-level threshold is `1.578` (§10). A saddle-point bound for `[x^r]𝒜` would lower the constant; this is not needed for `C* < C_*`.

### 5.3 Numbers (Part 1 of the `.out`) and uniformity in `N`

* **Levels (I1).** All 152 cells are certified. The worst certified rate is `−0.01464`, in the cell `v ∈ [0.07, 0.08]` with circle `(0.6607, 1.5633)`. Hence `max_i M_i ≤ e^{−56.75} = 2.3·10^{−25}` at `N = N_A`.
* **Level-0 constants.** `q̄ = 0.04137`, `p̄ = 0.03972`, `g_1 = 0.923707`, `θ_1 = 0.98608`.
* **Bottom (I2).** On the circle `(0.4013, 1.2387)` the rate is `≤ −0.03446`, so `ε_D ≤ e^{−140.28}`. The margin is `≥ 33.6`.
* **Top (I3).** On the circle `(0.0457, 1.5633)` the rate is `≤ −0.03451`, so `ε_T ≤ e^{−138.28} = 8.8·10^{−61}`. The margin is `≥ 2415`.
* **Middle (I4).** `√(2R_0)·Σ ≤ e^{−9.57} = 7.0·10^{−5}`. The worst cell exponent (rate + pairing) is `−0.00499`.
* **Uniformity.** Every certified bound has the form `(value at N_A) × poly(N) × exp((N − N_A)·r)`, with a certified `r < 0` for the pieces. `θ_1` is computed with `1/N_A`.
  * For the levels, `NΨ + E ≤ (N_AΨ + E) + (N − N_A)Ψ`.
  * In the middle regime, `√(2(c−1)N)(N·dτ + 1)e^{(N−N_A)·expo}` is decreasing for `N ≥ 1.5/|expo| = 301`.
  * So (I1)–(I4) hold for every `N ≥ 4096`. ∎

## 6. Finite certificates (`16 ≤ N ≤ 2048`, Part 2 of the `.out`)

For `N = 32, …, 2048` the exact fold gives a margin certificate: `max_{i≥1} M_i ≤ 1/2`, and at level 0 an exact Krawtchouk packet with margin `≥ 5`, `e_{0,0} = 1` and `|e_{0,R_0}| ≤ 1`. That is (R3), so Lemma R applies. For `N ≤ 256` the integral block is also constructed explicitly and verified.

| N | used? | max T(L)/L on [N,4N) | h | max_{i≥1} M_i | level-0 min margin | max_r \|e_{0,r}\|/C | \|e_{0,R_0}\| | margin certificate | explicit Lemma-R rounding (sites; max interior \|f−u\|) |
|---|---|---|---|---|---|---|---|---|---|
| 16 | yes | 28/17 = 1.6471 | 24 | 0.132 | 4.87 (r=1) | 0.459 | 1.2e-2 | FAIL (4.87 < 5) | VERIFIED (702; 4.36) |
| 32 | – | 55/34 = 1.6176 | 48 | 0.0928 | 7.09 | 0.606 | 6.6e-2 | PASS | VERIFIED (2 772; 4.08) |
| 64 | yes | 109/68 = 1.6029 | 97 | 0.0807 | 14.96 | 0.596 | 2.0e-3 | PASS | VERIFIED (11 018; 4.39) |
| 128 | – | 222/139 = 1.5971 | 194 | 0.0176 | 31.00 | 0.587 | 1.2e-3 | PASS | VERIFIED (43 938; 4.65) |
| 256 | yes | 51/32 = 1.5938 | 388 | 1.5e-3 | 63.00 | 0.583 | 4.9e-6 | PASS | VERIFIED (175 482; 4.77) |
| 512 | – | 823/517 = 1.5919 | 776 | 7.7e-6 | 125.00 | 0.586 | 1.2e-10 | PASS | — |
| 1024 | yes | 1653/1039 = 1.5910 | 1552 | 3.8e-10 | 249.00 | 0.588 | 2.6e-19 | PASS | — |
| 2048 | – | 3305/2078 = 1.5905 | 3104 | 1.0e-18 | 497.00 | 0.589 | 9.3e-37 | PASS | — |

The column "max T(L)/L" is `max_{N≤L<4N} ⌈1.59L⌉/L`, which tends to `1.59`. All level-0 minimum margins occur at `r = 1`.

`N = 16` fails the margin test narrowly: its level-0 margin is `4.87 < 5` at `r = 1`. It is covered by the explicit rounding, which is verified exactly. The needed sizes `N = 16, 64, 256, 1024` are therefore all certified. The unneeded sizes `32, 128, 512, 2048` are shown as a consistency check.

## 7. Proof of Theorem 1 (assembly)

* **Large sizes.** For `N = 16·4^k ≥ 4096`, §5 gives (R3) for the real center of Lemma F0. §2 gives (R1) and (R2). Lemma R gives integral packets, and Lemma S together with realizability turns them into a balanced super-block `[N, 4N)` with `Φ_block = w^N S_N(w)`.
* **Small sizes.** For `N = 16, 64, 256, 1024`, §6 gives the same.
* **Prefix.** `L = 1` and Long's three blocks are balanced (F5).
* **Fairness.** The blocks tile `[1, ∞)`, so the extractor is fair.
* **Deadlines.** `T(L) = t_i ≤ ⌈cL⌉` on every super-block, and `T(L) ≤ 2N ≤ 2L` on the prefix blocks. ∎

**Remark.** The margins are large. For `N ≥ 4096` the level masses are `≤ 10^{−25}`, and the level-0 bottom margin is at least `33 ≫ 5`. The constant is limited by (I2) alone. At block level the explicit integer super-blocks of the first note (`157/100` at `N = 256`) and the fold tables suggest `≈ 1.57` for this family.

## 8. THM-4467 at `γ = 377/1000`: `C* ≥ 1.377`

**Method.** This is the method of THM-4467 §3.4 (subordination), with one numerical repair: the inverse Riemann map of `W_γ` is computed by a *damped* Newton iteration that never leaves `W_γ`. The undamped version converged to spurious roots beyond the reentrant corner.

**The certificate** is `amm12592_procgen_20260923_hyp9128_gamma377_certificate.json` (SHA-256 `4fb1a3b3…2e9e58`).
* It is `P(z) = Σ_{k≤32767} b_k z^k` with real float64 coefficients, `P(0) = 0` and `b_1 = 1.00016088 > 1`.
* `Σ k|b_k| = 1.9522`.

**The check** (`--verify-only`, Part 4):
* `P(e^{it}) ∈ W_{0.377}` is checked on two interleaved grids of `2^22` nodes each, for the upper half circle (real coefficients).
* The explicit perturbation allowance is `δ = 7.32·10^{−7}` (node spacing plus the evaluation allowance `10^{−9}`), with `min|1−2p| = 1.3178`.
* The maximum on the grid is `0.9999882`, and the maximum certified bound is `0.9999902 < 1`.

**Consequence.** Schwarz's lemma gives `R(W_{0.377}, 0) ≥ b_1 > 1`, so `Λ(0.377) > 0`. THM-4467 then gives `C* ≥ 1.377`.

**Rigor caveat.** This is the same caveat as THM-4467's `3/8` certificate: a floating-point FFT evaluation, with an allowance `10^{−9}` that dominates the a-priori FFT error.

The remaining gap to the method's reach is `γ_1 − 0.377 = 5.2·10^{−4}`. Closing it needs a certificate radius closer to 1, i.e. more coefficients; that is not worth doing.

## 9. Structured majorant (the owner's question)

**The step in question.** THM-4467 continues `F` from the termwise bound `|p^m q W_m(p)| ≤ |p|^m |q| S(p)^{d_m}`. This involves two triangle inequalities:
* **(i)** inside `W_m`, over the Bernstein index `k`;
* **(ii)** across `m`, implicitly: the series is bounded term by term.

In the owner's words, the gap in each is "calculable from the values averaged". We compute it.

**9.1 The sum over `k` is sharp up to `1/π` (box).**
* For any `z_k ∈ C`, choosing `w_k ∈ {0, C(d,k)}` by the sign of `Re(e^{−iφ}z_k)` gives `|Σ_{chosen} z_k| ≥ Σ_k max(0, Re(e^{−iφ}z_k))`.
* The average over `φ` of the right side is `Σ|z_k|/π`.
* Hence `max_box |W(p)| ≥ S(p)^d/π`.

Part 5 [1] gives ratios `0.5000` at real `p < 0` (terms of alternating sign) and `0.3186–0.3216` at three complex points. So there is no exponential loss here.

**9.2 Parity forces only the Lucas skeleton.**
* `E_m = 2W_m − 1` has Bernstein coefficients `2w_k − C(d,k) ≡ C(d,k) (mod 2)`. These are odd exactly at the Lucas submasks `k ⊆ d`. This is Lemma P's `φ ≡ Σ w^{2^j}` seen blockwise, the Artin–Schreier input.
* The largest the forced part can be is `Σ_{k⊆d}|p|^{d−k}|q|^k = Π_{j∈bits(d)}(|p|^{2^j} + |q|^{2^j})`.
* Relative to `S^d`, this has rate `−0.32` at `p = −0.62` and `−0.54` at `p = 0.3+0.6i` (Part 5 [2]).
* Parity therefore constrains an exponentially thin skeleton and leaves the bulk coefficients free. It yields no upper-bound improvement.

**9.3 Fairness does not either. The cross-level loss is the fairness kernel.**
* The fairness constraints of a block are the linear equations `A f = Q`. Their kernel is spanned by the `K_{i,r}` of Lemma R.
* A move `f → f + Σ_r m_r K_{i,r}` with even `m_r` touches level `i` and level `i+1` only. It preserves:
  * the class polynomial, hence `Φ_block = w^N S_N(w)`, hence `F`;
  * parity;
  * every deadline.
* It changes the individual `E_i(p0)`.

**Proposition 9.3 (kernel saturation).**
* **Setting.** Take a super-block of Theorem 1. Let `h < i ≤ n−2` be a capped level with `R := R_i ≥ 3`, so `d_i = 0` and `|f| ≤ 2` at interior sites. Fix `0 < α ≤ 1/3` and `p0 ∈ C \ [0,1]`.
* **The move.** Put `m_r := ±2⌊α C(R−1, r−1)/2⌋` for `2 ≤ r ≤ R−1`. Choose the signs to align the terms `τ_r := (−1)^{i+r} p0^r (1−p0)^{R−r}` with a phase `φ` maximising `Σ_r |m_r| |Re(e^{−iφ}τ_r)|`. Then choose a global sign.
* **Conclusion.** The result is a valid fair super-block with the same deadlines, and

  ```text
  |E_i'(p0)| ≥ (2/π) Σ_r |m_r||τ_r| = (2α/π) (|p0|/S(p0)) S(p0)^{R} (1 − o(1))    (R → ∞).
  ```

  For real `p0` the factor `2/π` is `1`.

*Proof.*
* **Box.** At level `i`, `|f'_{i,r}| ≤ 2 + αC(R−1,r−1) ≤ C(R−1,r−1) ≤ C(R,r)`, because `m_r ≠ 0` forces `C(R−1,r−1) ≥ 6`. At level `i+1`, whose sites `r−1 ∈ [1, R−2]` are interior, `|f'| ≤ 2 + αC(R−1, r−1) ≤ C(R_{i+1}, r−1)`.
* **Boundary and parity.** Boundary sites are untouched (`m_1 = m_R = 0`), and `m_r` is even.
* **Size.** `max(|E+X|, |E−X|) ≥ |X|`. The `φ`-average of `|Re(e^{−iφ}τ)|` is `2|τ|/π`. The floors cost `O((R+1) max(|p0|,|q0|)^R) = o(S^R)`. ∎

**Exact examples** (Part 5 [4], every modified block re-verified exactly). For `α = 1/64` and `N = 64, 128`, the ratio `|E_i(p0)|/S(p0)^{d}` is:

| level | d/m | p0 (on ∂Ω_0(γ)) | before | after |
|---|---|---|---|---|
| capped, `N=64` (`m=185`, `d=70`) | 0.378 | `−0.7152` (`γ=0.3775`) | `3.1·10^{−11}` | `4.60·10^{−3}` |
| capped, `N=128` (`m=371`, `d=140`) | 0.377 | `−0.7152` | `6.7·10^{−22}` | `4.60·10^{−3}` |
| capped, `N=64` | 0.378 | `−0.3684+0.6380i` | `8.8·10^{−13}` | `3.31·10^{−3}` |
| capped, `N=128` | 0.377 | `−0.3684+0.6380i` | `6.6·10^{−25}` | `3.26·10^{−3}` |
| uncapped, `N=64` (`m=112`, `d=66`) | 0.589 | `−0.6210` (`γ=0.59`) | `4.5·10^{−3}` | `2.01·10^{−2}` |
| uncapped, `N=128` (`m=225`, `d=132`) | 0.587 | `−0.6210` | `5.0·10^{−3}` | `2.06·10^{−2}` |

The "after" column is independent of `N`, so the rate `(1/d)log(|E|/S^d)` is `O(1/d) → 0`. At `p0 = −0.7152` it equals `α|p0|/S(p0) = 4.598·10^{−3}` to the digit shown. Before the move, the capped levels sat at rate `≈ −0.35`.

**Corollary 9.4 (no termwise improvement).**
* **(a)** Let `Φ(d, p)` be any function with `|W_m(p)| ≤ Φ(d_m, p)` for every level of every exactly fair extractor. It suffices that this holds for the kernel modifications of Theorem 1's extractors. Then for every `p0 ∉ [0,1]` and all `d ≥ d_0(p0)`,

  ```text
  Φ(d, p0) ≥ (α/2π) (|p0|/S(p0)) S(p0)^d.
  ```

  *Why.* Every `d ≥ 3` occurs as `R_i` at a capped level of the super-block `[N, 4N)`, where `N = 16·4^k` and `1.48N ≥ d`. Also `|W_m| ≥ (|E_m| − 1)/2`, since `W_m = (E_m + 1)/2`.
* **(b)** The same holds for any `Φ(m, d, p)` along the realized pairs `(m, d) = (N+i, 3N−1−i)`. Their ratios `d/m` fill `(0, 0.589)` with spacing `O(1/N)`.
* **(c)** Consequently, for a hypothetical extractor with `d_m = ⌊γm⌋`, the termwise series `Σ|p|^m Φ(d_m, p)` diverges wherever `|p| S(p)^γ > 1`. So the convergence domain obtainable from `Φ` is contained in the closure of `Ω_0(γ)`.

So no termwise majorant built from box, parity, Bernstein positivity and fairness can enlarge `Ω_0(γ)`, and hence `U_γ` and `W_γ`. **γ_1 = 0.3775246 is the reach of the continuation-plus-Pólya method with any termwise input**, and no `γ' > 3/8` is certified by this route. The owner asked: "which extractors attain the triangle bound near `∂W_γ`?" They are the kernel-saturated super-block extractors above, with one move per super-block at the capped level where `d_m/m ≈ γ`, phase-aligned at `p0 ∈ ∂U_γ`. They are fair with `T ≤ ⌈1.59L⌉`, so they do not contradict THM-4467: their *global* profile is `0.59`, but their *individual* terms saturate the bound at `d_m/m ≈ γ`.

**9.5 Where the lever actually is.**
* The triangle inequality across `m` discards exactly the kernel cancellation between adjacent levels.
* The invariant object is the block polynomial (`Φ_block = w^N S_N(w)`), not its terms.
* For block-structured extractors, bounding the block polynomial is Theorem B (`S_N = ±1` gives the golden bound) and the H2/H4 potential theory of §10.
* For a general fair extractor, the analogous invariant is the antisymmetrised partial sum `D_M(p) = Φ_{<M}(p) − Φ_{<M}(q)`, a polynomial of degree `≤ max_{L<M} T(L)`.
  * By fairness, `D_M(p) = Φ_{≥M}(q) − Φ_{≥M}(p)`, the difference of the two tails, which are `O(q^M)` and `O(p^M)` respectively.
  * It is unchanged by any fairness-preserving move inside `[1, M)`.

A domain larger than `U_γ` for all extractors would have to come from a bound on these handoff polynomials, not on the `W_m` (hypothesis candidate H8).

**Defect of the unmodified extremal objects** (Part 5 [3]). The rate is `max_i (1/R_i) log(|E_i(p)|/S^{R_i})`:
* Long's blocks: `−0.070` (`N=32`) and `−0.038` (`N=64`) at the golden point.
* Our super-blocks: `−0.073` (`N=32`) and `−0.055` (`N=64`).
* Medians are `−0.09` to `−0.34`.

The constructions themselves are exponentially below the triangle bound, but that is a property of the constructions, not of the class.

## 10. HYP-9129 numerics (NUMERICAL)

**10.1 Model** (`…_hyp9128_h2.py`).
* A ratio-`B` handoff state `S_N` has a normalised zero measure `μ` of mass `≤ M = (B−2)/2`.
* `S_N(0) = ±1` and an integral leading coefficient force the Mahler constraint `∫log|ζ|dμ ≤ 0`.
* `(1/N)log|S_N| → U_μ(w) = ∫log|1 − w/ζ|dμ`.
* The **fold-sufficient** threshold is the least `c` for which the inequalities (I1)–(I4) hold at exponential order, on circles, with `N·G` replaced by `N(a + U_μ)`. Two extra conditions are added: the `r = 1` condition `|∫Re(−1/ζ)dμ − (2−c)| ≤ 0.95(c−1)`, and the bottom regime via a Cauchy bound in `z`.
* **H2** minimises this threshold over `μ` by alternating LPs: fix the extremal circles, solve an LP in the dictionary weights, re-select.
* **H4** is the necessary side. It imposes the single-block evaluation inequality `U_μ(w) ≤ min_{x ↦ w} [max_v(v log|x| + α(v) log|1+x| + ρ(v) log(1+|x|)) − (B−2)log|1+x|]` at test points, together with mass and Mahler. (This is Long's Thm 11.1 generalised to palindromic `P_+ = S(w)(1+x)^{(B−2)N}`.)

**10.2 Results** (Part 6 of the `.out`)

| quantity | B | threshold | comment |
|---|---|---|---|
| golden-zero family, `κ = 0` (uniform circle only) | 4 | 1.6154 | no golden zero |
| golden-zero family, `κ = 1/32` | 4 | 1.590 | |
| golden-zero family, `κ = 1/16` (the proved family) | 4 | 1.5777 | rigorous certificate: `1.59` (§5) |
| golden-zero family, `κ = 0.09` | 4 | 1.5738 | best `κ` on this grid |
| golden-zero family, `κ = 1/8` | 4 | 1.5797 | |
| golden-zero family, best `κ` (`0.09`) | 8 | 1.5988 | |
| H2, 2-power cyclotomic states (realizable) | 4 | **1.5668** | weights `0.300 (1+w^16)`, `0.214 (1+w^8)`, `0.178 (1+w^32)`, `0.145 (1+w^4)`, `0.104 (1+w)`, `0.058 (1+w²)` |
| H2, integer factors of degree ≤ 4 plus uniform (realizable) | 4 | 1.5685 | `0.556` uniform, `0.243 (1+w)²(1+w²)`, `0.175 (1+w⁴)`, `0.020 (1+w²)` |
| H2, real atoms (relaxation; **not realizable**) | 4 | 1.5017 | Mahler active; mass on the arc `\|ζ\| = 1`, `arg ζ ∈ [1.4, 2.9]` |
| H2, real atoms (relaxation) | 8 | 1.5397 | `0.427` uniform plus arc atoms |
| H4 necessary LP | 2 | 1.5933 | Long's golden bound, up to grid error (`C_* = 1.5980`) |
| H4 necessary LP | 4 | 1.4352 | heaviest zeros at `−0.23+0.86i` and along the arc near `−1` |
| H4 necessary LP | 8 | 1.4135 | |

An independent scratch discretization of the relaxation (full `161`-point grid and a log-spaced atom grid) gave `1.4996` (`B = 4`) and `1.553` (`B = 8`). It is not part of the `.out`.

**10.3 Reading**

* **The explicit family.** Its measure-level threshold is `1.578` (`κ = 1/16`). The optimum over `κ` is `≈ 1.574` near `κ ≈ 0.09`, which matches the first note's contour asymptotics (`1.570`, deformed contours) up to the circles-only loss. The rigorous `1.59` of §5 is consistent.
* **Realizable versus relaxed.** Realizable states reach `1.567`. The best cyclotomic mixture puts its weight on `1 + w^{16}`, `1 + w^8`, `1 + w^{32}`, `1 + w^4`, `1 + w` and `1 + w²`. The integer-factor LP picks `(1+w)²(1+w²)` (the exact `N = 4` CP-SAT optimum of the first note), `1 + w⁴`, and the uniform unit circle. The unconstrained real measure reaches `≈ 1.50` at `B = 4`. It does so by spreading atoms along the arc `|ζ| = 1`, `arg ζ ∈ [1.4, 2.9]`, with the Mahler constraint active.
* **Why the relaxation optimum cannot be realized (argument; the ingredients are rigorous, the limit passage is sketched).** Write `S_N = ±Π_j P_j^{e_j}` with `P_j ∈ Z[w]` irreducible, `P_j(0) = ±1` and `D_j = deg P_j`. For each factor, `disc(P_j)` is a nonzero integer and `|lead(P_j)| = 1/Π|ζ|`. So its root-counting probability measure `ν_j` satisfies

  ```text
  Σ_{k≠l} log|ζ_k − ζ_l| ≥ 2(D_j − 1) Σ_k log|ζ_k|,
  i.e.  ∫∫ log|z−w| dν_j dν_j ≥ 2∫ log|ζ| dν_j − O(1/D_j)    (off the diagonal).
  ```

  Mahler gives `∫ log|ζ| dν_j ≤ 0` for each factor. Now suppose the limit measure `μ = lim Σ_j (e_j D_j/N) ν_j` is carried by `|ζ| = 1`, with `∫ log|ζ| dμ = 0`. This is the configuration of the atoms-LP optimum: Mahler active, mass on the unit circle.
  * The factors carrying weight then have their roots asymptotically on the circle.
  * **Factors with `D_j → ∞`.** The inequality forces `∫∫ log|z−w| dν_j dν_j ≥ −o(1)`. On the circle this forces `ν_j → uniform`: the uniform measure is the unique maximiser, with value `0 = log cap`. An arc, which has capacity `< 1`, is excluded.
  * **Factors of bounded degree.** A factor carrying positive weight produces atoms at all of its roots, so all of them lie on the circle. Its reciprocal roots are then algebraic integers with every conjugate on the unit circle, so by Kronecker's theorem (classical) they are roots of unity.
  * **Parity.** Every super-block packet vector has `f ≡ c (mod 2)`. So `Q ≡ A c`, and hence `S_N ≡ 1 + w^N = (1+w)^N (mod 2)`. By unique factorisation in `F_2[w]`, every factor is a power of `1 + w` mod 2 (up to degree drop). Among cyclotomics this leaves only `Φ_{2^k} = 1 + w^{2^{k−1}}`. For odd `m > 1`, `Φ_m(1)` is odd, so `Φ_m` is coprime to `1 + w` over `F_2`; and `Φ_{2^k m} ≡ Φ_m^{2^{k−1}} (mod 2)`.

  So on the circle the realizable measures are exactly the `cyclo` mixtures (optimum `1.567`). The arc measure found by the atoms LP (`≈ 1.50`) is excluded. This is the quantitative form of the classical Fekete / Fekete–Szegő rigidity (Math. Z. 63 (1955) 158–172; CITED, not re-read; only the discriminant inequality above is used). Off the circle, where the Mahler constraint has slack, realizability is governed by the same per-factor energy inequality. Whether such measures can approach `1.50` is H6.
* **`B = 8`.** The relaxation is worse at `B = 8` (`1.54`). The circles-only sufficient criterion degrades on the longer level range, even though the necessary thresholds decrease with `B`.
* **H4.** The necessary-condition LP stays feasible down to `≈ 1.435` (`B = 4`) and `≈ 1.41` (`B = 8`). At `B = 2` it returns `1.593 ≈ C_*`, reproducing Long's golden bound; this is the sanity check. So evaluation inequality + mass + Mahler cannot prove `C ≥ 3/2`, even for ratio-4 super-blocks. Any proof of `C* ≥ 3/2` must use the arithmetic of realizable states (the per-factor discriminant, Kronecker and parity rigidity above) or the full fold/rounding structure.
* **Rigor.** The LP places atoms on a grid and constraints at test points. Zeros sitting exactly on test points are not excluded, so this is a heuristic and not a certificate.

## 11. Hypothesis candidates (not filed; for the orchestrator)

* **H5 (mechanical sharpening).** Replace the bottom-regime majorant (b) by a saddle-point bound. The same interval certificates should then reach `c = 1.58` for `κ = 1/16`, and about `1.575` for `κ = 0.09`. This is routine.
* **H6 (realizable handoff frontier, HYP-9129 refinement).** Conjecture: the infimum over *realizable* ratio-`B` handoff measures of the fold-sufficient threshold is `≥ 1.56` for every `B`, so fold-plus-rounding super-blocks cannot reach `3/2`.
  * By §10.3 (per-factor discriminant inequality, Kronecker, parity), the realizable measures carried by `|w| = 1` are exactly the `cyclo` mixtures.
  * Off the circle, each irreducible factor's root measure satisfies `∫∫ log|z−w| dν_j dν_j ≥ 2∫ log|ζ| dν_j`.
  * A convex program over mixtures of such `ν_j` would decide H6 at the measure level.
* **H7 (H4 rigidity).** A lower bound `C* ≥ 3/2` for super-block extractors would need an arithmetic input beyond the evaluation inequality. The candidate is the discriminant-energy inequality plus Mahler on the parts of `μ` near the golden arc. The necessary LP without it is feasible at `1.435`, and its optimal measures are arc-concentrated just like the non-realizable relaxation.
* **H8 (handoff majorant).** Is there an extractor-independent bound on the antisymmetrised partial sums `D_M(p)` (§9.5) that beats the termwise sum on a region beyond `U_γ`? By Corollary 9.4, this is the only form in which the owner's "structured majorant" can still move `γ_1`.

## 12. Files and reproduction

* Scripts, all in `04-computation/experiments/`:
  * `amm12592_procgen_20260923_hyp9128_contours.py` — §5 (interval certificates, `--NA 4096`).
  * `amm12592_procgen_20260923_hyp9128_finite.py` — §6 margins and the independent Lemma R rounding (`--nmin 16 --nmax 2048 --round-max 256`).
  * `amm12592_procgen_20260923_hyp9128_lemmas.py` — exact checks F1–F6.
  * `amm12592_procgen_20260923_hyp9128_gamma377.py` — §8 build and `--verify-only`.
  * `amm12592_procgen_20260923_hyp9128_majorant.py` — §9.
  * `amm12592_procgen_20260923_hyp9128_h2.py` — §10.
  * `amm12592_procgen_20260923_hyp9128_run.sh` — runs everything, about 12 min, one process at a time under an 880 MB watchdog.
* Output: `05-knowledge/results/amm12592_procgen_20260923_hyp9128.out`. It records script SHA-256 digests and the peak RSS of every step.
* Certificate: `05-knowledge/results/amm12592_procgen_20260923_hyp9128_gamma377_certificate.json`.
* External input: C. D. Long's distributed certificate `finite_blocks.json` (SHA-256 `1be015fd…84cb2a9`, as in the first note). It is used for the three prefix blocks (F5) and the Long-block defect table. It is fetched with the generic user agent if absent.

SHA-256 of the deliverables (raw LF bytes; the `.out` header records the script digests of the run):

```text
8ef7d42def100d0d4b95983a7c8864a96706fa95747fea57cc3ad5a5ba8e744c  amm12592_procgen_20260923_hyp9128_contours.py
206f5f7e3ed90ad307c6fc65935a13cf4ba3e5118de20dc8e908ef5fb12a439f  amm12592_procgen_20260923_hyp9128_finite.py
847ec0d74db26531d60a20a7eb2df6bc0f702e24ad174d1591a6e0530ad8b020  amm12592_procgen_20260923_hyp9128_gamma377.py
994ade6671af656591a7a3257f2293826c6c291cd78a4e894f4fc4145dacede1  amm12592_procgen_20260923_hyp9128_h2.py
c29070c39ee45499f6523e25e2db8dbff8e5748c9488077dbf291c5964488eac  amm12592_procgen_20260923_hyp9128_lemmas.py
b176afa27b189439c472de7e7fc0191f5be538b63f25f5b62f5e1c84c2e4590c  amm12592_procgen_20260923_hyp9128_majorant.py
f5403b08c87ee62fb3cc8e03fd093a15b38c2886ff5dcc68bdb2350d44f2bd30  amm12592_procgen_20260923_hyp9128_run.sh
4fb1a3b3ab4cab79ff8bc514ecdc9aecc2c57de96de063a61a152653e32e9e58  amm12592_procgen_20260923_hyp9128_gamma377_certificate.json
8edace740e77533687d28004eae776473d437d3612f8fb671c6653892091afde  amm12592_procgen_20260923_hyp9128.out
```

The run of 2026-09-23T20:37Z took about 12 minutes. The peak RSS per step was at most 805 MB (Part 6e), with all steps exiting 0.
