# Where 21/20 and 1/700 go: an exponent ledger for the ordinary two-point Chowla paper, and what it touches in the corpus

2026-10-14, mac-mini-2026-10-14-wilson. Audited by audit D (corrections recorded as MISTAKE-612); audit D read pp. 1–85 of the PDF.

**Source.** "Ordinary two-point correlations of multiplicative functions" (OpenAI, 24 September 2026; 87 pp.; supplied by the owner as `final.pdf`, trusted per the owner's instruction, not re-verified here).

**Its Theorem 1.1.** There is an absolute `c > 0` with

`|Σ_(n<=X) λ(a₁n + b₁) λ(a₂n + b₂)| <= C X/(log X)^c`

for every fixed nonproportional pair of affine forms, i.e. ordinary two-point Chowla with a power-of-log saving. Theorem 1.2 is the binary corrected Elliott conjecture for ordinary averages.

The owner asked how the numbers **21/20** (Lemma 4.4: `≪ L^(−21/20)`) and **700** (Theorem 4.3 = MRT15 Theorem 1.3: `(log Z)^(−1/700)`) are used, and how they relate to the corpus's past work.

**Typing.**
- ARITHMETIC: exact bookkeeping of the paper's displayed inequalities.
- STRUCTURAL: the same mechanism appears in both places.
- ANALOGY.
- COINCIDENCE.

## 1. The ledger (ARITHMETIC)

**Parameters (§2, eq. (2.1)).**
- `L = (log X)^(1/A)` with `A >= 1000`.
- `δ = 1/200`, so `1 − δ = 0.995`.
- `J = ⌊δ log L/(6W)⌋`, `η = e^(−J)`, `K = e^(4J)`.
- Prime supplies `P_i ⊂ {p ≡ 1 (mod 5)}` with `log p ∈ [Y_(i−1), Y_i)`, and padding primes `Q = {p ≢ 1 (mod 5), p <= e^L}`.

**Rough-number Fourier bounds (Lemma 4.2).** Use integers in `[D, 2D)` free of primes up to `y = exp(L^0.99)`:
- one form gives `‖B‖_∞ ≪ L^(−0.99)`;
- the four forms `x₁, x₂, x₃, x₁ + x₂ − x₃` give `∫|B|⁴ ≪ D^(−1) L^(−3.96)`, with `3.96 = 4·0.99`.

**Where 700 goes (Lemma 4.4, proof of (4.6)).**
- Theorem 4.3 is applied with `H = D >= ½ exp(L^0.995)` and `log Z >= L^A/2`.
- The first term is `≪ L^(−0.995) log L`.
- The second is `(log Z)^(−1/700) <= 2^(1/700) L^(−A/700)`.
- The lemma needs both `≪ L^(−0.8)`. That requires only **`A >= 0.8·700 = 560`**.
- The paper takes `A >= 1000`, giving `L^(−10/7)`. The 1000 is the paper's own floor (§2, §8.4, Lemma 4.4's hypothesis). Lemma 3.1 itself needs `A > C₁ + 20`, and `D/Y <= exp(2L − L^A/2)` needs only `A > 1`.
- So **the MRT exponent is not binding**. Any short-interval exponent `1/θ` with `θ <= 1250` would serve unchanged.
- The final exponent is `c = δ/(6WA)` (§8.4). It sees 700 only through `A`, and only if `0.8·700` exceeded `max(1000, C₁ + 20)`, which it does not.
- With `W >= 10`: `c <= 1/(1.2·10⁷)`, and smaller still, since `W` must make `e⁵C₃/√W <= e^(−1)`.

**Where 21/20 goes (Lemma 4.4, after (4.7)).** Split the `θ`-integral at `|B| = L^(−τ)`.
- **Small `|B|`** (Parseval and Cauchy–Schwarz): `≪ L^(−τ)`.
- **Large `|B|`.** Chebyshev gives `|E| <= L^(4τ) ∫|B|⁴ ≪ D^(−1) L^(4τ − 3.96)`. Combined with `‖B‖_∞ ≪ L^(−0.99)` and `Σ_v |F_v| ≪ Y D L^(−0.8)`, this part is `≪ L^(4τ − 5.75)`.
- With `τ = 21/20` these are `L^(−21/20)` and `L^(−31/20) = L^(−1.55)`, the paper's displayed values.
- The balanced choice would be **`τ* = 23/20`**: `4τ − 5.75 = −τ`, giving `L^(−23/20)` for both. So 21/20 is a convenient, not an optimal, cut.

**Why the cut is `1 + 1/20`.**
- There are `O(L/η)` short multiplicative bins.
- Lemma 4.1 sums the per-bin error `L^(−21/20)` over them, leaving `η^(−1) 2^J L^(−1/20)`.
- In §8.4, `e^J 2^J L^(−1/20) <= e^(−J)` iff `(2 + log 2)J <= (1/20) log L`. That is guaranteed by `J <= δ log L/(6W)` together with the paper's condition `δ(2 + log 2)/(6W) < 1/20`.
- At `W = 10`, the left side is `0.000224`. **The 1/20 margin is about 220 times what is needed.** The lemma's exponent is `min(τ, 5.75 − 4τ)`, so any `τ ∈ (1.000224, 1.1874)` would do.

**The other "just above one" exponent: 6/5 in (4.13).** `φ(t) ≪ (1 + L|t|)^(−6/5)` with `6/5 = (8/5)·(3/4)`.
- The factor `8/5` comes from the two tilted padding divisors `q₁, q₂`: each available prime is selected with probability `4/5`, so a prime `r` contributes `+log r` or `−log r` to `log(q₁/q₂)` with probability `4/(5(r + 4))` each, a total of `(8/5)/(r + 4)`. This tilt comes from the weight `u(q) = 4^ω(q)`.
- The factor `3/4` is the density of primes `≢ 1 (mod 5)`.
- Integrability over `|t| <= 2` needs `(8/5)α_Q > 1`, i.e. `α_Q > 5/8`.
  - Modulus 3 gives `α_Q = 1/2` and **fails**.
  - **5 is the smallest modulus that works** (`φ(m) >= 3`).
- The tilt `P*(r | n) = 5/(r + 4)` comes from `w(n) = 5^(ω_Q(n))`.
- With weight `k^ω` in place of `4^ω`, integrability needs `φ(m) > 2k/(k − 1)`. So 5 is the least modulus for every `k >= 3`, and it would be 7 for `k = 2`. "The 5 is forced by the 4" overstates (audit D).
- The second 5, in `w(n) = Σ_(q|n) 4^ω(q) = 5^(ω_Q(n))`, is the binomial theorem. They meet as `φ(5) = 4` and `5 = 1 + 4`: an internal COINCIDENCE of the parameter choice.

**Part II: the other 20s and the other 700.**
- **MRT15 Theorem 1.7, corrected (the paper's Theorem 12.2).** For a 1-bounded multiplicative `b`, with `Q = min((log X)^(1/125), (log D)^5)`:

  `sup_α ∫_0^X |Σ_(y<m<=y+D) b(m)e(αm)| dy ≪ DX (e^(−𝓜(b;X,Q)/20) + log log D/log D + (log X)^(−1/700))`.

  So the **20 and the 700 both originate in MRT15**: the pretentious distance enters as `e^(−𝓜/20)` and the frequency-uniform saving as `(log X)^(−1/700)`.
- **Part II uses the 20 and the 700 only qualitatively.** In (12.7) the graph parameters `B` and `D` are fixed and then `Y → ∞`. So `(log X)^(−1/700) → 0`, and `e^(−𝓜/20) → 0` by uniform nonpretentiousness, for any positive constants in place of 1/20 and 1/700. (The middle term `log log D/log D` is used quantitatively, in (12.17).)
- **Part II's own `ρ = 1/20`, in (10.1) and (13.1).** It sets the forbidden-path length `L = B^(1−ρ) = B^(19/20)`. It is used throughout §§13–16, not only in §13 (corrected after audit D):
  - `O(B^ρ)` blocks and `B^(1−2ρ)` thresholds (Prop. 14.1, Lemma 14.3);
  - Lemma 14.5's `B^((1−4ρ)/8)`;
  - Lemma 14.6's `1 − 8ρ = 3/5 > 2ε`, which needs **`ρ < (1 − 2ε)/8 = 0.124975`**;
  - (15.6)'s `M₀ = O(B^(1−ρ+4ε) log B)`;
  - Lemma 16.3's `N_exc = B^(1−ρ/2)` with `1 − ρ + 4ε < 1 − ρ/2`, which needs **`ρ > 8ε = 8·10^(−4)`**.

  So the text needs `0.0008 < ρ < 0.125`: margins of 62.5 times below and only 2.5 times above. (The pre-audit text said "only `ρ > 4·10^(−4)`, margin 125"; MISTAKE-612.)
- **Summary (ARITHMETIC).** No 20 or 700 anywhere in the paper is sharp.
  - 1/700 needs `A >= 560` in Part I and is irrelevant in Part II.
  - 21/20 needs only `> 1.000224`.
  - `ρ = 1/20` needs only `0.0008 < ρ < 0.125`.
  - MRT's `1/20` needs only some positive constant.
  - Other appearances, none sharp either (audit D):
    - Lemma 3.1's circuit depth 20 and `A > C₁ + 20` (its (3.4) needs only `A > max(C₁ + 1, 11)`);
    - Theorem 5.4's "combined depth at most 20";
    - **Lemma 6.5's `exp(−L^(21/20))`**, a second 21/20 backed by an `L^(647/600)` saving;
    - §15's `t₀/20`.

  The final exponent `c = δ/(6WA)` is governed by `δ = 1/200`, the Braverman constant `C₁` (through `A`), and the graph constant `C₃` (through `W`). It is not governed by these numbers.

## 2. Exact refinement of (4.3) (PROVED)

For the four forms `x₁, x₂, x₃, x₁ + x₂ − x₃` (the uniform matroid `U_(3,4)`: every three are independent), Whitney's formula gives the exact avoiding proportion

`χ(p)/p³`, with `χ(t) = t³ − 4t² + 6t − 3 = (t − 1)(t² − 3t + 3)`.

So `ν₄(p) = 4/p − 6/p² + 3/p³` exactly (checked by direct count at `p = 2, 3, 5, 7, 11, 13`). The paper uses the vector `(1,1,1)` as a witness that `ν₄(p) < 1` for every `p`. Over `F₂` it is the unique avoider (`χ(2) = 1`).

This is the finite-field method of Crapo–Rota and Athanasiadis (found independently by Björner–Ekedahl). It is the same inclusion–exclusion over a hyperplane arrangement that computes injective Paley pattern sums in THM-4652 and THM-4653, with quadratic characters inserted.

## 3. Relations to the corpus

1. **Nonbacktracking operators and the K6 − P6 terminal (STRUCTURAL in skeleton only; the mechanisms differ).**
   - Theorem 5.4's `H_(d,d') = 1_(d≠d') A_(d')` records walks with unequal consecutive tuples. Lemma 8.1 transfers its spectral radius to the full sum `Σ B_i`.
   - The Paley terminal `V_q = I_(K6−P6)(q)` ([THM-4653](../../01-canon/theorems/THM-4653-paley-k6-minus-p6-terminal-reduction-to-a-legendre-bilinear-form.md)) has the same skeleton: a sum with "unequal consecutive labels" compared with the full sum. Every edge forces distinctness through `χ(0) = 0`, so the only surviving distinctness conditions are the five **consecutive** pairs of the deleted path.
   - The mechanisms differ. Lemma 8.1 is a resolvent inequality, while THM-4653 is an exact inclusion–exclusion.
   - THM-4653 is the exact, closed-form version of the nonbacktracking-to-full transfer for one pattern: `W_q = −εM_q + 7 + 13ε − 2εq`.
   - HYP-9276 then identifies the full form `M_q`, conjecturally, with the symmetric square of a weight-3 newform.
2. **Singleton cancellation and isolated-edge vanishing (STRUCTURAL, but generic to all centred moment methods).**
   - Lemma 5.5 reads `E[(1_(Z=a) − 1/p) F(Z)] = (1/p) E[F(a) − F(Z)]`: a prime that occurs once is integrated against a mean-zero factor.
   - THM-4652 (3): a pattern component that is a single edge has `Σ* = 0`.
   - Both moment methods are driven by the labels that occur at least twice. Compare the paper's "at most `Jk + S₁/2` distinct labels" with the Paley run expansion, where only paired cherries survive at first order (HYP-9274: `c₁ = 3/2`).
3. **Prohibited words and the safe Mahler language (ANALOGY; retyped after audit D).**
   - Definition 5.1 deletes the sites where a minimal forward-prohibited word starts. That word is defined by a **suffix divisibility** condition `p | Σ_(i=a)^m D_i`, and the deletion costs density (Lemma 5.3) but buys structure on the surviving paths (Lemma 7.1).
   - [THM-4654](../../01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md) deletes nothing. Its safe language is defined by a suffix-sum **inequality** (the `β = 3/2` shift), which Z-number codes satisfy by hypothesis. Its minimal forbidden factors (`11`, `1010^k1` for `1 <= k <= 4`, `(2,6)` followed by a 2, …) force:
     - the block law;
     - odd density at most 1/3;
     - real contraction of the Collatz maps.
   - The only common thread is "suffix conditions define forbidden factors".
4. **Ordinary versus logarithmic averages and Collatz density (ANALOGY; a target).**
   - The paper's advance is control of the terminal scale. It covers every cutoff, rather than the logarithmic averages of Tao 2016, Helfgott–Radziwiłł and Pilatte.
   - Tao's Collatz theorem ("almost all orbits attain almost bounded values", Forum Math. Pi 10 (2022) e12, Theorem 1.3) is stated in logarithmic density. Tao's Remarks 1.4 and 1.16 call a natural-density version plausible but leave it open.
   - The paper's device is short multiplicative bins `T_j = X e^(jη)`, with a per-bin saving `L^(−21/20)` that beats the `L/η` bins. It is the visible template for such an upgrade.
   - The corpus's own merging theorems (THM-4581, THM-4627) are natural-density statements for specific events.
5. **"Saving beats counting" (ANALOGY).** `τ > 1` (a per-bin saving beating the number of bins) has the same shape as the Collatz cycle budget of THM-4633, where the convergent denominators of `log₂3` must beat the number of cycle candidates.
6. **20 and 21 in the K6 − P6 constants (COINCIDENCE).**
   - For `q ≡ 1 (mod 4)`, THM-4653 and HYP-9276 give `W_q = −M_q − 2q + 20` and `W_q = |b(q)|² − q² + 21`. These are the values of `7 + 13ε` and `7 + 14ε` at `ε = 1`.
   - No mechanism links them to the paper's 21/20. They are recorded because the owner asked for every appearance.
7. **Mod 5 (COINCIDENCE).**
   - The paper's modulus 5 is forced by `α_Q > 5/8` (§1).
   - The corpus's mod-5 phenomena are unrelated mechanisms:
     - the Mahler mixed-power maximum 2/5 at the mod-five phases;
     - the Collatz cycle 1/5 with code `100`;
     - THM-4654's equality family `A_t = (2^(18t+8) − 9)/19`.

## 4. What to take from it

- **Paley.** Recast run expansions as nonbacktracking operators on pattern labels. Item 1 says the terminal-by-terminal comparison is exactly computable. The next terminals (HYP-9275's twelve-step cores) could be organised the same way.
- **Collatz.** The bin-and-budget device of §1 is the cleanest available model for converting logarithmic-density statements into every-scale statements. It is worth testing on THM-4627's density-one merging before attempting Tao's theorem.
- **Mahler.** THM-4654 is only an analogue of the prohibited-word deletion (§3, item 3). Its remaining questions are a super-linear factor-complexity bound for Mahler codes of positive integers (Dubickas gives `>= 1.71ℓ`) and a Mahler-side argument for low-complexity codes.
