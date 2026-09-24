# Fixed points, the owner's triangle-inequality sandwich, and arXiv:2502.20642: the claimed proof fails at a swapped quantifier, no contraction in |x−y| can work, and the sign lives at the centre 0

**Status.**
* **REFUTED:**
  1. **The claimed proof of the Collatz conjecture** in T. Kawasaki, *A proof of the Collatz conjecture*, arXiv:2502.20642v1 [math.GM].
     * Theorems 2.1(5), 2.2(5) and 2.3(5) are **false** under the pairwise reading of condition (5), which is the reading the paper's own §3 uses (Theorem B). The counterexamples live on the paper's own space `(N, |x−y|)` and use its own constants `λ ≡ 0, A = 1/2, B = 2, M = 2`. The least counterexample has 3 points (Proposition C).
     * Under the global reading the three theorems are true, but §3's coefficients do not satisfy it, and **no** coefficients can (Theorem D).
     * The error is a swapped quantifier in the proof of Theorem 2.1(5) (Theorem A). For the Collatz coefficients it is hit at **every** odd step `p ≥ 3`.
  2. The paper's intermediate claim `d(T^n x, T^(n+1) x)^2 ≤ A^n d(x,Tx)^2`. It fails at `x = 5, n = 1` (`16 > 9/2`), and at `n = 1` for every odd `x ≥ 3`. The orbit-radius bound it implies fails at `x = 27` (bound `74.80`, actual maximum `4616`).
  3. The literal reading of "the sign decides the Sharkovskii type". Every continuous extension of `T` has points of every period on both half-lines (Theorem K(c)). What the sign does decide is stated in Theorem K(b).

  The Collatz conjecture itself is untouched and remains **OPEN**.
* **PROVED** (hand proofs below; every finite ingredient is re-checked by the scripts):
  * the audit: Theorems A, B, D, E, F, Propositions C and U;
  * the paper's Theorem 3.1 is **true**, with one harmless algebra slip (§2.1);
  * the owner's inequality: Propositions G and H, Theorem I, Propositions J and L;
  * fixed points: Theorem K(a)–(c), and the contraction-metric trichotomy of Theorem F(5).
* **FINITE-EXACT** (script output):
  * Theorem 3.1 at every pair `x, y < 400` (exact) and `< 3000` (vectorized int64), for `3x+1` and verbatim for `3x−1`;
  * condition (5) at every pair `< 400`;
  * the swap-gap census for `p < 5000` (`2499` uncovered steps = the odd `p ≥ 3`); a scratch audit extends it to `p ≤ 10^5`;
  * the counts of the fixed-point chain for periods `p ≤ 24`: the in-`Bad` counts of the mod-192 note are reproduced for `p ≤ 18` and extended;
  * the integral periodic points for `p ≤ 24` are exactly the five known cycles;
  * Rédei's theorem for all `33867` labelled tournaments on `n ≤ 6` vertices.
* **CITED:**
  * Chamberland 1996, via Lagarias's annotated bibliography (arXiv math/0309224v13, entry 35, read) and via Chamberland's survey *An update on the 3x+1 problem* (2003, §6.5, read). The primary paper was not read. Its numbers μ1, μ3 and A2 are re-derived here.
  * Letherman–Schleicher–Wood 1999, via Lagarias's entry 115 and the same survey (§6.6). The primary is **UNVERIFIED**: the publisher page returned a bot-check, which was not bypassed.
  * Borovkov–Pfeifer 2000, via the survey only (**UNVERIFIED**; its claim is re-proved here as Theorem K(c)).
  * Bessaga 1959, via Wikipedia's "Banach fixed-point theorem", section *Converses* (read; primary **UNVERIFIED**).
  * Standard theorems used by name only, with primary sources not read this session (bibliographic details **UNVERIFIED**): Sharkovskii 1964, Li–Yorke 1975, Sperner 1928, Knaster–Kuratowski–Mazurkiewicz 1929, Rédei 1934, Caristi 1976, and Block–Coppel turbulence. Every consequence this note needs from them is either re-proved (the IVT chains of §4.2) or used only as context.
* **OPEN:**
  * Collatz;
  * its cycle half, equivalently: the Banach fixed-point chain meets `Z_{>0}` only in `{1,2}` (§4.1);
  * the integer-cycle classification, i.e. the cycle part of the five-cycle Z-Collatz statement (the chain meets `Z` only in 18 points).
* **UNRESOLVED:** Lagarias's summary of LSW puts a real critical point of `f_0` at `−1/2`. For the formula as printed, `f_0'(−1/2) = 2`, and the non-integer real critical points are `−1.6444, −1.1791, −0.1014` (§3.3(ii)).
* **Verdicts on the owner's items** (details in §5):
  * **REAL:** Lemma 2.1 = AM–QM; the sign law as the two equality cases with `0` in the middle; the `b = 0` central sheet and the Mahler bridge; `{−1,0}` (definitional).
  * **ANALOGY:** orthogonality in `C`; "|x| as the precursor of polynomials" (the underlying identities are REAL); Sperner/Brouwer/Rédei.
  * **REFUTED:** "the sign decides the Sharkovskii type" as a statement about `f`. The stability version is PROVED; the integer-cycle version is FINITE-EXACT for `p ≤ 24` and OPEN in general.

Session `collatz-procgen-20260922` (mac-mini), lane "fixed points / the owner's inequality / arXiv:2502.20642", 2026-09-24. No HYP or THM file was created.
* **Scripts:** `04-computation/experiments/collatz_procgen_20260924_fixedpt_{kawasaki,inequality,chains}.py`, run by `..._fixedpt_run.sh`.
* **Output:** [`collatz_procgen_20260924_fixedpt.out`](collatz_procgen_20260924_fixedpt.out), with sections K1–K12, I1–I7 and F1–F5. Every check is an explicit `raise`.
* **Cost:** about 34 s in total. The peak is 428 MB (the `p = 24` level of F2), with one process at a time. A rerun is byte-identical.

## 0. The request and the short answer

The owner asked to "think Brouwer's fixed point theorem ... connect Rédei, Hamiltonian path, fixed point chain growth" and the paper, and to hone in on the sandwich

`|d(x,z) − d(z,y)| ≤ d(x,y) ≤ d(x,z) + d(z,y)`,  squared: `d(x,z)^2 − 2d(x,z)d(z,y) + d(z,y)^2 ≤ d(x,y)^2 ≤ d(x,z)^2 + 2d(x,z)d(z,y) + d(z,y)^2`.

In the owner's reading the two sides are the positives and the negatives, and the central part is `0`. The sandwich is the proof of the paper's Lemma 2.1.

**Short answer.**
1. **The paper's Lemma 2.1 is correct.** It is the triangle inequality followed by AM–QM, and it throws away exactly the cross term `±2ab`.
2. **The paper's error is elsewhere.** Its proof of the contraction step needs one alternative of condition (5) at the pair `(p, Tp)` or the other alternative at the *swapped* pair `(Tp, p)`. Condition (5) only offers an alternative at each pair separately.
   * Collatz's coefficients fall into the gap at every up-step.
   * The gap is fatal: the theorem as stated is false, already for `x ↦ x+1` on `N` with the paper's own constants.
   * Once repaired, the theorem says only "consecutive distances shrink geometrically". Collatz violates that at every odd `p ≥ 3`, where the distance grows by `3/2` or by `(3p+1)/(2(p+1)) ≥ 4/3`.
   * The paper's coefficient table, copied verbatim, "proves" the false statement that every positive integer reaches 1 under `3n−1`.
3. **The owner's sandwich is where the sign law lives.**
   * For every parity word, `2^p T^p(x) = 3^a x + c_w` with `c_w ≥ 0`. The positive half-line realizes the upper equality of the sandwich and the negative half-line the lower one, with the middle point `z` equal to `0`.
   * The metric `|x−y|` cannot see which case occurs, because `x ↦ −x` is an isometry.
4. **The correct fixed-point theorem for Collatz is Banach's, in `Z_2`, for the inverse branches.**
   * Every parity word has exactly one fixed point, the cycle gate `c_w/(2^p − 3^a)`.
   * This fixed-point chain grows like `2^p`, and its part in `Bad` like `2^(~0.95p)`.
   * It meets `Z` in only 18 points for `p ≤ 24`.
5. **Brouwer's contribution is one-dimensional (the IVT).** Every continuous extension of `T` has points of all periods on both sides. In Chamberland's extension the sign decides stability, not the period set.

## 1. The paper, read exactly

Source: the arXiv HTML of 2502.20642v1, read in full (fetched 2026-09-24; the HTML header reads "arXiv:2502.20642v1 [math.GM] 28 Feb 2025"; sha256 of the fetched file `e48466fe…0f524c9`).

* **§2, the WGP inequality.** Coefficient *functions* `α, β, γ, δ, ε, ζ : X × X → R`. `T` is an `(α,…,ζ)`-weighted generalized pseudocontraction (WGP) if
  `α d(Tx,Ty)^2 + β d(x,Ty)^2 + γ d(Tx,y)^2 + δ d(x,y)^2 + ε d(x,Tx)^2 + ζ d(y,Ty)^2 ≤ 0` for all `x, y`.
* **Lemma 2.1:** `2min(θ,0)(d(x,z)^2 + d(z,y)^2) ≤ θ d(x,y)^2`.
* **Lemma 2.2:** the `λ`-mixtures of the coefficients at `(x,y)` and `(y,x)` again satisfy the WGP inequality.
* **Theorem 2.1** (orbits are Cauchy), **Theorem 2.2** (orbits converge; `X` complete) and **Theorem 2.3** (`T` has a fixed point `u = lim T^n x`) each have five alternative hypotheses. Only (5) is non-degenerate. Write
  `s1 = α+ζ+2min(β,0)`, `r1 = δ+ε+2min(β,0)`, `s2 = α+ε+2min(γ,0)`, `r2 = δ+ζ+2min(γ,0)`
  (all `λ`-mixed), and
  * **alt1(x,y):** `s1 > 0` and `−r1 ≤ A s1`;
  * **alt2(x,y):** `s2 > 0` and `−r2 ≤ A s2`.

  Condition (5) asks for "`A ∈ (0,1)` such that for any `x,y`, alt1 or alt2". Theorem 2.3 adds `α+β+ζ ≥ B` to alt1, `α+γ+ε ≥ B` to alt2, and `|coefficients| ≤ M`.
* **§3.** `T(1) := 1`, `T(x) = x/2` (x even), `(3x+1)/2` (x odd, `x ≥ 3`) on `N = {1,2,…}` with `d = |x−y|`. The piecewise table of Theorem 3.1 has entries in `{−2,…,2}`; on odd/odd pairs they depend on `k−l`, `11k−10l+1` and `−10k+11l+1` (`x = 2k+1`, `y = 2l+1`). Theorem 3.2 takes `λ = 0, A = 1/2, B = 2, M = 2` and concludes `T^n x → 1`, hence Collatz.
* **Typos (harmless).**
  * Theorem 2.1(4) prints `δ+ε+2min(β,0) > 0` where the proof needs `δ+ζ+2min(γ,0) > 0`.
  * The second displayed inequality in the proof of Theorem 2.1 prints `min{γ d(x,y), 0}`.
  * Theorem 2.3's proof prints `α_λ(u, T^n u)` and "for any `x ∈ N`".

## 2. The audit

### 2.1 What is correct (PROVED; K1, K2)

* **Lemma 2.1 and Lemma 2.2 are correct.** Lemma 2.1 is analysed in §3.1.
* **Theorem 3.1 is true.** The table satisfies the WGP inequality at every pair (K1: all pairs `< 400` exactly, all pairs `< 3000` vectorized). Nine case polynomials were re-derived symbolically; see the table below.
  * The odd/odd general form `(18+4δ0)(k−l)^2 − 5β0(k−l)(k+l+2) + ε0(k+1)^2 + ζ0(l+1)^2` agrees with the paper.
  * So do all five sub-case factorizations: `2(k+1)(11k−10l+1)`, `4(k−l)(6k−l+5)`, `2(l+1)(−10k+11l+1)`, `4(k−l)(k−6l−5)` and `(k−l)^2(4−5(k+l))`.

  | pair type | derived | printed in the paper |
  |---|---|---|
  | `(1, 2l)` | `−2l^2+2l` | same |
  | `(1, 2l+1)` | `−6l^2+4l+2` | same |
  | `(2k, 1)` | `−2k^2+2k` | same |
  | `(2k, 2l)` | `−k^2+2kl−2l^2` | same |
  | `(2k, 2l+1)` | **`−2l^2+1`** | `−2l^2−4kl−2k−1` (slip) |
  | `(2k+1, 1)` | `−6k^2+4k+2` | same |
  | `(2k+1, 2l)` | `−2k^2+1` | same |

  The slip is harmless: `−2l^2+1 ≤ −1`. At `(2,3)` the left side is `−1`, not the printed `−9`.
* **Condition (5), pairwise reading, holds for the table at every pair `< 400`**, with `|coefficients| ≤ 2`.
  * Census: 99302 pairs satisfy alt1 only, 59700 alt2 only and 199 both.
  * The odd/odd coefficient tuples are exactly the paper's lists: three for `x > y`, four for `x ≤ y`.
* **Theorem 2.1(1)–(4) are correct but degenerate.** They force `T^2 = T` or `T = id`.

So the Collatz-specific algebra is right. The failure is in §2 of the paper.

### 2.2 Theorem A: the swap gap (PROVED; K3)

**Theorem A.** Let `T` be a WGP with coefficient functions and `λ`, and let `A ∈ (0,1)`.
1. For every `p ∈ X`:
   * `s1(p,Tp) d(Tp,T²p)^2 + r1(p,Tp) d(p,Tp)^2 ≤ 0`;
   * `s2(Tp,p) d(Tp,T²p)^2 + r2(Tp,p) d(p,Tp)^2 ≤ 0`.

   Hence **P(p) := alt1(p,Tp) or alt2(Tp,p)** implies `d(Tp,T²p)^2 ≤ A d(p,Tp)^2`.
2. These are the only uses of (5) in the proof of Theorem 2.1(5). The proof therefore needs P(p) at every `p` on the orbit. Hypothesis (5) gives "alt1 or alt2" at `(p,Tp)` and, separately, at `(Tp,p)`. The case "alt2 only at `(p,Tp)` and alt1 only at `(Tp,p)`" is not covered.
3. For the paper's table, P(p) fails exactly at the odd `p ≥ 3`, although (5) holds at both pairs of every step.

*Proof.*
1. At `(x,y) = (p,Tp)`, Lemma 2.1 with `z = y` gives `β d(x,Ty)^2 ≥ 2min(β,0)(d(x,y)^2 + d(y,Ty)^2)`. Here `d(Tx,y) = 0`, `d(x,y) = d(x,Tx) = d(p,Tp)` and `d(Tx,Ty) = d(y,Ty) = d(Tp,T²p)`, so the WGP inequality becomes the first line.
2. At `(x,y) = (Tp,p)`, Lemma 2.1 with `z = x` bounds `γ d(Tx,y)^2` in the same way, and `d(x,Ty) = 0`; this gives the second line.
3. *Item 3.* Let `p` be odd with `p ≥ 3`.
   * If `Tp` is even, the pairs have types `(odd, even)` and `(even, odd)`. The table gives `s1 = 0−2−4 = −6` at the first, so alt1 fails, and `s2 = 0−2−4 = −6` at the second, so alt2 fails.
   * If `Tp` is odd, then `Tp > p`. At `(p,Tp)` we have `β0 ∈ {−1,−2}` and `ζ0 = 0`, so `s1 = 2+2β0 ≤ 0`. At `(Tp,p)` we have `γ = −β0` with `β0 ∈ {1,2}` and `ε0 = 0`, so `s2 = 2−2β0 ≤ 0`.
   * Condition (5) holds via alt2 at `(odd, even)`: `s2 = 2`, `−r2 = 1 = A s2`, `α+γ+ε = 2 = B`. It holds via alt1 at `(even, odd)`. On odd/odd pairs it holds via alt2 when `x < y` and via alt1 when `x > y` (the paper's own lists, §2.1).
   * Even steps are covered by alt1 at `(p, p/2)`.
   * K3 confirms this for every `p < 5000`: 2499 uncovered steps, which are exactly the odd `p ≥ 3`. ∎

### 2.3 Theorem B: counterexamples, and Proposition C: minimality (PROVED; K6, K7)

**Theorem B.** Under the pairwise reading, Theorems 2.1(5), 2.2(5) and 2.3(5) are false.
* **(i) The paper's own space and constants.**
  * `X = N`, `d = |x−y|`, `T(x) = x+1`, `λ ≡ 0`, `A = 1/2`, `B = 2`, `M = 2`.
  * Coefficients: `V1 = (1,1,−2,0,0,0)` if `x ≥ y`, and `V2 = (1,−2,1,0,0,0)` if `x < y`.
  * With `e = x−y`, the WGP left side is `e^2 + (e−1)^2 − 2(e+1)^2 = −6e−1` for `e ≥ 0`, and `e^2 − 2(e−1)^2 + (e+1)^2 = 6e−1` for `e < 0`. Both are `< 0`.
  * `V1` satisfies alt1: `s1 = 1`, `r1 = 0`, `α+β+ζ = 2`. `V2` satisfies alt2: `s2 = 1`, `r2 = 0`, `α+γ+ε = 2`.
  * Yet `T^n x = x+n` is not Cauchy, and `T` has no fixed point.
  * The swap appears at every step: `(p, p+1)` is alt2-only and `(p+1, p)` is alt1-only.
* **(ii) Three points inside N, with the paper's constants.** `X = {5,7,10}`, and `T` is the `3n−1` map `5 → 7 → 10 → 5`. The table (K6b; at most two non-zero entries per pair):

  | pair | role | `(α,β,γ,δ,ε,ζ)` | alternative |
  |---|---|---|---|
  | `(5,5), (7,7), (10,10)` | diagonal | `(2,0,0,0,0,0)` | alt1 |
  | `(5,7)` | `(p,Tp)` | `(0,−1,0,0,2,0)` | alt2 |
  | `(7,10)` | `(p,Tp)` | `(0,0,0,0,2,−1)` | alt2 |
  | `(10,5)` | `(p,Tp)` | `(0,0,0,−1,0,2)` | alt1 |
  | `(7,5)` | `(Tp,p)` | `(0,0,−1,0,0,2)` | alt1 |
  | `(5,10)` | `(Tp,p)` | `(0,0,0,−1,2,0)` | alt2 |
  | `(10,7)` | `(Tp,p)` | `(0,0,0,0,−1,2)` | alt1 |

  The same table works on `{0,1,2}` with `0 → 1 → 2 → 0`. The coordinator's table `(M,0,−M,0,0,0) / (1,−M,M,0,0,0) / (1,M,−M,0,0,0)` also works, but needs `M ≥ 4` (at the pair `(1,2)`).
* **(iii) The fixed-point step has the same gap.**
  * `X = {0} ∪ {1/n : n ≥ 1}` (complete), with `T(1/n) = 1/(n+1)` and `T(0) = 1`.
  * Coefficients: `(1,1,−5,0,0,0)` when `d(Tx,y) ≥ d(x,Ty)`, else `(1,−5,1,0,0,0)`; `λ ≡ 0, A = 1/2, B = 2, M = 5`.
  * Every orbit converges to `0`, so the conclusion of Theorem 2.2 holds. But `T(0) = 1`, and `T` has no fixed point.
  * The fixed-point step of Theorem 2.3's proof needs alt1 at `(T^n x, u)` or alt2 at `(u, T^n x)` with `u = 0`. For `n ≥ 2` the table has alt2 at `(1/n, 0)` and alt1 at `(0, 1/n)`.
  * *Proof of the hypotheses:* Proposition U with `C = 2`.
    * For `m > n`, `d(x,Ty) = 1/n − 1/(m+1)` exceeds `|x−y|`, which is at least `|Tx−Ty|`, because `T` is 1-Lipschitz on `{1/k}`.
    * Symmetrically for `m < n`.
    * For `(1/n, 0)` and `(0, 1/m)` with `n, m ≥ 2`, the larger of `d(Tx,y)` and `d(x,Ty)` is at least `1/2`, while `d(Tx,Ty) < 1`.
    * The pairs involving `1` have ratio 1.
  * K6d checks this exactly on `{0} ∪ {1/n : n ≤ 300}`. ∎

**Proposition C (2-cycles; minimality).**
1. Let `{a,b}` be a 2-cycle of `T` in any metric space. For every `λ`, every `A ∈ (0,1)` and every `B, M`, no coefficient values at `(a,b)` satisfy both the WGP inequality and "alt1 or alt2".
2. Hence every counterexample to Theorems 2.1(5)–2.3(5) has at least 3 points, and Theorem B(ii) shows that 3 suffice.
3. The shortcut map's own 2-cycle `{1,2}` makes (5) unsatisfiable for the true `T`. **The paper's redefinition `T(1) := 1` is therefore forced by its hypothesis**, not a normalization.

*Proof.*
1. At `(a,b)` the distance vector is `D^2 (1,0,0,1,1,1)`, so WGP reads `α+δ+ε+ζ ≤ 0`. But alt1 gives `α+δ+ε+ζ = s1 + r1 − 4min(β,0) ≥ (1−A)s1 − 4min(β,0) > 0`, and alt2 gives the same with `(s2, γ)`. By Lemma 2.2 the `λ`-mixed coefficients satisfy the same WGP inequality at the same pair, so no `λ` helps. A grid sanity check over `{−3,…,3}^6` finds no vector for `A = 1/2` or `A = 99/100` (K7).
2. A 1-point space has a fixed point. A 2-point map without a fixed point is a swap. A 2-point map with a fixed point has eventually constant orbits.
3. With `T(1) = 2` the pair `(1,2)` has the vector `(1,0,0,1,1,1)`. The paper's own entry there gives left side `+1`. ∎

**Two readings of (5).**
* **Pairwise** (`∀(x,y)`: alt1 or alt2). This is the reading §3 uses: its verification assigns alt1 to the types `(1,1), (1,even), (1,odd), (even,even), (even,odd), (odd,odd; x>y)` and alt2 to the rest. Under it, Theorem B applies.
* **Global** (alt1 at every pair, or alt2 at every pair). Under it Theorems 2.1–2.3(5) are **true**; the proofs go through as written. But the table violates both global alternatives: alt1 fails at `(3,2)` and alt2 at `(2,3)` (K2). Theorem D shows that no table can satisfy either for Collatz.

### 2.4 Theorem D: the corrected theorem is one-step decay, and Collatz fails it for every coefficient choice (PROVED; K10)

**Theorem D.**
1. **(Corrected Theorem 2.1.)** Suppose P(p) holds at every `p` of the forward orbit of `x`. Then:
   * `d(T^n x, T^(n+1) x)^2 ≤ A^n d(x,Tx)^2`;
   * `{T^n x}` is Cauchy, with `d(T^n x, T^m x) ≤ A^(n/2)(1 − A^(1/2))^(−1) d(x,Tx)`.

   In a uniformly discrete space the orbit is therefore eventually fixed.
2. **(Nothing more.)** The inequality `d(Tp,T²p)^2 ≤ A d(p,Tp)^2` holds iff the single vector `(B,0,0,−AB,0,0)` placed at `(p,Tp)` satisfies WGP there. That vector also satisfies alt1, with the `B`-clause. Such placements conflict only on 2-cycles. **So the corrected hypothesis is exactly one-step geometric decay of consecutive distances along the orbit; the six-coefficient machinery adds nothing to it.**
3. **(Collatz.)** For odd `p ≥ 3`:
   * `d(Tp,T²p)/d(p,Tp) = 3/2` if `p ≡ 3 (mod 4)`;
   * `d(Tp,T²p)/d(p,Tp) = (3p+1)/(2(p+1)) ≥ 4/3` if `p ≡ 1 (mod 4)`.

   So P(p) fails for **every** `λ`, every coefficient choice and every `A < 1`. At `p = 5` it would force `16 ≤ 9A < 9`. Every `x` that is not a power of 2 has an odd number `≥ 3` in its orbit (its odd part), so the corrected theorem applies to no nontrivial orbit.

*Proof.*
1. Item 1 is Theorem A(1), iterated.
2. Item 2: the WGP inequality at `(p,Tp)` with this vector reads `B d(Tp,T²p)^2 − AB d(p,Tp)^2 ≤ 0`; alt1 holds because `s1 = B`, `−r1 = AB` and `α+β+ζ = B`.
3. Item 3:
   * If `p = 4m+1`, then `Tp = 6m+2` is even and `T²p = (3p+1)/4`.
   * If `p = 4m+3`, then `Tp = 6m+5` is odd and `d(Tp,T²p) = (Tp+1)/2 = (3p+3)/4`.
   * In both cases `d(p,Tp) = (p+1)/2`. ∎

K10 confirms, for every `p < 5000`, that the steps covered by the table are exactly the decay steps, and that the failing steps have squared ratio `≥ 16/9`. For `3n−1` every up-step stretches by at least `3/2` (checked for odd `p < 20001`).

### 2.5 Proposition U: the hypothesis carries no dynamics (PROVED; K8)

**Proposition U.** Suppose `T` has no 2-cycle and
`C := sup d(Tx,Ty) / max(d(Tx,y), d(x,Ty)) < ∞`, the sup taken over the pairs where the maximum is positive.

Put `K = (B/2)(C^2+1)`. Use `(B/2, B/2, −K, 0,0,0)` [alt1] when `d(Tx,y) ≥ d(x,Ty)`, and `(B/2, −K, B/2, 0,0,0)` [alt2] otherwise. These satisfy WGP and (5) with `λ ≡ 0`, any `A ∈ (0,1)` and `M = max(B/2, K)`.

*Proof.*
* Write `u = d(Tx,y)`, `v = d(x,Ty)` and `w = d(Tx,Ty)`.
* If `u ≥ v` and `u > 0`, then `(B/2)(w^2 + v^2) ≤ (B/2)(C^2+1)u^2 = K u^2`. The case `v > u` is symmetric.
* If `u = v = 0`, then `Tx = y` and `Ty = x`, so (no 2-cycle) `x = y` is a fixed point and all six distances vanish. ∎

For the paper's `T` and for `T_−` on `N`, `C ≤ 11`.
* Write `Tx = λ_x x + μ_x` with `λ ∈ {1/2, 3/2}` (and `λ = 1` at the paper's `T(1)`), `|μ| ≤ 1/2`, and `m = max(u,v) ≥ 1`.
* Then `|1 − λ_x λ_y| ≥ 1/4` unless `x = y = 1`. This gives `x ≤ 5 + 10m`, hence `|x−y| ≤ 3 + 6m` and `w ≤ 3 + 8m ≤ 11m`.
* Numerically, `C = 5` for the paper's `T` (attained at `(12,7)`) and `C < 4` for `T_−`, over `x, y < 2500`.

**So the hypothesis of Theorem 2.3(5) holds, with bounded coefficients, for essentially any map without (near-)2-cycles: `3x+1`, `3x−1`, `x+1`.** It contains no dynamical information.

### 2.6 Theorem E: the SHEET control (PROVED; K9)

**Theorem E.** Let `T_−(x) = x/2` (x even) and `(3x−1)/2` (x odd) on `N`; then `T_−(1) = 1` with no redefinition needed. The paper's table, **verbatim**, makes `T_−` a WGP with (5).
* The case polynomials, with `k, l ≥ 1`:
  * `(1,2l)`: `−2l^2+2l`;
  * `(1,2l+1)`: `−6l^2`;
  * `(2k,1)`: `−2k^2+2k`;
  * `(2k,2l)`: `−(k−l)^2 − l^2`;
  * `(2k,2l+1)`: `−2l^2−4l−1`;
  * `(2k+1,1)`: `−6k^2`;
  * `(2k+1,2l)`: `−2k^2−4k−1`.
* The odd/odd case is `(18+4δ0)(k−l)^2 − 5β0(k−l)(k+l) + ε0k^2 + ζ0l^2`, with sub-cases:
  * `2k(11k−10l)`, which is `< 0` because `11k−10l ≤ −1`;
  * `4(k−l)(6k−l)`, which is `< 0` because `11k ≥ 10l` gives `6k > l`;
  * `2l(11l−10k)` and `4(k−l)(k−6l)`, symmetrically;
  * `(k−l)^2(14−5(k+l))`, which is `≤ 0` because `k = l` or `k+l ≥ 3`.
* Condition (5) depends only on the table.
* The only fixed point of `T_−` in `N` is 1.

So **the argument of Theorem 3.2, verbatim, "proves" that every positive integer reaches 1 under `3n−1`**. That is false: `5 → 7 → 10 → 5`, and the 11-cycle through 17.

K9 checks every pair `< 400` exactly and every pair `< 3000` in vectorized form, and re-derives all the polynomials symbolically.

**Reading.** In the foundry's typing the paper's method is SHEET-blind, now as a theorem: its hypothesis holds on both sheets, so it cannot imply a statement that is true on one sheet and false on the other. This is Corollary 7 of the [mod-192 note](collatz_procgen_20260924_inverse_tree_mod192.md), instantiated. Proposition U shows more: the hypothesis is satisfiable for any map without near-2-cycles.

### 2.7 Theorem F: no contraction certificate in |x−y| (PROVED; K4, K5, K11, K12)

**Theorem F.** Let `T` be the Collatz map on `N`, with `T(1) = 1` or `T(1) = 2`.
1. `sup_p d(Tp,T²p)/d(p,Tp) = 3/2`, attained at every `p ≡ 3 (mod 4)`. Every odd `p ≥ 3` stretches the distance by at least `4/3`; even steps never stretch (checked to `2·10^5`).
2. For `x = 2^m − 1` one has `T^j x = 3^j 2^(m−j) − 1` for `j ≤ m`. Hence:
   * `d(T^(m−1)x, T^m x) = (3/2)^(m−1) d(x,Tx)`;
   * `d(x, T^m x) = ((3^m − 2^m)/2^(m−1)) d(x,Tx)`.
3. **Consequences.**
   * (a) Any hypothesis implying `d(Tp,T²p) ≤ d(p,Tp)` for every `p` fails at every odd `p ≥ 3`. This covers Banach (`k < 1`), Kannan (`d(Tx,T²x) ≤ k(d(x,Tx) + d(Tx,T²x))`, `k < 1/2`), Chatterjea (`d(Tx,T²x) ≤ k d(x,T²x)`, `k < 1/2`), Reich (ratio `≤ (a+b)/(1−c) < 1`) and the corrected Theorem 2.1. Each reduction is one line with `y = Tx`.
   * (b) Any uniform orbit-radius bound `d(x,T^n x) ≤ K d(x,Tx)` fails; the paper's bound `K = (1−√A)^(−1)` is one such.
   * (c) Any uniform decay `d(T^n x,T^(n+1)x) ≤ g(n) d(x,Tx)` needs `g(n) ≥ (3/2)^n`.
4. **Caristi.** Take the paper's `T`. A function `φ ≥ 0` with `d(x,Tx) ≤ φ(x) − φ(Tx)` exists iff every orbit reaches 1. Then `φ(x) = Σ_{n<τ(x)} ceil(T^n x/2)`, the orbit's total variation, satisfies it with equality. For example `φ(27) = 20328`.
5. **Contraction metrics (the trichotomy).**
   * `|x−y|`: no contraction (items 1–3).
   * *Some* complete metric makes `T` a Banach contraction iff `N` has no `T`-cycle other than `{1}`. This is Bessaga's converse (CITED), and it is the **cycle half only**: a divergent orbit would still converge to 1 in such a metric.
   * A **uniformly discrete** complete contracting metric exists iff Collatz holds.
     * (⇒) `ρ(T^n x, 1) ≤ q^n ρ(x,1) < inf ρ` forces `T^n x = 1`.
     * (⇐) `ρ(x,y) = K^(max(h(x),h(y)))` (`x ≠ y`, `h` = steps to 1, `K = 1/q`) is an ultrametric, `≥ 1` off the diagonal, and `ρ(Tx,Ty) ∈ {0, qρ(x,y)}`. K12 checks this for `x, y ≤ 1500`.

*Proof.*
* Items 1 and 2 are the computations of Theorem D(3) and induction on the up-run from `2^m − 1`.
* Item 4:
  * (⇐) telescoping.
  * (⇒) the partial sums of the step lengths are bounded by `φ(x)`. Steps are integers, so the orbit is eventually constant. The only fixed point is 1.
  * The formula uses `|T(n) − n| = ceil(n/2)` for `n ≥ 2` (THM-4470), checked to `10^5`.
* Item 5 is as displayed. ∎

**Conclusion of the audit (precise).**
1. The paper does not prove Collatz. Its Theorems 2.1(5)–2.3(5) are false as used (Theorem B), and their correct form cannot apply (Theorem D).
2. Its §3 hypotheses are equally satisfied by `3n−1` (Theorem E) and by `x+1` (Theorem B(i)).
3. No contraction-type argument in `|x−y|` can prove Collatz. Here contraction-type means that the argument's engine is consecutive-distance monotonicity, a uniform orbit-radius bound or uniform decay; up-steps stretch consecutive distances by up to `3/2`, and Mersenne starts defeat every uniform bound (Theorem F(1)–(3)).
4. The metric principles that survive (Caristi; contraction metrics) are reformulations: Caristi and the uniformly discrete version are equivalent to Collatz, and Bessaga's version to its cycle half (Theorem F(4)–(5)).
5. The paper's last step ("the distance between elements in `N` is 1 or more") is exactly where the divergence half would have to enter.

## 3. The owner's inequality, made exact

### 3.1 (a) Lemma 2.1 is AM–QM, and it discards the cross term (PROVED; I1, I6)

**Proposition G.** Write `a = d(x,z)`, `b = d(z,y)`, `c = d(x,y)`.
* For `θ ≥ 0`, Lemma 2.1 reads `0 ≤ θc^2`: the lower half of the sandwich, `(a−b)^2 ≤ c^2`, is discarded.
* For `θ < 0` it reads `c^2 ≤ 2(a^2 + b^2)`, i.e. the triangle inequality `c^2 ≤ (a+b)^2` followed by AM–QM, `(a+b)^2 ≤ 2(a^2+b^2) ⟺ 2ab ≤ a^2+b^2 ⟺ (a−b)^2 ≥ 0`.
* Equality (`θ < 0`) holds iff `c = a+b` and `a = b`, i.e. `z` is a metric midpoint.

**The cross term `±2ab` is replaced by its sign-blind bound.**

Along Collatz triples `(p, Tp, T²p)`:
* every triple is an equality case of the owner's sandwich: upper for monotone triples, lower for turning triples;
* Lemma 2.1 loses exactly `(d1−d2)^2` on monotone triples and `(d1+d2)^2` on turning ones;
* it is never sharp, since `d1 = d2` occurs only at the turning triple `2, 1, 2` (checked for `p ≤ 10^5`).

**The owner's thesis, made exact.** Cauchy–Schwarz and AM–GM are "places where structure was abstracted away". Here the abstracted structure is the orientation, i.e. the sign of the cross term.
* The THM-4470 form is exact. On consecutive pairs `{2i−1, 2i} → {3i−1, i}`:
  * `(a+b)` is fixed: a **linear** invariant;
  * `ab` is scaled by `3/4 + 1/(8i−4)`: a **quadratic** contraction;
  * `Δ(2ab) = −2i(i−1) = −Δ(a^2+b^2)` (checked to `i = 10^6`).
* Equivalently `QM^2 + GM^2 = 2AM^2`: with the AM fixed, `T` lowers the GM and raises the QM by the same amount. Lemma 2.1 keeps only `QM ≥ AM`. The Collatz drift `log(2/√3)`, the AM–GM gap of the step factors `1/2` and `3/2`, lives in the discarded cross term.

The paper's error is not in this lemma (§2); the lemma is lossy but sound.

### 3.2 (b) The sign law is the pair of equality cases, with `0` in the middle (PROVED; I2, I3)

**Theorem I.** Let `w` have length `p` and `a` ones, and let `c_() = 0`, `c_(z e) = 3^e c_z + e 2^|z|` (THM-4469's carry `R`), so `c_w ≥ 0`, and `c_w > 0` if `a ≥ 1`. For every real `x`, the affine branch `L_w(x) = (3^a x + c_w)/2^p` satisfies `2^p|L_w(x)| = d(3^a x, −c_w)`, and with the owner's middle point `z = 0`:
* `x ≥ 0` ⟺ `0` lies between `3^a x` and `−c_w` ⟺ **upper equality**: `2^p|L_w(x)| = d(3^a x, 0) + d(0, −c_w) = 3^a|x| + c_w`. The cross term is `+2·3^a|x|c_w`.
* `x ≤ 0` ⟺ **lower equality**: `2^p|L_w(x)| = |3^a|x| − c_w|`. The cross term is `−2·3^a|x|c_w`.

For `x` an integer (or a 2-adic rational) in the cylinder of `w`, `L_w(x) = T^p(x)`. For negative integers `3^a|x| ≥ c_w` (orbits stay negative).

**Transport** (mod-192 note, Theorem 6): `T_+^p(−x) = −T_−^p(x)`. So the lower side for `3n+1` on `−N` *is* `3n−1` on `N`.

Checked (I3): all words of length `p ≤ 12`, with lifts of both signs, 57318 cases.

**Proposition H (the cross-term sign on orbits is side-blind).** On the line every consecutive orbit triple is an equality case of the sandwich.
* `d(x_(n−1), x_(n+1)) = d1 + d2` iff `x_(n−1)` and `x_n` have the same parity, and `|d1 − d2|` iff they do not. This holds on **both** half-lines (I2: 4.24 million triples from `[−2·10^4, 2·10^4]`).
* A zero step, the degenerate "centre", occurs only at the fixed points `0` and `−1`.
* A strictly interior triangle, with the cross term `0` and both steps non-zero, needs a second dimension, where it is orthogonality (§3.3(ii)).

**Reading.** The owner's instinct is exact: the sign lives at the centre `0` of the sandwich. Two consequences follow.
* Restoring the cross term that Lemma 2.1 discards would *not* make the metric method side-aware: on orbits the cross term's sign is the parity-agreement bit, a function of the parity word alone.
* The side-aware datum is the position of the points relative to `0`. It is an odd observable, which `|x−y|` never sees, because `x ↦ −x` is an isometry.

### 3.3 (c) "The central part representing 0": three candidates

**(i) The `b = 0` sheet — REAL (Proposition J; I4).** Write `T_b^w(x) = (3^a x + b c_w)/2^p`.
1. **Midpoint and odd part.** `T_0^w = (T_+^w + T_−^w)/2` is the odd part of `T_+^w`, and `c_w/2^p` is its even part. In `(T_b^w x)^2 = (9^a x^2 + 2b·3^a c_w x + b^2c_w^2)/4^p` the cross term has sign `b·sgn(x)`. So `b = −1, 0, +1` are the owner's left side, centre and right side (for `x > 0`). Checked for all 2046 words with `p ≤ 10`.
2. **Recentring.** On odd `x`, `T_b(x) + b = (3/2)(x + b)`: each sheet's odd branch is the central map `x ↦ 3x/2` (Mahler's), recentred at its fixed point `−b`. Along an up-run `x_k + b = (3/2)^k(x_0+b)`. The even branch is centred at `0`. The `b = 0` sheet is where the two centres coincide.
3. **Shadowing.** `T_b^n(x) = (3^(a_n)/2^n) ξ_n` with `ξ_n = x + b c_n/3^(a_n)` (exact; `|x| ≤ 500`, `n ≤ 60`). Every orbit point is a **central-sheet image of a shifted start**; `ξ_n` increases on the plus sheet and decreases on the minus sheet.
4. **The Mahler bridge (THM-4469).** For the blocks `0111101110` and `1101100111` (`R = 4726, 4727`; cylinder classes `990` and `187 mod 1024`, re-derived in I4) and 40 random block sequences of length 6, realized by actual integers:
   * `x_j = ξα^j − t_j`, with `t_j ∈ [4726, 4727]/1163 · (1 − ρ^(6−j))` exactly;
   * the `3n−1` orbit of `−x` is `(−ξ)α^j + t_j`.

   So THM-4469 is precisely the statement that the `±1` orbits shadow the central (Mahler) orbit `ξ(2187/1024)^j` with a *confined carry*, below it on the plus sheet and above it on the minus sheet. Mahler's own `3/2` problem is the all-odd word, which THM-4469 already records as not transferring.

**(ii) Orthogonality in `C` — ANALOGY (exact complex analysis; I5).**
* **Pythagorean middle.** `|z+w|^2 = |z|^2 + |w|^2 + 2Re(z w̄)`: the middle case is `Re(z w̄) = 0`. For `L_w(z)` the cross term is `2·3^a c_w Re z`:
  * the right half-plane is `+2ab`, the left half-plane `−2ab`;
  * the imaginary axis is the orthogonal middle.
* **Chamberland's `f(x) = x + 1/4 − ((2x+1)/4)cos(πx)`.** `f = T` on `Z` and `f'(n) = 1 − (−1)^n/2`, i.e. `1/2` or `3/2`. It keeps the branch slopes and is conformal at the integers.
* **LSW's `f_h`.** `f_0(z) = z/2 + (1−cos πz)(z+1/2)/2 + (1/2 − cos πz) sin(πz)/π`, as printed in Chamberland's survey §6.6. It satisfies `f_0 = T` on `Z` and `f_0'(n) = 0`: every integer is critical, as Lagarias's entry 115 states for every `f_h`.
  * Here `f_0'(x) = sin(πx)(2 sin(πx) + (π/2)(x+1/2))`. At each critical integer `f_0(n+y) − T(n) ≈ c_n y^2` and `f_0(n+iy) − T(n) ≈ −c_n y^2` (checked at `n = −5, −2, 1, 2, 7`).
  * **The fold sends the orthogonal direction onto the opposite real ray.** This is the exact content of "the sharp corner at 0 introducing orthogonality": squaring, the smooth version of `|x|`, identifies `±x` and turns the orthogonal axis into the sign.
* **Collatz content: stability, and only in Chamberland's extension.**
  * In Chamberland's `f` an integer cycle has multiplier `3^a/2^p`, so it is attracting iff its points are `≥ 0` (sign law; Theorem K(b)).
  * In the LSW maps every integer cycle is superattracting on both sides: stability is sign-blind there.
  * Nothing else transfers.
* **Discrepancy (UNRESOLVED).** Lagarias's summary of LSW lists `−1/2` as a real critical point of `f_0`, but for the printed formula `f_0'(−1/2) = 2`, and the non-integer real critical points are `−1.644421083, −1.179060746, −0.1013635846`. Lagarias's own printed formula has `1 + cos πz`, which does not interpolate `T` at even integers, so it is likely a typo. The primary paper is not read, and nothing here depends on this point.

**(iii) The fixed points `0` and `−1` — REAL but definitional (I7).**
* `Fix(T_b) ∩ Z = {0, −b}`. These are the Banach fixed points of `D` and `E_b` (the one-letter words, §4.1) and the centres of the two affine branches; `E_b(x) + b = (2/3)(x+b)`.
* The only consecutive pairs mapped onto themselves (`|i| ≤ 10^5`) are:
  * for `3n+1`: `{−1,0}` (pointwise) and `{1,2}` (swapped);
  * for `3n−1`: `{0,1}` and `{−2,−1}`, the mirror images.
* `{−1,0}` is THM-4470's pair 0, whose obstruction is Applegate–Lagarias's `−1`.
* The "central pair straddling the sign" is where `b` places the odd branch's centre: `−1` or `+1`. This is the same fact as (i).2, and it gives no new leverage.

### 3.4 (d) "|x| as the precursor of polynomials, with its corner at 0" (PROVED identities; I6)

**Proposition L.**
1. `|T_+(x)| = T_(sgn x)(|x|)` for every integer `x` (checked for `|x| ≤ 10^6`, and iterated for `p ≤ 50` on `|x| ≤ 10^4`), with `T_0(0) = 0`. On odd `x`, `|T(x)| = (3|x| + sgn x)/2`: the correction is `sgn x = d|x|/dx` (`x ≠ 0`).
   * On the size coordinate `|x|` the single map `T_+` becomes the pair `T_(+1), T_(−1)`, selected by the derivative of `|x|`.
   * The corner of `|x|` at `0` is exactly where the selector jumps.
   * Iterated: `|T^p(x)| = (3^a|x| + sgn(x) c_w)/2^p` (Theorem I).
2. **Linear versus quadratic.** Pair sums are conserved (linear) and pair products contract by `3/4 + 1/(8i−4)` (quadratic), with the exact transfer `Δ(2ab) = −Δ(a^2+b^2)` (§3.1).
3. **Smoothing the corner erases the sign.** Even observables (`|x|`, `x^2`, `|x−y|`, any polynomial in `x^2`) are invariant under `x ↦ −x`, hence side-blind (mod-192 note, Theorem 6). The sign law enters only through odd observables (`x`, `sgn x`, the cross term `2·3^a c_w x`).
   * Passing from `|x|` to the "polynomial" `x^2`, or to the holomorphic fold of LSW, keeps the side-blindness and loses the selector.

Verdict: the identities are REAL; "precursor of polynomials" is an ANALOGY. Its exact content is that the corner, i.e. the odd part `sgn x`, is the only place the sheet is recorded.

## 4. Brouwer, Banach and fixed-point chains

### 4.1 (a) Banach in `Z_2`: the fixed-point chain and its growth (PROVED; FINITE-EXACT counts; F1, F2, F3)

**Theorem K(a).**
1. `D(x) = 2x` and `E(x) = (2x−1)/3` are contractions of `Z_2` with factor `1/2` (`3` is a 2-adic unit), with fixed points `0` and `−1`.
2. For a forward word `w` of length `p`, the inverse word `W = G_(w_0) ∘ … ∘ G_(w_(p−1))` (`G_0 = D`, `G_1 = E`) contracts by `2^(−p)`. So it has **exactly one** fixed point, `x_w = c_w/(2^p − 3^a)`, the cycle gate.
3. `T^p(x_w) = x_w` with parity word `w`.
4. `Per_p(T) ∩ Z_2 = {x_w : |w| = p}`: `2^p` points, primitive count `Σ_(d|p) μ(d) 2^(p/d)`.
5. `x_w ∈ Z` iff `(2^p − 3^a) | c_w`.

**The fixed-point chain is the set of Banach fixed points of words.** F1 checks all 8190 words with `p ≤ 12`: Banach iteration from 0 reaches `x_w` modulo `2^64`, `T^p(x_w) = x_w` exactly, and distinct words give distinct points.

Growth by period (F2; the counts are computed for all words and Möbius-inverted, since every statistic used is invariant under `u ↦ u^k`):

| `p` | points `2^p` | primitive | cycles | primitive `x_w > 0` | primitive `x_w < 0` | in `Bad` | integral |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 2 | 2 | 0 | 1 | 1 | 2 |
| 2 | 4 | 2 | 1 | 2 | 0 | 0 | 2 |
| 3 | 8 | 6 | 2 | 3 | 3 | 1 | 3 |
| 4 | 16 | 12 | 3 | 8 | 4 | 2 | 0 |
| 5 | 32 | 30 | 6 | 25 | 5 | 3 | 0 |
| 6 | 64 | 54 | 9 | 36 | 18 | 6 | 0 |
| 8 | 256 | 240 | 30 | 208 | 32 | 16 | 0 |
| 11 | 2048 | 2046 | 186 | 1485 | 561 | 127 | 11 |
| 12 | 4096 | 4020 | 335 | 3252 | 768 | 216 | 0 |
| 18 | 262144 | 261576 | 14532 | 230544 | 31032 | 7451 | 0 |
| 24 | 16777216 | 16772880 | 698870 | 15502080 | 1270800 | 286339 | 0 |

(Every `p ≤ 24` is in F2.)

* **Comparison with the mod-192 note.** The in-`Bad` column is exactly the note's K5 sequence `1, 0, 1, 2, 3, 6, 12, 16, 36, 60, 127, 216, 366, 721, 1290, 2095, 4227, 7451` for `p = 1..18`, recomputed by independent code. It extends to `14989, 27262, 46597, 93094, 168806, 286339` for `p = 19..24`.
  * These are the chain points whose period lies above the line: the rational skeleton of the incongruent set.
  * `log2(count)/p` is `0.690, 0.737, 0.755` at `p = 16, 20, 24`, creeping toward the PROVED exponent `h(log_3 2) = 0.950` of the word count `|Bad_p|` (choice ladder).
* **Integral points, all `p ≤ 24`.** Exactly the 18 points of the five known cycles, `{0}, {1,2}, {−1}, {−5,−7,−10}, {−17,…,−136}`, of periods 1, 2, 1, 3, 11. The positive ones sit on contracting words and the negative ones on expanding words (sign law).
* **So:**
  * the chain grows like `2^p`;
  * its `Bad` part grows like `2^(~0.95p)`;
  * its integral part is 18 points and does not grow;
  * Collatz's cycle half is exactly the statement that the chain meets `Z_{>0}` only in `{1,2}` (OPEN; FINITE-EXACT for `p ≤ 24`);
  * the cycle part of the five-cycle Z-Collatz statement (no integer cycles beyond the five known ones) is that the chain meets `Z` only in these 18 points.
* **Banach in `R` and in `Z_3` (F3).**
  * In `R` the inverse word has slope `2^p/3^a`, so it is a real contraction iff `3^a > 2^p` iff `x_w < 0` (sign law). Of the 2036 words of length `≤ 10` with `a ≥ 1`, the 411 contractions are exactly the negative points.
  * Backward iteration converges to `−1` (slope `2/3`), `−5` (`8/9`) and `−17` (`2048/2187`), and is repelled by `1`.
  * In `Z_3`, `D` is an isometry and `E` expands by 3. There is no contraction principle at the prime where integrality is decided.

**This is the correct fixed-point theorem for Collatz.** Banach in `Z_2` gives existence and uniqueness of one periodic point per word, for free. The whole problem is which of these rational points are integers (the gate congruence mod `2^p − 3^a`), and the orbit-level statements on top of that.

### 4.2 (b) Brouwer (the IVT), Sharkovskii and Chamberland's extension (PROVED unless marked; F4)

**Theorem K(b) (stability is decided by the sign).**
* Chamberland's `f(x) = (x/2)cos^2(πx/2) + ((3x+1)/2)sin^2(πx/2) = x + 1/4 − ((2x+1)/4)cos(πx)` satisfies `f(n) = T(n)` exactly (`cos πn = (−1)^n`; checked for `|n| ≤ 10^4`) and `f'(n) ∈ {1/2, 3/2}`.
* An integer `T`-cycle is an `f`-cycle with multiplier `3^a/2^p`. By the sign law it is **attracting iff its points are `≥ 0`**.
* The five integer cycles have multipliers `1/2` and `3/4` (attracting), and `3/2`, `9/8` and `2187/2048` (repelling).
* This agrees with Chamberland's claim (CITED) that a nontrivial positive-integer cycle would be attracting, and adds the converse on the negative side.

**Chamberland's cited numbers, re-derived.**
* Fixed points `μ1 = 0.277733766172`, `μ2 = 1.57737324444`, `μ3 = 2.44570769366`.
* Second attracting 2-cycle `A2 = {1.19253190705, 2.13865633552}`, multiplier `−0.2308`.
* `f([μ1, μ3]) ⊂ [μ1, μ3]` on a `4·10^5`-point grid.

These match Lagarias's entry 35 (`μ1 = 0.27773…`, `μ3 = 2.44570…`, `A2 = {1.19253…, 2.13865…}`). The negative-Schwarzian and homoclinic-point claims were **not** re-checked.

**Theorem K(c) (all periods on both sides, for every continuous extension).** Let `F: R → R` be continuous with `F = T` on `Z`.
* **Positive side.** `F(2) = 1`, `F(3) = 5`, `F(4) = 2`, so `F([2,3]) ⊃ [1,5]` and `F([3,4]) ⊃ [2,5]`, and both contain `[2,4]`: a turbulent pair.
  * For every word `u ∈ {0,1}^n` the IVT chain gives a periodic point with itinerary `u`. Choose nested compact pullbacks `K_k ⊂ J_(u_k)` with `F(K_k) = K_(k+1)` and `K_n = J_(u_0)`; then `F^n − id` changes sign on `K_0`.
  * No periodic orbit inside `[2,4]` contains `3`, since `F(3) = 5`. So the itinerary is determined by the orbit, and primitive words give exact period `n`.
  * Hence **every period occurs, and `#Fix(F^n) ∩ [2,4] ≥ 2^n`: the same `2^n` as the Banach count in `Z_2`.**
* **Negative side.** The same holds on `[−4,−3] ∪ [−3,−2]` (`F(−4) = −2`, `F(−3) = −4`, `F(−2) = −1`; no periodic orbit inside `[−4,−2]` contains `−3`). Also, the integer 3-cycle `{−5,−7,−10}` is a period-3 orbit of `F`, so Sharkovskii / Li–Yorke (CITED) give all periods and Li–Yorke chaos.
* F4 constructs one periodic orbit of Chamberland's `f` per admissible Lyndon itinerary, for periods 1–7:
  * `[2,3] ∪ [3,4]`: `2, 1, 2, 3, 6, 9, 18` orbits;
  * `[−10,−7] ∪ [−7,−5]` (golden mean): `1, 1, 1, 1, 2, 2, 4`;
  * `[−4,−3] ∪ [−3,−2]`: `2, 1, 2, 3, 6, 9, 18`.

  Each orbit is verified periodic to `10^−30`, with the correct itinerary and exact period. This re-proves the Borovkov–Pfeifer statement as reported in Chamberland's survey (any continuous extension of `T` has periodic orbits of every period).
* **Monotone divergent orbit** (Chamberland; LSW: a Cantor set of them). `F([2j+2, 2j+3]) ⊃ [2j+4, 2j+5]` for every `j ≥ 0`, so nested pullbacks give points whose orbit climbs through these intervals. F4 exhibits `x0 = 2.88652508471819` for `j ≤ 13`; compactness gives the infinite orbit. Divergence of a real extension happens between the integers and is not Collatz divergence.

**"The sign decides the Sharkovskii type", made precise.**
* As a statement about continuous extensions it is **REFUTED**: both half-lines carry every period (K(c)).
* What the sign does decide:
  * **the stability** of the integer cycles in Chamberland's `f`: PROVED (K(b)). In LSW's maps every integer cycle is superattracting, so there the sign decides nothing;
  * **which integer cycles exist.** The positive side has periods `{1 (the point 0), 2}` and the negative side `{1, 3, 11}`. This is FINITE-EXACT for `p ≤ 24` (F2). On the positive side it is, in general, the cycle half of Collatz (OPEN; CITED bounds in the [atlas](procgen_atlas_20260924_collatz_implication_atlas.md)).
* The negative side's integer period 3 forces all periods by Sharkovskii. The positive side's integer periods do not, but the positive side's turbulent integer triple `T(2), T(3), T(4)` does.

### 4.3 (c) Sperner, Brouwer, Rédei — ANALOGY (F5; details belong to the tournament lane)

* **Sperner** (1928; CITED by name). A Sperner labelling of a triangulated simplex has an *odd* number of fully labelled cells. The door-to-door argument pairs cells along paths in a graph of maximum degree 2. Knaster–Kuratowski–Mazurkiewicz (1929) derive Brouwer from it. F5 checks random 1D and 2D labellings.
* **Rédei** (1934; repo THM-001). Every tournament has an *odd* number of Hamiltonian paths. The repo's routes include involution proofs and the OCF, `H(T) = I(Ω(T), 2) ≡ 1 (mod 2)` (THM-002). F5 checks all 33867 labelled tournaments on `n ≤ 6`.
* **Common core:** `|S| ≡ |Fix(ι)| (mod 2)` for an involution `ι`, i.e. existence from parity.
* **Collatz.**
  * The per-word count is not a parity but exactly **one** (Banach, §4.1).
  * The real count per itinerary is at least one (IVT, §4.2).
  * The difficulty is integrality at the gate `(2^p − 3^a) | c_w`: a congruence modulo an odd number prime to 6, invisible to every mod-2 count.
  * The "fixed point chain growth" is the same `2^p` in `Z_2` (exact), in `R` (lower bound per turbulent pair) and for the itineraries. The integers see 18 points.

## 5. Verdicts on the owner's items

| item | verdict | exact content |
|---|---|---|
| the paper (arXiv:2502.20642v1) | **REFUTED** (proof); Collatz OPEN | swap gap (Theorem A); false Theorems 2.1(5)–2.3(5) (Theorem B); corrected version = one-step decay, violated at every odd `p ≥ 3` (Theorem D); SHEET-blind (Theorem E) |
| minimal counterexample | **PROVED**: 3 points | 2-cycles violate (5) (Proposition C); `{5,7,10}` under `3n−1` with the paper's constants |
| Lemma 2.1 = AM–QM | **REAL** | Proposition G; discards the cross term; never sharp on orbits |
| sign law = the two sides of the sandwich, `0` = centre | **REAL** (PROVED) | Theorem I; the cross-term sign on orbits is side-blind (Proposition H); the sign of the *points* is the side-aware datum |
| `b = 0` sheet = centre | **REAL** | midpoint/odd part; recentring `T_b(x)+b = (3/2)(x+b)`; shadowing identity; THM-4469 = confined carry around the central (Mahler) orbit |
| orthogonality in `C` | **ANALOGY** | Pythagorean locus `Re z = 0`; LSW's critical integers fold orthogonal into opposite; the Collatz content is only stability (Chamberland) |
| `{−1, 0}` | **REAL**, definitional | centres of the branches = Banach fixed points of `E`, `D`; pair 0 of THM-4470 |
| `\|x\|` as the precursor of polynomials | identities **REAL**; the phrase an **ANALOGY** | `\|T(x)\| = T_(sgn x)(\|x\|)`; linear sums, quadratic products; even observables are side-blind |
| Brouwer / fixed-point chain growth | **REAL** | Banach in `Z_2` (`2^p` per period, 18 integral points for `p ≤ 24`); the IVT gives `≥ 2^n` per turbulent pair on both sides |
| sign decides the Sharkovskii type | **REFUTED** for `f`; stability version **PROVED**; integer-cycle version FINITE-EXACT (`p ≤ 24`) and otherwise the OPEN cycle half | Theorem K(b), K(c) |
| Sperner / Brouwer / Rédei | **ANALOGY** | parity-of-fixed-points, `\|S\| ≡ \|Fix ι\|`; no mod-2 principle reaches the gate congruence |

## 6. What this adds to the session's picture

* **Foundry card: metric contraction, disposition REFUTED as a route.**
  * It is blind to SHEET, by Theorem E (its hypotheses hold verbatim for `3n−1`) and Proposition U.
  * It is broken locally by the up-step stretch: the odd branch's multiplier `3/2`, seen in consecutive distances.
  * Its unrefuted forms, Caristi and the discrete contraction metric, are equivalent to the conjecture. Bessaga's metric is equivalent to the cycle half.
* **The owner's sandwich locates the sign law exactly.**
  * The sign law is the choice between the two equality cases, with `0` as the middle point (Theorem I).
  * The metric sees only which equality case each orbit triple realizes (Proposition H), and that is side-blind. This is Corollary 7 of the mod-192 note in the owner's own coordinates.
* **The correct fixed-point theorem is Banach in `Z_2` for the inverse tree.** It is the mod-192 note's `D/E` tree read as a contraction system: every word has exactly one periodic point, and integrality at the gates is the whole remaining content of the cycle half.

## 7. Sources read, and reproduction

**Read this session.**
* arXiv:2502.20642v1, HTML, in full.
* Lagarias, *The 3x+1 problem: an annotated bibliography (1963–1999)*, arXiv math/0309224v13, entries 35 and 115. A local cached copy from an earlier lane of this session was used: `scratch/procgen_transversality/lit_A/lag_bib1.txt`.
* Chamberland, *An update on the 3x+1 problem* (2003), §§6.4–6.6 (local cached copy, `…/lit_A/chamb.txt`).
* Wikipedia, "Banach fixed-point theorem", section *Converses* (Bessaga 1959).

**Not read.** The primaries of Chamberland 1996, LSW 1999 (publisher bot-check, not bypassed) and Borovkov–Pfeifer 2000, and the classical theorems named in the status. All requests used the generic User-Agent `Mozilla/5.0 (research; math-repo)` and no personal data.

**Checks.**
* A scratch audit (`scratch/procgen_fixedpt/indep_audit.py`, not a deliverable) re-typed the table from the paper independently and re-checked:
  * the WGP inequality at 300000 random large pairs for each sheet;
  * the swap census to `p ≤ 10^5`;
  * the `x+1` counterexample on random pairs up to `10^12`;
  * the THM-4469 cylinder classes;
  * the `Bad` counts for `p ≤ 16`.
* All agree.

```bash
cd <worktree>
sh 04-computation/experiments/collatz_procgen_20260924_fixedpt_run.sh      # writes the .out (about 34 s, peak 428 MB)
# or individually (stdout = sections, stderr = timing and memory):
python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_kawasaki.py    # K1-K12
python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_inequality.py  # I1-I7
python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_chains.py      # F1-F5
```
