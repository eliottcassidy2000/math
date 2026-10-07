# Audit A — THM-4581 (Haar coalescence) and its corollaries

Independent adversarial audit, 2026-10-07, of worktree `math-wt-chessboard-20261006` (read-only; no repository file edited, no state-changing git command).

All code is my own: Python with exact integers, plus a GMP C simulator. Code and outputs are in `audit/A_coalescence_scripts/` (`a1`–`a12`, `.out` files).

**Overall verdict.**
* **Statement 3 (almost-sure coalescence) is CORRECT**, and so are statements 1, 2 and 5.
* Statement 4 is sound as a sketch. It should stay typed "PROVED at sketch level".
* The corollaries 6(a)–(f) hold. (b), (e) and (f) need wording or gap fixes, and (c) a typing fix; none changes a conclusion.
* The errors are in surrounding claims:
  * "equal odd-step counts";
  * "Γ_C merges" read as equal-time merging;
  * "template total = Terras merge time";
  * E[J] ≈ 14.8 in statement 7;
  * the rate typed PROVED (not sketch) in several hypothesis files;
  * a false-looking any-lag lower bound in THM-4556 Update 2.

---

## 1. Chain table, parity, measurability, integrality — CONFIRMED WITH CORRECTION

**Derivation.** I re-derived the four branches. For example, `(σ, β) = (1, 1)`: `u` is even and `v` is odd, so `u' = u/2` and `3^k v = 3^(k−1)(2v' − 1)`, which gives `(k − 1, (e − 3^(k−1))/2)`.

**Integrality.** In `N = 3^max(0,−k) e` coordinates every branch has an even numerator, so `N` stays an integer. At `k = 0`, `e = N`.

**Measurability.** `(k_n, e_n)` is a function of `β_0, …, β_(n−1)`. The coin `c_n` is `β_n` XOR an `F_n`-measurable bit, so it is fair and independent of `F_n`.

**Machine check (`a1`).**
* Setup: 15 starts, including `(−2, 5/9)`, `(−3, −13/27)` and `(5, 242)`, against direct 2-adic orbits (1200-bit precision). In total 600 trials and 660,000 steps.
* Mismatches: 0 in each of the following.
  * `u_n ≡ 3^(k_n) v_n + e_n (mod 2^(P−n))`.
  * Parity `β ⊕ σ`.
  * Integrality.
  * `(0, 0) ⇒ u_n = v_n`.
  * The theorem's rational table equals my integer table.
  * **`k_n − k_0 = O_u(n) − O_v(n)`** (the difference of odd-step counts).

**Correction.** By that last identity, at absorption `O_u − O_v = −k_0`. So "a merge at equal Terras time **with equal odd-step counts**" (title and Setting, bullet 3) is false when `k_0 ≠ 0`.
* Witness (`a9`): start `(1, 0)` with `n = 15`. `T^6(45) = T^6(15) = 20`, after 3 and 4 odd steps respectively.
* Replace with: "at equal Terras time; `u` has made exactly `k_0` fewer odd steps than `v` (equal counts iff `k_0 = 0`)". Statement 5(iii) and the proof of 6(d) already use this correct form.

## 2. Lemma 1 (a)–(d) — CONFIRMED (one wording fix)

**Machine check (`a2`).** 700,048 (state, bit) pairs: all `|k| ≤ 12`, `|N| ≤ 3000`, plus 200,000 random large states (`|k| ≤ 60`, `|N|` up to `2^212`, and run-regime edge states `N = (3^h − 1)·r`). Failures: 0 for each of:
* the size bound (a), including the departure exception `|f|/2 + 1/6`;
* the sharper away-bound `|f'| ≤ |f|/2 + 1/18`;
* direction (b);
* the exact `k = 0` factors (c);
* every valuation transition of (d), including `e = 0` and `k < 0`;
* the one-step drift at `θ = 1/2`: `E[V'] ≤ V + ε_h` at flips, and `≤ 0.966V + ε_h/2` at runs. The worst relative slack is `4·10^(−16)`, i.e. exact equality in the multiplicative part.

In every state the map `β ↦ c` is a bijection, so `c` is fair.

**Wording fix.** "At odd `h` it is exactly `R ~ Geom(1/2)`" holds only for nonempty runs (Geom on `{1, 2, …}`). A run is empty when the arriving `e` is odd, and that event is past-measurable: at `h = 1` there are 2.03 steps per visit, so about half the runs are empty. The proof uses only `R ≤ m_h − 1 + G`, which is correct.

## 3. Statement 2 (return drift) — CONFIRMED

**The constants.**
* `s` solves `½((3/2)^θ s^(−1) + 2^(−θ) s) = 1`, so `s = 2^θ(1 − √(1 − (3/4)^θ))`.
* `s(θ) < 1` on all of `(0, 1)`, tending to 1 only at the endpoints (0.99997 at `θ = 0.9999`), and `ρ(1/2) = 0.633975`.
* The additive terms are `½·2^(−θ) s^(h−1) + ½·18^(−θ) s^(h+1) ≤ ε_h`.
* `C_exc = 58.35` and `C = 58.71` at `θ = 1/2` (`a3`).

**The skeleton.** Coins at flip times are fresh (optional skipping). For simple random walk killed at 0, `G(1, h) = 2`.
* Simulated over 20,000 excursions: 1.97–2.00 visits per level for `h = 1, …, 32`.
* Steps per visit: 1.90–2.03, within `m_h + 2`.

**Optional stopping.** The phrase "supermartingale up to additive errors" is rigorous via `M_n = V_(n∧τ) − Σ_(m<n∧τ) ε_(h_m)`, which is a true supermartingale. Then Fatou applies.
* `E[#steps at h] ≤ (m_h + 2)·E[#visits to h]` is legitimate, because the level of a visit is known at its start and the run bound is conditional on that start.

**No adversarial loophole.** Stays at `k = 0` carry no additive term (factor `0.966^(v_2(e))`), and each run continuation beyond `m_h − 1` steps costs a fresh coin. All bounds are one-step conditional, so correlations between run length and size are harmless.

**Numerics.** `E|e_ret|^(1/2) / |e|^(1/2)` is 0.578–0.631 for `|e|` from `10^3` to `2^200` (4000 excursions each), below 0.634. For small `|e|` the excess over `ρ|e|^(1/2)` is at most 1.32.

**Minor.** For `k_0 ≠ 0`, `E|e_(τ_0)|^θ < ∞` needs the same argument with `G(|k_0|, h) = 2 min(|k_0|, h)`. Statement 3(ii) should then read `max(E|e_(τ_0)|^θ, C/(1 − ρ))`.

## 4. Statement 3 (almost-sure absorption) — CONFIRMED

* **(i)** `k` is an `F_n`-martingale with increments `σ(1 − 2β) ∈ {−1, 0, 1}`, so Durrett's dichotomy applies.
  * If `k` converges, then eventually `σ ≡ 0`, the parity vectors agree, and `u_N = v_N`.
  * If moreover `k_N ≠ 0`, then `v_N` equals an `F_N`-measurable value. That has conditional probability 0, because `v_N` is Haar given `F_N`. A countable union over `N` finishes (i).
* **(ii)** Correct, using sup over `j` of `E Z_j` and Fatou.
* **(iii)** Machine check (`a4`): from every `(0, e)` with `0 < |e| ≤ 10^4`, an explicit coin string reaches `(0, 0)` using only the generic step.
  * Longest string: 33 steps, at most 12 excursions.
  * For `e > 0` the landing value is `(U(e) − 1)/2 ≤ (3e − 1)/4`, as claimed.
  * For `e < 0` the mirror path works: depart with `c = 1` to `k = −1`, halve with `c = 0`, return with `c = 1`. This gives `e ↦ −(U(|e|) − 1)/2`, a direct chain proof that can replace the translation argument (which is also valid).
* **(iv)** Lévy's 0–1 law plus the strong Markov property at `τ_j`.

## 5. Statement 4 (rate) — CONFIRMED as a sketch; keep typed PROVED at sketch level

**Lower bound.** Valid for every non-absorbed admissible start. Runs at `k = 0` never absorb, so a first passage of the skeleton from `|k| ≥ 1` to 0 is required.

**Upper bound.** The gaps are routine:
* **(R1)** The restart step needs `E[Z·1_fail] ≤ ρ^D K + C/(1 − ρ)`.
* **(R2)** `Σ_(i≤F) G_i` must be dominated by an i.i.d. geometric sum (true, since each continuation costs a fresh coin), then a union over `{F > t/12}`.
* **(R3)** For `k_0 ≠ 0`, add the tail of `τ_0`. The arithmetic `r^(3/2)(log T)^(1/2) T^(−1/2) = O((log T)² T^(−1/2))` checks.

**Numerics (`a5`, GMP, 400,000 paths from `(0, 1)`).**
* `P(J > r)` is geometric, with ratio 0.924 per excursion for `r ≤ 36`.
* `√t·P(X > t) = 1.10–1.12`, approaching `√(4/π) = 1.128`.
* `√T·q(T) = 10.98, 11.06, 11.11, 11.18` at `T = 2.56·10^4, …, 1.64·10^6`. No log growth is visible.

## 6. Corollaries

**(a) CONFIRMED.**

**(b) CONFIRMED WITH CORRECTION.**
* "Every element of `Γ_C` merges a.e." is true only in the grand-orbit sense, `T^m g(y) = T^n y`.
  * Under equal-time comparison, the 2-power of the multiplier `2^a 3^k` is invariant. For `a ≠ 0`, `T^n u = T^n v` would pin `v_n` to a single value, which has probability 0.
  * So `y ↦ 2y` a.s. never merges at equal time.
* THM-4569 has **no written proof** of (7). Supply one here. For `g(y) = 2^a 3^b y + c_0/(2^m 3^l)` with `M ≥ m` and `M + a ≥ 0`, set `y' = 2^(M+a) y`. Then
  * `y ~ y'`;
  * `g(y) ~ 2^M g(y) ~ 3^l·2^M g(y) = 3^(b+l) y' + 2^(M−m) c_0`, which is admissible;
  * apply statement 3 at `(l, 0)` and at `(b + l, 2^(M−m) c_0)`. Affine maps are nonsingular, so the a.e. statements compose.

**(c) CONFIRMED.** Density one is PROVED from statement 3.
* Own residue census (`a8`), `K = 8, 12, 16`, for both `n` vs `n + 1` and `n` vs `3n`: absorption by `K` is a function of `n mod 2^K` and equals the first integer equal-time merge, with 0 disagreements.
* Unmerged at `K = 16`: `40295/65536` (`n + 1`) and `50911/65536` (`3n`). The first matches the session.
* The bound `q(K) ≤ C K^(−1/2)(log K)²` comes from statement 4, so it is sketch level.

**(d) CONFIRMED.**
* `p = 3^D q_D + 3^D − 1`. Absorption implies `U^i(p) = U^(i+D)(q_D)`, because `O_q − O_p = D` at absorption. The converse holds off a Haar-null set (unequal totals pin `q_D`).
* Integer check (`a6`):

  | lag `D` | odd `a` | absorbed before 1 | disagreements with S19's switch |
  |---|---|---|---|
  | 1 | ≤ 1201 | 453 of 599 | 0 |
  | 3 | ≤ 801 | 264 of 398 | 0 |
  | 5 | ≤ 801 | 246 of 397 | 0 |

* `σ(M_a) = σ(M_(a−1))` in all 453 lag-1 cases.
* `q_1 ≡ 5 (mod 16)`, and the forced bits 1, 0, 0, 0 lead to the state `(2, 1)`.
* Lag 1 alone gives `μ_2(S) = 1`, since `S` contains the lag-1 set.
* Proposition 6 uses exactly `μ_2(S) = 1` and THM-4556 (iv), plus (ii) for even `a`.
* Suggest adding one line: since `p + 1 = 3^D(q_D + 1)`, an equal-total coincidence is exactly a THM-4555 collision `f_u(−1) = f_u'(−1)`. That is the definition of `S` in Proposition 6.

**(e) CONFIRMED WITH CORRECTION.**
* The chain identification `(n, n − 1)` from `(0, 1)` is right.
* The reset-2 condition matches the sampler: `n`'s reset exponent is 2 iff `v_2(3^r t − 1) ≥ 2`, checked for all `n < 2·10^5`.
* **Gap.** The sampler counts a merge only when the common odd value is not 1. "The merge value exceeds 1" is not enough.
  * Witness (`a7`): `n = 3465223915` (`B = 32`). The chain is absorbed at Terras time 77 at value 128; the sampler reports no merge. One 64-bit source behaves the same way.
* **Fix.** Use absorption by `K < B − 2`, and discard the sources where `T^t(n)` is a power of 2. That proportion is at most `(B + K)·2^(K+3−B)`, because `T^t` is injective on each class mod `2^t`. Also, the reset-2 set is a countable union of residue classes: split off `v_2(n + 1) > K − 3` (mass at most `2^(3−K)`).
* Agreement otherwise: 940/941, 1254/1255, 1783/1783 and 611/611 at `B = 32, 64, 256, 1024`. The conclusion `P_B → 1` stands.

**(f) CONFIRMED WITH CORRECTION.**
* "The template total is the Terras time of the merge" is false. In fact template total = Terras merge time + `v_2(merge value)`.
* Witness: `a = 13`, where `T^44(p) = T^44(q_1) = 5690` but the template total is 45.
* For odd `a ≤ 1201`, 246 of the 453 absorbed cases differ. The overshoot counts are 207 (0), 119 (1), 70 (2), 38 (3), and so on.
* The overshoot is Geom(1/2) and independent of the σ-field at the merge, so both bounds transfer. The upper bound is sketch level.

## 7. Statement 5 — CONFIRMED

* At `n = 40` (`a10`, 100,000 paths), `corr(β, β ⊕ σ) = 0.4427 = 1 − 2P(σ = 1)`.
* Correlations at distinct times are at most 0.004 in absolute value (standard error 0.003).

## 8. Statement 7 (one big jump) — CONFIRMED WITH CORRECTION

**Start `(0, 1)`.** `E[J]` counted up to `T` is 9.42, 9.69, 9.81 at `T = 10^5, 4·10^5, 1.6·10^6`. The increments halve with each ×4 in `T`, which extrapolates to **`E[J] ≈ 9.93`**. The session's 9.83 is biased low by truncation. Then `9.93·√(4/π) = 11.21`, against 11.06–11.18 measured. The heuristic is supported.

**Lag-1 Mersenne start** (post-prefix state `(2, 1)`, 150,000 paths).
* `√T·q = 16.41, 16.53, 16.35` at `T = 2.6·10^4, 10^5, 4·10^5`, consistent with S19's 16.7.
* But **`E[J]` (excursions from 0) ≈ 12.8, not 14.8.**
* For a start at level `|k_0|`, the first passage carries `|k_0|` times the tail of one excursion. So the constant is `(|k_0| + E[J])·√(4/π) = (2 + 12.8)·1.128 = 16.7`.

## 9. Prior art — nothing found

Six searches, plus a text check of Akin's paper, found no statement of almost-sure merging of `y` and `y + 1` in `Z_2`, or of density-one coalescence of `n` and `n + 1`. Related work:
* Garner (1985), and [arXiv:1511.09141](https://arxiv.org/abs/1511.09141) (counterexamples to Garner's conjecture).
* [Kontorovich–Lagarias](https://arxiv.org/abs/0910.1944): `8n+4` and `8n+5` coalesce.
* Burson: [arXiv:2005.09456](https://arxiv.org/abs/2005.09456) (withdrawn) and [arXiv:1906.10566](https://arxiv.org/abs/1906.10566).
* [Bernstein–Lagarias](https://dept.math.lsa.umich.edu/~lagarias/doc/bernstein.pdf), [Monks–Yazinski](https://monks.scranton.edu/files/pubs/autoconjabs.txt), [Akin](https://math.sci.ccny.cuny.edu/document/Akin2004_WhyThreeXisHard.pdf), [Lagarias bibliography II](https://arxiv.org/abs/math/0608208).

This is not a novelty certificate; the theorem's "not a claim of novelty" is right.

## 10. Typing and overclaims

* **Title.** It allows `e ∈ Z[1/3]`, but the proof assumes `3^max(0,−k) e ∈ Z`.
  * The general case is true. The 3-adic excess `max(0, −v_3(e) − max(0, −k))` never increases, and it drops with probability at least 1/2 at each step while positive, because one of the two coins always lowers it. Checked (`a11`): 6 inadmissible starts, 1800 paths, all admissible within 20 steps.
  * Fix: add this argument, or restrict the title.
* **THM-4581 status** types "statements 1–3, 5, 6" PROVED. But 6(c)'s `q(K)` bound and 6(f)(1)–(2) rest on statement 4, so they are sketch level.
* **HYP-9220** status ("rate PROVED up to logarithms") and its update's residue bound: should be sketch level. The constant "10.8" should read about 11.1–11.2.
* **HYP-9217** status and Update 2 call the exponent 1/2 PROVED; it is sketch level. The phrase "(= Terras time of the merge)" is wrong (see 6(f)).
* **HYP-9213** status ("α ≥ 1/2"): sketch level.
* **INDEX.md** (HYP-9217 and HYP-9220 lines), **CURRENT-FRONTIER**, the **THM-4569 update** ("sharp up to logarithms") and the **THM-4564 update** ("`T^(−1/2+o(1))` deficit"): each should add "at sketch level".
* **THM-4556 Update 2** says "the density of odd `a` without a switch of template total ≤ `K` is between `c K^(−1/2)` and …", and calls the certificates of (v) finite-`K` instances.
  * Statement (v) certifies switches at any lag. For any lag, only the upper bound follows from THM-4581.
  * The lower bound looks false (`a12`; 8000 Haar samples; odd `D ≤ 31`):

    | `T` | 200 | 6400 |
    |---|---|---|
    | any-lag `√T·q` | 6.1 | 3.4 |
    | lag-1 `√T·q₁` | 8.9 | 15.7 |

    The local exponent of the any-lag share is about 0.69, matching HYP-9217 (2).
  * Restrict the sentence to lag 1.

---

## Required corrections

1. Title and Setting: "odd-step counts differing by exactly `k_0`", not "equal" (§1).
2. Title: add admissibility, or the 3-adic-excess argument (§10).
3. 6(b): grand-orbit merging only; replace the THM-4569 (7) citation with the reduction in §6(b).
4. 6(e): common odd value ≠ 1, and the infinite residue union.
5. 6(f) and HYP-9217: template total = Terras merge time + `v_2(merge value)`.
6. Statement 7: `E[J] ≈ 12.8` for the lag-1 start, constant `(|k_0| + E[J])·√(4/π)`; `E[J](0, 1) ≈ 9.93`.
7. "At sketch level" on every downstream use of the rate (the files listed in §10).
8. THM-4556 Update 2: the `c K^(−1/2)` lower bound holds for lag 1 only.

## Optional improvements

Write out the compensated supermartingale, the `τ_0` moment, the mirror descent, and (R1)–(R3); doing so would upgrade statement 4 to PROVED. Also say "nonempty runs" in 1(d).

**Final verdict:** THM-4581 statement 3 is **correct**. Every admissible affine pair `u = 3^k y + e` merges at equal Terras time for Haar-almost every `y`, and so does every pair with `e ∈ Z[1/3]` by the reduction in §10. HYP-9220, HYP-9213 and HYP-9214, and `[R_A : R_C] = 1`, follow. They need the corrections above, none of which affects the conclusions.
