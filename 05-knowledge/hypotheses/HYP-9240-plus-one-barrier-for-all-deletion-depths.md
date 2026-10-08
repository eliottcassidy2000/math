---
id: HYP-9240
title: "+1 barrier for every deletion depth: for every D >= 1 the 2-adic limit chain at the trivial cycle (source x* = 1, child y* = (2 - 3^D)/3^D, pair-chain state (D, 3^D - 1), Terras shift D) never absorbs; equivalently no (generalized) deletion child's class-decided collision with a residual first-reset-2 source completes by Terras time j after the run end, so every K-uniform collision certificate on unrefined 2-adic classes has depth >= the two-run length"
status: >
  OPEN. FINITE-EXACT for D <= 8000 (THM-4601 (ii); D <= 3000 in twoanchor_core.py, 3001-8000 in barrier_extend.py, both reproduced by
  audit A): D = 1, 2 go to 0; D >= 3 land on positive integers N_D (<= 880 for D <= 3000, <= 2527 for D <= 8000) and enter the trivial
  cycle with debt k_inf(D) >= 3 (in phase / out of phase 2288 / 710 for D <= 3000, 3610 / 1390 for 3001..8000). NUMERICAL: in-phase
  k_inf/D in [0.75, 1.47] for 3 <= D <= 3000, [0.90, 1.10] for 163 <= D <= 3000. HEURISTIC for all D: the landing debt is about D
  (s_0 in [1.79D, 2.22D] for D >= 100) while the odd-step excess of N <= 880 is at most 9 (N = 871). Corrected after audit A (MISTAKE-586).
  ADDENDUM 2026-10-08 (mac-mini-2026-10-08-reframes): the children satisfy y*_(D-1) = 3 y*_D + 2 at equal times (the Mersenne
  lag-1 relation), so they coalesce in rivers on which the landing (N_D, s_0) and the entry debt are constant (PROVED for merged
  pairs); FINITE-EXACT: 195 rivers for 3 <= D < 4500, 3929 of 4496 consecutive pairs share the landing, largest river 144 depths;
  per river the landing has the stationary Kesten-Goldie tail P(N > x) ~ 2.87/x of the level recursion (NUMERICAL), so the
  apparent scarcity of large N_D is a sample-size effect of the rivers, and the barrier is decided river by river; consecutive
  depths share their landing for a density-one set of D (PROVED from THM-4581 3 by equidistribution of 3^-D in Z_2).
source: mac-mini-2026-10-07-twoanchor, 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 3)
related:
  - 01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md
  - 01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md (the sign barrier at -1; this is its +1 twin)
scripts:
  - 04-computation/experiments/twoanchor_20261007/twoanchor_core.py (part B; argument DMAX)
  - 04-computation/experiments/reframes_20261007/landing_map.py, landing_windows.py, landing_rivers.py, landing_growth.py (+ .out) (2026-10-08 addendum)
---

# HYP-9240 — the +1 barrier for every deletion depth

**Statement.**
* For every D ≥ 1, the rational Terras orbits of `x* = 1` and `y* = (2 − 3^D)/3^D` never coincide at equal time with zero debt.
* Here the debt is `k_s = D + ⌈s/2⌉ − #odd steps of y* before time s`.

**What it would give.**
* For every D, the deletion rule `n ⇝ (n+1)/2^D − 1`, and every generalized deletion `3^a(n+1)/2^b − 1` with b − a = D, completes no earlier than Terras time j + 1 after the run end.
* Hence every K-uniform collision certificate on unrefined 2-adic classes (no condition on t mod 3) has depth at least the two-run length.
* 3-adic refinements escape the barrier. For example, `(2n−1)/3` is a depth-0 certificate when 3 | t.
* The two-anchor rule of THM-4601 attains that depth plus O(1).

**Reduction.**
* The orbit of y* makes exactly D denominator-clearing odd steps (`a/3^d ↦ (a + 3^(d−1))/(2·3^(d−1))` when odd). It lands, at time s_0 ≥ D, on an integer `N_D ≥ 0`. This holds because `y* + 1 = 2/3^D > 0` and `w ↦ 3w/2`, `w ↦ (w+1)/2` preserve positivity of `w = y + 1`.
* The landing debt is `⌈s_0/2⌉`. Absorption would then need the orbit of N_D to reach 1 in phase with odd-step excess exactly that debt.
* The exact entry debt is `k_entry = ⌈(s_0 + σ_T(N_D))/2⌉ − odd(N_D)`.
* So HYP-9240 follows from `excess(N_D) < s_0(D)/2`, equivalently `odd(N_D) < ⌈(s_0 + σ_T(N_D))/2⌉`. The excess of N is odd(N) − σ_T(N)/2.
* The weaker bound `excess < ⌈s_0/2⌉` is not sufficient: when s_0 and σ_T are both odd, excess = s_0/2 gives zero debt (audit A).

**Nature of the problem (remark, 2026-10-07 continuation; corrected after audit C).**
* The clearing word of `y* = −1 + 2·3^(−D)` is the Terras parity vector of `−1 + 2·(3^(−D) mod 2^(s_0+1))`. So the landing `N_D = (2 − 3^D + B_word)/2^(s_0)` and the entry debt are functions of D and of the first s_0 binary digits of 3^(−D) (checked for D ≤ 3000).
* Real-size bounds give only `N_D + 1 < (3/2)^D`: with `w = y + 1`, `max(w, 1)` never grows under `w ↦ (w+1)/2` and grows by at most 3/2 per odd step. Observed landings are far smaller: `log(N_D + 1)/log((3/2)^D) ≤ 0.455`.
* An all-D proof needs 2-adic control of the orbit of `−1 + 2·3^(−D)` up to its clearing time s_0 (1.79D–2.22D for 100 ≤ D ≤ 3000), i.e. of a nonlinear function of the first s_0 binary digits of 3^(−D). We know no method. The zeroless-power problems (THM-4580) are an analogy, not a reduction.
* The finite range D ≤ 8000 is settled exactly.

**Evidence.**
* For D ≤ 3000 every landing satisfies `N_D ≤ 880`, and for D ≤ 8000 `N_D ≤ 2527`.
* The odd-step excess of every N ≤ 880 is at most 9 (attained at N = 871), while the required excess is about D.
* D = 3, 4: debt 3.
* The in-phase debts lie within [0.75D, 1.47D] for D ≤ 3000, and within [0.90D, 1.10D] for 163 ≤ D ≤ 3000.

## Addendum 2026-10-08: the children form their own Mersenne line

**Level recursion (PROVED, elementary).**
* Write the child at 3-adic level m as `a_m/3^m`, with `a_D = 2 − 3^D`.
* One level is: `v_2(a_m)` halvings, then one denominator-clearing odd step. This gives

      a_(m−1) = (3^(m−1) + odd(a_m))/2,      N_D = a_0,

  where `odd(a)` is the odd part of `a`, with its sign.
* In real terms, `b_m = a_m/3^m` obeys `b_(m−1) = 1/2 + (3/2) 2^(−v_2(a_m)) b_m`. This is a perpetuity `X = B + AX` with `A = (3/2)2^(−v)`.
* If `v` is fair-geometric, `E[A] = 1`, so `κ = 1`, and the Kesten–Goldie constant is `C = E[B]/E[A ln A] = 0.5/0.1744 = 2.87`.
  * This is 3/2 times the per-step constant 1.91, because a landing is observed right after an odd step.

**Rivers (PROVED structure + FINITE-EXACT counts).**
* `y*_(D−1) + 1 = 3(y*_D + 1)`, i.e. `y*_(D−1) = 3y*_D + 2`, at equal times. This is the pair-chain state `(1, 2)` of the Mersenne line.
* When two children merge at equal time, they land together: same `N_D` and same `s_0`. Their barrier debts also coincide, because the merge has zero relative debt.
* So the landing data and the entry debt `⌈(s_0 + σ_T(N_D))/2⌉ − odd(N_D)` are constant on the river.
* **Density one (PROVED, corollary of THM-4581 3).**
  * The pair `(y*_(D−1), y*_D)` is the pair chain from `(1, 2)`, driven by the parities of `v = 2·3^(−D) − 1`.
  * As D runs over a period of `2^(m−2)`, `3^(−D) mod 2^m` runs once over `⟨3⟩`. So `v mod 2^(m+1)` is uniform on the coset `2⟨3⟩ − 1`.
  * Absorption by time τ depends only on `v mod 2^(τ+1)`. Hence the density of D absorbed by τ equals the conditional Haar probability, which tends to 1 (THM-4581 3); this is the transfer of THM-4581 6(d).
  * Since `s_0 ≥ D → ∞`, a merge by a fixed τ happens before landing.
  * So consecutive depths share `(N_D, s_0)` for a set of D of natural density 1. The number of rivers below X is `o(X)`.
* For `3 ≤ D < 4500` there are 195 distinct landings `(N_D, s_0)`.
  * 3929 of the 4496 consecutive pairs share one.
  * The largest rivers have 144, 126, 124, 114, 114 depths.
  * The record landings form rivers: `N = 880` for D in [253, 268] (8 depths), `732` for [815, 852] (22 depths), `412` for {43, 44}, `327` for [129, 134], `304` for [641, 648], `217` for [3287, 3338].
  * The key must be `(N_D, s_0)`, not `(N_D, s_0 − D)`: the relation holds at equal times, unlike the Mersenne line's shifted runs.

**Tail per river (NUMERICAL).**
* Counted per river, `P(N > 4, 16, 64, 256) = 0.467, 0.221, 0.067, 0.026`, against `2.87/x = 0.72, 0.18, 0.045, 0.011`.
  * The agreement is loose: 5 rivers lie above 256 against 2.2 expected (audit F).
  * There are only 40 distinct `N_D` values but 195 distinct landings `(N_D, s_0)`. Audit F reproduced every count exactly by a different method (2-adic truncation plus the parity-vector formula).
* Counted per depth, the tail looked truncated: none above 1024 for `D ≤ 4500`, and windows of 500 depths with maxima 4–217. That is a sample-size effect, since depths in one river are one sample.
* The exhaustive landing from level m over all `c = 3^m(v+1) ∈ [1, 2·3^m]` has maximum `≈ 2(3/2)^m`. The ratios are 1.50 ± 0.02 for `m = 3..11`, the all-expansion path. Its tail is close to the i.i.d. model's (landing_map.out).

**Consequence.**
* HYP-9240 is decided river by river: one check per river (195 rivers for `D < 4500`).
* An all-D proof would have to control how often the children's rivers break and where they land. That is the Mersenne-line coalescence problem (THM-4581, HYP-9242/9243) transported to the rational line `2/3^D − 1`.

