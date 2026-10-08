---
id: HYP-9240
title: "+1 barrier for every deletion depth: for every D >= 1 the 2-adic limit chain at the trivial cycle (source x* = 1, child y* = (2 - 3^D)/3^D, pair-chain state (D, 3^D - 1), Terras shift D) never absorbs; equivalently no (generalized) deletion child's class-decided collision with a residual first-reset-2 source completes by Terras time j after the run end, so every K-uniform collision certificate on unrefined 2-adic classes has depth >= the two-run length"
status: >
  OPEN. FINITE-EXACT for D <= 8000 (THM-4601 (ii); D <= 3000 in twoanchor_core.py, 3001-8000 in barrier_extend.py, both reproduced by
  audit A): D = 1, 2 go to 0; D >= 3 land on positive integers N_D (<= 880 for D <= 3000, <= 2527 for D <= 8000) and enter the trivial
  cycle with debt k_inf(D) >= 3 (in phase / out of phase 2288 / 710 for D <= 3000, 3610 / 1390 for 3001..8000). NUMERICAL: in-phase
  k_inf/D in [0.75, 1.47] for 3 <= D <= 3000, [0.90, 1.10] for 163 <= D <= 3000. HEURISTIC for all D: the landing debt is about D
  (s_0 in [1.79D, 2.22D] for D >= 100) while the odd-step excess of N <= 880 is at most 9 (N = 871). Corrected after audit A (MISTAKE-586).
source: mac-mini-2026-10-07-twoanchor, 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 3)
related:
  - 01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md
  - 01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md (the sign barrier at -1; this is its +1 twin)
scripts:
  - 04-computation/experiments/twoanchor_20261007/twoanchor_core.py (part B; argument DMAX)
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
