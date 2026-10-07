---
id: HYP-9240
title: "+1 barrier for every deletion depth: for every D >= 1 the 2-adic limit chain at the trivial cycle (source x* = 1, child y* = (2 - 3^D)/3^D, pair-chain state (D, 3^D - 1), Terras shift D) never absorbs; equivalently no deletion child's class-decided collision with a residual first-reset-2 source completes inside the source's two-run, and K-uniform certificates have depth >= the two-run length"
status: >
  OPEN. FINITE-EXACT for D <= 3000 (THM-4601 (ii)): D = 1, 2 go to 0; D >= 3 land on positive integers N_D <= 880 and enter the trivial cycle
  with debt k_inf(D) >= 3, k_inf/D -> 1 (2288 in phase, 710 out of phase). HEURISTIC for all D: the landing debt ceil(s_0/2) >= D/2 would have to
  be repaid by the odd-step excess of the orbit of a small positive integer N_D.
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
* For every D, the deletion rule `n ⇝ (n+1)/2^D − 1` is useless inside the two-run of a residual source.
* Hence every K-uniform certificate for that branch has depth at least the two-run length.
* The two-anchor rule of THM-4601 attains that depth plus O(1).

**Reduction.**
* The orbit of y* makes exactly D denominator-clearing odd steps (`a/3^d ↦ (a + 3^(d−1))/(2·3^(d−1))` when odd). It lands, at time s_0 ≥ D, on an integer `N_D ≥ 0`. This holds because `y* + 1 = 2/3^D > 0` and `w ↦ 3w/2`, `w ↦ (w+1)/2` preserve positivity of `w = y + 1`.
* The landing debt is `⌈s_0/2⌉`. Absorption would then need the orbit of N_D to reach 1 in phase with odd-step excess exactly that debt.
* So HYP-9240 follows from a bound such as "`excess(N_D) < ⌈s_0(D)/2⌉`", e.g. from a bound on `N_D`. The excess of N is #odd steps − σ_T(N)/2.

**Evidence.**
* For D ≤ 3000 every landing satisfies `N_D ≤ 880`, so the odd-step excess is at most about 20, while the required excess grows like D/2.
* D = 3, 4: debt 3.
* The in-phase debts are within 0.75D–1.1D.
