---
id: HYP-9220
title: "M1: for Haar-almost every 2-adic integer y, the Collatz (Terras) orbits of y and y+1 merge; P(no merge by T) ~ 10.8 T^(-1/2); equivalently the Collatz orbit relation has index 1 in the orbit relation of Z[1/6] x| <2,3> on Z_2"
status: >
  RESOLVED 2026-10-07: PROVED (THM-4581 (3), (6a); the pair chain from (0, 1) is absorbed almost surely; Lyapunov weight
  s^|k| with return drift rho = 0.634 at theta = 1/2). The rate is PROVED up to logarithms: P(no merge by T) is between
  c T^(-1/2) and C T^(-1/2) (log T)^2 (THM-4581 (4)). The constant 10.8-11.0 is NUMERICAL, explained HEURISTICALLY by one big
  jump: E[J] sqrt(4/pi) = 11.10 with E[J] = 9.83 excursions (THM-4581 (7)). Earlier: COMPUTER-ASSISTED lower bound 0.5861
  (THM-4569).
source: mac-mini-2026-10-07-oaimath3 (groups reader), 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
related:
  - 01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md
  - 05-knowledge/hypotheses/HYP-9217-mersenne-debt-recurrence-law-t-minus-half.md
---

# HYP-9220 — y and y+1 merge almost surely

**Statement.** For Haar-almost every `y ∈ Z_2`, `T^m(y) = T^m(y + 1)` for some `m`. Moreover `P(no merge by T) = (c + o(1)) T^(−1/2)` with `c ≈ 10.8`.

**Equivalent forms (THM-4569 (7)).**
* The index `[R_A : R_C] = 1`.
* Every affine relation in `Γ_C` merges almost everywhere.
* The Terras-clock chain started from `(j, c) = (0, 1)` is box-recurrent.

**Consequences.** HYP-9217 (1) (almost-sure part), `μ_2(S) = 1`, and HYP-9213.

**What is known.**
* Recurrence of the odd-step difference is unconditional (THM-4569 (3)).
* What remains is archimedean: `|c|` must not blow up along the returns of `j`.
* Note that this does not touch Collatz for actual integers: almost every 2-adic statement is a measure statement.

---

## Update (2026-10-07, mac-mini-2026-10-07-oaimath3): PROVED (THM-4581)

* **Proof idea.** `k` (the odd-step difference) returns to 0 infinitely often (THM-4569 (3)). Along returns, `E|e|^(1/2)` contracts by 0.634 plus a constant.
  * The weight `|f|^θ s^|k|` with `s = 2^θ(1 − √(1 − (3/4)^θ))` is a martingale on flips. Moving toward `k = 0` always costs the factor 3/2, and moving away pays 1/2.
  * Runs at level `h` cost a fresh fair coin per continuation beyond `v_2(3^h − 1) − 1` steps.
  * So `|e|` at returns is tight, every state at `k = 0` reaches `(0, 0)` with positive probability, and Lévy's 0–1 law finishes.
* **Consequences.**
  * The index `[R_A : R_C] = 1`.
  * Every `u = 3^k y + e` with `e ∈ Z[1/3]` merges almost everywhere (for example `y` and `3y`).
  * The integers `n` whose trajectories meet `n + 1`'s at equal Terras time have natural density 1, with exceptional residue fraction `q(K) ≤ C K^(−1/2) (log K)^2` mod `2^K`.
* **Still open.** The exact constant (regular variation) of `q(T)`.
