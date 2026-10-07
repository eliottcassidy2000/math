---
id: HYP-9220
title: "M1: for Haar-almost every 2-adic integer y, the Collatz (Terras) orbits of y and y+1 merge; P(no merge by T) ~ 10.8 T^(-1/2); equivalently the Collatz orbit relation has index 1 in the orbit relation of Z[1/6] x| <2,3> on Z_2"
status: >
  OPEN. COMPUTER-ASSISTED lower bound P(y ~ y+1) >= 0.5861 (THM-4569, box value iteration). NUMERICAL: sqrt(T) q(T) ~ 10.8
  on [1.28e4, 2e5], tail exponent 0.485, median merge time 70 T-steps. PROVED equivalences (THM-4569 (7)).
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
