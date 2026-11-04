---
id: HYP-9454
title: "After THM-4862: (1) sweep conjecture - for every odd v >= 307 there is a strong tournament T with a sole-exit vertex of digit sigma and 3 h(T) + 2 sigma = v, of order at most log_2 v + c; with the inherited C2 law this implies every odd number >= 307 is a Hamilton count, i.e. the Chan-Pak completeness; equivalently the ready pairs (h(U) - e_r(U) + [P A_r](U), e_r(U)) over strong (U, r) meet every line 3H + 2 sigma = v; (2) the two-branch spectra R1, R2 have density tending to zero within their order windows (per-order density 1, 1, 0.64, 0.38, 0.33, 0.24, 0.22 for sigma = 1 at orders 4..10), so the lineage's two-digit closure criterion (every ternary descendant of 1429 in R1 n R2) fails at some depth unless the descendants are special; (3) the name-retaining all-strong tree dies at depth 4 (order 19) if the sigma = 2 ready floor keeps doubling per vertex (floor_2(n) ~ 0.29 * 2^n), so any unbounded selector must choose targets of order about log_2 v, below the birth count; (4) the +4 rearming of 12875 and 12877 at order 15 is OPEN: the +2 branches are closed (THM-4862 (e)), the +4 branches need a strong order-15 tournament with a sole-exit vertex of digit 2 and H = 12875 or 12877, and the nearest digit-2 states found lie at H = 12879 and 12879"
status: >
  OPEN (klein-2026-11-04). (1) is FINITE-EXACT on [307, 26551] with bases of order <= 10 (THM-4862 (d)); it is a
  density statement about a two-parameter family and no proof route is known; (2) is read off complete censuses through
  order 10 and the Fibonacci count F(2n - 6) of the sigma = 1 classes (Garey 1972 via THM-4862 (b)): the number of sigma = 1-ready
  classes with floor <= X is about (X/0.4)^1.39, enough in count, but their H-values spread over [0.4 * 2^n, 3^n] and
  the fraction of odd numbers in [floor, max] attained falls with the order; (3) uses the cited floors 4841, 9545,
  18905 at orders 14-16 and the observed doubling; a proof would need floor_2(n) >= c * 2^n, which is open (the
  elementary bound from (a) is only H >= 2 f(n - 1) - 2 with Moon's strong floor f); (4) searched by constrained annealing, two- and three-flip neighbourhoods of the digit-1 witnesses, and hard-digit-2 walks; not found.
session: klein-2026-11-04 (codex/session-spectrum-20261104)
source: 05-knowledge/results/two_branch_rearming_sigma_sweep_20261104/README.md
related:
  - 01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md
  - 05-knowledge/results/nonzero_preparation_20261103/algebra/PREPARATION_REDUCTION.md
  - 05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md
---

# HYP-9454 — the full-digit C2 sweep as the unbounded target selector

## 1. The sweep conjecture (OPEN; implies Chan–Pak)

For every odd `v ≥ 307` there is a strong tournament `T` with a sole-exit vertex `a → r` such that `3h(T) + 2e_r(T − a) = v`. By the inherited C2 law the child is a strong tournament with `H = v`, so the conjecture plus the finite base (every odd number below 307 other than 7 and 21 is attained, THM-1370/THM-4744) gives the full odd spectrum. The selector is defined at every scale by one rule: pick a ready pair on the line `3H + 2σ = v`. Exact through 26,551 with bases of order ≤ 10; the digit sets grow like `0.29 H` in size, so the expected number of pairs on a line of height `v` grows linearly in `v`.

## 2. Sparsity of the two-digit spectra (OPEN as a limit statement)

`R₁` (digit 1 available) has per-order densities 1, 1, 0.64, 0.38, 0.33, 0.24, 0.22 at orders 4–10 inside `[floor, max]`, and `R₂` 1, 0.38, 0.50, 0.40, 0.39, 0.34 at orders 5–10. If these tend to zero, the lineage's criterion "every ternary descendant of the seed lies in `R₁ ∩ R₂`" fails at some depth for a generic seed, and the two-branch grammar cannot be unbounded without the full digit set. The particular descendants of 1429 are not shown to fail: the level-3 values have digit-2 realisers (lineage `surgery/TARGETS.md`), and the digit-1 question at 38627 … 38635 is open.

## 3. The name-retention wall (OPEN as an asymptotic law)

With `floor₂(n) ≈ 0.29 · 2^n` for `n ≥ 14` and tree nodes `H ≈ 1429 · 3^ℓ` on `11 + 2ℓ` names, the inequality `H ≥ floor₂` fails from depth 4; so an unbounded selector must take targets of order near `log₂ v` and discard names, which the integer spectrum permits. A proof needs an exponential lower bound on ready floors.

## 4. The +4 rearmings at order 15

The +2 branches are closed by THM-4862 (e). For the +4 branches, digit-2 states at order 15 exist one and two arc flips away from the 12875 digit-1 witness (`H = 14925, 15273, 12467`), and hard-digit-2 walks move freely among digit-2 states, so the obstruction, if any, is arithmetic (which H occur among digit-2 strong tournaments of order 15 near 12,875) rather than structural; the weak-allowed alternative (a smaller-order digit-2 witness padded by a transitive prefix, PREPARATION_REDUCTION) is equally open, since digit 2 forces near-uniquely-Hamiltonian cores (THM-4862 (b)) whose counts are sparse.
