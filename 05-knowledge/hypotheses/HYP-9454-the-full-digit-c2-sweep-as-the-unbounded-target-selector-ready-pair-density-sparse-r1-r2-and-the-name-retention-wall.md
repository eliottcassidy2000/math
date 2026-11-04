---
id: HYP-9454
title: "After THM-4862: (1) sweep conjecture - for every odd v >= 307 there is a strong tournament T with a sole-exit vertex of digit sigma and 3 h(T) + 2 sigma = v, of order at most log_2 v + c; with the inherited C2 law this implies every odd number >= 307 is a Hamilton count, i.e. the Chan-Pak completeness; equivalently the ready pairs (h(U) - e_r(U) + [P A_r](U), e_r(U)) over strong (U, r) meet every line 3H + 2 sigma = v; (2) [CORRECTED the same night] the two-branch spectra R1, R2 have dyadic-window densities 0.65-0.77 through [1024, 2048) (exact through orders 13 and 12) and do not decay; the lineage's two-digit closure criterion (every ternary descendant of 1429 in R1 n R2) is a supercritical-branching question, open either way; (3) the name-retaining all-strong tree dies at depth 4 (order 19) if the sigma = 2 ready floor keeps doubling per vertex (floor_2(n) ~ 0.29 * 2^n), so any unbounded selector must choose targets of order about log_2 v, below the birth count; (4) the all-strong +4 rearming of 12875 and 12877 at order 15 is OPEN: the +2 branches are closed (THM-4862 (e)), 12877 has a strong digit-2 realiser at order 12 (THM-4862 (h); closes its +4 branch in the weak-allowed language of PREPARATION_REDUCTION), 12875 is in neither computed spectrum, and the nearest order-15 digit-2 states found lie at H = 12871 and 12879; (5) the digit-2 ready classes are counted by A377857(n - 2) = (n - 2) F(2n - 9) - L(2n - 10) (eight terms through order 12)"
status: >
  OPEN (klein-2026-11-04). (1) is FINITE-EXACT on [307, 26551] with bases of order <= 10 (THM-4862 (d)); it is a
  density statement about a two-parameter family and no proof route is known; (2) is read off the exact spectra through orders
  13 (digit 1) and 12 (digit 2), THM-4862 (f), computed by the self-certifying block-lemma generators; (3) uses the cited floors 4841, 9545,
  18905 at orders 14-16 and the observed doubling; a proof would need floor_2(n) >= c * 2^n, which is open (the
  elementary bound from (a) is only H >= 2 f(n - 1) - 2 with Moon's strong floor f); (4) the order-15 all-strong version searched by constrained annealing, two- and three-flip neighbourhoods of the
  digit-1 witnesses, and hard-digit-2 walks; not found; the order-12 realiser of 12877 comes from the exact digit-2 spectrum;
  (5) is a conjecture from eight exact terms, the last predicted before it was computed.
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

## 2. The two-digit spectra: density about two thirds, not decaying (CORRECTED; limit OPEN)

`R₁` (digit 1 available) has dyadic-window densities `0.72, 0.66, 0.65, 0.77` on `[128,256), [256,512), [512,1024), [1024,2048)` and `R₂` `0.52, 0.65, 0.67, 0.76` (exact, THM-4862 (f)). The per-order densities quoted in the first version of this entry (1, 1, 0.64, 0.38, 0.33, 0.24, 0.22) measured one order at a time inside `[floor, max]` and were misread as decay; the title's "density tending to zero" is withdrawn. Under independence a two-digit tree is a supercritical branching process (mean offspring `ρ₁ + ρ₂ ≈ 1.4`), so the lineage's criterion "every ternary descendant of the seed lies in `R₁ ∩ R₂`" is plausible for some seeds and false for others; for 1429 the seed and its first children are ready for both digits at orders 10–13 (THM-4862 (h)). OPEN: do the window densities converge, and to what?

## 3. The name-retention wall (OPEN as an asymptotic law)

With `floor₂(n) ≈ 0.29 · 2^n` for `n ≥ 14` and tree nodes `H ≈ 1429 · 3^ℓ` on `11 + 2ℓ` names, the inequality `H ≥ floor₂` fails from depth 4; so an unbounded selector must take targets of order near `log₂ v` and discard names, which the integer spectrum permits. A proof needs an exponential lower bound on ready floors.

## 4. The +4 rearmings at order 15

The +2 branches are closed by THM-4862 (e), and the +4 branch of 12877 has a strong digit-2 realiser at order 12 (THM-4862 (h)), which closes it in the weak-allowed language of PREPARATION_REDUCTION. For the all-strong +4 branches at order 15, digit-2 states at order 15 exist one and two arc flips away from the 12875 digit-1 witness (`H = 14925, 15273, 12467`), and hard-digit-2 walks move freely among digit-2 states, so the obstruction, if any, is arithmetic (which H occur among digit-2 strong tournaments of order 15 near 12,875) rather than structural; the weak-allowed alternative (a smaller-order digit-2 witness padded by a transitive prefix, PREPARATION_REDUCTION) is equally open, since digit 2 forces near-uniquely-Hamiltonian cores (THM-4862 (b)) whose counts are sparse.

## 5. A second Fibonacci law (CONJECTURE)

The digit-2 ready classes are counted by `A377857(n − 2) = (n − 2) F(2n − 9) − L(2n − 10)`: `1, 5, 18, 60, 191, 589, 1775, 5257` at orders 5–12, the last term predicted from the first seven and then computed (THM-4862 (g)). Through the bijection of THM-4862 (b) this counts strong sole-exit tournaments with exactly two Hamiltonian cycles. A proof would run through the block lemma (PROOFS.md §6) the way Garey's proof of `F(2n − 6)` runs through Douglas's characterisation of the unique-cycle tournaments; the generating function `x³(1 − x)²(1 + x)/(1 − 3x + x²)²` suggests a two-block decomposition with a Fibonacci-convolution count.
