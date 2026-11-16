---
id: HYP-9832
title: "After THM-5282 and THM-5292: (a) is there a synchronised bank at least as long as its block count h (through block order 8 the longest are 4 at h = 9 and 6 at h = 45)? By the composition law two coprime banks with L >= 2h reproduce themselves at every scale; (b) for which N is K(A, B) = N solvable, K = 432AB + 174A + 303B + 45 - 30(A mod 2) - 15(B mod 2) - 90(J(A) + J(B))? The solvable N have about (X ln X)/432 presentations below X; (c) classify the polynomial r(s, c) and inversions tau for which Psi_r o tau is a polynomial map; the known Keller counterexamples should be exactly the solutions whose parameter c vanishes on the inverted hyperplane"
status: >
  OPEN (QUESTIONS), mac-mini-2026-11-15-tapes. (a) FINITE-EXACT data only (1,823,696 seeds on blocks of order
  3..8). (b) NUMERICAL: 16.9 per cent of the odd numbers below 10^7 are reached by the one bank, growing by about
  2.6 points per decade; no solvability criterion. (c) one instance known (THM-1300 = r = (c s^3 - 2s^2)/2); the
  planar analogue with constant coefficients is excluded (THM-5292 (E)).
  Chan-Pak and JC_2 remain OPEN.
session: mac-mini-2026-11-15-tapes
date: 2026-11-15
source: 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_two_tape_banks_jacobian_anatomy_20261115/scripts/cross_bank.c, bank_coverage.py,
    jacobian_levels.py
related:
  - THM-5282, THM-5292, HYP-9822, HYP-9343, THM-4625
---

# HYP-9832 - self-covering banks, the two-tape selector, fold families

Sources: [README](../results/chanpak_two_tape_banks_jacobian_anatomy_20261115/README.md) sections 4-5 and
[JACOBIAN](../results/chanpak_two_tape_banks_jacobian_anatomy_20261115/JACOBIAN.md) section 6.

**(a)** Find seeds `(C_i, X_i)`, `i = 0..L-1`, with equal `(h, G_e, G_o, a_o, a_e, b_o, b_e)`, consecutive cut
values and `L >= h`. Falsifier of usefulness: a proof that `L < h` always.

**(b)** Describe `{K(A, B)}`. A solution of `K(A, B) in {N - 94, N - 92, N - 90}` realises the odd number `N`.

**(c)** For which `(r, tau)` is `Psi_r o tau` polynomial, `Psi_r = (sw + r - s r_s, 2w - 2 r_s, c)`,
`tau = (a(x, y, z), 1/x, c(x, y, z))` with `c` vanishing on `x = 0`?
