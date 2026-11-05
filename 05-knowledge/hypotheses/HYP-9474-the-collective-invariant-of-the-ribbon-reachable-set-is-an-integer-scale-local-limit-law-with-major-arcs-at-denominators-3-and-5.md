---
id: HYP-9474
title: "The collective invariant of the ribbon reachable set (after THM-4880 and THM-4881). (1) Local limit form: with r_n(N) = #{u : |u| = n, H(u) = N} over the 4^(n-1) span-2 ribbon words, r_n(N) = 4^(n-1) rho_n(N) (1 + o(1)) uniformly for N in a bulk window around exp(n lambda), where rho_n is the density of a log-normal law (Lyapunov exponent lambda > 0.44, variance sigma^2) times an arithmetic factor in N modulo powers of 3 and 5; the ribbon tiling law (HYP-9473 (3)) is the statement that this density is positive at integer scale on [c_n, d_n]. (2) Minor-arc decay: |S_n(theta)| / 4^(n-1) -> 0 geometrically for every theta not of the form a / (3^i 5^j), with the major-arc masses converging to the limit measure of H modulo 3^i 5^j. (3) What a proof needs: a spectral gap for the twisted transfer operators of the linear random walk y -> A_a y on the 13-dimensional integer states, restricted to the coset of THM-4881 (A), away from the major arcs. (4) Closed routes: quotient numerations and digit blocks (THM-4880 (G), (H)), arithmetic progressions ((H2)), dense response sets (THM-4881 (C)), dense extension sets."
status: >
  OPEN (mac-mini-2026-11-04-invariant). (1) NUMERICAL: multiplicities and residue frequencies at levels 9-13
  (collective.py); the log-normal shape is the standard random-matrix-product picture (Guivarc'h-Le Page;
  Benoist-Quint, CITED as context), not verified here beyond (1/n) E log H = 0.386 ... 0.443 for n = 8..13. (2)
  NUMERICAL to n = 13 at twelve frequencies (fourier.py): geometric decay at 1/7, 1/8, 1/16, 3/11, 5/13 and four
  irrationals; growth at 1/3, 1/5, 1/9. (3) is a reformulation (the standing hypotheses - irreducibility THM-4880
  (G), proximality THM-4881 (E) - hold; strong irreducibility is unverified). (4) PROVED or FINITE-EXACT as cited.
  Chan-Pak remains OPEN (MISTAKE-559).
session: mac-mini-2026-11-04-invariant
date: 2026-11-04
source: 05-knowledge/results/ribbon_reachable_set_collective_invariant_20261104.md
related: [THM-4881, THM-4880, HYP-9473, HYP-9453, THM-4094, THM-1370]
scripts: 04-computation/experiments/ribbon_spectrum_20261104/collective.py; fourier.py; bsets.py
---

# HYP-9474 - the collective invariant as an integer-scale local limit law

**Why distributional.** THM-4880 showed that no quotient of the 13-dimensional representation exists (G) and that no
short block acts as a digit or inserts an arithmetic progression (H). THM-4881 (C) shows that the response sets at a
fixed count are sparse and saturated, so the near-join induction of HYP-9453 cannot be the mechanism inside the
family, and THM-4880's extension data show the same for the three letters. Yet every odd number `<= 100,000` except
7, 21, 47 is a count (THM-4880). The covering is therefore collective: for a target `N` the representations come from
many different `(H(u), extension)` pairs, with multiplicity growing like `4^n / exp(n lambda)` on the bulk.

**(1) Local limit form.** Let `lambda = lim (1/n) E log H(u)` over uniformly random words and `sigma^2` the variance
rate. Conjecture: uniformly for `N` odd with `|log N - n lambda| <= C sqrt(n)`,
`r_n(N) = 4^(n-1) (2 / (N sigma sqrt(2 pi n))) exp(-(log N - n lambda)^2 / (2 sigma^2 n)) kappa(N) (1 + o(1))`, where
`kappa(N) > 0` depends only on `N` modulo `3^i 5^j` for bounded `i, j` (the arithmetic factor: at `n = 13` the class
frequencies deviate from uniform by 0.21 modulo 3 and 0.16 modulo 5, and the deviation grows). The ribbon tiling law is
the positivity of the main term at integer scale throughout `[c_n, d_n]`, with `c_n` and `d_n` outside the Gaussian
bulk handled by the nested levels.

**(2) Minor arcs.** `S_n(theta) = sum_{|u| = n} e(theta H(u))`. Conjecture: `|S_n(theta)| / 4^(n-1) <= C q(theta)^n` with
`q(theta) < 1` for every `theta` not of the form `a / (3^i 5^j)`, and `S_n(a / 3^i 5^j) / 4^(n-1)` converges to the
Fourier coefficient of the limit distribution of `H` modulo `3^i 5^j`. Observed at `n = 9, 11, 13`: `1/8`: 0.063,
0.039, 0.020; `0.1234...`: 0.062, 0.031, 0.003; `1/3`: 0.190, 0.265, 0.324.

**(3) What a proof needs.** A spectral gap for the twisted operators `(P_theta f)(y) = sum_a e(theta phi(A_a y))
f(A_a y)` on functions of the integer state (restricted to the coset of THM-4881 (A)), uniform away from the major
arcs - the integer-lattice analogue of the Guivarc'h-Le Page / Benoist-Quint local limit theorem for the norm cocycle
of products of random matrices, whose standing hypotheses (irreducibility, proximality) the family satisfies. The
major arcs at `3^i 5^j` come from the multiplicative closure (ordinal sums of 3-cycles and strong 4-tournaments).

**(4) Numerical programme.** (a) decay rates `q(theta)` on a frequency grid for `n <= 16` with the deduplicated level
enumerator (states carry multiplicities); (b) `lambda` and `sigma^2`; (c) the limit measures of `H` modulo 9, 27, 25;
(d) the same for span 3 (rank 62), where 47 is attained.
