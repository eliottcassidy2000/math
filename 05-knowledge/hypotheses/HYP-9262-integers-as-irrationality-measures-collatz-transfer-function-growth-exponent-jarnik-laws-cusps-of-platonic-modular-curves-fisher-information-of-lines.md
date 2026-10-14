---
id: HYP-9262
title: "The positive integers as irrationality measures: the boundary growth exponent of the Collatz transfer function is 2(mu(log2 3) - 1) (THM-4632), so the integer n labels the laws of exponent 2n - 2 (n = 1 the rational parameters with poles on the line, n = 2 the generic class where the Collatz law is expected to sit, n >= 3 the Jarnik continua of dimension 2/n); the trivial cycle 1 -> 2 -> 1 is the only pole on the line; the characters of THM-4629 modulo m up to the flip are the cusps of the principal modular curves X(m), Platonic for m = 2, 3, 4, 5 (3, 4, 6, 12 cusps), with the twelve cusps of X(5) mapping two-to-one onto P^1(F_5), the carrier of the exotic S_5 of Out(S_6); the Fisher information of a line is (2 - theta)^6/(8 pi^2 theta_f^4) and the Dirac limit is the trivial cycle; the owner's sources typed: mu(log r) = 2 (preprint; does not reach the ratio log 3/log 2), Kaczynski's boundary functions (no bridge), the Ramanujan-inequality manuscript (effective substitutes for conjectures: 7/8 for RH there, Rhin's 26.6 for the exponent 2 here; quasi-RH => Artin => THM-4631 densities), the quasi-Riemann half-planes"
status: >
  PROVED: THM-4632 statements 1-5 (growth exponent identity for every irrational parameter in (1, 2); rational parameters
  iff poles on the line; Jarnik continua; Fisher information at leading order); the cusp identification of the Farey
  colours (classical). NUMERICAL: the exact peaks (thirty convergents), the local exponents of six parameters to 2600
  partial quotients, the Jarnik laws for mu = 3, 4, 5, the cusp census on X(5), the Artin density to 10^5. CITED: Jarnik,
  Besicovitch, Roth, Lagrange, Khinchin; Rhin (via Simons-de Weger, from memory); Rhin-Viola 2001; Marcovecchio 2009;
  Salikhov (from memory); Hooley 1967; Klein; the exotic S_5 in S_6; Becher-Reimann-Slaman arXiv:1601.00153. PREPRINT
  CLAIM (abstracts fetched, not checked): Liu-Jiang-Zhang arXiv:2610.10192; Cureton-Thorner hexagon:2610.00200; the
  pasted Weil-group 7/8 abstract; the OpenAI primitive-roots preprint (unread); the owner's Ramanujan-inequality
  manuscript (README read). ANALOGY: the Kaczynski theme; the Dirac reading. OPEN: mu(log2 3); the per-character exponents;
  whether the quasi-RH preprints reach Hooley's hypothesis for the Kummer fields.
session: opus-2026-10-14-S31
source: 05-knowledge/results/collatz_integers_as_irrationality_measures_20261014.md
scripts:
  - 04-computation/experiments/collatz_irrationality_exponent_growth_20261014.py (+ .out: local exponents of log2 3; exact peaks; Jarnik laws; cusp census; Artin density; Fisher information)
  - 04-computation/experiments/collatz_irrationality_exponent_zeta_20261014.py (+ .out: local exponents of zeta(3), zeta(5), zeta(7), 1 + log 2, log 3, log2 3 to 2600 terms)
related:
  - THM-4632 (the theorem), THM-4629 (the law), HYP-9260 (line bound and peak law), HYP-9261 (Farey colouring; record lines are the convergents), THM-4631 (fast modes; Artin), THM-4600 (cycle points)
  - 05-knowledge/results/collatz_integers_as_irrationality_measures_20261014.md
  - 05-knowledge/results/odd_zeta_parallels_petersen_20261001.md (the failed zeta(5), zeta(7), zeta(9) claim), collatz_expense_diophantine_20261005.md (the record table of near-returns)
---

# HYP-9262 — integers as irrationality measures

## Statements

1. **(The ladder; PROVED, THM-4632.)** For every irrational `x` in `(1, 2)` the formal law `R_x(s) = 1/(1 - e^(-s)/2 - e^(s(x-1))/2)` has `limsup log|R_x(2 pi i f)|/log f = 2(mu(x) - 1)` on the line; the integer `n >= 2` labels the laws of exponent `2n - 2`, a continuum of Hausdorff dimension `2/n` (Jarnik), with explicit members `a_(k+1) = ceil(q_k^(n-2))` for `n >= 3`; `n = 1` is the rational class, exactly the laws with poles on the line (HYP-9260's rational approximants); the Collatz law `x = log2 3` has exponent `2(mu(log2 3) - 1)`, conjecturally `2`, at most `26.6` by Rhin's effective measure.
2. **(The only abnormality; PROVED for the law.)** The trivial cycle `1 -> 2 -> 1` is the pair `(Q, P) = (0, 0)`, the only exact return `2^P = 3^Q` and the only pole of the Collatz law on its line; every other near-return is a strictly off-line pole seen as a Lorentzian peak. For the integers this is the parabolic version of the cycle statement, not the cycle statement: cycles need `2^P - 3^Q | c_w`, and their exclusions (Steiner, Simons-de Weger) use the same Diophantine invariant, an effective `mu(theta)`.
3. **(The sources; typed.)** Liu-Jiang-Zhang's `mu(log r) = 2` (PREPRINT CLAIM) would put `log 2` and `log 3` in the generic class and leaves `mu(log 3/log 2)` open (a ratio of two numbers of exponent `2` can have any exponent); the odd zeta values are in the same position (`mu(zeta(3)) <= 5.513891`; `zeta(5)` irrational iff the law with parameter `zeta(5)` has no pole on its line). Kaczynski's boundary-function theorems (Lean-formalised) do not apply to a function continuous up to its line: no bridge. The owner's Ramanujan-inequality manuscript is the same move as THM-4632's Rhin bound (an effective substitute for a conjecture in an explicit exponent) and, through Hooley, links the quasi-Riemann half-planes to THM-4631's Artin density (`0.37566` of the odd primes below `10^5`; conditional on the unread preprints).
4. **(Cusps; PROVED, classical.)** The Farey colours of HYP-9261 up to sign are the cusps of `X(m)` (`J_2(m)/2` of them for `m >= 3`): `3, 4, 6, 12` for `m = 2, 3, 4, 5`, the vertices of the dihedron, tetrahedron, octahedron and icosahedron; the convergents of `log2 3` visit the twelve cusps of `X(5)` with counts `195-238` and the six points of `P^1(F_5)` (the exotic `S_5 = PGL(2, 5)` of the outer automorphism of `S_6`) with counts `422-453`.
5. **(Fisher information; PROVED at leading order.)** The information a line carries about its own frequency is `1/(2 gamma^2) = (2 - theta)^6/(8 pi^2 theta_f^4)` (`6.47343 x 10^(-5) theta_f^(-4)`, to six digits from the ninth line), growing like `q_(n+1)^4` along the convergents; the Dirac limit is the trivial cycle.
6. **(Local exponents; NUMERICAL.)** `log2 3`: the largest `e_n = log Q_(n+1)/log Q_n` is `1.0262` beyond `n = 100`, `1.0045` beyond `n = 1000`, `1.0024` beyond `n = 2000` (mean `1.0004`); `zeta(3), zeta(5), zeta(7), 1 + log 2, log 3`: largest `e_n` beyond `n = 100` / `1000` of `1.0435 / 1.0037`, `1.0335 / 1.0056`, `1.0244 / 1.0059`, `1.0362 / 1.0054`, `1.0304 / 1.0047`; the largest partial quotients in 2600 terms are `9283` (`zeta(7)`) and `170644` (`zeta(5)`); none moves the local `mu` above `2.01` beyond `n = 1000`, the generic behaviour and no evidence beyond that.

## Cheapest decisive tests

- Part 1 to `10^4` and `10^5` partial quotients of `log2 3`; the per-colour exponents for `m = 2, 3, 5`; a Jarnik law with `mu = 2.5`; the raw OpenAI PDF; the weakest effective measure of `theta` that yields Simons-de Weger's `m <= 68` (the manuscript's Theorem-B pattern).

## Why

The owner asked to see the positive integers as irrationality measures, with `1 <-> 2` the only abnormality, and to read the odd-zeta/logarithm hypothesis, Fisher information, `Out(S_6)`, the Ramanujan-inequality repository and the quasi-Riemann preprints into Collatz. THM-4632 is the exact form of the first request: the irrationality exponent of `log2 3` is a growth exponent of the Collatz transfer function, the integers label the Jarnik classes of formal laws, the rational class is the class with poles on the line, and the trivial cycle is the only pole of the Collatz law. The sources are typed to what they would and would not give.
