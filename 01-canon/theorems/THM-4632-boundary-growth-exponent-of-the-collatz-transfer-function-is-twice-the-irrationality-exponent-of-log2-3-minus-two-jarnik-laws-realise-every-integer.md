---
id: THM-4632
title: "The boundary growth exponent of the Collatz transfer function is twice the irrationality exponent of log2 3 minus two: for the trivial character, R(s) = 1/(1 - e^(-s)/2 - e^(s(theta-1))/2) with theta = log2 3 satisfies limsup_(f -> infinity) log |R(2 pi i f)| / log f = 2 (mu(theta) - 1), where mu(theta) is the irrationality exponent of theta; the same holds for the formal law with any irrational parameter x in (1, 2), so by Jarnik every integer n >= 2 is the half-exponent-plus-one of a continuum of laws (Hausdorff dimension 2/n) with exponent 2n - 2, the rational parameters are exactly the laws with poles on the line (exponent infinity), and mu(log2 3) = 2 (open; the generic value) is equivalent to the Collatz transfer function growing like f^(2 + o(1)) on its line; Rhin's effective measure gives the exponent at most 26.6, and the lines carry Fisher information (2 - theta)^6/(8 pi^2 theta_f^4) about their frequency"
status: >
  PROVED: statements 1 and 2 (the growth exponent identity, for theta = log2 3 and for every irrational x in (1, 2)),
  from HYP-9260 statement 4 (the line bound |R(2 pi i f)| <= x^2/(4 theta_f(f)^2), Jordan's inequality) and the peak
  law at the convergents (THM-4629 statement 3's root expansion and the residue computation of HYP-9260: |R| at the
  drifted frequency equals 2/((2 - x)|Re w|)(1 + O(theta_f^2)) = (2 - x)^2/(2 pi^2 theta_f^2)(1 + O(theta_f))), together
  with Lagrange's characterisation of the convergents and mu(x) = 1 + limsup log q_(n+1)/log q_n. PROVED: statement 3
  (rational parameters have exact poles on the line, HYP-9260 statement 2) and statement 4 (Jarnik's continua, CITED:
  Jarnik 1931 and Besicovitch 1934 for dim {mu(x) >= mu} = 2/mu, Jarnik for the exact-exponent set; Becher, Reimann and
  Slaman, arXiv:1601.00153, for the effective version). CITED: Rhin's effective measure as used by Simons and de Weger
  (from memory; the hypothesis of the bound unverified, S29); Marcovecchio 2009 (mu(log 2) <= 3.5746) and Salikhov
  (mu(log 3) <= 5.125) for the single logarithms; the irrationality exponent of theta = log 3/log 2 itself is OPEN and is
  not a consequence of mu(log 2) = mu(log 3) = 2. PREPRINT CLAIM (not checked): Liu, Jiang and Zhang, arXiv:2610.10192
  (October 2026), mu(log r) = 2 for every positive rational r != 1; it concerns single logarithms and rational linear
  combinations of them and does not cover the ratio theta. NUMERICAL: the exact peaks at the first thirty convergents
  equal the peak law to all printed digits from n = 8; the local exponents 2 log Q_(n+1)/log Q_n of log2 3 fall below
  2.05 beyond the hundredth convergent and below 2.005 beyond the two-thousandth; the Jarnik laws with mu = 3, 4, 5 show
  local exponents 3.92, 5.98, 7.98 at their largest computed convergents (targets 4, 6, 8); the Fisher-information law of
  statement 5 to six digits. OPEN: mu(log2 3); whether the per-character exponent equals the global one (it does along
  the convergents of the character's colour, HYP-9261, which are a third of all convergents for m = 2 to n = 2600).
session: opus-2026-10-14-S31
source: 05-knowledge/results/collatz_integers_as_irrationality_measures_20261014.md
scripts:
  - 04-computation/experiments/collatz_irrationality_exponent_growth_20261014.py (+ .out: Part 1 local exponents; Part 2 exact peaks against the law; Part 3 Jarnik laws; Part 6 Fisher information)
related:
  - THM-4629 (the law; statement 3 the root expansion), HYP-9260 (the line bound and the peak law; the rational approximants with exact circle poles), HYP-9261 (the Farey colouring; the record lines are the convergents)
  - HYP-9262 (integers as irrationality measures: the dictionary and the sources), THM-4631 (the fast modes; Artin)
  - 05-knowledge/results/collatz_integers_as_irrationality_measures_20261014.md
---

# THM-4632 — the boundary growth exponent is 2(mu - 1)

## Setting

For an irrational `x` in `(1, 2)` let `R_x(s) = 1/(1 - e^(-s)/2 - e^(s(x-1))/2)`, the scalar law of THM-4629 for the trivial character with `x` in place of `theta = log2 3` (for `x = theta` it is the Collatz transfer function of the untwisted step-count law). Write `delta = 2 - x`, `||t||` for the distance from `t` to the nearest integer, `p_n/q_n` for the convergents of `x`, and `mu(x) = 1 + limsup_n log q_(n+1)/log q_n` for the irrationality exponent (the infimum of the `mu` for which `|x - p/q| < q^(-mu)` has finitely many solutions; `mu = 2` for almost every `x` and for every algebraic irrational, Roth; `mu = 1` for rationals; `mu = infinity` for Liouville numbers). For `f > 0` let `Q = Q(f)` be the integer nearest `f`, `P` the integer nearest `Q x`, and `theta_f(f) = Q x - P`, so `|theta_f(f)| = ||Q x||`.

## Statements

1. **(Growth exponent; PROVED.)** `limsup_(f -> infinity) log |R_x(2 pi i f)| / log f = 2 (mu(x) - 1)` (both sides may be `+infinity`). More precisely: (a) for every `f >= 1`, `|R_x(2 pi i f)| <= x^2/(4 ||Q(f) x||^2)`; (b) at the convergents, with `theta_f = q_n x - p_n`, the value of `|R_x|` on the line at the drifted frequency `q_n + Im(w_n)/(2 pi)` is `2/(delta |Re w_n|) (1 + O(theta_f^2)) = (delta^2/(2 pi^2 theta_f^2)) (1 + O(theta_f))`, where `w_n` is the root of `e^(-w) + e^(2 pi i theta_f) e^(w(x-1)) = 2` continued from `0`; (c) `|q_n x - p_n| = 1/(q_(n+1) + q_n t_(n+2))` with `t_(n+2) = [0; a_(n+2), a_(n+3), ...]` in `(0, 1)`, so `|R_x| = (delta^2/(2 pi^2)) (q_(n+1) + q_n t_(n+2))^2 (1 + O(1/q_(n+1)))` at the `n`-th line.
2. **(Collatz; PROVED.)** For `x = theta = log2 3` the Collatz transfer function has growth exponent `2(mu(theta) - 1)` on its line: it grows like `f^(2 + o(1))` iff `mu(theta) = 2`, and like `f^(2n - 2 + o(1))` iff `mu(theta) = n`. Rhin's effective measure (`|P ln 2 - Q ln 3| > Q^(-13.3)`, CITED as in S29) gives `|R(2 pi i f)| <= C f^(26.6)`; the local exponents `2 log Q_(n+1)/log Q_n` along the convergents are below `2.05` from the hundredth convergent and below `2.005` from the two-thousandth (Part 1), consistent with `mu(theta) = 2`, which is open.
3. **(The rational parameters; PROVED.)** For `x = P/Q` rational the law is `2 pi i Q`-periodic in `s` and has poles on the line (HYP-9260 statement 2: `R_(P/Q)(w) = w^(P-Q)/p(w)` with the circle root `w = 1` for the trivial character); the growth exponent is `+infinity` in the sense that `|R|` is unbounded on every period. These are the only parameters with a pole on the line: for irrational `x`, `R_x` is finite at every point of the line, and for `x = theta` the single pole at `s = 0` (the pair `(Q, P) = (0, 0)`, the trivial cycle `1 -> 2 -> 1`) is the only exact resonance because `2^P != 3^Q` for `(Q, P) != (0, 0)` (HYP-9260 statement 5: no parity word is parabolic).
4. **(Jarnik laws; PROVED and CITED.)** For every real `mu >= 2` the set of `x` in `(1, 2)` with `mu(x) = mu` has Hausdorff dimension `2/mu` (Jarnik; Besicovitch), and for `mu > 2` the continued fraction with `a_(n+1) = ceil(q_n^(mu - 2))` is an explicit member. Hence for every integer `n >= 2` there is a continuum of laws `R_x` of exponent exactly `2n - 2`: `n = 2` is the generic class (almost every `x`, every algebraic `x`, conjecturally `theta`), `n = 1` is the rational class of statement 3, and `n >= 3` are the constructed classes. Part 3: the laws built from `a_(n+1) = ceil(q_n^(mu-2))` for `mu = 3, 4, 5` have local exponents `3.92, 5.98, 7.98` at their largest computed convergents (`q_n` with `23, 78, 79` digits).
5. **(Fisher information of a line; PROVED at leading order.)** Near the `n`-th line `|R_x(2 pi i f)|^2` is a Lorentzian of half-width `gamma_n = |Re w_n|/(2 pi) = 2 pi theta_f^2/delta^3 (1 + O(theta_f))` in `f`; the Fisher information of the corresponding Cauchy location family is `1/(2 gamma_n^2) = delta^6/(8 pi^2 theta_f^4) (1 + O(theta_f))` (Part 6: the ratio `I theta_f^4` equals `delta^6/(8 pi^2) = 6.47343 x 10^(-5)` to six digits from the ninth line on, `6.96 x 10^(-5)` at the sixth). Along the convergents the information grows like `q_(n+1)^4`, i.e. like `q_n^(4(mu - 1))`; the Dirac limit (`gamma = 0`, infinite information) is the pole on the line, attained only by the trivial cycle.

## Proof

**1(a).** HYP-9260 statement 4 for `m = 1`, `alpha = beta = 0`: with `phi, psi` the principal arguments of `e^(-s)` and `e^(s(x-1))` at `s = 2 pi i f`, `Re(1 - e^(-s)/2 - e^(s(x-1))/2) = (1 - cos phi)/2 + (1 - cos psi)/2 >= (phi^2 + psi^2)/pi^2` (Jordan), and `psi + (x - 1) phi = 2 pi (Q x - P)` with `Q` the integer nearest `f` and `P` the integer making `psi` principal, so `max(|phi|, |psi|) >= 2 pi |Q x - P|/x` and `|R_x| <= pi^2/(phi^2 + psi^2) <= x^2/(4 (Q x - P)^2)`; `|Q x - P| >= ||Q x||`. (The proof in the S29 note is written for `theta`; nothing uses more than `1 < x < 2`.)

**1(b).** THM-4629 statement 3 with `theta_f = q_n x - p_n` and `delta = 2 - x`: the root `w_n` continued from `0` is simple and `w_n = 2 pi i theta_f/delta - 4 pi^2 theta_f^2/delta^3 + O(theta_f^3)`, with constants depending only on `x` because the equation depends on `theta_f` alone. Let `D(s) = 1 - e^(-s)/2 - e^(s(x-1))/2` and `s_n = 2 pi i q_n + w_n` the pole; `D(s_n) = 0` and `D'(s_n) = e^(-s_n)/2 - (x-1) e^(s_n(x-1))/2 = delta/2 + O(theta_f)`. At the point `s_n^* = s_n - Re w_n` of the line, `|s_n^* - s_n| = |Re w_n|` and `D(s_n^*) = -D'(s_n) Re w_n + O(|Re w_n|^2)`, so `|R_x(s_n^*)| = 2/(delta |Re w_n|) (1 + O(theta_f))`, and with `|Re w_n| = 4 pi^2 theta_f^2/delta^3 (1 + O(theta_f))` this is `delta^2/(2 pi^2 theta_f^2) (1 + O(theta_f))` (the audit of S29 recorded that the first correction is in fact `O(theta_f^2)` with a large coefficient, HYP-9260 statement 4).

**1(c).** `q_n x - p_n = (-1)^n/(q_(n+1) + q_n t_(n+2))` is the standard identity for continued fractions.

**1.** Upper bound: given `epsilon > 0`, by the definition of `mu = mu(x)` there is `Q_0` with `||Q x|| >= Q^(1 - mu - epsilon)` for all `Q >= Q_0`; for `f >= Q_0 + 1/2`, `Q(f) >= Q_0` and `Q(f) <= f + 1/2`, so by 1(a) `|R_x(2 pi i f)| <= (x^2/4) (f + 1/2)^(2(mu - 1) + 2 epsilon)`, and the `limsup` is at most `2(mu - 1)` (if `mu = infinity` there is nothing to prove). Lower bound: by 1(b), 1(c), at `f_n = q_n + Im(w_n)/(2 pi)` (so `f_n = q_n (1 + o(1))`), `|R_x(2 pi i f_n)| >= (delta^2/(4 pi^2)) q_(n+1)^2` for `n` large, hence `log |R_x| / log f_n >= 2 log q_(n+1)/log q_n - O(1/log q_n)`, whose `limsup` is `2(mu - 1)` by the definition of `mu`. Both bounds hold for every irrational `x` in `(1, 2)`; `theta` is one of them. ∎

**2.** Statement 1 with `x = theta`; the local exponents are Part 1; the Rhin bound is S29 (HYP-9260 statement 4). ∎

**3.** HYP-9260 statement 2 (the rational approximants) for the poles; for irrational `x`, `D(2 pi i f) = 0` forces `e^(-2 pi i f) = e^(2 pi i f(x-1)) = 1`, i.e. `f` and `f x` both integers, i.e. `f = 0`; for `x = theta` the pole at `s = 0` is the pair `(0, 0)` and `2^P = 3^Q` has no other solution. ∎

**4.** Jarnik and Besicovitch for the dimension (CITED); for the explicit members, `mu(x) = 1 + limsup log q_(n+1)/log q_n = 1 + limsup log(a_(n+1) q_n)/log q_n = 1 + (mu - 2) + 1 = mu` when `a_(n+1) = ceil(q_n^(mu-2))`; statement 1 gives the exponent `2(mu - 1)`; the dimension of the set of parameters is `2/mu` because the shift `x -> x - 1` and the restriction to `(1, 2)` do not change Hausdorff dimension. Part 3 is the computation. ∎

**5.** Near `s_n`, `R_x(s) = (2/delta)/(s - s_n) (1 + O(theta_f))`, so on the line `|R_x(2 pi i f)|^2 = (4/delta^2)/((2 pi (f - f_n))^2 + (Re w_n)^2) (1 + O(theta_f))`, a Lorentzian in `f` of half-width `gamma_n = |Re w_n|/(2 pi)`; the Fisher information of the Cauchy location family with scale `gamma` is `1/(2 gamma^2)` (classical). Substituting `Re w_n` from 1(b) gives the law; Part 6 checks it. ∎

## Scope and boundary

- The theorem is about the model law (the scalar law of THM-4629), not about the integers: the step from the law to the integer statistics is THM-4629's open 3-adic closure and is untouched. What is new is the identification of a classical Diophantine invariant of `log2 3` with a growth exponent of the Collatz transfer function, and the realisation of every integer `n >= 2` as the exponent class of a continuum of formal laws (the owner's "integers as irrationality measures", HYP-9262).
- The exponent is a `limsup` along the convergents; between the lines `|R|` is small, and the per-character exponents for nontrivial characters depend on the convergents of the character's colour (HYP-9261); the theorem is stated for the trivial character, where every convergent resonates.
- `mu(theta) = 2` would follow from a Roth-type theorem for `log 3/log 2`, which is not known; Liu, Jiang and Zhang's preprint (if correct) gives `mu(log 2) = mu(log 3) = 2` but says nothing about their ratio, and the effective bounds for the ratio remain of Baker type (Rhin). The owner's hypothesis "all logarithms and all odd zeta values have `mu = 2`" is the generic expectation; the Collatz-relevant instance is the ratio.
