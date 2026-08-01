---
id: HYP-9076
title: "The Factorial Conjecture is a moment-determinacy question, and the difficulty is indexed by DEGREE"
status: >
  FRAMING + CENSUS. L is the moment functional of the PRODUCT EXPONENTIAL
  measure, so FC(n) says: if every moment (m>=1) of the pushforward
  nu = f_* mu vanishes then f = 0 -- exactly a determinacy question. The
  Krein/Carleman criterion then locates the difficulty by DEGREE, not
  dimension: a Stieltjes density decaying like exp(-t^rho) is determinate iff
  rho >= 1/2, and the pushforward of e^{-|x|} under a degree-d polynomial has
  tail exp(-t^{1/d}), so the argument closes for deg f <= 2 in EVERY number of
  variables and fails from cubics on. Census: zero survivors among all small
  integer f of degree <= 2 in 2 and 3 variables with the first 3-4 moments
  vanishing. Degree >= 3 untested and is the real regime.
source: opus-2026-07-31-amm12592-writeup
related:
  - THM-2922  # first-window SFC(4), Macaulay-Newton atlas
  - THM-2891
  - HYP-9075  # the shell-imbalance module, same moment-vanishing shape
external:
  - "van den Essen, Wright, Zhao, On the image conjecture (Factorial Conjecture)."
script: 04-computation/fc_factorial_conjecture_census.py
output: 05-knowledge/results/fc_factorial_conjecture_census.out
---

# HYP-9076 -- FC as determinacy

## 1. The framing

With `L(x^alpha) = alpha!`,

```text
L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx,
```

so `L` is the moment functional of the product exponential measure `mu`, and
`L(f^m)` is the `m`-th moment of the pushforward `nu = f_* mu`. FC(n) reads:

```text
all moments (m >= 1) of nu vanish   =>   f = 0.
```

If `nu` is DETERMINATE this is immediate: vanishing moments plus total mass
force `nu = delta_0`, hence `f = 0` `mu`-a.e., hence `f = 0` as a polynomial.
**So FC is exactly a moment-determinacy statement.**

## 2. The difficulty is indexed by degree, not dimension

Krein/Carleman: a Stieltjes density decaying like `exp(-t^rho)` is determinate
iff `rho >= 1/2`. The pushforward of `e^{-|x|}` under a degree-`d` polynomial
has tail `exp(-t^{1/d})`, so

```text
deg f <= 2 :  rho = 1/d >= 1/2  ->  DETERMINATE  ->  FC closes
deg f >= 3 :  rho < 1/2         ->  INDETERMINATE -> the argument fails
```

in EVERY number of variables. This is worth stating because FC is indexed by
`n` in the literature and by supports/slots in the repo's SFC work, while the
determinacy obstruction is indexed by `d`. Caveat: for complex `f` the
pushforward is a complex measure and the criterion needs the real/imaginary
split; the clean statement is for real `f`.

## 3. Census

Small integer `f` with `L(f^m) = 0` for `m = 1..M`:

```text
n=2, deg<=2, coeffs in [-2,2], M=3 :  0 survivors / 3124
n=2, deg<=2, coeffs in [-3,3], M=4 :  0 survivors / 16806
n=3, deg<=1, coeffs in [-3,3], M=4 :  0 survivors / 342
n=3, deg<=2, coeffs in [-1,1], M=3 :  0 survivors / 19682
```

Nothing survives even the first three moments in degree `<= 2`, consistent
with (and much stronger than) the determinacy argument, which only says the
infinite family of moments suffices.

## 4. Open

Degree `>= 3` is untested here and is exactly where section 2 says the moment
argument stops working. That is the regime the repo's SFC machinery
(THM-2922's four-slot supports of diameter five) actually lives in, and the
census should be pushed there rather than to higher `n`.

## 5. Relation to HYP-9075

The shell-imbalance module has the same shape -- a signed object whose moments
against a PRODUCT measure vanish -- with Bernoulli in place of exponential and
the family indexed by `p` rather than by powers. The difference is that our
Bernstein basis makes the determinacy step trivial (`D_m = 0` iff the
composition vector vanishes), so all our difficulty is in the integer box,
whereas FC has no box and all its difficulty is in the determinacy. The two
problems are complementary in exactly the way HYP-9075 sec 8 found for
THM-2922's machinery.
