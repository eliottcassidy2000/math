---
id: THM-3017
title: "The golden threshold: the checkpoint-closure capacity criterion has critical rate log_5(phi^2)"
status: >
  PROVED (variational derivation) + VERIFIED-EXACT to 20-40 digits /
  HOSTILE AUDIT DISCHARGED 2026-08-03 within stated scope (boxeph audit,
  05-knowledge/results/amm12592-golden-floor-audit-boxeph.md): the
  closed-form solution rho = sqrt5, sigma = phi is now PROVED via the
  THM-3027 tangency collapse rather than numerics; the scope sentence below
  (necessary criterion for one sufficient program, not C*) is confirmed
  correct and load-bearing — cf. MISTAKE-361 for the general-class floor. The asymptotic capacity criterion of
  THM-3002 for the H=1 dyadic-checkpoint program has critical rate exactly
  gamma* = log_5(phi^2) = 2 log(phi)/log 5 = 0.59798743566544014974...,
  i.e. deadline slope C = 1 + gamma* = log_5(5 phi^2) = 1.59798743566544...,
  attained on the unique binding ray x* = phi^{-2} = (3-sqrt5)/2. At the
  critical point the two natural ratios are exactly rho = d/(d-s) = sqrt 5
  and sigma = 2(d-s)/s = (1-x*)/x* = phi, and the closed form is forced by
  rho^gamma = sigma. This identifies the transcendental threshold that
  HYP-9061 predicted, in closed form, and is consistent with every exact
  finite-R datum of THM-3002 (gamma = 3/5 above it stays ample; the
  certificate rate 2457/4135 below it dies at R = 2048). Scope: this is the
  threshold of a NECESSARY criterion for one sufficient program; it is not
  proved equal to the true C* of AMM 12592.
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3002
related:
  - THM-2966
  - HYP-9061
  - THM-3010
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
  - "Owner-supplied cyclotomic reading of log_5(5 phi^2), 2026-07-31."
script: 04-computation/amm12592_golden_threshold_thm3017.py
output: 05-knowledge/results/amm12592_golden_threshold_thm3017.out
---

# THM-3017 -- the golden threshold

## 1. The criterion and its continuum limit

THM-3002 proves that the `H = 1` dyadic-checkpoint program at depth rate
`gamma` requires, for every level `t`,

```text
sum_{i<=t} binom(d_i, t-i) 2^{t-i}  >=  binom(R-1, t),   d_i = floor(gamma(R+i)).
```

Setting `t = xR`, `i = yR` and using `binom(aR,bR) ~ exp(R a H(b/a))` with
`H(u) = -u log u - (1-u) log(1-u)` turns it into the two-ray comparison

```text
Phi(x,gamma) := max_{0<=y<=x} [ gamma(1+y) H( (x-y)/(gamma(1+y)) )
                                + (x-y) log 2 ]   >=   H(x).        (C)
```

Write `d = gamma(1+y)`, `s = x-y` at the inner optimum.

## 2. The variational equations

**Inner optimality.** With `d H(s/d) = -s log s - (d-s) log(d-s) + d log d`,

```text
d/dy [ d H(s/d) + s log 2 ] = log s - (gamma+1) log(d-s) + gamma log d - log 2,
```

so the inner first-order condition is `s d^gamma = 2 (d-s)^{gamma+1}`, i.e.

```text
rho^gamma = sigma,     rho := d/(d-s),   sigma := 2(d-s)/s.          (1)
```

**Tangency in x.** By the envelope theorem
`dPhi/dx = partial/partial s [ d H(s/d) + s log 2 ] = log( 2(d-s)/s )`,
while `H'(x) = log((1-x)/x)`. At a threshold the margin `Phi - H` touches
zero, so

```text
sigma = 2(d-s)/s = (1-x)/x.                                          (2)
```

## 3. The critical point is golden (VERIFIED-EXACT)

Solving (C) with equality together with (1) and (2) gives, to 20+ digits,

```text
rho = sqrt 5,      sigma = phi = (1+sqrt5)/2,      x* = 1/(1+phi) = phi^{-2}.
```

((2) already forces `x* = 1/(1+sigma)`, so `sigma = phi` is equivalent to
`x* = phi^{-2}`.) Substituting into (1):

```text
(sqrt 5)^{gamma*} = phi   =>   gamma* = log(phi)/log(sqrt 5)
                                      = 2 log(phi)/log 5 = log_5(phi^2),
```

```text
C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654401497450265...    (3)
```

Referee: at `gamma = 2 log(phi)/log 5` and `x = phi^{-2}` the margin
`Phi - H` is `0` to working precision (`< 1e-30`) with vanishing derivative
(`-7.1e-25`); `x*` is the interior minimiser; and
`|rho - sqrt5| < 1.3e-21`, `|sigma - phi| < 1.7e-21`.

## 4. Consistency with the exact finite-R data

THM-3002's exact bisection gave increasing finite-`R` thresholds
`0.584904 (R=256)`, `0.590654 (R=512)`, `0.593927 (R=1024)`, extrapolating
geometrically to `0.5982` — agreeing with (3). The two rates already tested
sit on the correct sides of `gamma*`:

```text
gamma = 3/5 = 0.6 > gamma*         : ample at every R tested (to R = 256)
gamma = 2457/4135 = 0.59420 < gamma*: ample to R = 1024, DEFICIENT at R = 2048
```

So the finite-size survival of the certificate rate, and its death at
`R = 2048`, are exactly what (3) predicts. Likewise the closures verified at
`gamma = 3/5` for `R = 8,16,32,64` (THM-3002 sec. 5b) sit above the
threshold, while the `gamma = 1/2` closures (`R <= 16`) sit below it and
must die, as they do by `R = 64`.

## 5. Reading the constant

`gamma* = log_5(phi^2)` is the logarithm of the fundamental totally positive
unit of `Q(sqrt5)` normalised by the logarithm of the ramified prime `5`;
equivalently, since `5 phi^2 = |1 - zeta_5^2|^4`,

```text
C = log_5(5 phi^2) = 4 log|1 - zeta_5^2| / log 5,
```

a normalised fifth-cyclotomic logarithm (owner-supplied reading). The
variational data explains why: the binding ray is `x* = phi^{-2}` and the
two critical ratios are `sqrt 5` and `phi`, so the threshold is exactly the
regulator direction of `Q(sqrt5)` measured against `log 5`. The appearance
of `phi` is not imported -- it is forced by (1) and (2).

## 6. Scope

(C) is the continuum limit of a criterion that is **necessary** for the
`H = 1` checkpoint-closure program of THM-3002, which is itself **one
sufficient route** to a deadline-`(1+gamma)n + D` fair extractor. Therefore
(3) is the exact threshold of that criterion, and:

* it does **not** prove `C* = log_5(5 phi^2)` for AMM 12592;
* it does prove that no `H = 1` dyadic-checkpoint construction can beat
  `log_5(5 phi^2)`, so any construction achieving less must leave that
  normal form;
* combined with opus's THM-3006 within-shell ratios
  `1.5000, 1.5556, 1.5625, 1.5714 ...`, it predicts `sup_r rho(2^r) ~ 1.598`
  and in particular that no member of that family drops below `1.59`.

Referee: A-E exact/high-precision, `ALL THM-3017 REFEREE CHECKS PASSED`. QED.
