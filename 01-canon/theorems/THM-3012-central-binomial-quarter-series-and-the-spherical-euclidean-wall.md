---
id: THM-3012
title: "The central-binomial quarter series S(k) and the spherical/Euclidean wall at k=4"
status: >
  PROVED + VERIFIED-EXACT. For S(k) = sum_{n>=0} C(2n,n)C(4n,2n)/((kn+1)64^n):
  (i) C(2n,n)C(4n,2n)/64^n = (1/4)_n(3/4)_n/(n!)^2 exactly, so the generating
  function is Ramanujan's signature-4 series 2F1(1/4,3/4;1;x) and
  S(k) = 3F2(1/4,3/4,1/k; 1,1+1/k; 1); (ii) the b-a=1/2 quadratic
  transformation gives the elementary kernel
  2F1(1/4,3/4;1;z) = (1/pi) int_0^pi dphi/sqrt(1+sqrt z cos phi), hence a
  K-moment form and the uniform representation
  S(k) = (2/pi) int_0^{pi/2} 2F1(1,4/k;1+2/k;(1-sqrt2 sin th)/2) dth;
  (iii) THOMAE NORMAL FORM S(k) = (8 sqrt2/(3 pi k)) 3F2(1-1/k,1,1;5/4,7/4;1),
  in which only ONE parameter carries k; (iv) THE WALL: the Schwarz angle sum
  of the controlling 2F1(1,4/k;1+2/k;.) equals 4/k + |1-4/k|, which is
  8/k - 1 > 1 (SPHERICAL, finite monodromy) exactly for k = 1,2,3 and is
  EXACTLY 1 (EUCLIDEAN) for every k >= 4. So the three classical closed forms
  are exactly the spherical cases, k=4 is the degeneration point, and no
  closed form of that type exists for k >= 4; a PSLQ battery at 120-150
  digits over the natural fields finds none for k = 4,5,6,12 while the same
  sweep recovers k=3 immediately. Also arctan(sqrt2/5) = pi - 3 arctan(sqrt2),
  so pi S(3) = 2 sqrt3 log(sqrt3+sqrt2) - 2 pi + 6 arctan(sqrt2).
source: death-star-2026-07-31-coinC2
depends_on: []
related:
  - THM-2000
external:
  - "Owner-supplied problem, 2026-07-31 (series over products of central binomial coefficients)."
script: 04-computation/central_binomial_quarter_series_thm3012.py
output: 05-knowledge/results/central_binomial_quarter_series_thm3012.out
---

# THM-3012 -- the central-binomial quarter series and the k = 4 wall

For an integer `k >= 1`,

```text
S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ( (kn+1) 64^n ).
```

## 1. The summand is a signature-4 coefficient (PROVED)

```text
C(2n,n) C(4n,2n)/64^n = (1/4)_n (3/4)_n / (n!)^2.                     (1)
```

*Proof.* `C(2n,n) = 4^n(1/2)_n/n!`, `C(4n,2n) = 16^n(1/2)_{2n}/(2n)!`, and
Legendre duplication gives `(1/2)_{2n} = 4^n(1/4)_n(3/4)_n`,
`(2n)! = 4^n n!(1/2)_n`; the `(1/2)_n` cancel. (Referee C1: exact rationals,
`n < 60`.) Hence the generating function is **Ramanujan's theory of
signature 4**, and with `1/(kn+1) = int_0^1 t^{kn}dt`,

```text
S(k) = int_0^1 2F1(1/4,3/4;1;t^k) dt = 3F2(1/4,3/4,1/k; 1,1+1/k; 1).   (2)
```

The `3F2` is balanced but non-terminating; Watson and Whipple fit the
numerator pair `(1/4,3/4)` but force `k = 1`, which is exactly Gauss:
`S(1) = 2F1(1/4,3/4;2;1) = 16/(3 pi sqrt2) = 8 sqrt2/(3 pi)`.

## 2. An elementary kernel, and two normal forms (PROVED)

Since `b - a = 1/2`, a quadratic transformation applies, and reduces the
signature-4 function to a completely elementary integral:

```text
2F1(1/4,3/4;1;z) = (1+sqrt z)^{-1/2} 2F1(1/2,1/2;1;2 sqrt z/(1+sqrt z))
                 = (1/pi) int_0^pi dphi / sqrt(1 + sqrt z cos phi).     (3)
```

Consequently (referee C4, C6; verified `k = 1..6`)

```text
S(k) = (16/(k pi sqrt2)) int_0^1 mu^{4/k-1}(2-mu^2)^{-2/k-1/2} K(mu) dmu, (4)
S(k) = (2/pi) int_0^{pi/2} 2F1(1, 4/k; 1+2/k; (1 - sqrt2 sin th)/2) dth.  (5)
```

**Thomae normal form** (referee C5; verified `k = 1,...,12`):

```text
S(k) = (8 sqrt2/(3 pi k)) * 3F2(1 - 1/k, 1, 1; 5/4, 7/4; 1).            (6)
```

In (6) the whole `k`-dependence sits in the single numerator parameter
`1 - 1/k`; at `k = 1` that parameter is `0` and the series terminates at
`n = 0`, returning `8 sqrt2/(3 pi)` immediately.

## 3. The wall (PROVED)

The controlling function in (5) is `2F1(a,b;c;.)` with
`a = 1, b = 4/k, c = 1 + 2/k`. Its Schwarz exponent differences are

```text
lambda = |1-c| = 2/k,  mu = |c-a-b| = 2/k,  nu = |a-b| = |1 - 4/k|,
lambda + mu + nu = 4/k + |1 - 4/k|
                 = 8/k - 1   (k <= 3),      = 1   (k >= 4).             (7)
```

So the angle sum is `> 1` — **spherical**, i.e. the monodromy is finite and
the function is algebraic — exactly for `k = 1, 2, 3` (values `7, 3, 5/3`),
and is **exactly 1 — Euclidean — for every `k >= 4`** (referee C7, checked
for `k < 40`). `k = 4` is precisely the degeneration point.

This explains the classical data completely: the three closed forms are the
three spherical cases, and no fourth one of that type can exist. For the
record, the `k = 3` value simplifies, since
`arctan(sqrt2/5) = pi - 3 arctan(sqrt2)` (referee C8):

```text
pi S(3) = 2 sqrt3 log(sqrt3 + sqrt2) - 2 pi + 6 arctan(sqrt2).          (8)
```

## 4. The negative search, with a positive control (VERIFIED-NUMERIC)

PSLQ at 120-150 digits over the natural constant pool -- `pi`,
`log(5+2sqrt6)`, `sqrt3 log(5+2sqrt6)`, `arctan(sqrt2/5)`, `log(1+sqrt2)`,
`log(2+sqrt3)`, `arctan(sqrt2)`, `arctan(sqrt2/3)`, `log 2`, Catalan's `G`,
the lemniscate constant `varpi = Gamma(1/4)^2/(2 sqrt(2 pi))` -- subset sizes
2-3, coefficient bound 800:

* **positive control:** the sweep applied to `pi S(3)` recovers
  `pi S(3) - sqrt3 log(5+2 sqrt6) + 2 arctan(sqrt2/5) = 0` immediately
  (14 relations found);
* **no relation** is found for `k = 4, 5, 6, 12` (referee C9).

Wider sweeps (fields `Q(sqrt2)`, `Q(2^{1/4})`, `Q(sqrt2,sqrt5)` with the
golden ratio; lemniscatic `varpi, varpi^2, varpi pi, varpi log(1+sqrt2)`;
weight-two constants `pi^2, G, log^2 2, log^2(1+sqrt2), pi log(1+sqrt2)`;
subset sizes up to 4; coefficient bounds to `10^6`) likewise find nothing.
For `k = 4` the explicit reductions

```text
S(4) = (2/pi) int_0^1 [arcsin v + arcsinh v]/sqrt(1-v^4) dv
     = (4/(pi sqrt2)) int_0^1 K(mu)/(2-mu^2) dmu
```

hold exactly (the integrand pair is the `v -> iv` symmetry of
`sqrt(1-v^4)`, as `arcsin(iv) = i arcsinh v`).

## 5. Scope

(7) is a statement about the *monodromy class* of the controlling `2F1`. It
proves that the mechanism producing the `k <= 3` evaluations is unavailable
for `k >= 4`, and together with the validated PSLQ battery it is strong
evidence that `S(4)`, `S(5)` have no closed form in logarithms of algebraic
numbers, `pi`, and the standard weight-<=2 constants. It is **not** a proof
of irrationality or of non-existence in a wider constant class. What *is*
proved for all `k` is the closed form (6): a single balanced `3F2` at 1 with
one `k`-dependent parameter, equivalently the `K`-moment (4).

Referee: C1-C9 exact/high-precision, `ALL THM-3012 REFEREE CHECKS PASSED`. QED.
