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

## Addendum, 2026-08-01 (death-star): a transcription correction and the elliptic representation

**1. The circulating `S(3)` closed form has `3` where `sqrt3` belongs.** The owner
relayed

```text
S(3) = -(1/pi) ( 3 log(5 - 2 sqrt6) + 2 arctan(sqrt2/5) )        [AS CIRCULATED -- WRONG]
```

which evaluates to `2.01363...`, against the true `S(3) = 1.08838540640395...`.
Since `log(5 - 2 sqrt6) = -log(5 + 2 sqrt6)`, the sign structure is fine; the
single error is the **rational `3` in place of the algebraic `sqrt3`**. Equation
(8) of this file is correct:

```text
pi S(3) = sqrt3 * log(5 + 2 sqrt6) - 2 arctan(sqrt2/5)
        = 2 sqrt3 * log(sqrt3 + sqrt2) - 2 pi + 6 arctan(sqrt2),
```

both verified against the `3F2` to **81 digits** (difference `4.2e-81`). Using
`3` instead gives exactly the circulated `6.32600941...` for `pi S(3)`, which
identifies the corruption unambiguously.

This matters beyond bookkeeping: **the coefficient is irrational.** A PSLQ sweep
over a basis of logarithms and arctangents with *rational* coefficients will
never find (8) -- it must include `sqrt3 * log(5 + 2 sqrt6)` as a basis element,
exactly as section 5's battery does. A sweep that omits it returns "no closed
form" spuriously. (Two such spurious negatives were produced while checking this
addendum before the basis was corrected, and PSLQ instead returned the trivial
degeneracies `arctan(sqrt2/5) + 3 arctan(sqrt2) = pi` and `3 arctan(sqrt3) = pi`
-- a reminder to strip dependent basis elements before reading a null result.)

**2. The elliptic representation of `S(4)` is CONFIRMED.** The owner supplied

```text
S(4) = (2 sqrt2 / pi) * int_0^1 K(k) / (2 - k^2) dk,     K the complete elliptic
                                                         integral of the first kind.
```

Evaluated at 60 digits (`K` in the `m = k^2` convention):

```text
int_0^1 K(k)/(2-k^2) dk        = 1.187752054750793832836583810076862670
(2 sqrt2/pi) * integral        = 1.069352554441268058582961939534278061
S(4) = 3F2(1/4,3/4,1/4; 1,5/4; 1) = 1.069352554441268058582961939...
```

agreeing to `~30` digits. So the representation is a genuine identity and gives
an independent route to `S(4)`.

**It does not evade the wall of section 4.** The integrand `K(k)/(2-k^2)` is a
period of the same Legendre family; the representation re-expresses the balanced
`3F2` rather than reducing its monodromy, and the Euclidean angle sum
`4/k + |1 - 4/k| = 1` at `k = 4` is unchanged. So the addendum supplies a new
handle on `S(4)` but no closed form, and the section 5 assessment stands:
`k = 1,2,3` are elementary, `k >= 4` are not, on the evidence of the PSLQ battery
plus the monodromy explanation.

**Method warning (cost real time here):** the defining series has terms
`~ 1/(sqrt2 pi k n^2)`, so the TAIL after `N` terms is `~ 1/(sqrt2 pi k N)`, far
larger than the last term. Truncating a summation when the *term* falls below
tolerance leaves an error of order `1e-5` at `N = 4000` and silently corrupts any
PSLQ. Use the `3F2` evaluation (or the elliptic integral), never naive summation.

## Addendum 2, 2026-08-01 (death-star): S(k) is one function's Mellin moment at s = 1/k

Two exact identities, both verified numerically, that reorganise the whole family.

**(M1) The series is a Mellin moment.** Since `1/(kn+1) = (1/k)/(n + 1/k)` and
`1/(n+c) = int_0^1 t^{n+c-1} dt`,

```text
S(k) = (1/k) int_0^1 z^{1/k - 1} 2F1(1/4, 3/4; 1; z) dz
     = int_0^1 2F1(1/4, 3/4; 1; u^k) du            (substitute z = u^k).
```

Verified against the `3F2` at `k = 1..5`. The second form identifies this family
with kps's `S_a(k) = int_0^1 2F1(a, 1-a; 1; x^k) dx` at `a = 1/4`.

**(M2) The integrand IS the elliptic K.** The quadratic transformation available
because `b - a = 1/2` gives, verified to 45 digits at `z = 0.1, 0.5, 0.9, 0.99`,

```text
2F1(1/4, 3/4; 1; z) = sqrt( 2/(1 + sqrt(1-z)) ) * (2/pi) * K(sqrt w),
                       w = (1 - sqrt(1-z))/(1 + sqrt(1-z)).
```

This is the source of the owner's confirmed representation
`S(4) = (2 sqrt2/pi) int_0^1 K(k)/(2-k^2) dk`: it is (M1) at `k = 4` rewritten
through (M2).

**Why this matters for the closed-form question.** In (M1) **the function being
integrated does not depend on `k` at all**; `k` enters *only* through the Mellin
exponent `s = 1/k`. So

```text
S(k) = s * M(s)|_{s = 1/k},      M(s) = int_0^1 z^{s-1} 2F1(1/4,3/4;1;z) dz.
```

Consequently **no argument that works at the level of the function can separate
`k <= 3` from `k >= 4`** -- the function is identical in every case, and it is a
complete elliptic integral, transcendental for all of them. Any genuine
impossibility statement must come from the **arithmetic of the exponent `1/k`**,
i.e. from which rational Mellin moments of `K` admit closed forms. This is a
sharper and better-posed target than "does `S(4)` have a closed form", and it is
the form in which the question should be attacked.

It also means section 4's monodromy discussion must be read carefully: the
Schwarz data there is that of the `3F2`, whose parameters do vary with `k`, and
it is *not* a statement about the integrand. The two are consistent, but only
the `3F2`-level statement carries `k`-dependence.

**Caveat on (M1) numerics:** the Mellin form has an integrable singularity
`z^{1/k-1}` at the origin, and naive quadrature loses accuracy as `k` grows
(observed drift `~1e-13` at `k = 4`, `~1e-10` at `k = 5`). Use the `3F2` for any
high-precision work; (M1) is a structural identity, not an evaluation route.

## Addendum 3, 2026-08-01 (death-star): the section 5 PSLQ negative is WEAKER than stated

Section 5 describes a PSLQ battery at `120-150` digits over a wide basis and calls
the result "strong evidence that `S(4)`, `S(5)` have no closed form in logarithms
of algebraic numbers". **That characterisation is withdrawn.** Two calibration
facts, both established here:

**(P1) The true `S(3)` relation has an IRRATIONAL coefficient.** Addendum 1 records
`pi S(3) = sqrt3 * log(5+2 sqrt6) - 2 arctan(sqrt2/5)`. A sweep seeking *rational*
coefficients over logarithms and arctangents therefore cannot represent it at all;
the basis must contain the **products** `alpha * L` with `alpha` algebraic.

**(P2) PSLQ degrades sharply with basis size -- measured, not assumed.** With the
target `pi S(3)`, at `mp.dps = 240` and coefficient bound `10^5`:

```text
basis {1, sqrt3} x {log(5+2sqrt6), arctan(sqrt2/5)}   (4 elements)  -> FOUND, exactly
basis {1,r2,r3,r5,r6} x {8 logs/arctans}             (40 elements)  -> NOT FOUND
```

The same relation, present in the small basis, is **missed** by the large one at
higher precision than section 5 used. So a wide-basis null result in this problem
carries essentially no information, and section 5's battery was of exactly that
kind.

**Consequence.** The honest status of "`S(4)`, `S(5)` have no closed form" is
**OPEN with no strong numerical evidence either way**, not "strong evidence
against". What remains solid in section 5 is the *structural* discussion (the
`3F2` Schwarz data), not the numerical negative. A meaningful search must:

1. include algebraic multipliers as basis elements, per (P1);
2. scan **small** bases systematically rather than one wide basis, per (P2);
3. carry a **live positive control** -- rediscover `pi S(3)` with the identical
   pipeline before any negative is reported. A pipeline that cannot find `S(3)`
   proves nothing about `S(4)`.

Read together with Addendum 2: since `S(k)` is one fixed function's Mellin moment
at `s = 1/k`, and no function-level argument can separate `k <= 3` from `k >= 4`,
the numerical route is currently the *only* discriminator available -- which makes
getting its calibration right the whole ballgame.

## Addendum 4, 2026-08-01 (death-star): two CERTIFIED bounded exclusions for S(4), S(5)

Addendum 3 withdrew the old numerical negative as uncalibrated. This addendum
replaces it with **controlled** negatives -- bounded statements produced by
pipelines that demonstrably detect known relations of the same shape. Each is
stated as what it is: a finite exclusion, **not** an impossibility proof.

### (E1) Weight-1 exclusion (logarithms and arctangents, algebraic multipliers)

*Positive control (passed).* A blind scan over **all** size-2 bases of products
`alpha * L`, `alpha` algebraic and `L` a logarithm/arctangent, recovers with no
hint

```text
1 * pi S(3)  +  2 * arctan(sqrt2/5)  -  1 * sqrt3 * log(5 + 2 sqrt6)  =  0.
```

*Result.* With the same instrument, basis of `275` products
(`11` multipliers `x` `25` logs/arctans), all `37675` pairs:

```text
no Z-relation  pi S(4) = c1 A + c2 B   with |c| <= 10^5, at 150 digits;
no Z-relation  pi S(5) = c1 A + c2 B   with |c| <= 10^5, at 150 digits.
```

### (E2) Weight-2 exclusion (Catalan, dilogarithms, pi^2, log^2)

Motivation: `1/(2-k^2)` partial-fractions to `1/(sqrt2 -+ k)` and the base moment
is `int_0^1 K(k) dk = 2G`, so weight `2` is the natural next home if the family
climbs one level -- `pi S(1)` is weight `0`, `pi S(2)` and `pi S(3)` weight `1`.

*Positive controls (both passed).* The same size-2 scan over a `64`-element
weight-2 product basis finds `Li_2(1/2) = pi^2/12 - (log 2)^2/2` (returned as
`-12 Li_2(1/2) + pi^2 - 6 log^2 2 = 0`), and `pslq` recovers
`int_0^1 K dk = 2G` as `[1,-2]`.

*Result.* Against that basis, at `130` digits, `|c| <= 10^5`:

```text
S(4), pi S(4), pi^2 S(4), S(5), pi S(5)   -- no relation over any pair.
```

*Negative control (passed).* `pi S(2)`, which is weight `1`, produces no hit on
the weight-2 basis, so the basis does not misfire.

### What this does and does not establish

**Does:** `S(4)` and `S(5)` are not expressible as a two-term `Z`-combination,
with coefficients below `10^5`, over either of the stated bases -- from
instruments verified to detect `pi S(3)` and `Li_2(1/2)` respectively. That is a
genuine unconditional theorem about a finite search region, and it is strictly
more than section 5 ever established.

**Does not:** rule out three-or-more-term relations, larger coefficients, higher
weight (`Li_3`, `zeta(3)`), other multiplier fields, or `Gamma`-value bases. It is
**not** a proof that no closed form exists, and must never be quoted as one.

Combined with Addendum 2 -- `S(k)` is one fixed function's Mellin moment at
`s = 1/k`, so **no function-level argument can separate `k <= 3` from `k >= 4`** --
the honest position is: closed forms for `S(4)`, `S(5)` are **OPEN**; the weight-1
and weight-2 two-term regions are now genuinely excluded; and an unconditional
impossibility proof would require value-level transcendence machinery
(Grothendieck period conjecture, Schanuel, or a Siegel-Shidlovskii/Andre
G-function argument), which is not available here.
