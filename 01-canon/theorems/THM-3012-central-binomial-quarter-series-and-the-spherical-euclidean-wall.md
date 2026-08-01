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
script:
  - 04-computation/central_binomial_quarter_series_thm3012.py
  - 04-computation/sk_S4_lemniscatic_eisenstein_reduction_thm3012.py
  - 04-computation/sk_S4_lambda_bounded_exclusion_thm3012.py
output:
  - 05-knowledge/results/central_binomial_quarter_series_thm3012.out
  - 05-knowledge/results/sk_S4_lemniscatic_eisenstein_reduction_thm3012.out
  - 05-knowledge/results/sk_S4_lambda_bounded_exclusion_thm3012.out
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

## Addendum 5, 2026-08-01 (death-star, S4-gamma lane): S(4) IS lemniscatic -- the exact reduction

Addenda 3-4 left `S(4)` OPEN with two bounded exclusions and no structure. This
addendum supplies the structure. **`S(4)` does contain the lemniscatic constant
`Gamma(1/4)^2`, exactly, with an algebraic coefficient.** What hid it is that the
factor multiplying `Gamma(1/4)^2` is *Catalan's `G` plus a geometrically
convergent Eisenstein tail* -- a shape absent from every basis used above.

### (R1) The exact reduction (PROVED; refereed to 200 digits)

```text
S(4) = (2 sqrt2/pi^2) K(1/sqrt2) Lambda = sqrt2 Gamma(1/4)^2 Lambda / (2 pi^{5/2}),
```

equivalently

```text
pi S(4) = (sqrt2 Gamma(1/4)^2/pi^{3/2}) [ G + 2 sum_{n>=1} chi_{-4}(n)/(n^2(e^{pi n}-1)) ],
```

with `K(1/sqrt2) = Gamma(1/4)^2/(4 sqrt pi)` the lemniscatic singular value
(`tau = i`) and

```text
Lambda = 2G + D,      D = int_{1/sqrt2}^1 [K(k) - K'(k)] dk/k
       = 2 sum_{n odd} (-1)^{(n-1)/2} coth(n pi/2)/n^2
       = 5 pi^2/24 - (1/2) sum_{m>=1} sech(pi m)/m^2
       = 2.01255816455161437222338092610049050178560412519436787901575898621705...
D      = 4 sum_{m>=1} Ti_2(e^{-pi m}) = 0.180626976197176342114173896235722280237...
```

Both tails are geometric in `e^{-pi} = 0.0432...`: the `D`-series (odd `n` only)
gains `~2.7` digits per term, so `~37` terms give `100` digits, and the `sech`
series `~74`. Contrast the defining series' `~1/N` tail (the trap of addendum 1):
`N = 4000` there leaves `~1e-5`.

### (R2) Proof chain

1. **Moment recurrence (exact).** For `M_{2m} = int_0^1 k^{2m}K(k)dk` the Legendre
   ODE `(k(1-k^2)K')' = kK` plus two integrations by parts (boundary term
   `E(1) = 1`) gives
   ```text
   4 m^2 M_{2m} = (2m-1)^2 M_{2m-2} + 1,     M_0 = 2G.
   ```
   Hence `M_{2m} = 2 c_m G + b_m`, `c_m = [(1/2)_m/m!]^2`, `b_m` **rational** with
   `4m^2 b_m = (2m-1)^2 b_{m-1} + 1`, `b_0 = 0`.
2. **Where the lemniscatic constant comes from.** `1/(2-k^2) = (1/2)sum_m (k^2/2)^m`
   and `sum_m c_m 2^{-m} = 2F1(1/2,1/2;1;1/2) = (2/pi)K(1/sqrt2)` give
   ```text
   I := int_0^1 K(k)/(2-k^2) dk = (2/pi) K(1/sqrt2) * G + beta,
   beta = (1/2) sum_m b_m 2^{-m}.
   ```
   The coefficient of `G` is the *generating function of the Catalan parts of the
   moments evaluated at the singular point* `x = 1/2`. That is the whole trick.
3. **The Legendre relation supplies the rest.** `B(x) = sum_m b_m x^m` solves
   ```text
   x(1-x)B'' + (1-2x)B' - B/4 = 1/(4(1-x)),
   ```
   an inhomogeneous hypergeometric equation with homogeneous solutions
   `y1 = (2/pi)K(sqrt x)`, `y2 = (2/pi)K(sqrt(1-x))` and Wronskian
   `W = -1/(pi x(1-x))` -- **this Wronskian is exactly the Legendre relation
   `E K' + E' K - K K' = pi/2`**. Because `f/W = -pi/(4(1-x))` is elementary,
   variation of parameters closes in one line:
   ```text
   B(x) = (pi/4)[ y1 int_0^x y2 dt/(1-t) - y2 int_0^x y1 dt/(1-t) ].
   ```
   At `x = 1/2` the two solutions coincide and `k -> k'` folds the two inner
   integrals into one: `beta = (K(1/sqrt2)/pi) D`.
4. **Modular normal form.** With `x = lambda(tau)`, `tau = i K'/K`, one has
   `K dk/k = (i pi^2/4) theta_3^2 theta_4^4 d tau` and `x = 1/2 <-> tau = i`, so
   `D = (pi^2/4) int_1^infty (s-1) theta_2^4 theta_3^2(is) ds`. Both weight-3
   forms are **purely Eisenstein** (checked exactly to `q^159`, `q = e^{i pi tau}`):
   ```text
   theta_2^4 theta_3^2 = 16 sum_n B_n q^n,  B_n = sum_{d|n} chi_{-4}(n/d) d^2,
   theta_3^2 theta_4^4 = 1 - 4 sum_n A_n q^n, A_n = sum_{d|n} chi_{-4}(d) d^2.
   ```
   Therefore `D = 4 sum_{n odd} chi_{-4}(n)/(n^2(e^{pi n}-1))`.
5. **The Fricke relation.** Residues of `pi cot(pi z) sech(pi z)/z^2` (poles at
   `z in Z`, at `z = i(k+1/2)`, and `Res_0 = -5 pi^2/6`) give
   `Lambda = 5 pi^2/24 - (1/2) sum_m sech(pi m)/m^2`.

### (R3) Why every earlier sweep in this file missed it

The leading term of `pi S(4)` is `sqrt2 Gamma(1/4)^2 G/pi^{3/2}` (about `91%` of
the value): a **product of the lemniscatic constant with Catalan's constant**.
Addendum 4's weight-2 basis contained `G`, `pi^2`, `log^2`, `Li_2` but no
`Gamma(1/4)`-times-`G` product; the weight-1 basis contained neither. So (E1) and
(E2) were excluding the wrong region -- correct as stated, but blind to the
answer. This is the addendum-1 lesson one level up: the missing basis element is
again a **product**, now of two transcendentals rather than an algebraic times a
log. It also explains the sibling result in
`04-computation/amm_sk_S4_oneD_elliptic_moment_opus_S4.py` ("`pi S(4)` independent
of `{K(i), lemniscate, Catalan, pi}` at 40 digits"): true, and consistent, because
`Lambda` is not a rational combination of those.

### (R4) What remains, and the structural obstruction (EVIDENCE, not proof)

`Lambda` is a **weight-3 Eisenstein Eichler value at the CM point `tau = i`**.
In `z = tau/2` that point is `z = i/2`, the fixed point of the Fricke involution
`W_4: z -> -1/(4z)`, which supplies **exactly one** linear relation between the
two Eichler values `xi_B(i) = D/4` and `xi_A(i) = U/2`,
`U = sum_m sech(pi m)/m^2`:

```text
4 xi_B(i) + xi_A(i) = 5 pi^2/24 - 2G.
```

`W_16` does not fix `z = i/2`, so no second relation is produced by the group.
Two unknowns, one relation: on the available modular mechanism `Lambda` is not
determined. Contrast the weight-2 analogue one level down,
`sum_{n>=1} 1/(n(e^{2 pi n}-1)) = log 2 + (3/4)log pi - log Gamma(1/4) - pi/12`,
which *is* closed form -- because a **first** Eichler integral of a weight-2
Eisenstein series is `log eta`, and `eta` at a CM point is Chowla-Selberg. A
**second** Eichler integral has no such handle. This is an obstruction argument
about available mechanisms, **not** an impossibility proof.

### (R5) Bounded exclusion for `Lambda` (CERTIFIED, live controls)

*Positive controls, all passed by the identical instrument:* `pi S(3)`
rediscovered blind over `154` products `alpha*L`, all `11781` pairs, at `130`
digits; `sum 1/(n(e^{2 pi n}-1))` -- the weight-2 Eisenstein Eichler value at
`tau = i`, the direct analogue one weight down -- recovered as a size-4 relation
over `{pi, log2, log pi, log Gamma(1/4)}`; `sum coth(pi n)/n^3 = 7 pi^3/180`;
`sum 1/(n^3(e^{2 pi n}-1)) = 7 pi^3/360 - zeta(3)/2`;
`int_0^1 K'^3 dk = Gamma(1/4)^8/(128 pi^2)`;
`int_0^1 K^3 dk = 3 Gamma(1/4)^8/(1280 pi^2)`; `K(sqrt2-1)` in the
`Gamma(1/8)Gamma(3/8)` sector; `K((sqrt3-1)/(2 sqrt2)) = 3^{1/4}Gamma(1/3)^3/(2^{7/3} pi)`.

*Result.* Three tiers, script
`04-computation/sk_S4_lambda_bounded_exclusion_thm3012.py`, output
`05-knowledge/results/sk_S4_lambda_bounded_exclusion_thm3012.out`. `K` below is
`K(1/sqrt2)`, `Ls = log(1+sqrt2)`.

```text
TIER 1  targets Lambda, U, S(4).  size <= 2 over 18 simple atoms
        {1, pi, pi^2, pi^3, G, piG, K, K^2/pi, K^2/pi^2, pi^2/K^2, pi^3/K^2,
         KG/pi, zeta(3), log2, Ls, K^4/pi^3, sqrt2, sqrt2 pi^2}
        |c| <= 10^8, 10^12, 10^20;  300 dps, 190-digit tol   -> NO relation
TIER 2  targets Lambda, S(4).  ALL 5456 size-3 subsets of a 33-atom core
        (1, pi, pi^2, pi^3, G, piG, G/pi, K, K/pi, K^2/pi, K^2/pi^2, pi/K,
         pi^2/K^2, pi^3/K^2, KG/pi, K^2G/pi^2, K^4/pi^3, pi^4/K^4, log2, Ls,
         log2^2, Ls^2, pi log2, pi Ls, zeta(3), zeta(3)/pi, Gamma(1/3)^6/pi^4,
         log Gamma(1/4), log pi, sqrt2, sqrt2 K^2/pi, sqrt2 pi^2, sqrt2 G)
        |c| <= 10^5;  210 dps, 135-digit tol                 -> NO relation
TIER 3  target Lambda.  ALL 19110 pairs of the 196 products alpha*atom over the
        49-atom pool (adds Gamma(1/8)Gamma(3/8) and its square/pi,
        Gamma(1/3)^3/pi, zeta(3)/K, more K^j/pi^i, G-products, log products),
        alpha in {1, sqrt2, sqrt3, 2^{1/4}}
        |c| <= 10^5;  210 dps, 135-digit tol       -> EXTENSION, see output file
```

**Certified region.** Tiers 1 and 2 are complete scans and are what this addendum
claims. Tier 1 completed inside the recorded run; tier 2's exact region was also
completed independently at `260` dps / `160`-digit tolerance (all `5456` subsets,
targets `Lambda` and `S(4)`, no relation) before the recorded run. Tier 3 is a
strictly larger region whose completed record is written by the same script into
`05-knowledge/results/sk_S4_lambda_bounded_exclusion_thm3012.out`; read the
verdict there, not here.

(`S(4)` and `U` are exact rescalings of `Lambda` by `2 sqrt2 K/pi^2` and
`5 pi^2/12 - 2(.)`, so a scan of `Lambda` covers them for any basis closed under
those factors; tiers 1-2 scan the targets directly anyway.)

**Does:** exclude a finite region for a genuinely new and better-posed target.
**Does not:** prove `Lambda` has no closed form. Never quote it as one.

### (R5b) Independent concurrent confirmation

A sibling lane (`04-computation/sk_S4_lemniscatic_lambert_thm3012.py`) reached the
same lemniscatic prefactor by a different route -- Jacobi `sn`/`dn` moments rather
than the `K`-moment recurrence -- in the form
`S(4) = varpi/4 + (2 varpi/pi^2)(2V - T)` with `varpi = sqrt2 K(1/sqrt2)`,
`T = sum_{n odd} 1/(n^2 cosh(pi n))`, `V = sum_{n>=0}(-1)^n/((2n+1)^2 sinh((2n+1)pi/2))`.
That is exactly this addendum's statement, since (verified to `50` digits, residual `0`)

```text
Lambda = pi^2/8 + 2V - T.
```

The two derivations agree and are independent; only the present one exhibits the
`G` (Catalan) content, i.e. the leading term `sqrt2 Gamma(1/4)^2 G/pi^{3/2}`.

### (R6) By-product for the whole `K`-moment family

**Route audit (recorded stopping reasons).** Two of the obvious attacks on
`int_0^1 K/(2-k^2)dk` are dead ends and should not be retried:

* *Partial fractions* `1/(2-k^2) = (1/(2 sqrt2))[1/(sqrt2-k) + 1/(sqrt2+k)]`
  expands into `log(1+sqrt2) * sum_n c_n 2^n - (correction)`, and
  `sum_n c_n 2^n` **diverges** (`c_n ~ 1/(pi n)`): the cancellation is essential,
  the split is formal only.
* *Fourier-Legendre.* With `K(sqrt x) = sum_n (2/(2n+1))P_n(2x-1)`,
  `I = sum_n (2/(2n+1)) J_n` with
  `J_n = int_0^1 P_n(2y^2-1)/(2-y^2)dy = rho_n + P_n(3) log(1+sqrt2)/sqrt2`,
  `rho_n` rational (`0, -2, -8, -118/3, -200, ...`; verified `n <= 8`). But
  `P_n(3) ~ (3+2 sqrt2)^n`, so separating the `log(1+sqrt2)` part **diverges**
  again; the FL series converges only as a whole, and slowly (`n <= 8` gives
  `1.18928` against `I = 1.18775`).

The route that works is the one above: the moment *recurrence* (which isolates
`G` with an exact coefficient) plus the Legendre relation as a *Wronskian*.

The same argument at general `x` gives, for `|x| < 1`,

```text
int_0^1 K(k)/(1 - x k^2) dk = (4/pi) K(sqrt x) * G + B(x),
```

with `B(x)` as in step 3. So **Catalan's constant enters every such moment with
coefficient `(4/pi)K(sqrt x)`**; `x = 1/2` is distinguished only by being the
lemniscatic singular point, where `K(sqrt x) = K(sqrt(1-x))` and the two-term
`B(x)` collapses to a single integral.

Referee: `04-computation/sk_S4_lemniscatic_eisenstein_reduction_thm3012.py`,
all checks passed. QED for (R1)-(R2); (R4) is an obstruction argument; (R5) is a
bounded exclusion.

## Addendum 6, 2026-08-01 (death-star, S4-contiguous lane): contiguity reduces at k=1 and nowhere else; S(4) in one balanced 3F2; a new closed form at k=2

**Headline.**

```text
(a)  A classical 3F2(1) summation or contiguous reduction applies to the family
     -- anywhere in its COMPLETE Thomae orbit -- iff k = 1.          [PROVED]
(b)  S(4) = varpi/2 - W/(2 varpi),  W = 3F2(3/4,3/4,1/2 ; 3/2,3/2 ; 1),
     varpi the lemniscate constant.               [VERIFIED-EXACT, 170 digits]
(c)  3F2(3/4, 3/4, 1/4 ; 3/2, 5/4 ; 1) = varpi (pi - 2 log(1+sqrt2))/pi
     -- a NEW closed form, and the product-of-transcendentals calibration
     object no CERTIFIED region of this file has ever covered.  It is now a
     required control for this lane.                       [VERIFIED-EXACT]
(d)  Contiguity's exact content is a FIRST-ORDER Mellin ladder
     (s+1/4)(s+3/4) M(s+1) = s^2 M(s) + 1/(pi sqrt2);  the integer ladder is
     Gauss-anchored, the ladder through 1/k (k>=2) is not.           [PROVED]
(e)  Bridge to the concurrent S4-quadratic lane, exact:
     W_q = varpi log(1+sqrt2)/2 + pi W/(4 varpi).  The two lanes isolate the
     SAME obstruction; W is the cleaner representative.    [VERIFIED-EXACT]
S(4) and S(5) remain OPEN.  Nothing here proves non-existence.
```


The configuration `3F2(1/4, 3/4, 1/k ; 1, 1+1/k ; 1)` is *contiguous*: the upper
parameter `1/k` sits one below the lower `1+1/k`; and at `k = 4` that upper
`1/k = 1/4` additionally *repeats* the upper `1/4`. This addendum settles what
the degeneracy buys, exactly.

Scripts and outputs (all four in `04-computation/` and `05-knowledge/results/`):

```text
sk_S4_contiguous_thomae_thm3012.py         structural, exact: orbit + rule scan
                                           + the Mellin ladder
sk_S4_contiguous_pslq_thm3012.py           170-digit verification + integralities
sk_S4_contiguous_pslq_weight1_thm3012.py   the B1/B3 bounded exclusions
sk_S4_contiguous_product_basis_thm3012.py  the B4 PRODUCT-basis exclusion, and the
                                           new closed form of (C7)
```

Each script writes its own `.out` of the same basename into
`05-knowledge/results/`; the two PSLQ scripts are long (the pair scans are
`~4.5 x 10^4` PSLQ calls per target at 170 dps), so their outputs are
regenerable rather than precious. Everything PROVED or VERIFIED-EXACT below is
in the first two scripts, which run in minutes.

### (C1) The complete Thomae orbit, symbolically in `x = 1/k` (PROVED)

Every parameter is an exact affine function of `x` over `Q`. Closing the family
under Thomae's transformation and all `3! * 2!` permutations gives **7** distinct
parameter multisets. The classical count is `120/(3! 2!) = 10`, and the generator
was validated on three generic seeds -- generic, generic *balanced*, and generic
`x`-dependent -- each of which returns exactly 10. So the drop to 7 is a genuine
extra degeneracy of this family (it is caused by the lower parameter `1`), not an
artefact of balancedness or of the implementation:

```text
[0] 3F2(1-x, 1/4, 1/4 ; 1, 5/4 ; 1)          s = x+3/4
[1] 3F2(1-x, 3/4, 3/4 ; 1, 7/4 ; 1)          s = x+1/4
[2] 3F2(1-x, 1, 1 ; 5/4, 7/4 ; 1)            s = x     <- eq (6) of this file
[3] 3F2(1/4, 3/4, x ; 1, x+1 ; 1)            s = 1     <- the defining form
[4] 3F2(1/4, 1, x+1/4 ; 5/4, x+1 ; 1)        s = 3/4
[5] 3F2(3/4, 1, x+3/4 ; 7/4, x+1 ; 1)        s = 1/4
[6] 3F2(x, x+1/4, x+3/4 ; x+1, x+1 ; 1)      s = 1-x
```

Each form, times its Gamma prefactor, was checked numerically against `S(k)` for
`k = 4, 5, 7` (max deviation `~1e-26` at `dps = 25`).

### (C2) THE ANSWER: a classical reduction exists iff k = 1 (PROVED, exhaustive)

Six rules were applied to **all 7 forms**, each solved as an *exact linear
equation in `x`* -- so the result is a statement about all real `k >= 1`, not a
scan:

* **G** Gauss cancellation (an upper equals a lower) -> `2F1(1)` = Gamma quotient;
* **R** contiguous reduction (upper - lower `= m in Z_{>=1}`) -> finite sum of `2F1(1)`;
* **T** termination (an upper in `Z_{<=0}`), which with balancedness is Saalschutz;
* **D** Dixon (well poised); **W** Watson; **P** Whipple.

```text
hits: k = 1 ONLY, and there FIVE of the six rules fire at once
      (G on forms [3],[4],[5];  T on [0],[1],[2];  D on [0],[1],[4],[5];
       W on [3],[4],[5],[6];    P on [2],[3])
no hit for any k >= 2, in particular none for k = 4 and none for k = 5.
```

*Positive control:* `k = 1` is detected, i.e. `S(1) = 2F1(1/4,3/4;2;1) = 8sqrt2/(3pi)`.
*Negative controls:* `k = 2, 3` are **not** detected -- **yet both have closed
forms.** That is the calibration that makes the `k = 4, 5` null meaningful and
that also caps it: *"no classical 3F2(1) theorem applies" is strictly weaker than
"no closed form exists."*

**Mechanism.** `(1/k)_n/(1+1/k)_n = (1/k)/(n+1/k)` is a rational function of `n`
*with a pole*. The reduction of a `3F2` to `2F1`'s needs this ratio to be a
*polynomial* in `n`, i.e. `1/k - 1 in Z_{>=0}`, i.e. `k = 1`. For `k >= 2` the
`3F2` collapses to a Mellin moment (addendum 2), never to a finite sum.

The same six-rule battery, run on the complete Thomae orbit (5 forms) of the
*residual* family `3F2(3/4, 3/4, 3/4-x ; 3/2, 7/4-x ; 1)` of (C4) below, also
hits at `k = 1` and at no other `k`.

### (C3) The exact content of contiguity: a first-order Mellin ladder (PROVED)

With `M(s) = int_0^1 z^{s-1} 2F1(1/4,3/4;1;z) dz = sum_n c_n/(n+s)` and
`S(k) = M(1/k)/k`,

```text
(s + 1/4)(s + 3/4) M(s+1)  =  s^2 M(s)  +  1/(pi sqrt2).            (L)
```

*Proof.* `c_{n+1}(n+1)^2 = c_n (n+1/4)(n+3/4)`; divide the numerator by `n+s+1`,
`(n+1/4)(n+3/4) = (n+s+1)(n-s) + (s^2+s+3/16)`. Summing `n = 0..N-1` on the left
and `m = 1..N` on the right, the two copies of the (individually divergent)
`sum c_n (n-s)` cancel; the surviving boundary term is `c_N(N-s) -> 1/(pi sqrt2)`
because `c_N ~ 1/(pi sqrt2 N)`. QED

(L) **is** the contiguous relation, in its sharpest form: first order,
inhomogeneous, rational coefficients, inhomogeneity `1/(pi sqrt2)`. Homogeneous
solution `h(s) = Gamma(s)^2/(Gamma(s+1/4)Gamma(s+3/4))`, so there is exactly
**one free constant per ladder `s + Z`**. Consequences:

* the **integer** ladder is Gauss-anchored at `M(1) = S(1) = 8sqrt2/(3pi)`, hence
  every integer moment is rational times `sqrt2/pi`:
  ```text
  r_m := M(m) pi/sqrt2,   r_1 = 8/3,   r_{m+1} = 16(m^2 r_m + 1/2)/((4m+1)(4m+3))
  r = 8/3, 152/105, 10568/10395, 178328/225225, 47453768/72747675,
      782539544/1405485081, 51331850888/105411381075,
      838519635608/1933976154825, ...          (each verified to ~50 digits)
  ```
* the ladder through `s = 1/k`, `k >= 2`, contains **no integer**, hence no Gauss
  point. (L) relates `S(k)` only to `M(1/k + j)`, `j >= 1`, which are not members
  of the `S`-family. **Contiguity propagates `S(k)` but can never compute it.**

(L) was verified numerically at `s = 1, 2, 1/2, 1/3, 1/4, 1/5, 1/7` and at a
generic real `s`, residual `<= 6e-62` in every case.

### (C4) Two new two-term reductions, valid for all k (PROVED)

Euler's integral for `2F1(1/4,3/4;1;z)` plus `M1` makes the inner `u`-integral
`int_0^1 (1-tu^k)^{-1/4} du = 2F1(1/4, x; 1+x; t)` -- itself contiguous, hence an
*incomplete Beta*, for **every** `k`. Two Fubinis (positive integrands) and one
interchange give, with `beta = 3/4 - x` (needs `k >= 2`) and `gam = 1/4 - x`
(needs `k >= 5`):

```text
(A)  (pi sqrt2 / x) S(k) = B(x, 3/4) B(3/4-x, 1/4)
                         - [Gamma(3/4)^2/(beta Gamma(3/2))] 3F2(3/4,3/4,beta; 3/2,1+beta; 1)

(B)  (pi sqrt2 / x) S(k) = B(x, 1/4) B(1/4-x, 3/4)
                         - [Gamma(1/4)^2/(gam Gamma(1/2))]  3F2(1/4,1/4,gam;  1/2,1+gam;  1)
```

Both are of **non-terminating-Saalschutz shape**: a balanced `3F2(1)` equals a
Gamma quotient *minus* a second balanced `3F2(1)`. The residual is **not** in the
Thomae orbit of (C1) (orbits verified disjoint at `k = 4, 5`) -- a two-term
relation is not a one-term one. Verified `k = 2..8` for (A) and `k = 5..9` for
(B) at `dps = 60`.

### (C5) `S(4)` in ONE balanced 3F2 (PROVED; VERIFIED-EXACT to 170 digits)

`k = 4` is exactly where (A) degenerates: `beta = 1/2`, so the lower parameter
`1+beta = 3/2` **equals** the other lower `3/2`. The residual acquires a repeated
upper *and* a repeated lower, and the leading Beta product collapses to the
lemniscate constant. Result:

```text
        S(4)  =  varpi/2  -  W/(2 varpi),
        W     =  3F2(3/4, 3/4, 1/2 ; 3/2, 3/2 ; 1),
        varpi =  Gamma(1/4)^2/(2 sqrt(2 pi))   (lemniscate constant)

W    = 1.267377930871768371750405432734219729951834196668985622433...
varpi= 2.622057554292119810464839589891119413682754951431623162817...
S(4) = 1.069352554441268058582961939534278061325932014506308037563...
```

Residual `0.0` at `dps = 170`, i.e. **VERIFIED-EXACT to 170 digits**, with `W`
computed two independent ways -- the defining `3F2` series and the iterated-Beta
integral of (A) -- agreeing to 131 digits.

**Name clash, resolved.** The concurrent `S4-quadratic` lane (next addendum)
also writes `W` for an object it isolates, the second-kind period
`W_q = int_0^1 V(v)[(1-v^2)^{-1/2}+(1+v^2)^{-1/2}]dv = 1.53513025889742...`,
`V` the lemniscatic arcsine. The two are DIFFERENT and are bridged exactly by

```text
   W_q  =  varpi log(1+sqrt2)/2  +  pi W /(4 varpi),          residual 1.1e-32 (limited by the direct quadrature of W_q; the S4-quadratic lane's own identity then closes to 3.1e-61)
```

which follows by eliminating `pi S(4)` between that lane's
`pi S(4) = varpi(pi/2 + log(1+sqrt2)) - 2 W_q` and (C5). Two consequences:
the two lanes isolate the *same* obstruction, the present `W` being the cleaner
one (it differs from `W_q` only by an elementary lemniscatic term); and `W_q`
becomes available to arbitrary precision through the `3F2`, since the bridge is
exact. The two lanes' decimals for `W_q` agree, independently:

```text
   from the bridge (via W)   1.535130258897420337752379135458933506929...
   S4-quadratic, own route   1.535130258897420337752379135458933506929...
```

Equivalently `W = varpi^2 - 2 varpi S(4)`, and this **rederives addendum 5**:
`W = varpi^2 (1 - 4 Lambda/pi^2)` to the full working precision. The gain is that
the entire transcendental content of `S(4)` is now a *single balanced `3F2` with
repeated parameters* -- **no Catalan constant, no Eichler integral, no Eisenstein
tail**.

### (C6) Why `k = 4` and no other integer (PROVED, four independent criteria)

```text
(i)   orbit collapse   : the 7 forms of (C1) stay distinct at every k in 1..12
                         except k = 1 (the Gauss point) and k = 4, where they
                         collapse to 4 distinct forms.
(ii)  ODE resonance    : the exponents at z = infinity are the upper parameters
                         (1/4, 3/4, 1/k); they are repeated iff 1/k in {1/4,3/4},
                         i.e. iff k = 4 among integers.  k = 4 is the unique
                         integer k at which the equation is RESONANT at infinity
                         (the configuration in which a logarithmic local
                         solution appears).
(iii) the (A) parameter: beta = 3/4 - 1/k equals 1/2 -- the coincident-lower
                         point -- iff k = 4.
(iv)  K-moment weight  : in eq (4) the weight mu^{4/k-1}(2-mu^2)^{-2/k-1/2} is a
                         RATIONAL function of mu iff BOTH exponents are integers:
                         4/k-1 in Z forces k | 4, i.e. k in {1,2,4}; and then
                         2/k+1/2 in Z forces k = 4.  This is exactly the
                         hypothesis under which addendum 5's moment recurrence
                         4m^2 M_{2m} = (2m-1)^2 M_{2m-2} + 1 can be applied.
                         (k = 1, 2 give mu^3 (2-mu^2)^{-5/2} and mu (2-mu^2)^{-3/2}:
                          half-integral, i.e. algebraic over a quadratic -- that is
                          the SPHERICAL mechanism of section 3, a different one.
                          k = 5 gives mu^{-1/5}(2-mu^2)^{-9/10}: denominator 10/5,
                          neither.)
```

`k = 5` fails all four. Its two residuals

```text
W5a = 3F2(3/4, 3/4, 11/20 ; 3/2, 31/20 ; 1) = 1.28783310672123019536336223730224486585...
W5b = 3F2(1/4, 1/4,  1/20 ; 1/2, 21/20 ; 1) = 1.01014347743357758443860712406668353805...
```

carry an irreducible denominator 20: no repeated parameter, no coincident lower,
no orbit collapse, no singular value. **The `k = 4` mechanism has no analogue at
`k = 5`** -- which is a *reason*, not merely a failed search.

### (C7) A NEW closed form -- and the calibration object the file has been missing

Identity (A) at `k = 2` is exactly solvable, because `S(2) = 4 log(1+sqrt2)/pi` is
known and, uniquely among `k >= 2`, the leading Beta product is ELEMENTARY:
`B(1/2,3/4) B(1/4,1/4) = 4 pi sqrt2`.
Solving for the residual gives (VERIFIED-EXACT, residual `0.0` at `dps = 50`,
against the defining `3F2` series and against (A) independently):

```text
   3F2(3/4, 3/4, 1/4 ; 3/2, 5/4 ; 1)  =  varpi (pi - 2 log(1 + sqrt2)) / pi
                                      =  varpi  -  2 varpi log(1+sqrt2)/pi
                                      =  1.150821447753979605424774878953724959491...
```

The same construction gives an explicit closed form for the `k = 3` residual too,
```text
W_3 = 3F2(3/4, 3/4, 5/12 ; 3/2, 17/12 ; 1)
    = [ B(1/3,3/4) B(5/12,1/4) - 3 pi sqrt2 S(3) ] (5/12) Gamma(3/2)/Gamma(3/4)^2
    = 1.231375424470671589432753496433802176703...      (checked to 50 digits)
```
but that one lives over `Gamma(1/3), Gamma(5/12), Gamma(13/12), Gamma(2/3)` and is
not a useful calibration object. `W_2` is, because it is clean.

Beyond being a new evaluation, **`W_2` is the calibration object this lane has
needed.** It is a *product of two transcendentals* -- the lemniscate constant
times a logarithm -- and it lies outside every **certified** region of this file:

* addendum 4's weight-1 basis `E1` (`algebraic x log`) cannot represent it;
* addendum 4's weight-2 basis `E2` (`G, pi^2, log^2, Li_2`) cannot represent it;
* addendum 5's tiers 1-2 contain `K(1/sqrt2)` and `log(1+sqrt2)` as SEPARATE atoms
  but not their product, and `varpi = sqrt2 K(1/sqrt2)`, so they cannot either.

(For the record: section 5's *wide* sweep did list `varpi log(1+sqrt2)` among its
atoms -- but that is precisely the battery addendum 3 withdrew as uncalibrated,
so it certifies nothing. No certified region has ever covered this shape.)

Yet it is an exact value of *exactly the object* those sweeps were probing. Any
future search on `S(4)`, `S(5)` or their residuals must carry a basis of
**products `Gamma`-type x elementary** and must demonstrate that it recovers this
`W_2` before its null is quoted. That is now a required control in this lane,
alongside `pi S(3)`.

### (C8) Bounded exclusions, each with its own live control (CERTIFIED)

Three instruments, three separately-stated regions. `H = 10^6`, `P = 170` dps.

```text
B1   {algebraic} x {1, pi, logs, arctans}          96 elements,  4560 pairs
     CONTROL 3/5: rediscovers pi S(3) BLIND, with its irrational sqrt3
     coefficient, from an unhinted scan of all pairs; also pi S(1), pi S(2).
     (It has no Gamma atoms, so varpi^2 and Gamma quotients are outside it.)
B3   B1 u ({algebraic} x {varpi, varpi^2, Catalan, zeta(3), Gamma(1/5)- and
     Gamma(1/8)-products, the Beta values of (A)/(B)})
                                                  282 elements, 39621 pairs
     CONTROL 5/5: pi S(3), pi S(1), pi S(2), varpi^2, a pure Gamma quotient.
     B1 and B3 contain algebraic x elementary and algebraic x Gamma products,
     but NOT Gamma x elementary products -- that class is B4's job.
B4   {algebraic} x {1, varpi, varpi^2, B(1/5,3/4)B(11/20,1/4), Gamma(1/5)-
     products, Gamma(1/8)Gamma(3/8)} x {1, 1/pi, 1/pi^2, log/pi, arctan/pi,
     G/pi^2}  -- the PRODUCT basis      308 elements (deduped), 47278 pairs
     CONTROL 5/5, and it is the decisive one: B4 rediscovers
     W_2 = varpi - 2 varpi log(1+sqrt2)/pi  BLIND -- a transcendental times a
     transcendental -- and also S(3), S(2), S(1), varpi^2.  (B4 normalises
     everything by 1/pi, so it is a basis for S(k), not for pi S(k).)
```

`B2` (a `Gamma`/weight-2 basis with no logarithm atoms) is deliberately reported
as **uninformative rather than negative**: it fails the `pi S(3)` control, i.e. it
is the exact failure mode addendum 3 diagnosed. Its scan is not quoted.

The excluded region is: *no* `T = c1 b_i + c2 b_j` with `c1, c2 in Z`,
`|c| <= 10^6`, `b_i, b_j` in the stated basis, **all pairs tested**, at 170
decimal digits with detection tolerance `1e-122` (`B4`: `1e-120`).

**Certified here** (scans that completed inside this session):

```text
B1, ALL EIGHT targets -- W, W_5a, W_5b, S(4), pi S(4), S(5), pi S(5), pi^2 S(5)
    -> no relation over any of the 4560 pairs, for any target.
B3, target W (the k = 4 residual)  -> no relation over any of the 39621 pairs.
B4, target W_5 (the k = 5 residual, the decisive one for S(5), and the basis
    whose control is the W_2 product) -> no relation over any of the 47278 pairs.
```

The remaining `B3` and `B4` targets form a strictly larger region; the completed
record is written by the two scripts into their `.out` files -- **read the
verdict there, not here.**

**This is a finite exclusion, not an impossibility proof, and must never be
quoted as one.**

### (C9) Scope

**PROVED:** (C1) the orbit; (C2) classical reducibility iff `k = 1`; (C3) the
Mellin ladder (L) and its two consequences; (C4) identities (A) and (B); (C6) the
four integrality criteria singling out `k = 4`.
**VERIFIED-EXACT:** (C5) `S(4) = varpi/2 - W/(2 varpi)` to 170 digits; (C7)
`3F2(3/4,3/4,1/4;3/2,5/4;1) = varpi(pi - 2 log(1+sqrt2))/pi`.
**Bounded exclusion:** (C8).
**OPEN:** closed forms for `S(4)` and `S(5)`. Nothing here proves non-existence.
The contribution to that question is (C2)+(C3): the *mechanism* of the `k <= 3`
evaluations -- reduction of a contiguous `3F2` to Gauss values -- is available at
`k = 1` and provably at no other `k`, and the contiguous ladder (L) can propagate
`S(k)` but never compute it.

## Addendum 7, 2026-08-01 (death-star, S4-quadratic lane): the quadratic-surd basis is closed off for `pi S(4)`, with a working positive control

*(Renumbered from "Addendum 2" by the S4-contiguous lane: that number was already
taken by the Mellin-moment addendum above. Content untouched.)*

The previous addendum warned that the coefficient in `pi S(3)` is the
*irrational* `sqrt3`, so that any PSLQ sweep seeking **rational** coefficients
over logarithms and arctangents returns "no closed form" spuriously. This
addendum redoes the `k = 4` search with that defect repaired: the basis consists
of **products `alpha * L`**, `alpha` a quadratic surd and `L` a
logarithm/arctangent/`pi`/Catalan, and the pipeline is required to rediscover
`pi S(3)` blind before any negative is reported. It does.

Script `04-computation/s4_quadratic_relation_lane.py`, output
`05-knowledge/results/s4_quadratic_relation_lane.out`.

### 1. Two corrections to the record

**(a) The `S(3)` decimal printed in Addendum 1 is wrong.** It states "the true
`S(3) = 1.08838540640395`". The correct value is

```text
S(3) = 1.0884041641172768712701774968372772011989808537277...
```

The *closed form* in Addendum 1 is right; only the decimal is mistyped. Equation
(8) reproduces `1.0884041641...` to **460 digits**, which also re-confirms (8).

**(b) The owner's elliptic representation is now confirmed to 460 digits**, not
`~30`: with `S(k)` evaluated by the elementary-kernel quadrature

```text
S(k) = (2/pi) int_0^1 (1+t^{k/2})^{-1/2} K( 2 t^{k/2}/(1+t^{k/2}) ) dt
```

(`K` in the `m = k^2` convention -- this is (3) rewritten, and it is the
evaluation route to use; it is fast and has only an integrable log singularity
at `t = 1`), the difference `S(4)_kernel - (2 sqrt2/pi) int_0^1 K(k)/(2-k^2) dk`
is `0.0` at 460 digits.

### 2. The basis, and the two elements that must be stripped

`alpha` ranges over `{1, sqrt2, sqrt3, sqrt5, sqrt6, sqrt10, sqrt15, sqrt30}` --
a `Q`-basis of `Q(sqrt2,sqrt3,sqrt5)` -- and `L` over `1, pi, log2, log3, log5,
log(1+sqrt2), log(2+sqrt3), log(5+2sqrt6), log(phi), atan(sqrt2), atan(1/sqrt2),
atan(sqrt2/5), Catalan`. Two of the thirteen are `Q`-dependent on the others,

```text
atan(1/sqrt2) = pi/2 - atan(sqrt2),      atan(sqrt2/5) = pi - 3 atan(sqrt2),
```

and were **stripped**. (Left in, the search returns exactly the degeneracy
`-alpha*atan(sqrt2) + 2 alpha*atan(1/sqrt2) - alpha*atan(sqrt2/5) = 0`; verified
-- PSLQ at 400 digits on the unstripped 104-element basis returns precisely that
and nothing else.) **88 basis elements remain.**

### 3. Method: exact integer LLL with a rigorous exclusion bound

Relations are sought by exact integer LLL (de Weger, `delta = 3/4`, all
arithmetic in Python `int`s) on the lattice spanned by
`r_i = (e_i | round(10^P x_i))`. Beyond reporting relations, each subset yields a
**rigorous** exclusion bound: `delta = 3/4` guarantees
`|b_1| <= 2^{(n-1)/2} lambda_1`, and a relation with `max|c_i| <= H` would put a
lattice vector of norm `<= sqrt(n) H (1 + sqrt(n)/2)` in the lattice, so no
relation exists with

```text
max|c_i|  <=  |b_1| / ( 2^{(n-1)/2} sqrt(n) (1 + sqrt(n)/2) )  =:  H_cert.
```

`H_cert` is independent of any acceptance cap, and is minimised over all subsets
scanned. **That minimum is the entire content of the null result.**

### 4. The positive control passes (blind)

The identical sweep, given `pi S(1)`, `pi S(2)`, `pi S(3)` and told nothing:

```text
size 1, P = 200:  +1*[pi S(2)] -4*[log(1+sqrt2)] = 0
                  +3*[pi S(1)] -8*[sqrt2]        = 0
size 3, P =  60:  -1*[pi S(3)] -2*[pi] +6*[atan(sqrt2)] +1*[sqrt3*log(5+2sqrt6)] = 0
```

all C(88,3) = 109736 size-3 subsets were scanned and **that is the only relation
reported, for either target**. The recovered `k = 3` form carries the irrational
multiplier `sqrt3`, i.e. precisely the feature a rational-coefficient basis
cannot see. So the negatives below are meaningful.

A second, sharper control comes from the size-4 sweep run deliberately with an
acceptance cap *above* its precision-supported level: of the 143 reported hits,
exactly **one** has `max|c| < 10^5` -- the genuine `pi S(3)` relation, with
`max|c| = 6` -- and the other 142 all have 8-digit coefficients in
`[2.7 x 10^7, 9.9 x 10^7]`, i.e. pure lattice noise. The true relation is
separated from noise by seven orders of magnitude, and **no `pi S(4)` hit has
`max|c|` below `2.7 x 10^7`.**

### 5. The negative for `pi S(4)` (BOUNDED, not absolute)

Over the 88-element basis `B` of section 2, with `P` digits of working
precision, no integer relation involving `pi S(4)` exists with

```text
 support size 1 :  P = 200 dig,       88 subsets,  max|c| <= 1.24 x 10^98
 support size 2 :  P = 200 dig,    3 828 subsets,  max|c| <= 3.42 x 10^64
 support size 3 :  P =  60 dig,  109 736 subsets,  max|c| <= 1.17 x 10^12
```

and over the 44-element sub-basis with `alpha in {1,sqrt2,sqrt3,sqrt6}` -- the
field `Q(sqrt2,sqrt3)` in which all three known closed forms live -- no relation
of support size 4 among C(44,4) = 135751 subsets at `P = 45` (see the output
file for the final bound).

**This is not a proof that no closed form exists.** It is exactly: *no
`Z`-linear relation over the stated basis, of the stated support size, with
coefficients below the stated bound, at the stated precision.*

### 6. The lemniscatic lane is closed off too

`k = 4` is the lemniscatic case: with `varpi = Gamma(1/4)^2/(2 sqrt(2 pi))`,

```text
K(m = 1/2) = varpi/sqrt2,     int_0^1 dv/sqrt(1-v^4) = varpi/2,
```

so the natural shape of a `k = 4` closed form is `varpi x (weight-one
constant)`, not a weight-one constant itself. Taking the target to be
`pi S(4)/varpi` against the *same* 88-element basis gives no relation with
`max|c| <= 9.57 x 10^11` over all 109736 size-3 subsets (`P = 60`), and likewise
nothing for `pi S(4)/varpi^2`, `A/varpi`, `B/varpi`, `A`, `B` at sizes 1-2 with
bounds `>= 4.3 x 10^17`. So the obvious lemniscatic normalisations do not help.

### 7. New exact structure for `pi S(4)` (VERIFIED-EXACT)

**(a) Split.** `pi S(4) = 2(A + B)` with

```text
A = int_0^1 arcsin(v) /sqrt(1-v^4) dv = int_0^{pi/2} a  da / sqrt(1+sin^2 a),
B = int_0^1 arcsinh(v)/sqrt(1-v^4) dv = int_0^{pi/2} arcsinh(sin a) da / sqrt(1+sin^2 a),
```

verified to 460 digits. The two right-hand integrands are **smooth on the closed
interval** (no endpoint singularity at all), so `A` and `B` are individually the
cheapest high-precision handle on `S(4)` this repo has.

**(b) By-parts normal form.** Let `V(v) = int_0^v dt/sqrt(1-t^4) =
v * 2F1(1/4,1/2;5/4;v^4)` be the lemniscatic arcsine, `V(1) = varpi/2`.
Integrating `A` and `B` by parts against `dV` gives, exactly,

```text
pi S(4) = varpi * ( arcsin 1 + arcsinh 1 ) - 2W
        = varpi * ( pi/2 + log(1+sqrt2) )  - 2W,
W = int_0^1 V(v) [ (1-v^2)^{-1/2} + (1+v^2)^{-1/2} ] dv
  = 1.535130258897420337752379135458933506929099268279853048340685809391999929598...
```

**verified to 101 digits** (split `W` at the two factors and substitute
`v = sin th` in the `(1-v^2)^{-1/2}` piece; both halves are then smooth, and
`2(A+B) - [varpi(pi/2+log(1+sqrt2)) - 2W] = 2.9e-101`). This isolates the *elementary
lemniscatic* part `varpi(arcsin 1 + arcsinh 1)` and shows the entire obstruction
sits in the single second-kind period `W`. It is the sharpest reformulation of
`S(4)` currently in the repo, and the natural next target: `W`, not `S(4)`.

### 8. Scope

Sections 4-6 upgrade the section-4/5 assessment of this file from a battery run
over a *rational*-coefficient pool to one run over a **quadratic-surd-coefficient
pool**, with the control that Addendum 1 showed to be indispensable. The
conclusion is unchanged in direction and much stronger in content: `k = 1,2,3`
are elementary, and `k = 4` is not reachable in that class within the recorded
bounds. The monodromy statement (7) remains the *explanation*; this addendum is
the *evidence*, now with the basis defect repaired. Nothing here bears on wider
constant classes -- in particular `W` of section 7(b) is unclassified.
