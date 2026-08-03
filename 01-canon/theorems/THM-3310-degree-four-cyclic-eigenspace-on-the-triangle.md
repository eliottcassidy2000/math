---
id: THM-3310
title: "Degree-four cyclic eigenspace on the triangle"
status: >
  PROVED + VERIFIED-EXACT (structure, basis, moment table) + FINITE-EXACT
  (support-<=3 exclusion); INDEPENDENT IMMUTABLE AUDIT PENDING.
  In the complex coordinate `z = s_1 + omega s_2 + omega^2 s_3` -- which sends
  `Delta_2` to the equilateral triangle on `1, omega, omega^2` and is exactly
  THM-3300's degree-one eigenvector -- the cyclic rotation acts by
  `z -> omega^2 z`, so `z^a zbar^b` is an `omega`-eigenvector precisely when
  `a - b = 2 (mod 3)`.  The eigenspace of degree `d` therefore has the monomial
  basis `{z^a zbar^b : a+b <= d, a-b = 2 mod 3}`, with dimensions
  `1,2,3,5,7,9,12` for `d = 1..7`; the degree-four basis is
  `{zbar, z^2, z zbar^2, z^3 zbar, zbar^4}` and the lower-degree eigenspaces
  are literally coordinate sub-loci of it.  In this basis the simplex average
  is a table lookup and every basis element is a SINGLE monomial.  All ten
  coordinate lines are excluded in exact `Q(omega)` and all ten coordinate
  planes by guarded modular resultants, so no `g` supported on at most three
  of the five monomials can satisfy `<g^m>=0` for all `m`.  Support four and
  five remain OPEN.
audit: >
  The exact companion verifies the rotation eigenvalues, `<z>=0`, `<|z|^2>=1/4`,
  and that the monomial model reproduces the projection-computed eigenspace
  dimension for every degree `1..7`.  It builds the mixed-moment table
  `<z^a zbar^b>` through total degree 36 in exact `Q(omega)`, checks the
  vanishing rule `a=b mod 3`, the named values `1/4, 1/10, 1/10, 1/10, 29/560`,
  and re-derives THM-3300's closed form `<z^m>=2/((m+1)(m+2))[3|m]` through
  `m=12`.  It excludes the ten coordinate lines by an exact `Q(omega)` gcd of
  `<g^3>` and `<g^6>` including the point at infinity, and the ten coordinate
  planes by interpolated modular resultants of degrees 18 and 27 with a
  degree-constancy check that makes the reduction sound.  Both gcd routines
  carry a positive control.  A MODULAR GUARD is enforced and demonstrated to be
  load-bearing.  Normal and `-O` replay are byte-identical.  Independent
  immutable audit is pending.
source: death-star-delta2-degree4-2026-08-03
depends_on: []
related:
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
script: 04-computation/degree_four_cyclic_eigenspace_thm3310.py
output: 05-knowledge/results/degree_four_cyclic_eigenspace_thm3310.out
script_sha256: 1a46fb4311a44d7807eccaf9339b291b8d930710d9f7c1fc8ebd03ac0c7acba6
output_sha256: 5d7292f2cc311c727e16d379c8134371ca3f3d1ac610ab2a1748fff8c8bb757c
hash_basis: LF-normalized bytes
---

# THM-3310 -- degree-four cyclic eigenspace on the triangle

**PROVED + VERIFIED-EXACT + FINITE-EXACT; INDEPENDENT IMMUTABLE AUDIT
PENDING.**

THM-3300 left the degree-four cyclic eigenspace on `Delta_2` open: it has
dimension five, and the resultant method used there handles two projective
parameters, not four.  The right first move is not a bigger elimination but a
change of coordinate.

## 1. The complex coordinate

Put

```text
z = s_1 + omega s_2 + omega^2 s_3,     zbar = s_1 + omega^2 s_2 + omega s_3.
```

This is exactly THM-3300's degree-one eigenvector.  It identifies `Delta_2`
with the equilateral triangle on the vertices `1, omega, omega^2`, sends the
barycentre to `0`, and satisfies `<z> = 0`, `<z zbar> = 1/4`.  On `Delta_2` the
functions `z` and `zbar` are algebraically independent, so every polynomial is
a polynomial in `(z, zbar)`.

The cyclic rotation `rho : s_1 -> s_2 -> s_3 -> s_1` acts by

```text
rho(z) = omega^2 z,      rho(zbar) = omega zbar,
```

so `rho(z^a zbar^b) = omega^(2a+b) z^a zbar^b`.  Hence:

**Theorem 1.**  The `omega`-eigenspace of functions of degree at most `d` on
`Delta_2` has the monomial basis

```text
{ z^a zbar^b : a + b <= d,  a - b = 2 (mod 3) }.                     (1)
```

Its dimensions for `d = 1,...,7` are `1, 2, 3, 5, 7, 9, 12`, each verified
against an independent group-projection computation.

For `d = 4` this gives the five monomials

```text
zbar,   z^2,   z zbar^2,   z^3 zbar,   zbar^4.                       (2)
```

Two structural facts fall out immediately and are worth stating.

- **The tower is a flag of coordinate sub-loci.**  `(1)` is nested in `d`, so
  the degree-1 eigenspace is the first coordinate of `(2)`, the degree-2
  eigenspace is the coordinate line `{zbar, z^2}`, and the degree-3 eigenspace
  is the coordinate plane `{zbar, z^2, z zbar^2}`.  The exclusions in section 3
  therefore re-derive THM-3300's degree-1, 2 and 3 results independently.
- **Each basis element is a single monomial**, so the simplex average becomes a
  table lookup and the whole computation collapses.  This is what makes the
  degree-four case tractable at all.

## 2. The mixed-moment table

Define `mu(a,b) = <z^a zbar^b>` against the uniform probability measure.  By
the `rho`-symmetry `mu(a,b) = 0` unless `a = b (mod 3)`.  The companion
computes the whole table exactly through total degree 36; the low entries are

```text
mu(0,0)=1,      mu(1,1)=1/4,    mu(3,0)=mu(0,3)=1/10,
mu(2,2)=1/10,   mu(1,4)=mu(4,1)=2/35,
mu(6,0)=mu(0,6)=1/28,           mu(3,3)=29/560.                      (3)
```

Along the pure-power axis this reproduces THM-3300's closed form
`mu(m,0) = 2/((m+1)(m+2))` when `3 | m` and `0` otherwise, checked to `m=12`.

## 3. Everything of support at most three is excluded

Write `g = alpha zbar + beta z^2 + gamma z zbar^2 + delta z^3 zbar + eps zbar^4`
and recall that a cyclic eigenvector needs `<g^m> = 0` only for `3 | m`
(THM-3300 Theorem 3), so the first obligations are `<g^3> = <g^6> = <g^9> = 0`.

**FINITE-EXACT.**

```text
all 10 coordinate LINES:  gcd_a( <g^3>, <g^6> ) is constant, in exact
                          Q(omega), and the remaining projective point of
                          each line is excluded separately;             (4)

all 10 coordinate PLANES: gcd_a( Res_b(<g^3>,<g^6>), Res_b(<g^3>,<g^9>) )
                          is constant, by interpolated modular resultants
                          of degrees 18 and 27.                         (5)
```

Since a plane's complement inside itself is a coordinate line, and a line's
complement is a coordinate point, `(4)`--`(5)` cover every `g` whose support is
at most three of the five monomials in `(2)`.  **No such `g` satisfies even the
first two surviving conditions.**

Soundness of `(5)`: a nonconstant gcd in characteristic zero reduces to a
nonconstant gcd mod `p` provided the degrees in the eliminated variable are
preserved, and the companion checks that those degrees are constant across all
interpolation nodes.

## 4. What remains open

Support four and support five.  The obstacle is not conceptual but the
elimination: with three or four projective parameters the resultant chain
requires bivariate interpolation at degrees in the thousands.  An exhaustive
rational scan of `P^4(F_q)` at the smallest valid prime `q=61` was implemented
and then dropped from the companion: it costs about `1.4 x 10^7` points and,
being a search for `F_q`-rational points only, would not be decisive anyway.

## 5. The modular guard, and why it is load-bearing

The uniform-simplex moment of a monomial of total degree `D` is
`2 a! b! c! / (D+2)!`.  **Reducing mod `p` is sound only when `p > D + 2`.**
Since `deg g <= 4`, testing `<g^m>` needs

```text
p > 4m + 2:   m=3 needs p>=19,  m=6 needs p>=31,  m=9 needs p>=43,
              m=12 needs p>=61, m=15 needs p>=67, m=18 needs p>=79.   (6)
```

Below that threshold the denominator is `0 mod p`, a naive modular inverse
returns `0`, and every moment is silently corrupted.  This is not hypothetical:
at `p=7`, sixteen distinct moment denominators through degree 12 vanish, and an
unguarded cascade reports `743` spurious solutions at the first condition
alone.  The companion enforces `p > MAX_TOTAL + 2` with an explicit `require`
and refuses to invert zero.

**Audit of prior work.**  THM-3300's degree-three modular certificate used
`g` of barycentric degree 3 with `m <= 9`, hence denominators at most `29!`,
against `p = 10^9 + 9`.  The guard is satisfied with enormous margin, so
THM-3300 is unaffected.

## 6. Scope

Proved: Theorem 1, the basis `(2)`, the nesting flag, and the moment table
`(3)`.  Finite-exact: `(4)` and `(5)`, hence the support-`<=3` exclusion.

**Not proved:** `FC(3)`, `HFC(3)`, or any statement about supports four and
five.  Nothing here bears on `FC(n)` for `n != 3`, on non-eigenvector
candidates, or on THM-3018's outstanding Laplace closure.  The exclusions
concern cyclic eigenvectors of degree at most four only.

Run

```text
python 04-computation/degree_four_cyclic_eigenspace_thm3310.py
python -O 04-computation/degree_four_cyclic_eigenspace_thm3310.py
```

and compare LF-normalized bytes with the declared output.  Exact cyclotomic and
guarded modular arithmetic only; no floating point, random sampling, imported
executable, or assertion-sensitive test.

**QED** for sections 1--2; sections 3 and 5 are finite-exact in their stated
universes.
