---
id: THM-3310
title: "Degree-four cyclic eigenspace on the triangle"
status: >
  PROVED + VERIFIED-EXACT (structure, basis, moment table) + FINITE-EXACT
  (support-<=3 exclusion) + INDEPENDENTLY VERIFIED-EXACT.
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
  of the five monomials can satisfy `<g^m>=0` for all `m`.  This theorem leaves
  supports four and five untreated; THM-3321 subsequently excludes support
  four.  Support five remains OPEN.
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
  load-bearing.  Normal and `-O` replay are byte-identical.  An independent
  characteristic-zero companion rebuilds the Fourier kernel over Q, excludes
  all ten lines by nonzero exact resultants and all ten planes by unit Groebner
  ideals, with no modular reduction.  A separate cubic-moment audit sharpens
  the shortest closed phase covering width to strictly greater than pi/3 and
  excludes five explicit
  lopsided axis boxes.
source: death-star-delta2-degree4-2026-08-03
depends_on: []
related:
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3304-fourier-dirichlet-kernel-and-alternating-quintic-hfc3-exclusion
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3321-hesse-moment-kernel-and-cyclic-quartic-support-four-exclusion
  - THM-3323-cyclic-quartic-support-five-exact-degree-21-rank
  - THM-3327-positive-moment-tensor-covering-arc-and-factorial-phase-dispersion
script: 04-computation/degree_four_cyclic_eigenspace_thm3310.py
output: 05-knowledge/results/degree_four_cyclic_eigenspace_thm3310.out
script_sha256: 1a46fb4311a44d7807eccaf9339b291b8d930710d9f7c1fc8ebd03ac0c7acba6
output_sha256: 5d7292f2cc311c727e16d379c8134371ca3f3d1ac610ab2a1748fff8c8bb757c
independent_script: 04-computation/factorial_hfc3_cyclic_quartic_sparse_strata_thm3310_independent.py
independent_output: 05-knowledge/results/factorial_hfc3_cyclic_quartic_sparse_strata_thm3310_independent.out
independent_script_sha256: d5e04a22a374ba430f948dc8ab3b840076bbe8d92c18ceabed464d90a4d0b69a
independent_output_sha256: a84049a4c4a313cb0bd61d0459d0bfba1f7906e40d578918f38831ead3b91b17
phase_script: 04-computation/factorial_hfc3_cyclic_quartic_lopsided_thm3310_audit.py
phase_output: 05-knowledge/results/factorial_hfc3_cyclic_quartic_lopsided_thm3310_audit.out
phase_script_sha256: 5aec456943b0e58ad7cfc1c02bb6993b7e4c0a3d74a33215016fe89b6fab78d0
phase_output_sha256: af05f1eb7517cb44e1d74d97b827b740bdaa39291f79da36eedea59beeacec36
hash_basis: LF-normalized bytes
---

# THM-3310 -- degree-four cyclic eigenspace on the triangle

**PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY VERIFIED-EXACT.**

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
at most three of the five monomials in `(2)`.  **No such `g` satisfies the first
three surviving conditions.**

Soundness of `(5)`: a nonconstant gcd in characteristic zero reduces to a
nonconstant gcd mod `p` provided the degrees in the eliminated variable are
preserved, and the companion checks that those degrees are constant across all
interpolation nodes.

### Independent characteristic-zero audit

The independent companion uses THM-3304's exact kernel

```text
<z^a zbar^b>
 = 2 a! b!/(a+b+2)! [u^a v^b](1-u^3-v^3-3uv)^(-1).       (6)
```

It does not import or execute the primary THM-3310 script and uses no modular
arithmetic.  On every coordinate line it computes a nonzero exact rational
resultant of `<g^3>,<g^6>`, including both projective endpoints.  On every
coordinate plane it dehomogenizes at the first selected coefficient and gets

```text
< <g^3>,<g^6>,<g^9> > = <1> in Q[y,z].                   (7)
```

The missing hyperplane in each affine chart is the coordinate line already
excluded.  All ten resultants and all ten unit-ideal tests are explicit runtime
truth gates.  Normal and optimized outputs are byte-identical.  This supplies
the previously pending independent audit and confirms that the support bound
is characteristic-zero, not an artifact of the guarded prime.

### Closed phase and lopsided-axis barriers

The cubic moment in the five coefficients has complete degree-three support:
all `35` coefficients are strictly positive.  Rotate any closed coefficient
phase arc of width `pi/3` to `[0,pi/3]`.  Every cubic term lies in the closed
upper half-plane.  Vanishing would force all terms to its boundary; pure cubes
put each occupied coefficient at an endpoint, while every ordered mixed term
`c_i^2c_j` then forces positive imaginary part if both endpoints occur.  If
only one endpoint occurs, all terms align and cannot cancel.  Hence every
cubic-null vector has shortest closed phase covering width **strictly greater**
than `pi/3`.

THM-3327 subsequently abstracts this argument: strict positivity of every
order-`m` tensor entry forces shortest closed coefficient covering width
greater than `pi/m`.  That statement remains basis-dependent.

A coefficientwise triangle inequality gives complementary magnitude barriers.
If coefficient `i` is largest and all others have modulus at most `q_i` times
its modulus, then the pure cube dominates for the following `q_i`.

In the basis order `(2)`, the radii are

```text
(q_0,q_1,q_2,q_3,q_4)=(1/10,1/16,1/18,1/22,1/23).        (8)
```

Thus at every projective cubic zero, the second-largest coefficient modulus is
strictly greater than `1/23` of the largest.  These phase/coamoeba exclusions
use only the necessary cubic moment; by themselves they close neither support
four nor support five.

## 4. What remains open

This theorem's resultant route leaves supports four and five untreated: with
three or four projective parameters the chain requires bivariate interpolation
at degrees in the thousands.  Subsequently, THM-3321 excludes support four by
five guarded homogeneous Macaulay certificates.  THM-3323 proves that the
support-five degree-21 map has exact rank `10980` and quotient dimension
`1670`; support five remains OPEN.  An exhaustive rational scan of
`P^4(F_q)` at the smallest valid prime `q=61` was dropped from this companion:
it costs about `1.4 x 10^7` points and would not be decisive anyway.

## 5. The modular guard, and why it is load-bearing

The uniform-simplex moment of a monomial of total degree `D` is
`2 a! b! c! / (D+2)!`.  **Reducing mod `p` is sound only when `p > D + 2`.**
Since `deg g <= 4`, testing `<g^m>` needs

```text
p > 4m + 2:   m=3 needs p>=19,  m=6 needs p>=31,  m=9 needs p>=43,
              m=12 needs p>=61, m=15 needs p>=67, m=18 needs p>=79.   (9)
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

**Not proved by this theorem:** `FC(3)`, `HFC(3)`, or any statement about
supports four and five.  Later theorems have the scoped consequences just
listed; support five is still open.  Nothing here bears on `FC(n)` for
`n != 3`, on non-eigenvector candidates, or on THM-3018's outstanding Laplace
closure.  The exclusions concern cyclic eigenvectors of degree at most four.

Run

```text
python 04-computation/degree_four_cyclic_eigenspace_thm3310.py
python -O 04-computation/degree_four_cyclic_eigenspace_thm3310.py
python 04-computation/factorial_hfc3_cyclic_quartic_sparse_strata_thm3310_independent.py
python -O 04-computation/factorial_hfc3_cyclic_quartic_sparse_strata_thm3310_independent.py
python 04-computation/factorial_hfc3_cyclic_quartic_lopsided_thm3310_audit.py
python -O 04-computation/factorial_hfc3_cyclic_quartic_lopsided_thm3310_audit.py
```

and compare LF-normalized bytes with the declared output.  Exact cyclotomic and
guarded modular arithmetic only; no floating point, random sampling, imported
executable, or assertion-sensitive test.

**QED** for sections 1--2; sections 3 and 5 are finite-exact in their stated
universes.
