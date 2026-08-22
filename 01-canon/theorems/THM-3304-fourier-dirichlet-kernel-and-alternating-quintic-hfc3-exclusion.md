---
id: THM-3304
title: "Fourier-Dirichlet kernel and alternating-quintic HFC(3) exclusion"
status: >
  PROVED + VERIFIED-EXACT (kernel, phase obstruction, and sign-quintic
  exclusion) + FINITE-EXACT MOD-103 SIDECAR; DEGREE-FOUR CELL OPEN.  In
  Fourier barycentric coordinates a,b on the triangle, the complete
  factorial-numerator moment kernel is 1/(1-u^3-v^3-3uv).  Thus <a^r b^s>
  vanishes exactly off the C3-balanced congruence and is strictly positive on
  it.  This gives an explicit phase obstruction for the five-dimensional
  degree-four C3 eigenspace: a third-moment null vector must have shortest
  closed coefficient-phase covering width at least pi/3.  Its complete moments M3,...,M15 were built
  modulo 103 but not projectively eliminated.  In the orthogonal S3-sign
  lane, every homogeneous alternating quintic is Vandermonde*(A e1^2+B e2),
  and its second/fourth moment numerators have resultant 846709600, with
  nonzero infinity chart.  Hence no nonzero sign quintic is an HFC(3)
  witness.  This does not prove HFC(3) or FC(3).
source: root/cross-frontier-35-81/2026-08-03
audit: >
  The Fourier companion proves the cyclotomic denominator identity, checks
  the positive coefficient formula on a complete 25-by-25 control box, builds
  all five degree-four moment forms through M15 modulo 103, and independently
  rehomogenizes M3.  The sign companion evaluates I2 and I4 by both sparse
  factorial-Dirichlet expansion and direct iterated triangle integration,
  then checks the affine resultant and projective infinity point.  Normal,
  optimized, and stored outputs are byte-identical for both companions.
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
related:
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3301-symmetry-vanishing-is-mathieu-compatible
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
  - THM-3321-hesse-moment-kernel-and-cyclic-quartic-support-four-exclusion
support_script: 04-computation/factorial_hfc3_symmetry_cells_support_thm3304.py
support_script_sha256: 1c7058df1bbd72b6f2eec7ebb9c06e49a09218df41a96bd4edf6d85e11d5394a
script: 04-computation/factorial_hfc3_fourier_kernel_degree4_probe_thm3304.py
output: 05-knowledge/results/factorial_hfc3_fourier_kernel_degree4_probe_thm3304.out
script_sha256: 6bb43cf838c50124677ddbb88bd034d23608d8a6e9aae1c5b4e079623a1fa817
output_sha256: 74c3fd06b528bb5e858b3db413788f058ea5d272484fe73e75354c7ae26c3ca5
independent_script: 04-computation/factorial_hfc3_sign_quintic_thm3304.py
independent_output: 05-knowledge/results/factorial_hfc3_sign_quintic_thm3304.out
independent_script_sha256: 8b94cc7792ce79df1b78813b60acb8a2ffc6e4994b0fb646319b59932de2f239
independent_output_sha256: 21935d028698ff34400e5652a2d81ee3cb17b3e7baab240a6cc4635cb2876c8f
hash_basis: LF-normalized bytes
---

# THM-3304 -- Fourier-Dirichlet kernel and alternating-quintic HFC(3) exclusion

**PROVED + VERIFIED-EXACT (kernel, phase obstruction, and sign-quintic
exclusion) + FINITE-EXACT MOD-103 SIDECAR; DEGREE-FOUR CELL OPEN.**

THM-3300 excludes the cyclic HFC(3) eigenspaces through degree three, while
THM-3301 explains why finite symmetry can only remove congruence classes of
moments.  This theorem opens the next two underreported symmetry cells.  One
produces a positive closed kernel and an explicit phase barrier but remains open;
the other closes completely after two surviving moments.

## 1. Fourier coordinates on the triangle

Use barycentric coordinates `x+y+z=1`, `x,y,z>=0`, and normalized area

```text
<h> = 2 int_(Delta_2) h dA,       so <1>=1.              (1)
```

Fix `omega^2+omega+1=0` and put

```text
a=x+omega*y+omega^2*z,
b=x+omega^2*y+omega*z.                                  (2)
```

The three-cycle `(x,y,z)->(y,z,x)` sends `a` to `omega^2 a` and `b` to
`omega b`.

For nonnegative `r,s`, expand

```text
a^r b^s = sum_(alpha) c_alpha x^alpha,
N_(r,s) = sum_(alpha) c_alpha alpha!.                   (3)
```

The simplex beta integral gives

```text
<a^r b^s> = 2 N_(r,s)/(r+s+2)!.                         (4)
```

## 2. The exact Fourier--Dirichlet kernel

**Theorem 1.**  The exponential generating function of the factorial
numerators `(3)` is

```text
sum_(r,s>=0) N_(r,s) u^r v^s/(r!s!)
 = 1/[(1-u-v)(1-omega*u-omega^2*v)(1-omega^2*u-omega*v)]
 = 1/(1-u^3-v^3-3uv).                                  (5)
```

*Proof.*  In the first expression, expand each geometric factor.  Choosing a
term from the `x`, `y`, and `z` factors gives exactly the coefficient of the
corresponding monomial in `a^r b^s`, multiplied by its multi-factorial.  The
second equality is the cubic identity

```text
(X+Y+Z)(X+omega Y+omega^2 Z)(X+omega^2 Y+omega Z)
 = X^3+Y^3+Z^3-3XYZ                                    (6)
```

at `(X,Y,Z)=(1,-u,-v)`.  QED

Expanding the last member of `(5)` gives the explicit positive formula

```text
[u^r v^s](5)
 = sum multinomial(i,j,ell) 3^ell,                     (7)
```

where the sum is over `i,j,ell>=0` with

```text
r=3i+ell,       s=3j+ell.                               (8)
```

There is a solution of `(8)` exactly when `r-s=0 mod 3`.  Therefore

```text
<a^r b^s> = 0    iff r-s != 0 mod 3,
<a^r b^s> > 0    if  r-s  = 0 mod 3.                    (9)
```

The strict positivity is stronger than the usual character-selection rule:
after the C3 congruence is satisfied, no further Fourier cancellation remains
inside a single monomial.

## 3. The degree-four C3 cell: a phase barrier, not a closure

On `x+y+z=1`, the degree-four `omega` eigenspace has basis

```text
b,       a^2,       a b^2,       b^4,       a^3 b.      (10)
```

These are the restrictions of the homogeneous degree-four basis obtained by
inserting the needed powers of `e1=x+y+z`.  Write

```text
P=c0*b+c1*a^2+c2*a*b^2+c3*b^4+c4*a^3*b.                (11)
```

Every parameter monomial in `<P^3>` has positive real coefficient.  Indeed,
the multinomial coefficient is positive and the total Fourier weight is
balanced, so `(9)` makes its triangle moment strictly positive.

Consequently, if all nonzero `c_i` lie in one open angular sector of width
less than `pi/3`, rotate `(11)` so their arguments lie in
`(-pi/6,pi/6)`.  Every cubic parameter monomial then has positive real part,
and

```text
<P^3> != 0.                                              (12)
```

Thus a degree-four cyclic null candidate must have shortest closed coefficient-
phase covering width at least `pi/3`.  This is a rigorous necessary geometry,
not a proof that no such candidate exists.  THM-3310 later sharpens the
inequality to strict by auditing the closed-arc endpoint case.

The exact companion constructs the complete forms `<P^m>` for
`m=3,6,9,12,15` over `F_103`; their term counts are

```text
35, 210, 715, 1820, 3874,                               (13)
```

and their serialized digest is

```text
42fb0c17ba4a3c7b3726181ce9941eee4fa11d01b0b2a1fead6823e1ef621658.
```

This is a `FINITE-EXACT MOD-103` sidecar only.  No projective elimination
certificate is claimed in characteristic zero or characteristic `103`.  The
degree-four C3 cell remains **OPEN**.

## 4. The alternating quintic cell closes at moments two and four

Let

```text
V=(x-y)(y-z)(z-x),
e1=x+y+z,       e2=xy+yz+zx.                            (14)
```

Every homogeneous alternating quintic is uniquely

```text
G=V(A e1^2+B e2).                                      (15)
```

This is because every alternating polynomial is divisible by `V`, and the
remaining symmetric homogeneous quadratic is a linear combination of
`e1^2,e2`.  A transposition reverses `G` and preserves triangle area, so all
odd moments vanish automatically.

Two exact evaluation routes give

```text
int_Delta G^2 dA
 = (165A^2+66AB+7B^2)/277200,                           (16)

int_Delta G^4 dA
 = (71060A^4+53295A^3B+15675A^2B^2
    +2123AB^3+111B^4)/29875045200.                      (17)
```

In the affine chart `A!=0`, set `t=B/A`.  The primitive numerators are

```text
q2=7t^2+66t+165,
q4=111t^4+2123t^3+15675t^2+53295t+71060,               (18)
```

and exact elimination gives

```text
Res_t(q2,q4)=846709600 != 0.                             (19)
```

The projective infinity point `A=0` has respective coefficients `7` and
`111`, both nonzero.  Hence `(16)` and `(17)` have no common nonzero
projective zero.  No nonzero alternating quintic can have all positive
moments zero.

## 5. Consequences and scope

Proved: the all-degree Fourier kernel `(5)--(9)`, the degree-four phase
obstruction `(12)`, and complete exclusion of the alternating quintic cell by
`(16)--(19)`.

Finite-exact only: the five modular degree-four forms and digest `(13)`.

Open: projective elimination of the degree-four C3 eigenspace, every other
degree-four/five representation cell, arbitrary HFC(3), and full FC(3).
MISTAKE-350 is load-bearing: THM-3018 identifies this compact triangle problem
with the **homogeneous** subclass HFC(3), not with all three-variable
polynomials.  Finite `S3` symmetry explains the automatic congruence gaps but,
in agreement with THM-3301, does not itself settle the surviving moments.

## 6. Reproduction

Run

```text
python 04-computation/factorial_hfc3_fourier_kernel_degree4_probe_thm3304.py
python -O 04-computation/factorial_hfc3_fourier_kernel_degree4_probe_thm3304.py
python 04-computation/factorial_hfc3_sign_quintic_thm3304.py
python -O 04-computation/factorial_hfc3_sign_quintic_thm3304.py
```

and compare LF-normalized bytes with the two declared outputs.  Both scripts
import only the declared support module and standard/SymPy libraries.  They
use exact integer, rational, algebraic, and finite-field arithmetic; neither
contains a Python assertion node or a floating literal.

**QED.**
