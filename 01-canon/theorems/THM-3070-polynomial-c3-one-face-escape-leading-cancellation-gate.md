---
id: THM-3070
title: "Polynomial C3 one-face escape leading-cancellation gate"
status: >
  PROVED + VERIFIED-EXACT (independent audit pending).  Let a polynomial
  planar Keller pair have a tame ramification-index-three valuation over the
  target coordinate line t=0, with residue degree one.  Suppose that in a
  coefficient-field uniformizer s exactly one affine source coordinate has a
  pole, x=A(u)s^(-p)+..., while y=B(u)s^r+... has r>=0.  If the two leading
  monomials contribute nontrivially to dx wedge dy, then necessarily r=p+3.
  In that case polynomiality of both target coordinates forces incompatible
  weight-zero and weight-three initial forms.  The proof reduces them to a
  rational autonomous ODE C'=lambda C H(C); every nonconstant rational
  solution has H(z)=az and C=constant/(u-u0), after which polynomiality of
  the other target would require u=G(C) for a polynomial G, impossible.
  Consequently every surviving one-coordinate C3 branch in this coordinate-
  line scope has r<=p+2, its leading wedge cancels, A^r B^p is constant, and
  the associated graded source has a primitive toric binomial relation.
  This is necessary, not sufficient.  Branches where both coordinates
  escape, non-coordinate Jelonek components, and the full A4/S4 and JC(2)
  problems remain open.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
related:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
script: 04-computation/jc_polynomial_c3_one_face_cancellation_thm3070.py
output: 05-knowledge/results/jc_polynomial_c3_one_face_cancellation_thm3070.out
script_sha256: cb8270e2a47df6e68faec00eb2c467976ee0efd916cc385dfc384e58567d04fa
output_sha256: 575b259ea6fab01736d8e3e44449b288ef22fd0b4c8f9f25d68058f98d791a16
hash_basis: LF-normalized bytes
---

# THM-3070 -- polynomiality forces a toric cancellation before C3 escape

**PROVED + VERIFIED-EXACT (independent audit pending).**

## 1. Inheritance and result

[THM-3064](THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md)
shows that a tame `C3` branch of a hypothetical polynomial Keller map needs
the exact inverse-different cofactor pole `-2`.
[THM-3068](THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift.md)
then shows why a bound on the affine-primitive cofactor cannot exclude it:
reciprocal coordinate change carries a five-step valuation invoice, and a
punctured rational map realizes the complete local pole and exactness ledger.

The next missing operation is therefore not another cofactor valuation.  It
is simultaneous polynomiality of the **two target coordinates**.  This
theorem proves that polynomiality already excludes the simplest one-face
local anatomy.  Any survivor must move into a resonant associated-graded
lane where the leading source monomials satisfy a binomial relation and the
Jacobian first appears in a later key form.

The target-coordinate-line hypothesis below is load bearing.  A general
Jelonek component need not be straightened by a constant-Jacobian polynomial
target change, so the theorem is not a full `C3` exclusion.

## 2. Precise local setup

Let

```text
P,Q in C[x,y],                    u=P, t=Q,
det partial(P,Q)/partial(x,y)=kappa in C*.               (1)
```

Let `w` be a divisorial valuation of `C(x,y)` over the target prime `(t)`
such that

```text
w(t)=3,                  w(C(u)*)=0,
residue_field(w)=C(u).                                  (2)
```

Thus `(2)` is a ramification-index-three, residue-degree-one branch over the
coordinate line `{t=0}`.  In the degree-four fixed-plus-`C3` situation the
three escaping sheets contribute exactly such a valuation.  We do not need
the rest of the quartic or assume a monodromy group.

Choose a coefficient-field uniformizer `s` in the completed valuation field.
Write

```text
t=tau(u)s^3+O(s^4),                 tau in C(u)*,         (3)

x=A(u)s^(-p)+O(s^(-p+1)),           p>=1,
y=B(u)s^r+O(s^(r+1)),                r>=0,
A,B in C(u)*.                                             (4)
```

The big-`O` terms are Laurent series in `s` with coefficients in `C(u)`.
The residue-degree-one condition in `(2)` is what places `A,B` in `C(u)`
rather than in an unspecified residue extension.  Exactly one displayed
source coordinate escapes: `x` has a pole and `y` does not.  Interchanging
the affine source coordinates covers the opposite ordering.

Call `(4)` **leading-wedge nondegenerate** when

```text
L_(p,r):=r A'(u)B(u)+p A(u)B'(u) !=0.                  (5)
```

The primary theorem is

```text
there is no polynomial Keller pair (1)--(4)
whose leading wedge is nondegenerate.                    (6)
```

Equivalently, every branch in this precise coordinate-line, one-coordinate
escape scope must obey the necessary conditions

```text
0<=r<=p+2,
r A'B+p AB'=0,
A^r B^p in C*.                                           (7)
```

No sufficiency or globalization is asserted in `(7)`.

## 3. The differential invoice fixes the nondegenerate face

Differentiate `(4)` coefficientwise in the independent local coordinates
`(u,s)`.  The first possible term is

```text
dx wedge dy
 =[r A'B+p AB']s^(r-p-1) du wedge ds
   +O(s^(r-p))du wedge ds.                              (8)
```

On the other hand, `(1)` and `(3)` give

```text
dx wedge dy=kappa^(-1)du wedge dt
 =3 kappa^(-1)tau(u)s^2 du wedge ds+O(s^3).             (9)
```

If `(5)` holds, comparison of exponent and coefficient in `(8)--(9)` gives

```text
r=p+3,
(p+3)A'B+p AB'=3 kappa^(-1)tau.                        (10)
```

If `(5)` fails, every later Laurent contribution has exponent strictly
larger than `r-p-1`.  Equality with `(9)` is impossible when `r-p-1>=2`.
Hence the degenerate case necessarily has `r<=p+2`, and `(5)=0` integrates
to

```text
(A^r B^p)'=0,                   A^r B^p in C*.          (11)
```

This proves `(7)` once the nondegenerate case `(10)` is excluded.

## 4. The nondegenerate associated graded is genuinely monomial

Assume `(10)`.  Put

```text
g=gcd(p,3)=gcd(p,p+3),
p_0=p/g,                    r_0=(p+3)/g,
C(u)=A(u)^r_0 B(u)^p_0.                                 (12)
```

Equation `(10)` gives

```text
C'/C
 =[(p+3)A'/A+pB'/B]/g
 =3 kappa^(-1)tau/(gAB),                                (13)
```

so `C` is nonconstant.  In particular it is transcendental over `C`.

Give source monomials the valuation weight

```text
weight(x^i y^j)=-pi+(p+3)j.                             (14)
```

Any two exponent pairs of the same weight differ by an integral multiple of
`(r_0,p_0)`.  Therefore the leading coefficient of a nonzero fixed-weight
polynomial is a nonzero monomial in `A,B` times a Laurent polynomial in
`C`.  Since `C` is nonconstant, that Laurent polynomial cannot vanish
identically.  Thus no same-weight cancellation is possible:

```text
w(R(x,y))=minimum weight of a monomial of R
for every nonzero R in C[x,y].                           (15)
```

This is the exact meaning of the **monomial initial-form** claim.  It is
proved here from `(13)`, not assumed.  Subleading terms in `(4)` cannot
cancel the smallest weight in `(15)`.

The argument is uniform when `3|p`.  Then `g=3`; for example

```text
p=3:  (r_0,p_0)=(2,1),
p=6:  (r_0,p_0)=(3,2).                                  (16)
```

The leading powers of `s` need not themselves display the inertia character;
only `(2)--(4)` and the weight lattice are used.

## 5. Polynomial target initials force one autonomous ODE

Because `Q(x,y)=t` has valuation three, `(15)` says that its smallest
monomial weight is exactly three.  The nonnegative integer solutions of

```text
-pi+(p+3)j=3                                            (17)
```

are precisely

```text
(i,j)=(1+r_0 k,1+p_0 k),                 k>=0.          (18)
```

Consequently the weight-three initial coefficient of the polynomial `Q` is

```text
AB H(C),                         H in C[Z] nonzero.      (19)
```

Comparison with `(3)` gives

```text
AB H(C)=tau.                                             (20)
```

Combining `(13)` and `(20)` cancels the arbitrary tame leading unit `tau`:

```text
C'=lambda C H(C),              lambda=3 kappa^(-1)/g.   (21)
```

This cancellation is why no choice of local uniformizer is hidden in the
argument.

Likewise `P(x,y)=u` has valuation zero.  The nonnegative solutions of

```text
-pi+(p+3)j=0                                            (22)
```

are exactly `(i,j)=(r_0k,p_0k)`, `k>=0`.  Hence the
weight-zero initial coefficient of the polynomial `P` is

```text
G(C),                            G in C[Z],
u=G(C(u)) in C(u).                                        (23)
```

Equations `(21)` and `(23)` are the two polynomial target gates.  Either one
alone can be realized by a Laurent hostile; together they contradict.

## 6. Rational autonomous ODE classification

We prove the elementary lemma used in that contradiction.

> **Lemma.**  Let `C in C(u)` be nonconstant, `H in C[Z]` nonzero, and
> `lambda in C*`.  If `C'=lambda C H(C)`, then
>
> ```text
> H(Z)=aZ,                  a in C*,
> C(u)=-1/[lambda a(u-u_0)] for some u_0 in C.           (24)
> ```

First, `C` cannot be a nonconstant polynomial.  If `deg C=n>=1` and
`deg H=d>=0`, the two sides would have degrees

```text
n-1                         and                         n(d+1),           (25)
```

which cannot agree.  Thus `C` has a finite pole.  At a finite pole of order
`m`, the left side has pole order `m+1`, while the right side has pole order
`m(d+1)`.  Therefore

```text
m+1=m(d+1),                  md=1,                       (26)
```

so `m=d=1`.  In particular

```text
H(Z)=aZ+b,                    a!=0,                      (27)
```

and every finite pole of `C` is simple.

The rational function `C` has no finite zero.  At a zero of order `n>=1`,
the derivative has order `n-1`; the right side has order `n` if `b!=0` and
order `2n` if `b=0`.  Neither equals `n-1`.

In lowest terms, a rational function with no finite zero has constant
numerator, so write `C=c/D(u)` with `D` a nonconstant polynomial of degree
`N`.  At `u=infinity`,

```text
C'=Theta(u^(-N-1)).                                     (28)
```

If `b!=0`, the right side of `(21)` is `Theta(u^(-N))`, impossible.  Thus
`b=0`.  It is then `Theta(u^(-2N))`, so `(28)` forces

```text
N+1=2N,                         N=1.                    (29)
```

Solving the remaining equation `C'=lambda a C^2` gives `(24)`.  This audits
all polynomial, finite-pole, finite-zero, and infinity cases.

## 7. The second polynomial target is impossible

By `(24)`, `C(u)` tends to zero as `u` tends to infinity.  Hence every
polynomial `G(C(u))` tends to the finite value `G(0)`, whereas the left side
of `(23)` is `u` and has a pole at infinity.  This contradiction proves
`(6)`.

Equivalently, if `deg G=d`, clearing the denominator `(u-u_0)^d` in
`G(C(u))-u` leaves a numerator of degree `d+1` with nonzero leading
coefficient.  The companion checks this exact algebra for `0<=d<=10`; the
degree argument proves it for every `d`.

## 8. What survives: a toric initial relation and delayed key form

Return to the only possible lane `(7)` and put

```text
h=gcd(p,r),
X_0=A(u)s^(-p),                  Y_0=B(u)s^r.            (30)
```

Equation `(11)` implies the primitive relation

```text
X_0^(r/h) Y_0^(p/h)=c in C*.                            (31)
```

Thus the associated-graded map from the polynomial source is no longer
monomially injective: its kernel contains the primitive weight-zero binomial

```text
X^(r/h)Y^(p/h)-c.                                       (32)
```

This is exactly the information destroyed in Section 4.  Powers and
multiples of `(32)` can cancel an apparent lowest face, allowing a later key
form to carry the target residue `u` and the order-three transverse term.
The gap

```text
delta=p+3-r>=1                                          (33)
```

is the minimum delay that future branch expansions must pay before the
Jacobian coefficient at `s^2` can appear.

The surviving proof obligation is therefore sharply smaller:

```text
classify the binomial-resonant key-form expansions (31)--(33),
including the first subleading coefficient that carries u and tau.       (34)
```

Condition `(31)` is necessary, not sufficient.  It does not construct a
polynomial pair, a quartic cover, or a nonproper component.

## 9. Three sharp controls

### 9.1 Polynomial `P`, Laurent `Q`

THM-3068's punctured cubic component is

```text
P=x^4y/3,                       Q=x^(-3),
det partial(P,Q)/partial(x,y)=1.                         (35)
```

Its inverse has

```text
x=s^(-1),                       y=3u s^4,
p=1, r=4,                       L_(1,4)=3.               (36)
```

It realizes the excluded nondegenerate face because `Q` is Laurent, not
polynomial.  This proves polynomiality of `Q` is load bearing.

### 9.2 Laurent `P`, polynomial `Q`

On the two-dimensional torus take

```text
P=9/(x^5y^2),                   Q=-x^6y^3/27.            (37)
```

Again the Jacobian is one.  Over `t=s^3`, the exact inverse is

```text
x=u^(-1)s^(-2),                 y=-3u^2s^5,
p=2, r=5,                       L_(2,5)=3.               (38)
```

Now `Q` is polynomial but `P` is Laurent.  Thus polynomiality of `P` is
independently load bearing.

### 9.3 The leading-cancellation equations are locally nonempty

The formal Laurent packet

```text
x=u s^(-1),
y=u^(-1)s+(3/4)s^4,
t=s^3                                                    (39)
```

satisfies exactly

```text
dx wedge dy=du wedge dt.                                (40)
```

Here `p=r=1`, the leading coefficient cancels, and `AB=1`.  This is only a
local symplectic control; no polynomial globalization is claimed.  It shows
why `(7)` must be treated as a live key-form lane rather than declared empty
from the differential equation alone.

## 10. Exact boundary and companion

```text
PROVED HERE:       exclusion of the nondegenerate one-face branch under
                   simultaneous polynomial target coordinates;
                   exact weight-zero and weight-three initial forms;
                   complete rational autonomous ODE classification;
                   uniform gcd treatment, including 3|p;
                   necessary toric cancellation and delay delta>=1;
                   two independent Laurent polynomiality hostiles.

NOT PROVED:        sufficiency or realizability of the cancellation lane;
                   a bound on its key-form depth;
                   branches where both affine coordinates escape;
                   straightening an arbitrary Jelonek component to t=0;
                   exclusion of every C3 component, A4, S4, G1,
                   JC(2), or DC(2).                                     (41)
```

Run

```text
python3 04-computation/jc_polynomial_c3_one_face_cancellation_thm3070.py
python3 -O 04-computation/jc_polynomial_c3_one_face_cancellation_thm3070.py
```

Both executions must LF-byte-match the stored transcript.  The companion
checks the leading wedge, all weight-zero/three lattice points for
`1<=p<=24` (including every multiple of three in that range), every integer
balance used in the ODE pole/zero/infinity proof, the unique rational normal
form, the polynomial-composition degree contradiction through degree ten,
both Laurent hostiles, and the exact cancellation control.  Every
truth-bearing executable check uses explicit runtime exceptions rather than
Python assertions.
