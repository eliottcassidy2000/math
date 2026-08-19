---
id: THM-3550
title: "Prime-degree exclusion and pencil-height eight floor for planar Keller counterexamples"
status: >
  PROVED / HOSTILE AUDIT IN PROGRESS.  No nonzero target-pencil member of a
  nonautomorphic complex planar Keller pair has prime total degree.  The
  first positive-weight Newton face leaving the forced pure-power top vertex
  would have a second vertex whose exponent is simultaneously divisible by
  the prime and of smaller total degree.  Together with THM-3544 and the
  homogeneous common-form lemma, this forces output-pencil height at least
  eight.  The height is total polynomial degree, not generic cover degree.
source: boxeph-2026-08-18-jacobian-dephasing
depends_on:
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-2740-polynomial-coordinate-first-target-triangularity-and-mixed-resolvent-shear-closure
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-3025-W-is-forced-to-vanish-for-every-jc2-counterexample-candidate
---

# THM-3550 -- prime-degree exclusion and pencil-height eight floor

**PROVED / HOSTILE AUDIT IN PROGRESS.**  Let

```text
F=(P,Q):C^2 -> C^2,             Jac(P,Q)=kappa in C*,  (1)
```

and suppose `F` is not a polynomial automorphism.  Then:

1. every nonzero `R=sP+tQ` has composite total degree; and
2. the output-pencil height

   ```text
   h(F)=max_((s,t)!=(0,0)) deg(sP+tQ)
       =max(deg P,deg Q)                                  (2)
   ```

   satisfies `h(F)>=8`.

Together with THM-3544, every pencil degree is therefore a composite integer
at least six.  In particular, if one pencil direction has degree six, every
independent complement has degree at least eight.

## 1. Prime-degree exclusion

Fix a nonzero pencil member `R` and complete it to a target `GL_2` change
`(R,S)`.  Then `Jac(R,S)` is a nonzero constant and the pair remains
nonautomorphic.  THM-2740 says `R` is not a coordinate.  Suppose for a
contradiction that

```text
deg R=p                                                       (3)
```

is prime.

At ordinary weight `(1,1)`, THM-2102 says the top homogeneous form is a
proper power.  Since its degree is prime,

```text
R_p=lambda L^p.                                               (4)
```

Make a linear source change taking `L=x`.  Because `R` belongs to a
nonautomorphic Keller pair, the source-fibre gate in THM-3544 (already
THM-2063 is enough here) says that `R` contains a monomial `x^i y^j` with
`j>0`.  Every such monomial lies below the pure degree-`p` top term, so

```text
i+j<=p-1.                                                     (5)
```

Over all support exponents with `j>0`, define

```text
rho=min (p-i)/j >0.                                           (6)
```

Choose positive integral weights proportional to `(1,rho)`.  Among the
exponents attaining the minimum in `(6)`, choose one with maximal `j` and
call it `(i_*,j_*)`.  The resulting top face is a segment whose two opposite
vertices include `(p,0)` and `(i_*,j_*)`; in particular `j_*>0`.
This endpoint choice is essential: an interior monomial of a power need not
have exponent divisible by the power.  Equation `(6)` is precisely the first
Newton slope at which a lower term ties `x^p`.

THM-2102 again says this weighted face is a proper power, say

```text
in_w(R)=mu H^m,                    m>=2.                       (7)
```

Newton polytopes satisfy

```text
Newt(H^m)=m Newt(H).                                         (8)
```

The vertex `(p,0)` in `(8)` forces `m|p`; hence primality and `m>=2` give
`m=p`.  Every other vertex of `(7)` consequently has both coordinates
divisible by `p`.  This is impossible for the opposite endpoint
`(i_*,j_*)`: it has `j_*>0` but `i_*+j_*<=p-1` by `(5)`.  Thus the first
departing face is not a proper power, contradicting THM-2102.  No prime
degree can occur.

The Newton-polytope step is insensitive to coefficient cancellation: an
extreme term of `H` contributes its nonzero `m`th power at the corresponding
extreme term of `H^m`.

## 2. Pencil height

Assume `h(F)<=7`.  THM-3544 gives degree at least six for every nonzero
pencil member, while Section 1 excludes degree seven.  Hence

```text
deg P=deg Q=6.                                               (9)
```

The top homogeneous part of `Jac(P,Q)=kappa` gives

```text
Jac(P_6,Q_6)=0.                                              (10)
```

The standard homogeneous common-form lemma in characteristic zero says
that two nonzero homogeneous binary forms with zero Jacobian are powers of
one common homogeneous form:

```text
P_6=a H^m,                    Q_6=b H^n.                     (11)
```

For completeness, `(10)` makes the ratio of suitable powers of `P_6` and
`Q_6` constant in the one-variable function field of `P^1`; unique
factorization of binary forms then gives `(11)`.  Equal degrees in `(9)`
force `m=n`, so `P_6` and `Q_6` are constant-proportional.  A nonzero target
combination cancels their degree-six terms and has degree at most five.
The combination is nonzero because a Keller pair is algebraically, hence
linearly, independent.  This contradicts THM-3544.  Therefore

```text
h(F)>=8.                                                     (12)
```

If `R` has degree six and `S` is an independent pencil member, `S` cannot
also have degree six by the same argument, cannot have degree seven by
Section 1, and cannot have smaller degree by THM-3544.  Thus `deg S>=8`.

## 3. Reduced pencil spectrum and typed leading cells

Every Keller pencil admits a target basis `(R,S)` with distinct degrees

```text
n=deg R < m=deg S.                                      (13)
```

Indeed, if a basis has equal degree `D`, its top forms have zero Jacobian;
the common-form lemma and equal degree make them proportional.  One target
combination cancels the top form and produces `(13)`.  In this reduced basis,
the projective pencil has exactly two degree values: the unique direction
`[R]` has degree `n`, while every other direction has degree `m`.

For a hypothetical counterexample, the preceding results say

```text
6<=n<m,             n and m composite,             m>=8. (14)
```

The top equation `Jac(R_n,S_m)=0` and the common-form lemma give, after
choosing the common base power-free,

```text
R_n=c H^a,          S_m=d H^b,
n=a deg H,          m=b deg H,          a,b>=2.         (15)
```

The inequalities `a,b>=2` are also forced separately by THM-2102.  Thus a
degree search can be organized by typed cells

```text
(n,m; deg H,a,b)
```

obeying `(14)`--`(15)`, rather than by arbitrary pairs of supports.

## 4. Exact anatomy at the minimal height

Suppose the lower bound is sharp, `h(F)=8`.  Choose a pencil member `S` of
degree eight and any independent member `P`.  Its degree is either six or
eight: the floor excludes smaller degrees and Section 1 excludes seven.  If
both have degree eight, the leading equation

```text
Jac(P_8,S_8)=0
```

and the common-form lemma make `P_8,S_8` proportional.  Subtracting the
leading term gives an independent pencil member of degree at most seven,
hence exactly six.  Thus every minimal-height pencil admits a target basis

```text
(deg R,deg S)=(6,8).                                      (16)
```

The top Jacobian equation now gives `Jac(R_6,S_8)=0`.  Choose the common
homogeneous base `H` power-free.  Then

```text
R_6=a H^r,                 S_8=b H^s,
r deg H=6,                 s deg H=8.                    (17)
```

Both exponents are at least two by THM-2102, and `deg H` divides
`gcd(6,8)=2`.  Exactly two leading architectures remain, up to source
`GL_2`:

```text
deg H=1:       (R_6,S_8)=(a L^6,b L^8),
deg H=2:       (R_6,S_8)=(a(xy)^3,b(xy)^4).             (18)
```

In the quadratic case, a power-free binary quadratic has two distinct roots
over `C` and is therefore equivalent to `xy`.  The source-fibre floor also
specifies the first repairs.  In the linear case, the `L=0` direction gets
all of its required fibre degree from rows below the top row.  In the
quadratic case, `R_6=(xy)^3` has fibre degree three along either root
direction, so a lower row must raise both to at least four.

A concrete support-positive pair of face skeletons for the first architecture
is

```text
R_0=(x^3+alpha y^2)^2,
S_0=(x^4+beta y^2)^2.                                  (19)
```

Their ordinary top forms are `x^6,x^8`, their rescue faces at weights
`(2,3)` and `(1,2)` are squares, both have `y`-fibre degree four, and every
positive-weight face of each displayed polynomial is a square.  They are
only face skeletons: their Jacobian contains the two displayed base factors,
so they are not a Keller pair.  The construction problem is to add lower
terms inside the Newton cages that remove the common gradient factors while
retaining all proper-power faces.  By contrast, the hostile sextic
`[xy(x-y)]^2` has a cubic base and therefore cannot be the degree-six member
of a height-eight common-leading-base pair; that skeleton first becomes
compatible with a complementary height divisible by three.

## 5. Failure boundary and non-consequences

The nonautomorphic hypothesis is essential.  For every prime `p`,

```text
(P,Q)=(x^p+y,x)
```

is a tame Keller automorphism with a prime-degree pencil member.  Its first
departing face is power-free, exactly as THM-2102 permits for a coordinate.

The prime argument genuinely stops at composite degrees.  The hostile
sextic

```text
[xy(x-y)]^2
```

has proper-power positive-weight faces and the directional fibre degrees
needed by the two gates used here, although its gradient has a common factor
and it is not a Keller component.  The lower bound `(12)` does not assert
that degree eight is realizable, bound generic cover degree, or prove
properness.  It only removes every prime degree and the entire pencil-height
six/seven box from a hypothetical counterexample.
