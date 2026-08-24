---
id: THM-3916
title: "Positive-genus collapsed valuation obstruction to a planar Keller model"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.  If a
  divisorial valuation of k(x,y) has positive-genus residue field and is
  positive on both entries of a polynomial pair, that pair cannot have
  scalar-unit Jacobian.  Applied to the rational chart of THM-3915, the
  common-zero divisor F=0 is irreducible of genus two: its cubic projection
  has eight simple finite ramifications and three unramified infinity
  places.  Hence the THM-3915 normal cubic completion contains no affine-
  plane open on which its target functions form a Keller pair, and that
  candidate is not a planar Jacobian counterexample.  Repeated-residual
  genus-zero designs and JC(2) remain open.
source: jc_zero_debt_lift / post-THM-3915 contracted-divisor lane, 2026-08-23
audit: >
  PROVISIONAL PROOF CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  The exact
  companion verifies the rational chart, both degree-three root descents,
  squarefree degree-eight discriminant, all three infinity Hensel branches,
  Riemann--Hurwitz ledger, and the rational F_z hostile.  Normal and optimized
  runs byte-match the frozen output and documentation passes.  The centre-
  trichotomy proof is written intrinsically on a common smooth model; no
  computation is being used to replace that birational-surface argument.
depends_on:
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
script: 04-computation/jc2_positive_genus_collapsed_valuation_thm3916.py
output: 05-knowledge/results/jc2_positive_genus_collapsed_valuation_thm3916.out
script_sha256: 1e658a4c13cb8cfe45d3a50b9f4e546abf283c62f40ccc654fc65686bded7078
output_sha256: 2b359fec32987cf934ba644ea458dee82220d11c877e1e832047604da2a07632
semantic_sha256: 75bbbaebe57df35e61244b7e5c28ee1b5d5a1f9660746c1308bc26b464b7868f
hash_basis: raw LF bytes
---

# THM-3916 -- a contracted genus-two divisor closes the rational decic candidate

**PROVED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**  Work over an
algebraically closed field `k` of characteristic zero.

## 1. Positive-genus common-zero valuations are not Keller

Let

```text
K=k(x,y),                    P,Q in k[x,y],                 (1)
```

and let `nu` be a divisorial valuation of `K/k`.  Denote by `C_nu` the
smooth projective curve with function field `kappa(nu)`.  If

```text
g(C_nu)>0,                   nu(P)>0, nu(Q)>0,              (2)
```

then

```text
Jac_(x,y)(P,Q) notin k*.                                  (3)
```

Indeed, realize `nu` by a prime divisor `E` on a smooth projective surface
`Y` dominating `P2_(x,y)`.  There are three possible centres on `P2`.

1. If `E` maps to a closed point, it is exceptional for a birational
   morphism between smooth surfaces.  Factoring that morphism into point
   blowups shows that `E` is rational.
2. If `E` maps onto the line at infinity, it is birational to that line and
   is again rational.
3. Since `(2)` excludes the first two cases, `E` maps birationally onto the
   closure of an irreducible affine plane curve `D=V(r)`.

In the third case the centre of `nu` on `A2` is the height-one prime `(r)`.
The inequalities in `(2)` therefore give

```text
P=r^a P_0,                 Q=r^b Q_0,       a,b>=1.        (4)
```

Taking exterior derivatives gives

```text
dP wedge dQ = r^(a+b-1) dr wedge
               (a P_0 dQ_0-b Q_0 dP_0)
             +r^(a+b)dP_0 wedge dQ_0.                     (5)
```

Thus `r` divides the Jacobian in `(3)`, a contradiction if it were a
nonzero scalar.  This proves the general valuation obstruction.  Notice
that positivity of the residue genus is load-bearing: a rational
exceptional divisor may be contracted by a birational plane map, for
example `(x,y)|->(x,xy)`, whose Jacobian pays the factor `x`.

## 2. The rational chart carries a hidden cubic curve

Use the rational coordinates of THM-3915, renaming its slope variable `u`:

```text
p=(u^2-4)(u^2+1)^2,
h=u^9-3u^7-9u^5+5u^3+30u,

F(u,z)=z^3-3p z+2h.                                      (6)
```

The target functions and the original cubic root are

```text
A=F/4,                 C=uF/4,                 Z=A^3z,    (7)
```

and conversely `u=C/A,z=Z/A^3`.  Hence the cubic field is `k(u,z)`, while

```text
Jac_(u,z)(A,C)=-F F_z/16.                                 (8)
```

The factor `F_z=0` is the ordinary ramification curve.  The other factor is
more subtle: the prime `F=0` is contracted by `(7)` to the target origin,
and

```text
ord_F(A)=ord_F(C)=1.                                      (9)
```

The next two sections prove that its residue curve has genus two.

## 3. Irreducibility of the collapsed cubic

The polynomial `F` is irreducible in `k[u,z]`.  Since it is monic of degree
three in `z`, reducibility over `k(u)` would give a root in `k(u)`.  The
integral closedness of `k[u]` makes that root a polynomial `r(u)`.

Put `n=deg r`.  If `n<=2`, the degree-nine term `2h` is uniquely highest in
`F(u,r)`.  If `n>=4`, the term `r^3`, of degree `3n`, is uniquely highest.
Only `n=3` remains.  Write

```text
r=a u^3+b u^2+c u+d.                                     (10)
```

The coefficient of `u^9` is

```text
(a-1)^2(a+2).                                             (11)
```

For `a=-2`, the coefficients of `u^8,u^7,u^6` successively force

```text
b=0,                      c=2,                       d=0,
```

but the coefficient of `u^5` is then `-72`.  For `a=1`, the coefficients
of `u^7,u^5,u^3` successively force

```text
b=0,                      c=-1,                      d=0,
```

but the coefficient of `u` is then `48`.  Both branches are impossible,
proving irreducibility.

## 4. Exact genus-two ledger

The discriminant of `(6)` in `z` is

```text
disc_z(F)=432 H(u),

H=24u^8-29u^6-246u^4-309u^2-16.                          (12)
```

Exact Euclidean division gives

```text
gcd(H,H')=1.                                              (13)
```

Thus the degree-three map from the normalization of `F=0` to `P1_u` has
eight finite simple branch points.  Equivalently, at each such point the
power discriminant has valuation one; an index correction has even
valuation, so the normalized cubic still has tame transposition inertia.

There is no hidden ramification at infinity.  Put

```text
w=1/u,                         v=z/u^3.                    (14)
```

Then

```text
w^9 F(1/w,v/w^3)
 =v^3-3v+2+6(v-1)w^2+(21v-18)w^4
   +(12v+10)w^6+60w^8.                                   (15)
```

At `w=0` the leading polynomial is

```text
(v-1)^2(v+2).                                             (16)
```

The root `v=-2` is simple.  For the double leading root substitute

```text
v=1-w^2-4w^4+w^5q.                                       (17)
```

After division by `w^10`, the residual polynomial at `w=0` is

```text
3(q^2-32).                                                (18)
```

It has two distinct roots.  Hensel therefore gives three distinct infinity
branches, all with uniformizer `w`; all three are unramified over infinity.
Riemann--Hurwitz now gives

```text
2g-2=3(-2)+8=2,                         g=2.               (19)
```

As a hostile control, the second factor in `(8)` is rational.  Indeed

```text
F_z/3=z^2-(u^2-4)(u^2+1)^2,
```

and after `z=(u^2+1)r` its normalization has function field

```text
r^2=u^2-4.                                                (20)
```

Thus the positive genus is carried specifically by the contracted
common-zero divisor `F=0`, not by every divisor of the natural Jacobian.

## 5. The THM-3915 candidate has no Keller plane model

Suppose the THM-3915 cubic completion were the finite normalization of a
degree-three planar Keller map.  By normalization-form Zariski Main, as in
THM-3801, the source `A2_(x,y)` would be an open subset of that completion,
with

```text
k(x,y)=Frac(B)=k(u,z),
A,C in k[x,y],                     Jac_(x,y)(A,C) in k*.   (21)
```

Take `nu=ord_F`.  Equations `(9),(19)` satisfy the two hypotheses `(2)` of
the general obstruction, contradicting `(21)`.  Consequently:

```text
the THM-3915 normal cubic completion contains no affine-plane open
on which A,C form a Keller pair.                            (22)
```

In particular THM-3915 is not a planar Jacobian counterexample.  This is
stronger than its Euler tariff: the maximal etale locus need not be a plane,
and now no proper plane open can realize the required degree-three Keller
field either.

## 6. Exact scope and replay

The general lemma applies to a positive-genus divisorial valuation that is
positive on both target functions.  It does not exclude a chart whose
collapsed valuations are all rational, a cubic residual with repeated roots
whose normalization loses the genus in `(19)`, a different cubic field, or
`JC(2)` in general.  In particular, the emerging repeated-residual
square-class designs are **OPEN** and are the correct successor rather than
an extrapolation from this one decic.

Reproduce the exact chart, irreducibility, discriminant, infinity, genus,
and hostile-control gates with

```bash
python3 04-computation/jc2_positive_genus_collapsed_valuation_thm3916.py
python3 -O 04-computation/jc2_positive_genus_collapsed_valuation_thm3916.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
