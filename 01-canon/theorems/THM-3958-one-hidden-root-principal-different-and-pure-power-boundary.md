---
id: THM-3958
title: "One hidden root forces a principal different and closes the natural cubic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. If the hidden cubic
  E+C h^2-2h^3 has exactly one root over k(t), that root is polynomial and
  the irreducible quadratic cofactor has the unique normal form
  h^2+a h+a r. The associated natural cubic surface is an integral normal
  finite-flat model. Its relative derivative has reduced divisor equal to
  the sum of the graph and residual ramification primes. Every Keller
  affine-plane open must delete both primes, making the derivative a
  forbidden nonconstant unit; equivalently, their principal class relation
  violates the affine-plane boundary-basis theorem. In the sharp
  one-address family r=c t^m, r+2a=d t^n, unequal exponents give residual
  genus |m-n|-1 and exactly two infinity places; equal exponents split.
  Common positive exponents also make the residual prime non-unibranch.
  The two remaining genus-zero hostiles have class groups Z and 0 and fail
  only at the principal-different gate.
source: jc-degree6-one-place / post-THM-3956 one-root residual audit, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-24). The audit
  independently reconstructed the monic integrality and unique-root normal
  form, the irreducibility criterion and finite singular locus proving
  normality, the reduced scheme quotient A/(J) and multiplicity-one graph
  and residual divisor packet, and the normalization-form Zariski Main
  forbidden-unit contradiction. It also checked the squarefree pure-power
  residual normalization, both infinity places, the common-positive
  two-address non-unibranch fibre, the exact localization
  A[L^-1]=k[t,L,L^-1], every Nagata valuation and class-group row, and the
  sharp genus-zero endpoints. Normal and optimized runs byte-match the
  frozen 161-gate output; raw and semantic hashes and documentation checks
  pass.
depends_on:
  - THM-3956-split-hidden-cubic-integrality-and-repeated-root-trichotomy
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
related:
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
script: 04-computation/jc2_one_hidden_root_principal_different_pure_power_thm3958.py
output: 05-knowledge/results/jc2_one_hidden_root_principal_different_pure_power_thm3958.out
script_sha256: f8f2d5415b966d03072438fadc6dc90f3b12d94ab4474e0ff8040041d6234f6e
output_sha256: 2e01d69d6e95a5c1744711d8ccea6bad6c9cf90a395b52206b41e0ab35efeb82
semantic_sha256: 74d85a5e93df9df20672bebdea8f6e8c4818accbdf72c5fee1905628e0a5c335
hash_basis: raw LF bytes
---

# THM-3958 -- one hidden root makes the different a forbidden unit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let

```text
C,E in k[t],
G(h,t)=E+C h^2-2h^3,                                      (1)
F(T,P,t)=T^3-3PT-(E+CP),
X=Spec A,       A=k[P,t,T]/(F),       pi:X -> A2_(P,t).    (2)
```

Assume that `G` has **exactly one** root over `k(t)`: it has one linear
factor and one irreducible quadratic factor. THM-3956 treats the complementary
case in which the quadratic also splits. The conclusion here concerns the
natural monogenic model `(2)` and same-function-field Keller charts, not an
arbitrary nonmonogenic overorder or JC(2).

## 1. The unique marked-root normal form

The root in `k(t)` is integral over `k[t]` because `-G/2` is monic. Since
`k[t]` is normal, it is a polynomial

```text
r in k[t].                                                   (3)
```

Polynomial division and the missing linear `h`-row give unique `a in k[t]`
such that

```text
G=-2(h-r)(h^2+a h+a r),
C=2(r-a),                    E=2a r^2.                    (4)
```

It is convenient to put

```text
s=r+2a,                       a=(s-r)/2.                  (5)
```

Then

```text
C=3r-s,                       E=(s-r)r^2,                 (6)
Q(h)=h^2+(s-r)h/2+(s-r)r/2,                               (7)
disc_h(Q)=(s-r)(s-9r)/4.                                  (8)
```

The hypothesis that `Q` is irreducible over `k(t)` is exactly that `(8)` is
not a square in `k(t)`. In particular

```text
r!=0,             s!=0,             s-9r!=0              (9)
```

as polynomial identities. The graph root meets the residual quadratic
precisely where

```text
Q(r)=r s=0.                                                (10)
```

Thus `r=0` and `s=0` are the two exact pair-incidence carriers: the first is
the smooth repeated-zero collision when `s!=0`, while the second is the
singular nonzero `2:-1` collision when `r!=0`. A common zero `r=s=0` is the
triple-root singular override.

## 2. The natural surface is integral and normal

For a general coefficient pair `(C,E)`, reducibility of the monic cubic `F`
would supply a root over `k(P,t)`. That root is integral over the UFD
`k[P,t]`, hence belongs to this normal ring. Comparison of `P`-degrees makes
it independent of `P`; its two coefficient rows then force

```text
T=-C/3,                         C^3+27E=0.                (11)
```

Conversely `(11)` gives that polynomial factor. In the marked-root normal
form, however,

```text
C^3+27E=s^2(9r-s).                                        (12)
```

Equation `(9)` makes `(12)` nonzero. Hence `F` is irreducible and `A` is a
domain.

Shift the graph root to the origin by

```text
X_1=T+r,                         Y_1=P-r^2.               (13)
```

The equation becomes

```text
H=X_1^2(X_1-3r)+Y_1(s-3X_1)=0.                           (14)
```

Its relevant partial derivatives are

```text
H_(Y_1)=s-3X_1,
H_(X_1)=3(X_1^2-2rX_1-Y_1),
H_t=-3r'X_1^2+s'Y_1.                                     (15)
```

At a singular point, the first two rows force

```text
X_1=s/3,                    Y_1=X_1^2-2rX_1.             (16)
```

Substitution in `H` then gives

```text
H=s^2(s-9r)/27.                                           (17)
```

By `(9)`, its zero set in `A1_t` is finite. For each such parameter `(16)`
gives at most one point, and the last row of `(15)` only removes candidates.
Thus the singular locus is finite. The hypersurface domain is `S2`; its
singular locus has codimension two, so it is `R1`. Serre's criterion proves

```text
X is normal.                                               (18)
```

Since `(2)` is monic of degree three, `pi` is finite free of rank three.
Normality also says that `X` is the actual finite normalization of the base
in its function field.

## 3. The different is exactly two reduced primes

In the shifted coordinates define

```text
J=X_1^2-2rX_1-Y_1=T^2-P.                                 (19)
```

Then

```text
F_T=3J.                                                    (20)
```

Modulo `J`, equation `(14)` factors exactly as

```text
H|_(J=0)=-X_1 R,
R=2X_1^2-(3r+s)X_1+2rs.                                  (21)
```

Under `h=r-X_1`, one has

```text
R=2Q(h).                                                   (22)
```

Therefore `R` is irreducible in `k[t,X_1]`. It is distinct from `X_1`
because `rs` is not the zero polynomial. The two height-one primes are

```text
E_G=(J,X_1),                         graph ramification,
E_R=(J,R),                           residual ramification. (23)
```

Equation `(21)` gives the reduced scheme identity

```text
A/(J)=k[t,X_1]/(X_1R),
div_X(J)=E_G+E_R.                                          (24)
```

There are no vertical components and both multiplicities are one. The
element `J` is nonconstant: in the free base-module basis `(1,T,T^2)`, its
`T^2` coefficient is one, so it cannot equal a scalar.

Suppose a same-function-field planar Keller map existed. Normalization-form
Zariski Main would identify its source with an open

```text
U isomorphic to A2 -> X                                   (25)
```

on which `pi` is etale. The differential presentation

```text
Omega_(A/k[P,t])=A dT/(3J dT)                             (26)
```

shows that the generic points of both primes in `(23)` are ramified. Hence
both primes lie in `X minus U`. By `(24)`, `J` has no zero on `U`, so

```text
J|_U in Gamma(U,O_U)^*.                                   (27)
```

But an affine plane has only scalar units, contradicting the nonconstancy of
`J`. Equivalently, `(24)` gives the forbidden positive boundary-class
relation

```text
[E_G]+[E_R]=0,                                             (28)
```

where THM-3922 requires the actual prime boundary classes to be a free
basis. This proves the universal no-atlas conclusion for the exactly-one-root
natural cubic.

## 4. Pure-power one-address classification

The universal proof above is already complete. The remaining calculation
explains the sharp geometry that survives every coarser gate. Translate the
only finite pair-incidence address to `t=0` and suppose

```text
r=c t^m,                 s=d t^n,
c,d in k*,               m,n>=0.                          (29)
```

This is exactly the case in which the zero divisors of the two incidence
carriers `r` and `s` have at most one finite support point.

If `m=n`, equation `(8)` is a scalar times `t^(2m)`, hence a square in
`k(t)`; the residual quadratic splits or repeats, routing back to THM-3956.
Assume `m!=n` and put

```text
kappa=|m-n|.                                               (30)
```

After removing the square `t^(2 min(m,n))`, the residual discriminant is

```text
(d t^kappa-c)(d t^kappa-9c),       if n>m,
(d-c t^kappa)(d-9c t^kappa),       if m>n.                (31)
```

The two degree-`kappa` binomials are separable and have disjoint zero sets.
Thus the squarefree degree is exactly `2kappa`. The smooth projective
normalization of `E_R` is the hyperelliptic curve defined by the square root
of `(31)`, so

```text
genus(E_R^bar)=kappa-1,          # infinity points=2.      (32)
```

Because the residual equation `(7)` is monic over `k[t]`, its affine
normalization is finite over `A1_t`; the two points over infinity are its
only projective punctures.

If `m,n>0`, both residual roots specialize to `h=0` at `t=0`, while the
square removed in `(31)` yields two distinct unramified normalization points
there. Both points map to

```text
t=0,                 T=P=0                               (33)
```

on the single prime `E_R`. Hence `E_R` is non-unibranch at `(33)` and cannot
be a boundary prime of an affine-plane open by THM-3920. If exactly one of
`m,n` is zero and `kappa>=2`, `(32)` instead makes `E_R` positive genus,
again impossible for an affine-plane boundary curve. These are independent
geometric obstructions, but neither covers the last equality cases.

## 5. The sharp genus-zero hostiles and exact class ledger

The only pure-power packets that survive the non-unibranch and genus gates
are

```text
(m,n)=(0,1) or (1,0).                                    (34)
```

Here `E_R^nu` is `P1` minus its two infinity points, hence `G_m`; the graph
and residual primes have exactly one finite incidence. Thus the
boundary-incidence forest, rationality, and unibranch tests alone do not
close `(34)`.

Their class groups also show why torsion is not the missing invariant. Put

```text
L=s-3X_1,                    M=s^2(s-9r).                 (35)
```

Localizing `(14)` at `L` gives

```text
A[L^(-1)]=k[t,L,L^(-1)].                                  (36)
```

If

```text
M=lambda product_j (t-alpha_j)^(e_j),                    (37)
```

then the primes above `L=0` are exactly

```text
Q_j=(L,t-alpha_j),                 ord_(Q_j)(L)=e_j.      (38)
```

Indeed at `L=0`, the non-`Y_1` term of `(14)` is `M/27`; at the generic
point of `Q_j`, `Y_1` is transcendental and the remaining coefficient is a
unit. Nagata localization gives

```text
A^*=k*,
Cl(A)=Z^(#Q_j)/<(e_j)>.                                   (39)
```

For `(29)` with unequal exponents, the exact multiplicity rows are

```text
n>m:   M=d^2 t^(2n+m)(d t^kappa-9c),
m>n:   M=d^2 t^(3n)(d-9c t^kappa).                       (40)
```

Every nonzero root in `(40)` is simple. Consequently all these class groups
are torsion-free. More explicitly,

```text
n>m:             Cl(A)=Z^kappa,
m>n>0:           Cl(A)=Z^kappa,
m>n=0:           Cl(A)=Z^(kappa-1).                       (41)
```

In the two sharp rows `(34)`, the Nagata relation vectors and class groups
are

```text
(m,n)=(0,1):        (2,1),             Cl(A)=Z,
(m,n)=(1,0):        (1),               Cl(A)=0.           (42)
```

Both have scalar units. The first even has the free positive rank expected
by the coarse torsion-free class-group gate, while the second has no class
group at all. What kills both is the finer, labelled relation `(28)` between
the **two mandatory ramification primes**: deleting both turns their
principal defining function `J` into a unit. In particular, the rank-one row
still cannot supply the rank-two boundary basis actually required here.
This is the exact structural lesson of the pure-power hostile.

## 6. Scope

Combining THM-3956 and the present theorem closes every natural hidden cubic
having at least one root over `k(t)`. A remaining natural cubic must have
irreducible hidden cubic ramification, and a viable nontrivial Keller
completion must in any event cease to be globally monogenic, in agreement
with THM-3862. The theorem does not transfer its principal different to a
different nonmonogenic overorder and does not prove JC(2). **QED.**
