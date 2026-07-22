---
id: THM-2132
title: "Factor-initial two-face closure and the local Newton-chord escape"
status: >
  PROVED. Let a proper-power positive-weight top face be h^m and let A be
  the first lower face. If rad(h) divides A and there are no further
  nonconstant faces, the component has no Jacobian mate: a radical
  factorization creates a critical point, with the only apparent exception
  reducing to a composite P(h) that has a critical fiber. Together with
  THM-2127 this closes every exact two-face proper-power component. With
  arbitrary later faces and noncentral A, any mate must either satisfy an
  explicit small-degree inequality or acquire a later face on or below a
  factorwise Newton chord before the terminal defect. The chord equality is
  sharp for x+(x^2+y)^k. This narrows but does not close the general
  factor-initial locus or planar JC.
source: codex-2026-07-22-JC2-factor-initial-newton-chord
depends_on:
  - THM-2127
related:
  - THM-2102
  - THM-2113
  - THM-2129
---

# THM-2132 -- factor-initial closure and Newton-chord escape

Let `w=(w_1,w_2)` be positive integral weights and put

```text
W=w_1+w_2.                                             (1)
```

Let `h in C[x,y]` be nonconstant, `w`-homogeneous, power-free, and of
weighted degree `d`. Power-free means that the gcd of the multiplicities in
the irreducible factorization of `h` is one.

## 1. The exact factor-initial two-face branch is impossible

Let `m>=2`, let `A!=0` be `w`-homogeneous of degree

```text
0<a<md,                                                (2)
```

and assume

```text
rad(h) divides A.                                     (3)
```

Then, for every constant `c_0`, the polynomial

```text
f=c_0+h^m+A                                           (4)
```

has no polynomial Jacobian mate.

### Radical factorization creates a critical point

Subtract `c_0`, which changes no Jacobian bracket, and put

```text
R=rad(h).
```

Factorization of a homogeneous element in this positive grading may be taken
through homogeneous irreducibles. Hence `R` is homogeneous, every one of its
nonconstant factors vanishes at the origin, and

```text
f=RQ,                 Q=h^m/R+A/R.                    (5)
```

Both `R` and `h^m/R` vanish at the origin. If

```text
deg_w(A/R)>0,                                          (6)
```

then `Q(0)=0`, and the product rule gives

```text
df(0)=Q(0)dR(0)+R(0)dQ(0)=0.                          (7)
```

This is incompatible with `{f,g}` being a nonzero constant.

It remains that `deg_w(A/R)=0`, so

```text
A=cR,                         c!=0.                    (8)
```

Now `Q(0)=c`. Avoiding (7) would require `dR(0)!=0`. If `R` had at least two
distinct irreducible factors, all of them would vanish at the origin and the
product rule would instead give `dR(0)=0`. Thus `R` is irreducible. Write

```text
h=u R^e,                         u in C^*.              (9)
```

Power-freeness forces `e=1`. After changing the nonzero scalar in (8), (4)
therefore has the form

```text
f=c_0+P(h),                 P(T)=T^m+c_1T,
                                      c_1!=0.          (10)
```

The nonconstant polynomial `P'(T)=mT^(m-1)+c_1` has a root `t in C`. The
nonconstant polynomial map `h:A^2->A^1` is surjective: choose a variable in
which `h` is nonconstant, specialize the other away from the zero set of the
leading coefficient, and apply the fundamental theorem of algebra. Hence
some point `p` has `h(p)=t`, and there

```text
df(p)=P'(h(p))dh(p)=0.                                (11)
```

This final contradiction proves the assertion.

### Complete exact two-face corollary

Suppose (4) has a Jacobian mate but impose no factor hypothesis. If (3)
holds, the theorem just proved gives a contradiction. If (3) fails,
THM-2127 applies with no later nonconstant faces: it makes `h,A` polynomial
coordinates and makes `f` a coordinate. Consequently

```text
every exact two-face proper-power Keller component is a coordinate;        (12)
```

the factor-initial half is empty and the complementary half is triangular.

## 2. Arbitrary tails must escape a local Newton chord

Let `f,g in C[x,y]` satisfy `{f,g}=1`, and use the full-face notation

```text
f=sum_(delta>=0) f_delta,
f_0=h^m,                 f_rho=A,
a=deg_w(A)=md-rho>0,                                  (13)
```

where `A` is the first lower face. Let a reduced mate have top face

```text
g_0=c h^n,           c!=0,            n=qm+r,
q>=0,                                  0<r<m,
D=(m+n)d-W>0.                                         (14)
```

Assume now

```text
rad(h) divides A,                 A notin C[h].        (15)
```

Because `A` is homogeneous, the second condition is equivalent to saying
that `A` is not a scalar power of `h`.

Factor

```text
h=u product_i pi_i^(e_i),        d_i=deg_w(pi_i),
b_i=v_(pi_i)(A).                                      (16)
```

Then some factor `pi=pi_i` satisfies

```text
Delta=e a-b d>0,                e=e_i, b=b_i.         (17)
```

For such a factor define

```text
mu=(me-b)/rho,
J=floor(en/(me-b))+1.                                 (18)
```

Here `me-b>0`. The following dichotomy is necessary for a mate:

```text
n(ea-bd)/(me-b)<=W-a,                                 (19)
```

or there is a later face `f_delta!=0` with

```text
rho<delta<=J rho<D,
v_pi(f_delta)<=me-mu delta.                           (20)
```

Thus `(delta,v_pi(f_delta))` lies on or below the line joining

```text
(0,me) and (rho,b).                                   (21)
```

In particular, if `a>=W`, the positive left side of (19) cannot be at most
`W-a`; every mate then forces the later Newton vertex (20).

### A noncentral factor exists

Write

```text
A=v C product_i pi_i^(b_i),              gcd(C,h)=1.  (22)
```

Weighted degrees give

```text
a=deg_w(C)+sum_i b_i d_i,
d=sum_i e_i d_i.                                      (23)
```

If (17) failed for every factor, then `b_i d>=e_i a`. Multiply these
inequalities by `d_i` and sum. Equation (23) forces equality everywhere,
`deg_w(C)=0`, and

```text
b_i=(a/d)e_i                         for every i.     (24)
```

Since `gcd_i(e_i)=1`, equation (24) makes `a/d` an integer. It follows that
`A` is a scalar power of `h`, contrary to (15). Hence a factor satisfying
(17) exists. Moreover `b/e<a/d<m`, so `me-b>0` as asserted.

For this factor,

```text
mu d-e
 =[d(me-b)-e rho]/rho
 =(ea-bd)/rho>0.                                      (25)
```

### The decisive tooth forces the chord

Assume first

```text
J rho<D.                                               (26)
```

Apply THM-2127's full formal train through defect `J rho`. In its seed-zero
train, selecting `A` exactly `J` times gives the nonzero term

```text
c binom(n/m,J) h^(n-mJ) A^J.                         (27)
```

Its coefficient is nonzero because `n/m` is not an integer. Its
`pi`-valuation is

```text
v_*=en-J(me-b)<0                                      (28)
```

by the definition of `J`.

Suppose, toward a contradiction, that every later face through this defect
lies strictly above the chord:

```text
v_pi(f_delta)>me-mu delta
                  for rho<delta<=J rho.               (29)
```

A general raw term at total defect `J rho` comes from a train seeded at
`ell d` and selects `N` lower faces of defects `delta_1,...,delta_N`, where

```text
ell d+sum_j delta_j=J rho.                            (30)
```

Its valuation is at least

```text
e(n-ell-mN)+sum_j v_pi(f_(delta_j))
 >=en-mu J rho+ell(mu d-e)
 =v_*+ell(mu d-e).                                   (31)
```

Equality is used for `A=f_rho`; inequality (29) is strict for every later
face. If `ell>0`, (25) makes (31) strictly larger than `v_*`. If `ell=0`,
equality can occur only for `N=J` with every selected face equal to `A`;
every other selection contains a later face and is strict by (29). Hence
(27) is the unique term of least `pi`-valuation.

Its valuation is negative, so it cannot cancel inside the polynomial face
`g_(J rho)`. This contradiction proves the existence of (20) under (26).

If (26) fails, then `D<=J rho`. Since

```text
J<=en/(me-b)+1,                                       (32)
```

we get

```text
(m+n)d-W<=rho[en/(me-b)+1].                           (33)
```

Use `rho=md-a` and rearrange (33) to obtain exactly (19). This completes the
dichotomy.

## 3. Equality is sharp and is a coarsened root

There are two genuinely different equality mechanisms.

### Central multiplicity equality

If `e_i a-b_i d=0` for every `i`, the same weighted-degree summation used in
Section 2 forces `A=lambda h^j`. The tame control

```text
f=y^4+y^2+x,                     g=y                 (34)
```

lies here: a resonant seed cancels the first pole. The Newton-chord theorem
intentionally excludes this central prefix.

### Noncentral chord equality

For every `k>=2`, with ordinary weights, put

```text
f=x+(x^2+y)^k,                    g=x^2+y.             (35)
```

Then `{f,g}=1`. In the notation above,

```text
h=x,       m=2k,       n=2,
A=k x^(2k-2)y,         rho=1,
e=1,       b=2k-2,     Delta=1,
mu=2,      J=2,        D=2k.                          (36)
```

The defect-two face

```text
binom(k,2)x^(2k-4)y^2                                (37)
```

has `x`-valuation `2k-4=me-2mu`, equality in (20). In Rees form,

```text
F(z)=(x^2+yz)^k+x z^(2k-1),
G(z)=x^2+yz,
{F,G}=z^(2k)=z^D.                                    (38)
```

Thus “on or below” in (20) cannot be strengthened to “strictly below.” The
collinear faces in this control assemble into an exact coarsened approximate
root.

### Coarsened approximate-root lemma

Let `Q` be the mate face at the same first defect `rho`, taking `Q=0` when
that face is absent, and suppose its synchronized first-defect class is zero:

```text
m h^(m-1)Q-c n h^(n-1)A=0.                           (39)
```

Let `s` divide `gcd(m,n)` and put `M=m/s`, `N=n/s`. The first faces arise
from a common coarsened root

```text
K=h^s+B,
f=K^M modulo defects >rho,
g=cK^N modulo defects >rho                           (40)
```

if and only if

```text
h^(m-s) divides A.                                    (41)
```

Here `B` is homogeneous of degree `sd-rho`. Necessarily

```text
B=A/[M h^(m-s)],
Q=cN h^(n-s)B.                                       (42)
```

Indeed, (42) follows from (39), while binomial expansion gives

```text
K^M=h^m+M h^(m-s)B+O(2rho),
cK^N=ch^n+cN h^(n-s)B+O(2rho).                       (43)
```

The notation `O(2rho)` means terms of defect at least `2rho`.

For `s=1`, this is ordinary approximate-root divisibility. Proper common
divisors `s>1` expose the coarsened factor-initial mechanism; in (35), `s=2`
and the full chord is the `k`th power of `x^2+yz`.

## 4. Exact scope and next object

The theorem closes the entire exact two-face locus. For arbitrary tails it
forces an explicit later Newton vertex whenever the noncentral first lower
face is too high for (19). It does not classify central prefixes
`A=lambda h^j`, the numerical branch (19) when `a<W`, or the chord initial
polynomial after equality occurs. The next honest object is the coarsened
approximate-root polygon together with its residual faces, not an isolated
first defect.

Tournament Analysis is not faithful here. Candidate vertices include factors
of `h`, weighted faces, resonant trains, valuation chords, and proof
obligations, but there is no intrinsic binary orientation preserving (31).
The faithful carrier is the lower Newton polygon labelled by factor
multiplicities and defects, together with the train seed and polynomiality
sidecars. A tournament would erase the metric slope `mu` and the equality
initial form that distinguish obstruction from the tame control (35).
