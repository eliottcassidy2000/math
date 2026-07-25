---
id: THM-2189
title: "The nonsplit quartic deck forces the remaining pole congruence"
status: >
  PROVED. In the reduced twice-odd planar Keller quartic branch, write
  P=V^2z^4+Vbz^3+gamma z^2+delta z+epsilon after THM-2180. If V is not a
  square in C(x), the genuine quadratic deck forces
  V|(4gamma-b^2). The proof first uses the boundary/flux triple to show that
  every hypothetical bad finite place has one exact square initial face.
  Deck parity removes all odd Faber seeds and makes the first flux identically
  zero. Completing the square at the next face then produces a unique
  noncancellable Faber tooth in either the first or second flux. Hence the
  canonical approximate root is polynomial throughout the nonsplit branch.
  The same congruence holds for reduced mate degree two even on a split deck.
  THM-2202 subsequently extends the same pole closure to every split
  twice-odd degree, so no quartic finite-pole survivor remains.
  Both polynomial-root branches still need the terminal nonmonic
  square-prefix descent.
source: codex-2026-07-24-JC-nonsplit-quartic-pole
depends_on:
  - THM-2129
  - THM-2158
  - THM-2180
related:
  - THM-2181
  - THM-2194
  - THM-2202
---

# THM-2189 -- the nonsplit deck kills the last finite pole

Let

```text
P=V^2z^4+beta z^3+gamma z^2+delta z+epsilon
```

belong to a planar polynomial Keller pair over `C`. Reduce the mate by
polynomial target shears, and assume its remaining fibre degree is

```text
n=4r-2,                         r>=1.                 (1)
```

THM-2180 gives `beta=Vb` with `b in C[x]`. Put

```text
D=4gamma-b^2.                                         (2)
```

Assume that

```text
V is not a square in C(x).                            (3)
```

Then

```text
V divides D.                                          (4)
```

By THM-2158, (4) says exactly that the canonical quadratic approximate root

```text
H_0=Vz^2+(b/2)z+D/(8V)                               (5)
```

is polynomial.

## 1. Deck parity and the two fluxes

Let `K=C(x)` and work in the genuine quadratic field

```text
L=K(U),                         U^2=V.                (6)
```

Its nontrivial deck involution is `sigma(U)=-U`. Monicize with `w=Uz` and
depress as in THM-2158:

```text
Z=w+b/(4U),
P=Z^4+pZ^2+qZ+s.                                     (7)
```

For every reduced Faber degree `m`, let

```text
E_m=Pol_Z(P^(m/4)).                                   (8)
```

THM-2180's triangular normal form is

```text
Q=H(P)+sum_(4 does not divide m, m<=n)c_mE_m,
H in C[T],                       c_m in C,            (9)
```

with `c_n!=0`.

Weighted parity gives

```text
sigma(E_m)=(-1)^m E_m.                               (10)
```

Indeed, `sigma(Z)=-Z`, the even depressed coefficients are fixed, the odd
coefficient is negated, and the normalized Laurent branch has leading term
`Z^m`. Since `Q` and `H(P)` are fixed and the Faber normal form is triangular
by fibre degree, uniqueness in (9) forces

```text
c_m=0                         for every odd m.        (11)
```

Thus every surviving reduced seed has degree `m=4j-2`.

Use THM-2129's notation

```text
P^(m/4)=E_m+c_(m,1)Z^-1+c_(m,2)Z^-2+...,
Phi_m=4c_(m,1),                 Psi_m=4c_(m,2).       (12)
```

Its exact Hamiltonian identity implies that the combinations

```text
Phi_Q=sum_m c_m Phi_m,
Psi_Q=sum_m c_m Psi_m                              (13)
```

have zero `x`-derivative. The constant field of the algebraic function field
`L` remains `C` because `C` is algebraically closed, so both belong to `C`.
For even `m`, Laurent parity gives

```text
sigma(Phi_m)=-Phi_m,             sigma(Psi_m)=Psi_m. (14)
```

Therefore

```text
Phi_Q=0,                         Psi_Q in C.          (15)
```

The exact zero in (15), rather than merely constancy, is the load-bearing
benefit of a nonsplit deck.

## 2. Every bad place enters the square cage

Fix an irreducible `pi|V`, extend its valuation `v` to `L`, and write

```text
e=v(V)>0,        f=v(b),        d=v(D),
g=v(gamma),      ell_0=v(delta).                     (16)
```

Suppose toward a contradiction that

```text
d<e.                                                  (17)
```

In the monic coordinate `w=Uz`, the reversed boundary polynomial at the
original fibre `w=0` is

```text
T(u)
 =1+(b/U)u+(gamma/V)u^2+(delta/U)u^3+epsilon u^4.    (18)
```

Put

```text
alpha_1=f-e/2,
alpha_2=g-e,
alpha_3=ell_0-e/2,
s_0=min(alpha_1,alpha_2/2,alpha_3/3).                (19)
```

One has `s_0<0`. If `alpha_1<0`, this is immediate. Otherwise
`2f>=e`; because `d<e<=2f`, the two terms in `D=4gamma-b^2` have unequal
valuations and `g=d<e`, so `alpha_2<0`.

Pass harmlessly to a valued extension containing `lambda` with
`v(lambda)=s_0`, and set

```text
S(X)=T(X/lambda).                                     (20)
```

Its residue is a nonconstant polynomial of degree at most three,

```text
S_0(X)=1+A X+B X^2+C X^3,                            (21)
```

because the `epsilon` slope is nonnegative while `s_0<0`.

For

```text
A_(m,k)=[X^k]S(X)^(m/4),                              (22)
```

translation from `Z=w+h`, `h=b/(4U)`, gives the exact boundary/flux
identities

```text
E_m(h)=lambda^m A_(m,m),

Phi_m=4lambda^(m+1)A_(m,m+1),

Psi_m=4lambda^(m+2)
  [A_(m,m+2)+(h/lambda)A_(m,m+1)].                   (23)
```

The last line uses residue invariance under translation:
the `w^-2` coefficient is `c_(m,2)-h c_(m,1)`.
Moreover

```text
v(h/lambda)=alpha_1-s_0>=0.                          (24)
```

Because `s_0<0`, the top degree `m=n` in (9) is uniquely deepest in each
line of (23). The polynomial boundary value `Q(x,0)-H(epsilon(x))` is
regular. Equations (15), used successively with (23)--(24), therefore force
the three leading residues

```text
[X^n]S_0^(n/4)
 =[X^(n+1)]S_0^(n/4)
 =[X^(n+2)]S_0^(n/4)=0.                              (25)
```

We now classify (25), including the cases where the linear face vanishes.

- If `A!=0`, rescale `X` so that `A=4`. THM-2129's all-degree
  boundary-triple classification, with `n=4r-2`, gives the unique square

  ```text
  S_0(X)=(1+(A/2)X)^2.                               (26)
  ```

- If `A=0` and `C!=0`, THM-2129's coefficient recurrence with three
  consecutive zeros descends to the contradiction that the constant
  coefficient is zero; `3n/4` is not an integer.

- If `A=C=0` and `B!=0`, then

  ```text
  [X^n](1+BX^2)^(n/4)
   =B^(n/2) binom(n/4,n/2)!=0,                       (27)
  ```

  because `n/4` is a half-integer.

The all-zero case contradicts the activity of the face. Thus (26) is the
only possibility.

Equality of the active linear and quadratic slopes in (26), strictness of
the cubic slope, and cancellation of the quadratic square give

```text
s_0=f-e/2,
g=2f,
d>2f,
ell_0+e>3f.                                          (28)
```

Together with (17),

```text
2f<d<e.                                               (29)
```

This derives the complete square cage; no separate genericity or balanced-
face assumption is being imported.

## 3. Complete the next square exactly

Put

```text
a=b/(2U),               A_0=-v(a)=e/2-f>0,
C_0=d-2f>0,             H=e-d>0.                     (30)
```

Then

```text
A_0=(C_0+H)/2.                                        (31)
```

Scale (18) canonically by `X=au`. The exact normalized polynomial is

```text
T_hat(X)
 =1+2X+(1+c)X^2+lX^3+mX^4,                          (32)

c=D/b^2,
l=8delta V/b^3,
m=16epsilon V^2/b^4.
```

Equations (28)--(31) give

```text
v(c)=C_0>0,            v(l)>0,            v(m)>=4A_0.
                                                               (33)
```

Complete the full quadratic prefix:

```text
K=1+X+(c/2)X^2,
T_hat=K^2+B,
B=Lambda X^3+Omega X^4,                              (34)

Lambda=(8delta V-bD)/b^3,
Omega=(64epsilon V^2-D^2)/(4b^4).
```

The two terms in `Omega` have valuations at least `4A_0=2(C_0+H)` and
exactly `2C_0`, respectively. Hence

```text
v(Omega)=2C_0,
rho=min(v(Lambda),v(Omega))
       satisfies 0<rho<=2C_0<4A_0.                  (35)
```

Write the first residue face as

```text
B_rho=X^3(L+MX),                  (L,M)!=(0,0).       (36)
```

## 4. The first rational tooth cannot cancel

For an even reduced degree `m=4j-2`,

```text
T_hat^(m/4)
 =sum_(q>=0) binom(j-1/2,q)
      K^(2j-1-2q) B^q.                               (37)
```

For `q<j`, the summand is polynomial of degree at most

```text
2(2j-1-2q)+4q=m,                                     (38)
```

so it contributes neither first nor second flux. The first rational tooth
is `q=j`. At the residue face it is a nonzero scalar multiple of

```text
X^(3j)(L+MX)^j/(1+X).                                (39)
```

For the top degree `n=4r-2`, put

```text
C_k=[X^k](L+MX)^r/(1+X).                             (40)
```

The first and second flux components of (39) are proportional to

```text
C_(r-1),                    C_r+(1/2)C_(r-1),         (41)
```

where `1/2=h/a`. They cannot both vanish. Indeed,

```text
C_k+C_(k-1)=[X^k](L+MX)^r.                           (42)
```

If `C_(r-1)=C_r=0`, equation (42) at `k=r` first gives `M=0`, and then

```text
C_(r-1)=(-1)^(r-1)L^r=0,
```

contradicting `(L,M)!=(0,0)`. This also covers `r=1`.

The valuations of the top first and second flux teeth are

```text
-(n+1)A_0+r rho,
-(n+2)A_0+r rho=r(rho-4A_0)<0.                       (43)
```

Only the second displayed valuation is asserted to be negative; the first
tooth is used through the exact identity `Phi_Q=0`. No lower Faber degree
can cancel either tooth. For `m'=4j-2<n` and any
nonpolynomial order `q>=j`, its `t`-th flux valuation is at least

```text
-(m'+t)A_0+q rho.
```

Subtracting this from the corresponding top valuation gives

```text
-4(r-j)A_0+(r-q)rho<0:                               (44)
```

if `q<=r`, use `q>=j` and `rho<4A_0`; if `q>r`, both terms are already
negative. Later orders in the top degree are also strictly shallower.

If the first component in (41) is nonzero, its unique depth contradicts
the exact identity `Phi_Q=0`. Otherwise the second component is nonzero,
and its negative valuation in (43) contradicts `Psi_Q in C`. This is the
desired contradiction to (17), proving (4).

## 5. Uniform degree-two addendum

The congruence (4) is also forced when the reduced mate degree is

```text
n=2,                                                   (45)
```

without the nonsplit assumption (3).

Repeat Section 2 before using deck parity. At its initial negative scale,
the top degree two is uniquely deepest against the sole lower degree one.
Boundary regularity and constancy of both fluxes therefore still force the
square cage (28)--(35).

Normalize the coefficient of `E_2` to one and write

```text
Q=H(P)+E_2+kE_1,                       k in C.         (46)
```

Let

```text
Delta=s-p^2/4
```

for the depressed coefficients in (7). Direct Laurent expansion gives

```text
E_1=Z,                     E_2=Z^2+p/2,
Phi_1=p,                   Psi_1=q,
Phi_2=2q,                  Psi_2=2Delta.              (47)
```

In the normalization of Section 3, `h=a/2`, and exact reverse square
completion gives

```text
p=a^2(c-1/2),
q=a^3 Lambda,
Delta=a^4(Omega-Lambda/2),                            (48)

E_1(h)=a/2,                    E_2(h)=a^2c/2.
```

Thus boundary regularity and flux constancy require the following three
expressions to be regular or constant:

```text
B_Q=(a^2c+k a)/2,

Phi_Q=2a^3 Lambda+k a^2(c-1/2),

Psi_Q=a^4(2Omega-Lambda)+k a^3 Lambda.                (49)
```

Write `A=A_0` and `C=C_0`. Since `0<C<2A`, regularity of `B_Q` forces
`C=A` and `k!=0`: its two polar terms must have equal valuation.
Now `c-1/2` is a unit because `v(c)=C>0`. Constancy of `Phi_Q` then forces

```text
v(Lambda)=A.                                         (50)
```

But in `Psi_Q`, the summand `-a^4Lambda` has valuation `-3A`, whereas

```text
2a^4Omega,                  k a^3Lambda
```

both have valuation `-2A`, using `v(Omega)=2C=2A`. The unique deeper pole
cannot cancel. This contradiction proves (4) uniformly for `n=2`.

THM-2071 already closes a Keller pair with a quadratic member; the point of
this addendum is the sharper local pole mechanism.

## 6. Exact surviving pole branch

We have proved the remaining pole congruence at every divisor of `V` when
the quadratic deck is nonsplit. If `V` is a square in `C(x)`, unique
factorization and algebraic closure of `C` give

```text
V=W^2                         for some W in C[x].     (51)
```

Thus every survivor of the quartic pole analysis has leading coefficient

```text
V^2=W^4.                                               (52)
```

THM-2202 subsequently closes every reduced degree `4r-2` in this split
fourth-power branch, so no quartic finite-pole survivor remains. Deck
anti-invariance did not force `Phi_Q=0` there; the later theorem replaces it
with an all-degree two-chamber filtration.

This theorem is only an exact **pole** result. After THM-2202, `H_0` is
polynomial on both branches but its leading term `Vz^2` is still nonmonic.
The terminal square-prefix/quadratic-member step remains. Neither theorem
proves general JC(2) or DC(2). QED.
