---
id: THM-3997
title: "Hasse repair forces residual diagonals and excludes the zero-residual reduced 2:3 cell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In THM-3992's
  normalized reduced 2:3 cusp-pole cell, the p-Taylor/Hasse transform gives
  an exact iff test for a finite Laurent tail to lie in k[p,y]. The first two
  source-normal diagonals force [p^2]R and [p^3]R as explicit functions of
  a,gamma and the first free scalar beta=[y]R; [py]R is next. In particular
  R=0 is impossible. An independent generic-fibre proof identifies the
  source and target as elliptic curves with incompatible potential-good
  reduction under the isogeny induced by a Keller pair. The full reduced
  2:3 cell and JC(2) remain open.
source: root + laurent_rows + bounded_residual + conductor_incidence / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (bounded_residual and conductor_incidence, 2026-08-24). The Laurent
  diagonal caps, two complete normal-row eliminations, every displayed
  coefficient, Hasse forward/reverse iff, and resonant R=0 contradiction were
  independently reconstructed. Two direct third-jet hostiles pass. The
  generic-fibre field identifications, smoothness, j-invariants, valuation
  signs, projective finiteness, and isogeny reduction argument were separately
  audited. Normal and optimized outputs are byte-identical.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
related:
  - THM-3979-two-color-formal-cusp-darboux-lifting
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-3998-reduced-two-three-three-by-at-most-three-source-weight-support-obstruction
external: >
  Jean-Pierre Serre and John Tate, "Good reduction of abelian varieties",
  Annals of Mathematics 88 (1968), 492-517, Section 1; imported only for the
  Neron--Ogg--Shafarevich criterion and potential-good reduction under isogeny.
script: 04-computation/jc2_reduced_23_hasse_residual_thm3997.py
output: 05-knowledge/results/jc2_reduced_23_hasse_residual_thm3997.out
script_sha256: 9a00d3d80c94591a8677b6c16f1306a747c00721a4e517147afd65ae38737d7e
output_sha256: 870b058333301f4d518faa41ecc4e05e02e50b432ff0e9244c51278c4c4407dd
hash_basis: raw LF bytes
---

# THM-3997 -- Hasse repair and the zero-residual no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem continues the normalized reduced `(2,3)` cell of THM-3989 and
THM-3992.  It proves necessary jet identities for a hypothetical pair and
closes the subcell `R=0`.  It does not construct a pair, show that the free
jets extend to elements of `B_2`, lower the pole depths in the remaining
cells, or prove `JC(2)`.

Work over an algebraically closed field `k` of characteristic zero.  Use

```text
s=xt,                       tau=t,
p=s^2+tau=t+x^2t^2,
y=s p=x t^2+x^3t^3,         u=s^2/tau=x^2t.
```

After THM-3992's determinant-one target normalization and centering, put

```text
h=gamma*s,                  gamma != 0,
I=3a^2/4,                   a != 0,
lambda=3a/(2gamma),

G=C^2-A^3+I A+a^3/4
 =gamma*u+lambda*p+R(p,y),  R in (p^2,y).                 (1)
```

The numerical coefficients below refer to this gauge.  The residual fifth
root acts by

```text
(gamma,a,G,R) -> (zeta*gamma,zeta^2*a,zeta*G,zeta*R),
zeta^5=1.                                                    (2)
```

Thus the value of a residual coefficient is gauge-weighted.  Its vanishing
is unchanged by `(2)`, but no invariance under arbitrary later target
automorphisms is claimed.

## 1. Laurent support gives sharp source-normal degree caps

For a finite Laurent polynomial

```text
F(s,tau)=sum_(i>=-d) f_i(s) tau^i,
```

substitution `s=xtau`, `tau=t` gives the exact diagonal formula

```text
[t^r]F(x,t)=sum_(i=-d)^r [s^(r-i)]f_i(s) x^(r-i).         (3)
```

In the normalized cell,

```text
a_-2=gamma^2 s^2,                    c_-3=gamma^3 s^3,
b=a_-1 in s^2 k[s],                 c_-2=(3/2)gamma*s*b.
```

The extreme rows contain exactly the displayed powers: `a_-2` has no term
above `s^2`, and `c_-3` has no term above `s^3`.  The next rows are
`a_-1=b` and `c_-2=(3/2)gamma*s*b`; they can supply the top powers in the
first two diagonals but cannot exceed the bounds obtained from `(3)`.
Consequently the first two positive source-normal diagonals satisfy

```text
A=A_b+tN+t^2M+O(t^3),     deg_x N<=2,  deg_x M<=3,
C=C_b+tK+t^2L+O(t^3),     deg_x K<=3,  deg_x L<=4.        (4)
```

This is only a necessary consequence of `B_2` membership, but that is the
direction used below.  No converse from degree caps to `B_2` membership is
assumed.

THM-3992 supplies the boundary rows

```text
X=gamma*x,
A_b=X^2+a,
C_b=X^3+(3a/2)X.                                         (5)
```

Let

```text
P(A,C)=C^2-A^3+I A+a^3/4,
q(x)=gamma*x^2+lambda=(X^2+3a/2)/gamma.                  (6)
```

Along `(5)`, direct differentiation gives the useful rotated-gradient law

```text
grad P(A_b,C_b)=q(x)*(-C_b',A_b').                       (7)
```

Indeed,

```text
P_A=-3(X^2+a/2)(X^2+3a/2),
P_C= 2X(X^2+3a/2),
A_b'=2gamma*X,
C_b'=3gamma*(X^2+a/2).
```

## 2. The first normal jet has one affine tangent polynomial

The constant term of `J_(x,t)(A,C)=1` is

```text
A_b' K-N C_b'=1.                                         (8)
```

One polynomial solution is

```text
N_*=-2/(3gamma*a),                 K_*=-x/a.             (9)
```

Because `a gamma !=0`, the two boundary derivatives `A_b'` and `C_b'` are
coprime in `k[x]`.  Every polynomial solution of `(8)` is consequently

```text
N=-2/(3gamma*a)+2gamma^2*x*theta(x),
K=-x/a+3gamma*(gamma^2*x^2+a/2)*theta(x).                (10)
```

The caps in `(4)` force

```text
theta(x)=theta_0+theta_1*x.                              (11)
```

This is the first information lost by retaining only the nodal boundary:
the arbitrary polynomial tangent reparameterization collapses to an affine
one once Laurent depth is retained.

## 3. Eliminating the second jet forces `[p^2]R`

Write `v=(N,K)` and `w=(M,L)`.  The coefficient of `t` in the Jacobian is

```text
2 det((A_b,C_b)',w)+det(v',v)=0.                         (12)
```

The coefficient of `t^2` in `P(A,C)` is

```text
g_2=grad(P)|_b dot w + (1/2) Hess(P)|_b(v,v).
```

Using `(7)` and `(12)` eliminates `w` completely:

```text
g_2=K^2-3A_b N^2-(q/2)(N'K-NK').                        (13)
```

Substitution of `(10)` simplifies before imposing `(11)` to

```text
g_2=2gamma*x*theta
    -(3a/(4gamma)+gamma*x^2/2)*theta'
    -5/(6a gamma^2).                                     (14)
```

For the affine row `(11)`, this is

```text
g_2=(3gamma theta_1/2)x^2+2gamma theta_0 x
    -3a theta_1/(4gamma)-5/(6a gamma^2).                 (15)
```

Now name the first residual coefficients

```text
alpha=[p^2]R,                    beta=[y]R.               (16)
```

Since

```text
p=t+x^2t^2,       y=xt^2+x^3t^3,       u=x^2t,
```

the ring side of `(1)` has

```text
[t]G=gamma*x^2+lambda,
[t^2]G=lambda*x^2+beta*x+alpha.                          (17)
```

Comparing `(15)` and `(17)` gives three exact identities:

```text
theta_1=a/gamma^2,
beta=2gamma*theta_0,
alpha=-3a^2/(4gamma^3)-5/(6a gamma^2)
     =-(9a^3+10gamma)/(12a gamma^3).                     (18)
```

Thus `[p^2]R` is not free.  It may vanish on the resonant divisor
`9a^3+10gamma=0`, so `(18)` is not a universal nonvanishing obstruction.
The scalar `beta=[y]R` is the first coefficient not fixed by this diagonal;
equivalently it is the remaining intercept `theta_0` in the fixed gauge.

### 3.1 Translation back to Laurent coefficients

Write ordinary coefficient expansions

```text
q_0=b/(2h)=q_01*s+q_02*s^2+...,
A_0=a+A_01*s+A_02*s^2+... .                              (19)
```

Formula `(3)` gives

```text
N=A_1(0)+A_01*x+2gamma*q_01*x^2.                         (20)
```

Comparison with `(10)` and `(18)` therefore yields

```text
A_1(0)=-2/(3gamma*a),
q_0'(0)=q_01=a/gamma,
[s^2]b=2a,
A_0'(0)=A_01=gamma*beta.                                 (21)
```

The first row of `(21)` recovers THM-3992's boundary Bezout value.  The
remaining three rows are the new divisibility/jet information.

## 4. One further normal diagonal forces `[p^3]R`

This section is optional for isolating the first free scalar, but it shows
that the mechanism repeats.  Put

```text
M=m_0+m_1x+m_2x^2+m_3x^3,
L=l_0+l_1x+l_2x^2+l_3x^3+l_4x^4.                        (22)
```

After substituting `theta_1=a/gamma^2`, the complete polynomial solution of
the Jacobian row `(12)` is

```text
m_0=(9a^3 theta_0^2 gamma^5+3a^3-2gamma)/(9a^3 gamma^3),

l_0=3a m_1/(4gamma)-theta_0(3a^3+2gamma)/(2a gamma),
l_1=3a m_2/(4gamma)
    +(-9a^6+36a^3 theta_0^2 gamma^6-6a^3 gamma-4gamma^2)
      /(12a^3 gamma^3),
l_2=3a theta_0 gamma+3a m_3/(4gamma)+3gamma m_1/2,
l_3=3a^2/(2gamma)+3gamma m_2/2,
l_4=3gamma m_3/2,                                       (23)
```

with `m_1,m_2,m_3` free at this stage.

Let `r` denote the third source-normal coefficient vector.  The `t^2`
Jacobian row is

```text
3 det((A_b,C_b)',r)+2det(v',w)+det(w',v)=0.              (24)
```

The cubic Taylor formula for `P` and `(7)`, `(24)` eliminate `r`:

```text
g_3=-q[2det(v',w)+det(w',v)]/3
    -6A_b N M+2K L-N^3.                                 (25)
```

Substitution of `(23)` and collection in `x` gives

```text
[x^3]g_3=2m_3/(3gamma),
[x^2]g_3=-(3a^2-5gamma^2m_2)/(6gamma^3),
[x]g_3=-(2a theta_0 gamma^2+a m_3-2gamma^2m_1)/(2gamma^3),

[1]g_3=(81a^6-27a^4gamma^2m_2+108a^3theta_0^2gamma^6
         +90a^3gamma-28gamma^2)/(108a^3gamma^5).         (26)
```

Put

```text
delta=[p y]R,                    epsilon=[p^3]R.          (27)
```

Direct expansion of `(1)` gives

```text
[t^3]G=beta*x^3+2alpha*x^2+delta*x+epsilon.              (28)
```

Comparing `(26)` with `(28)` and using `(18)` forces

```text
m_3=3gamma^2 theta_0,
m_2=2(-3a^3-5gamma)/(5a gamma^2),

delta=(2m_1-5a theta_0)/(2gamma),

epsilon=gamma theta_0^2
        +21a^3/(20gamma^5)+4/(3gamma^4)
        -7/(27a^3gamma^3)

       =beta^2/(4gamma)
        +21a^3/(20gamma^5)+4/(3gamma^4)
        -7/(27a^3gamma^3).                               (29)
```

Thus `[p^3]R` is also forced once `beta` is chosen, while `delta=[py]R` is
the next surviving scalar.

### 4.1 The zero-residual cell is empty

There is now an immediate exact contradiction if `R=0`.  That assumption
first gives

```text
beta=[y]R=0,                 alpha=[p^2]R=0.              (29a)
```

The forced value of `alpha` in `(18)` then gives

```text
gamma=-9a^3/10.                                             (29b)
```

Substitute `beta=0` and `(29b)` into the forced value of `epsilon` in
`(29)`.  Exact reduction gives

```text
[p^3]R=epsilon=4000/(6561 a^12) !=0,                     (29c)
```

because `a!=0` and the field has characteristic zero.  This contradicts
`R=0`.  Therefore

```text
no normalized reduced (2,3) Keller pair has R=0.          (29d)
```

The argument uses both normal diagonals.  The first forces the unique
possible resonance `(29b)`; the second excludes it.  It does not exclude a
nonzero coordinated residual.

Using the diagonal formula `(3)` once more,

```text
M=A_2(0)+A_1'(0)x+[s^2]A_0*x^2+[s^3]b*x^3.              (30)
```

Equations `(23)`, `(29)` translate into the Laurent jets

```text
q_0=(a/gamma)s+(3beta/4)s^2+O(s^3),
b=2a s^2+(3gamma beta/2)s^3+O(s^4),

[s^2]A_0=-6a^2/(5gamma^2)-2/(a gamma),
A_2(0)=beta^2/4+1/(3gamma^3)-2/(9a^3gamma^2),
A_1'(0)=gamma*delta+5a beta/(4gamma).                    (31)
```

### 4.2 Independent generic-fibre reduction proof

There is a second proof of `(29d)` which does not use either normal-diagonal
calculation. It is retained because it tests the local algebra against a
global invariant.

Assume `R=0`, put `lambda=3a/(2gamma)`, and use the target value itself as
a parameter:

```text
q=G(A,C)=gamma*u+lambda*p,                 K=k(q).       (29e)
```

Since `A,C` are algebraically independent, `q` is transcendental over
`k`. On the source generic fibre set

```text
s=y/p=xt,                  p=s^2+t,        u=s^2/t.
```

Equation `(29e)` is equivalent to

```text
s^2(q+gamma-lambda*p)=p(q-lambda*p).                    (29f)
```

With

```text
V=s(q+gamma-lambda*p)/lambda,
```

it becomes the smooth cubic

```text
E_S: V^2=p^3-(2q+gamma)p^2/lambda
             +q(q+gamma)p/lambda^2

       =p(lambda*p-q)(lambda*p-q-gamma)/lambda^2.        (29g)
```

Its roots `0,q/lambda,(q+gamma)/lambda` are distinct. Conversely,

```text
s=lambda*V/(q+gamma-lambda*p),    t=p-s^2,    x=s/t,     (29h)
```

so the smooth projective model of `(29g)` has function field

```text
K(E_S)=K(p,V)=k(x,t).
```

The target generic fibre is

```text
E_T: C^2=A^3-(3a^2/4)A+q-a^3/4,                         (29i)
```

and `K(E_T)=k(A,C)`. The Keller condition makes

```text
K(E_T) subset K(E_S)
```

a finite extension of function fields. Thus it gives a finite nonconstant
map between the smooth projective genus-one models. Translating the target
by the image of the source origin turns that map into an isogeny.

The exact invariants are

```text
j(E_T)=-216a^6/[q(2q-a^3)],

j(E_S)=256(q^2+gamma*q+gamma^2)^3
             /[gamma^2 q^2(q+gamma)^2].                 (29j)
```

At the place `q=infinity`, adjoin `q=rho^-6` and put
`A=rho^-2 X`, `C=rho^-3 Y`. Equation `(29i)` becomes

```text
Y^2=X^3+1-(3a^2/4)rho^4 X-(a^3/4)rho^6,
```

whose special fibre `Y^2=X^3+1` is smooth in characteristic zero. Hence
`E_T` has potential good reduction. But `v_infinity(j(E_S))=-2`; every
finite extension multiplies this negative valuation by its ramification
index, whereas good reduction has integral `j`. Hence `E_S` has no
potential good reduction.

The imported reduction fact is the Neron--Ogg--Shafarevich criterion in
Serre--Tate, *Good reduction of abelian varieties*, Ann. Math. 88 (1968),
Section 1: good reduction is equivalent to unramified rational
`ell`-adic Tate module, and potential good reduction to finite inertia.
For an isogeny, choose `ell` away from its degree; the rational Tate modules
are isomorphic. Potential good reduction is therefore isogeny-invariant.
This contradicts the isogeny above and independently proves `R!=0`.

## 5. Complete all-row Hasse repair transform

The preceding calculation used the first source-normal diagonals.  There is
also an exact all-row description of polynomiality in the cusp plane.

Remove the one negative row from `(1)` and write

```text
Delta=G-gamma*u=sum_(n=0)^N d_n(s) tau^n.                (32)
```

For every `m>=0`, define

```text
H_m(Delta)=sum_(n=m)^N binom(n,m)d_n(s)(-s^2)^(n-m).     (33)
```

Since `tau=p-s^2`, this is exactly

```text
H_m(Delta)=[p^m] Delta(s,p-s^2).                         (34)
```

The substitution map

```text
k[p,y] -> k[s,p],                    y -> s p            (35)
```

is injective and has image

```text
direct_sum_(m>=0) p^m * k[s]_(degree<=m).                (36)
```

For the forward direction, a monomial maps as

```text
p^c y^e -> s^e p^(c+e),                 e<=c+e.          (37)
```

For the reverse direction, every basis monomial in the right side of
`(36)` has the unique preimage

```text
s^e p^m <- p^(m-e)y^e,                  0<=e<=m.         (38)
```

Equations `(34)--(38)` prove the iff criterion

```text
Delta in k[p,y]
iff deg_s H_m(Delta)<=m for every m.                     (39)
```

The finite support of `Delta` makes this a finite test.  Moreover,

```text
Delta=lambda*p+R(p,y),       R in (p^2,y)
```

holds if and only if

```text
H_0=0,
H_1=lambda+beta*s for some beta in k,
deg_s H_m<=m for every m>=2.                             (40)
```

Indeed, `(37)` shows that the `m=1` layer consists exactly of `p` and `y`.
All layers `m>=2` reconstruct monomials `p^(m-e)y^e`, which lie in
`(p^2,y)`.  Conversely, reducing `lambda*p+R` modulo `p^2` gives precisely
`p(lambda+beta*s)`.  Thus `(40)` includes both directions and identifies

```text
beta=[y]R=H_1'(s).                                       (41)
```

The transform `(33)` is the exact repair quotient hidden across all positive
Laurent rows.  A scalar row-by-row pole count would not see the degree bound
in `(39)`.

## 6. Hostile controls and audit boundary

1. **Sharp image hostile.**  Adding `c*s^(m+1)p^m` to an arbitrary element
   of `k[s,p]` is still polynomial after `tau=p-s^2`, but `(33)` exposes an
   `H_m` of degree `m+1`; it is not in the image of `k[p,y]`.
2. **Geometry-only hostile.**  The boundary pair `(A_b,C_b)` lies on the
   nodal cubic but has Jacobian zero.  It does not satisfy `(8)` and cannot
   be used to infer `(18)`.
3. **Direct formal-jet hostiles.**  The companion script independently
   solves the next symplectic normal row and directly expands `J` and `P`
   for `(a,gamma,theta_0,m_1)=(1,1,0,0)` and `(2,3,1,5)`.  In both cases
   `J=1+O(t^3)` and the coefficients through `t^3` agree with `(17)`, `(28)`.
   These are controls of the identities, not elements of `B_2` and not
   Keller constructions.
4. **Resonant hostile.**  The forced coefficient `[p^2]R` vanishes when
   `9a^3+10gamma=0`; no argument may divide by it.  Instead, the next
   diagonal evaluates there to `(29c)` when `beta=0` and closes `R=0`.

Exact symbolic reproduction:

```text
python3 04-computation/jc2_reduced_23_hasse_residual_thm3997.py
python3 -O 04-computation/jc2_reduced_23_hasse_residual_thm3997.py
```

The normal and optimized outputs are byte-identical and end with
`ALL THM-3997 EXACT CHECKS PASSED`.

## 7. Honest scope

**Necessary and exact in the THM-3992 gauge:** `(18)`, `(21)`, `(29)`,
`(31)`, and the all-row iff transform `(39)--(40)`.

**Not proved:**

1. that an arbitrary value of `beta` or `delta` extends through every higher
   Laurent row;
2. that the displayed formal source jets lift to actual elements of `B_2`;
3. that `beta` is invariant under arbitrary target automorphisms outside the
   fixed normalization;
4. that the resonant divisor is empty;
5. that an arbitrary nonzero residual is impossible, the depths decrease,
   or the full `(2,3)` cell is empty;
6. `JC(2)`.

The precise gain is a new repair quotient and a much smaller residual
ledger: the pure coefficients `[p^2]R` and `[p^3]R` are forced, `[y]R` is
the first surviving scalar in the fixed gauge, `[py]R` is the next, and the
entire zero-residual subcell is impossible.
