---
id: THM-3253
title: "Affine binomial phases: cyclic residues exclude every linear perturbation of a pure power"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  Let ell be
  a nonconstant algebraic affine polynomial on the coordinate two-simplex.
  For every integer d>=3 and algebraic A,B,C,h with A nonzero, the phase
  `A(ell-h)^d+B(ell-h)+C` has exponential simplex period neither zero nor
  `1/2`.  For B nonzero, the endpoint primitives form a rank-(d-1) cyclic
  connection with residues `-j/d`.  Its spectrum is exactly the critical
  value set of `At^d+Bt`.  A weighted packet splits over endpoint
  exponentials iff every noncritical grouped source vanishes and every
  critical source lies on the one-dimensional derivative line
  `(1,2rho,...,(d-1)rho^(d-2))`.  Three-endpoint moment identities plus a complete
  1+1+1 / 2+1 / 3 collision audit keep the triangle packet off those lines.
  The scalar annihilators of its first two primitive coordinates have
  disjoint exact Frobenius exponent sets, so the two marked copies in the
  boundary period cannot cancel.  A doubled knot is recovered by the same
  coordinate coprimality.  B=0 is THM-3250/3251.  This closes an infinite
  non-pure family in every degree, not arbitrary one-variable polynomials of
  degree at least four, genuinely multivariate phases, or FC(3).
source: codex-2026-08-02-fc3-cubic-frontier
depends_on:
  - THM-3250-fc3-noncollinear-pure-power-turn-current-exclusion
  - THM-3251-fc3-collinear-pure-power-spline-residue-exclusion
related:
  - THM-3203-fc3-complex-affine-coordinate-quadratic-cauchy-green-cycle-exclusion
  - THM-3252-fc3-depressed-cubic-bessel-marked-copy-splitting-gate
external:
  - "Cauchy--Green formula."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc_affine_binomial_cyclic_residue_thm3253.py
output: 05-knowledge/results/fc_affine_binomial_cyclic_residue_thm3253.out
script_sha256: a31951101fd692993551ec2ffc27fe47daffcf19d505aabe5feee640895b103a
output_sha256: 3446b498a7aeb475e9c694e1d1f9eb9a528d93917ae92982a13db4845a39b608
---

# THM-3253 -- every affine binomial phase is excluded

## 1. Statement and novelty

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},          area(Delta)=1/2,
ell in Qbar[u,v] affine and nonconstant,
d>=3 an integer,
A,B,C,h in Qbar,                            A!=0,            (1)

q=A(ell-h)^d+B(ell-h)+C,
K(s)=integral_Delta exp(sq)du dv.
```

Then

```text
K(1)!=0,                         K(1)!=1/2.                 (2)
```

When `B=0`, this is the pure-power result THM-3250/3251.  The new branch is
`A*B!=0`.  It contains every depressed cubic, recovering THM-3252 at `d=3`,
and supplies a genuinely non-pure infinite family in every higher degree.
For `d>=4` it does not allow intermediate terms such as `t^2,...,t^(d-1)`;
those destroy the monomial cyclic connection used below.  Thus (2) is not
all affine-coordinate polynomials and not `FC(3)`.

Put

```text
t=ell-h,                 p(t)=At^d+Bt,                     (3)
r=d-1,
E_x(s)=exp(sp(x)),
J_(j,x)(s)=integral_0^x t^(j-1)exp(sp(t))dt,  1<=j<=r,
Y_x=(J_(1,x),...,J_(r,x))^T.                              (4)
```

All endpoint values `x` are algebraic.  Sections 2--3 identify the exact
rank-`r` object; sections 4--7 retain enough triangle geometry to forbid
every collision and marked-copy cancellation.

## 2. The rank-(d-1) cyclic primitive connection

For `1<=j<r`, Euclidean division gives

```text
t^(j-1)p=(t^j/d)p'+[B(d-1)/d]t^j.                         (5)
```

Integrating the first term by parts yields

```text
dsJ_(j,x)'+jJ_(j,x)-B(d-1)sJ_(j+1,x)=x^jE_x.              (6)
```

For the last row, use

```text
integral_0^x t^(d-1)exp(sp(t))dt
 =[E_x-1]/(dAs)-[B/(dA)]J_(1,x).                          (7)
```

With

```text
kappa=B(d-1)/(dA),                                        (8)
```

the final equation is

```text
dsJ_(r,x)'+rJ_(r,x)+[B^2(d-1)/(dA)]sJ_(1,x)
  =(x^r+kappa)E_x-kappa.                                  (9)
```

Define

```text
a=B(d-1)/d,                 c=-B^2(d-1)/(d^2A),
D=diag(-1/d,-2/d,...,-r/d),                                (10)

Ccyc_(j,j+1)=a,             1<=j<r,
Ccyc_(r,1)=c,
all other entries zero.                                   (11)
```

For the endpoint-source vector defined componentwise by

```text
(v_x)_j=x^j,                   1<=j<r,
(v_x)_r=x^r+kappa,                                          (12)
```

equations (6) and (9) become

```text
Y_x'=(Ccyc+D/s)Y_x+[v_xE_x-kappa e_r]/(ds).               (13)
```

The connection has a regular singularity at zero, pairwise noncongruent
residues `-j/d`, and a single nonzero cycle coupling all residue classes.

## 3. Critical spectrum and exact one-copy splitting

Let `X` be a finite endpoint multiset with algebraic weights `w_x`, and put

```text
Y_w=sum_x w_xY_x.                                         (14)
```

Group endpoints by `eta=p(x)`.  The source on `exp(eta s)` is

```text
b_eta=sum_(x:p(x)=eta) w_xv_x
      -1_(eta=0)(sum_x w_x)kappa e_r.                     (15)
```

Let `mathcal E_X` be the `Qbar(s)`-span of the distinct endpoint
exponentials and `1`.

### 3.1 The spectrum is the critical-value set

If `rho` is critical, then

```text
rho^r=-B/(dA),
eta=p(rho)=[B(d-1)/d]rho=a rho.                            (16)
```

The right and left eigenvectors of `Ccyc` at `eta` are

```text
q_rho=(1,rho,...,rho^(r-1))^T,
l_rho=(1,rho^-1,...,rho^(-(r-1)))^T.                      (17)
```

Indeed `c=a rho^r`, so the wraparound row also agrees.  Moreover

```text
eta^r=a^(r-1)c=-B^d(d-1)^(d-1)/(d^dA).                   (18)
```

The `r` critical values are distinct and satisfy the characteristic
equation (18); hence they are exactly `Spec(Ccyc)`.

### 3.2 Splitting criterion

The following are equivalent:

```text
(i)  Y_w belongs to mathcal E_X^r;

(ii) for every eta,
       b_eta=0,                                      eta noncritical,
       b_eta in L_rho,                               eta=p(rho) critical,

     where
       L_rho=span((1,2rho,...,r rho^(r-1))^T).              (19)
```

To prove necessity, distinct exponentials reduce a rational representation
of `Y_w` to one coefficient vector `R_eta in Qbar(s)^r` per phase class.
It must satisfy

```text
R_eta'=(Ccyc-eta I)R_eta+(D/s)R_eta+b_eta/(ds).            (20)
```

A rational solution has no finite pole.  At a nonzero pole, differentiation
raises the pole order.  At a pole of order `m>=1` at zero, its leading vector
would lie in the kernel of

```text
-mI-D=diag(-m+1/d,...,-m+r/d),                             (21)
```

which is invertible because `1<=j<=r<d`.  Thus `R_eta` is polynomial.

If `eta` is noncritical, its top coefficient is killed by the invertible
matrix `Ccyc-eta I`; downward induction gives `R_eta=0` and `b_eta=0`.
If `eta=p(rho)` is critical and the polynomial has positive degree `N`, its
top coefficient is a multiple of `q_rho`.  Solvability of the next equation
would require

```text
l_rho^T(NI-D)q_rho=0.                                     (22)
```

But every coordinate product of the vectors in (17) is one, so

```text
l_rho^T(NI-D)q_rho
 =sum_(j=1)^r(N+j/d)=r(N+1/2)!=0.                         (23)
```

There is no positive-degree exception.  A constant coefficient must lie in
`ker(Ccyc-eta I)` and (20) gives

```text
b_eta=-dDq_rho
       in span((1,2rho,...,r rho^(r-1))^T)=L_rho.          (24)
```

Conversely, every source in (19) has this constant rational particular
solution.  After subtracting all such solutions, the remainder is an entire
homogeneous solution.  Its first nonzero Taylor coefficient `h_n` would obey
`(nI-D)h_n=0`, impossible because every diagonal entry is `n+j/d`.  The
remainder is zero, proving sufficiency.

The same proof permits adjoining finitely many algebraic exponentials with
zero source: equation (20) then forces their rational coefficients to vanish.

## 4. The three-endpoint critical-fibre audit

For both three-distinct triangle geometries there are algebraic nonzero
weights `omega_j`, endpoints `x_j`, and a nonzero algebraic factor `Omega`
such that

```text
sum_j omega_j=0,
sum_j omega_jx_j=0,
sum_j omega_jx_j^2=-Omega!=0.                             (25)
```

For a noncollinear image, these are the vertex-turn weights and `Omega=W`
from THM-3250.  Every turn weight is nonzero because equality of two adjacent
edge slopes would make the image triangle collinear.  For three distinct
collinear knots, these are the B-spline weights and `Omega=1` from THM-3251;
their explicit divided-difference denominators are nonzero.

Put

```text
U=sum_j omega_jY_(x_j).                                   (26)
```

The constant correction in (15) vanishes.  Summing all grouped sources and
using (25), their first two coordinates are

```text
(0,-Omega).                                                (27)
```

Thus some grouped source is nonzero.  We now audit every partition of three
endpoints by phase value.

### 4.1 Singleton phase classes

A nonzero singleton source at a noncritical value violates (19).  Suppose
instead `p(x)=p(rho)` is critical.  For `d>=4`, if `v_x` lay on `L_rho`,
comparison of its first two coordinates would force

```text
x=2rho.                                                     (28)
```

Normalize `x=rho y`.  The critical fibre equation is

```text
Phi_d(y)=y^d-dy+d-1=0,                                    (29)
```

where `y=1` is the double critical root.  But

```text
Phi_d(2)=2^d-d-1!=0,                       d>=3.            (30)
```

Thus every critical singleton with `d>=4` is transverse.  For the boundary
`d=3`, the two distinct points of the critical fibre are
`rho,-2rho`; using `kappa=-2rho^2`, both endpoint sources are multiples of
`(1,-rho)`, whereas `L_rho=span((1,2rho))`.  They are transverse as well.

### 4.2 The complete partition audit

There are only three partition types.

1. In type `1+1+1`, every class is a nonsplitting singleton.
2. In type `2+1`, the singleton has nonzero weight and is nonsplitting,
   regardless of what happens inside the two-endpoint collision.
3. In type `3`, there is one grouped source.  It is nonzero by (27).  If its
   value is noncritical, (19) fails.  If critical, membership in `L_rho`
   and its zero first coordinate would force its line coefficient to be zero,
   contradicting the nonzero second coordinate in (27).

Hence

```text
U notin mathcal E_X^r                                      (31)
```

for every three-distinct noncollinear or collinear image.  The argument does
not assume distinct endpoint phase values and explicitly includes all roots
of a critical fibre, not only the double root.

## 5. The first two scalar annihilators are coprime

Modulo the derivative-stable exponential space, every packet obeys the
homogeneous connection

```text
Y'=(Ccyc+D/s)Y.                                            (32)
```

Because the residue eigenvalues are nonresonant, there is for each
`1<=k<=r` a unique Frobenius solution

```text
H_k=s^(-k/d)(e_k+O(s)).                                    (33)
```

The recurrence for its coefficient `h_n` is

```text
[(-k/d+n)I-D]h_n=Ccyc h_(n-1).                            (34)
```

The monomial cycle sends coordinate `k` to `k-1`, with `1` wrapping to `r`.
All its edges and all denominators in (34) are nonzero.  Therefore coordinate
`j` of `H_k` first occurs at exact exponent

```text
alpha_(j,k)=-k/d+delta_j(k),
delta_j(k)=(k-j mod r) in {0,...,r-1}.                    (35)
```

For fixed `j`, the `r` values have distinct fractional parts, so the
successive system rows recursively recover all coordinates from coordinate
`j`; equivalently that coordinate row is cyclic.  It has a scalar annihilator
`L_j` of order `r` with exact Frobenius exponent set

```text
S_j={alpha_(j,k):1<=k<=r}.                                 (36)
```

The first two sets are disjoint.  If
`alpha_(1,k)=alpha_(2,l)`, reduction modulo integers gives `k=l`; then
`delta_1(k)=delta_2(k)`, impossible because the two targets are distinct
vertices of one directed cycle.  Any common local solution decomposes into
the distinct monodromy characters indexed by `k`; a nonzero character
component would require exactly such an equality.  Hence `L_1` and `L_2`
have no common local solution and no common right factor over
`Qbar(s)[d/ds]`.  The Euclidean property of this Ore ring gives operators
`P,Q` with

```text
P L_1+Q L_2=1.                                             (37)
```

Consequently any class annihilated by both `L_1` and `L_2` is zero.  This is
the all-degree replacement for the order-`1/3` versus order-`2/3` Bessel
subtraction in THM-3252.  It also answers the common-higher-factor hostile:
the complete rank-`r` coupling introduces no shared scalar factor.

## 6. Three distinct endpoints: marked copies cannot cancel

Alongside `U` from (26), define

```text
V=sum_j omega_jx_jY_(x_j),
F=U_2-V_1.                                                  (38)
```

The Cauchy--Green formula in the noncollinear case and the B-spline formula
in the collinear case both read

```text
Omega K(s)=exp(Cs)F(s).                                    (39)
```

Suppose `F` belonged to the endpoint exponential space.  In the quotient,

```text
u_2=v_1=:f.                                                 (40)
```

As a second coordinate, `f` is killed by `L_2`; as a first coordinate it is
killed by `L_1`.  Equation (37) forces `f=0`.  Thus `u_2=0`.
Rows `2,3,...,r` of the cyclic homogeneous system successively force

```text
u_3=...=u_r=u_1=0,                                         (41)
```

because every cycle edge is nonzero.  This contradicts (31).  Therefore

```text
F notin mathcal E_X.                                       (42)
```

The boundary mixes two marked packets, but it asks them to identify scalar
coordinates with provably coprime annihilators.

## 7. One doubled knot

Suppose two vertex values are `a` and the remaining value is `b!=a`.  Put

```text
L=b-a,
T=Y_a-Y_b,
G=-bT_1+T_2.                                               (43)
```

The triangular marginal identity is

```text
L^2K(s)=exp(Cs)G(s).                                       (44)
```

First, `T` cannot split.  Its constant correction cancels.  If
`p(a)!=p(b)`, both phase classes are singletons and section 4.1 applies.  If
the common value is noncritical, the grouped source `v_a-v_b` has nonzero
first coordinate `a-b`, violating (19).

Suppose the common value is critical.  For `d=3`, the distinct points of the
normalized fibre (29) are `1,-2`, whose sum is `-1`, not the splitting value
`2`.  For `d>=4`, if `v_a-v_b` lay on `L_rho`, its first three coordinates
would give

```text
a+b=2rho,
a^2+ab+b^2=3rho^2.                                        (45)
```

Hence `ab=rho^2`, so `a,b` are both roots of
`z^2-2rho z+rho^2=(z-rho)^2`, contrary to `a!=b`.  The `kappa` in the last
coordinate cancels when `d=4`, so the third identity remains exactly (45).
Thus

```text
T notin mathcal E_{a,b}^r.                                (46)
```

If `G` were endpoint-elementary and `b!=0`, the quotient relation
`t_2=bt_1` would make a nonzero scalar function belong to both coordinate
solution spaces, contradicting (37).  Hence `t_1=t_2=0`, and the cyclic rows
force every component to vanish.  If `b=0`, the relation already says
`t_2=0`; rows `2` through `r` again force `t_3,...,t_r,t_1=0`.  Either case
contradicts (46), so

```text
G notin mathcal E_{a,b}.                                   (47)
```

This treats knot multiplicity uniformly in every degree.

## 8. E-functions, Beukers, and the two values

For `1<=j<=r`,

```text
J_(j,x)(s)=sum_(n>=0) s^n/n!
              integral_0^x t^(j-1)(At^d+Bt)^n dt.         (48)
```

The coefficients are algebraic with exponential conjugate and denominator
growth, and (6),(9) give a holonomic system.  Thus all primitives, weighted
packets, and algebraic endpoint exponentials are E-functions.  After
augmenting (13) by those exponentials and `1`, the system has denominator
only `s`, so `s=1` is ordinary.

By the zero-source extension following (24), adjoin `exp(-Cs)` to the
exponential spaces in (42) and (47) without changing either nonsplitting
conclusion.  Distinct exponentials are independent over `Qbar(s)`.  Beukers
Corollary 1.4 therefore transfers the functional independence to algebraic
linear independence of values at `s=1`.

In a three-distinct geometry, (39) says that `K(1)=0` or `1/2` would give

```text
F(1)=0,
F(1)-(Omega/2)exp(-C)=0,                                  (49)
```

respectively.  In the doubled-knot geometry, (44) would give

```text
G(1)=0,
G(1)-(L^2/2)exp(-C)=0.                                    (50)
```

All four relations are impossible.  This proves (2) for `B!=0`; the pure
case `B=0` is THM-3250/3251.

## 9. Exact boundary and next target

The proof closes exactly the translate-linear perturbations

```text
A(ell-h)^d+B(ell-h)+C.                                    (51)
```

Its preserved structure is the monomial cycle in (11).  Adding an
intermediate term changes the Euclidean remainders from one cycle edge to a
denser connection.  The residues `-j/d` survive after a suitable basis
choice, but three properties used here require a new proof:

1. the constant spectrum need not be the critical-value set through one
   monomial companion;
2. critical splitting sources need not have the derivative line (19); and
3. the first two scalar coordinate annihilators may acquire a common factor.

The next honest affine-coordinate target is therefore a trinomial with one
intermediate exponent, together with an exact common-factor and critical-
fibre audit.  Merely repeating the moment total is insufficient once a
two-endpoint critical collision has no nonsplitting singleton.

## 10. Reproduction contract

Run

```bash
python3 04-computation/fc_affine_binomial_cyclic_residue_thm3253.py
python3 -O 04-computation/fc_affine_binomial_cyclic_residue_thm3253.py
```

The exact verifier checks:

1. 120 Taylor coefficients of the rank-`d-1` system for `3<=d<=7`;
2. critical eigenvectors and the universal `N+1/2` obstruction through
   `d=9`;
3. critical/noncritical polynomial splitting ranks through degree three for
   `d=4,5,6`;
4. every `1+1+1`, `2+1`, and `3` critical-fibre partition mechanism in
   degrees `4,5,6`, plus the doubled-pair diagonal obstruction;
5. disjoint first/second Frobenius exponent sets for `3<=d<=10` and exact
   cyclic observability determinants for `d=4,5,6`; and
6. 45 direct Taylor-coefficient comparisons of the three period formulas
   for non-pure phases of degrees `4,5,6`.

Normal and optimized runs must be byte-identical to the archived output.
