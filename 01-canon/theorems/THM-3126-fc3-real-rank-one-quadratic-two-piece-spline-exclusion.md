---
id: THM-3126
title: "FC(3) real rank-one quadratics: the two-piece spline source cannot cancel"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  Let r be a
  real-valued affine polynomial on the coordinate two-simplex with algebraic
  coefficients, and let Q in Qbar[t] have degree at most two.  Then
  int_Delta exp(Q(r)) dA is nonzero, and it equals 1/2 iff Q(r) is identically
  zero.  Equivalently this covers an algebraic affine coordinate whose three
  vertex images lie on one real affine line, after an algebraic
  reparameterization of that line.  When the vertex values a<b<c are
  distinct, r_*dA has an exact two-piece linear B-spline density.  Removing
  explicit endpoint exponentials leaves one E-function whose first-order
  source vector is (d_a^2/(2Ap), -(p+q)d_b^2/(2Apq), d_c^2/(2Aq)).  No source
  vanishes after grouping any equal quadratic endpoint values; functional
  independence and Beukers Corollary 1.4 at z=1 exclude values 0 and 1/2.
  Repeated vertex values and affine Q reduce to THM-3116.  This does NOT cover
  affine coordinates with noncollinear complex vertex images, genuinely
  bivariate quadratic phases, nonflat leading forms, or FC(3) in general.
source: codex-2026-08-02-fc3-simplex
depends_on:
  - THM-3116
related:
  - THM-3039
external:
  - "Riesz--Markov--Kakutani representation theorem."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc3_rank_one_quadratic_spline_thm3126.py
output: 05-knowledge/results/fc3_rank_one_quadratic_spline_thm3126.out
---

# THM-3126 — real rank-one quadratic phases on the triangle

## 1. Statement and exact scope

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},
lambda_0=1-u-v, lambda_1=u, lambda_2=v,
mu_Delta=coordinate area on Delta,     mu_Delta(Delta)=1/2.     (1)
```

Suppose `r in (Qbar intersect R)[u,v]` is affine and
`Q in Qbar[t]` has degree at most two.  Then

```text
K(1):=int_Delta exp(Q(r(u,v))) du dv !=0,                       (2)

K(1)=1/2    iff    Q(r(u,v))=0 identically on Delta.            (3)
```

This is invariant under algebraic reparameterization of a real affine line.
Thus (2)--(3) also hold for `P(ell)`, where `ell` is algebraic affine and

```text
ell=x+eta r,       x,eta in Qbar, eta!=0,                       (4)
```

with `r` real-valued, and `P in Qbar[t]`, `deg P<=2`: replace `Q(t)` by
`P(x+eta t)`.  Equivalently, the three barycentric vertex values of `ell`
must lie on one real affine line with an algebraic parameter.  The theorem
does not order noncollinear complex vertex images and makes no claim there.

If `r` is constant, (2)--(3) are Hermite--Lindemann for
`exp(Q(r))/2`.  If exactly two vertex values agree, algebraically normalize
the other two values to `0,1`; then `Q(r)` is a degree-at-most-two polynomial
in one barycentric coordinate and THM-3116 applies.  If `Q` is affine, then
`Q(r)` is affine on `Delta`, also covered by THM-3116.  It remains to prove
the case

```text
r(vertex set)={a,b,c},       a<b<c in Qbar intersect R,
Q(t)=A t^2+B t+C,            A!=0.                              (5)
```

The result closes the full **real rank-one quadratic** branch of the flat
simplex exponential obstruction.  It does not address a genuinely
bivariate quadratic, a nonflat leading form in the Factorial Conjecture, or
FC(3) itself.

## 2. RMK with typed hypotheses: the two-piece spline

After permuting the barycentric vertices, write

```text
r=a lambda_0+b lambda_1+c lambda_2,
p=b-a>0, q=c-b>0, L=c-a=p+q.                                  (6)
```

Riesz--Markov--Kakutani represents the positive functional
`f -> int_Delta f dA` by `mu_Delta`.  Because `r` is real-valued, its
pushforward `nu=r_*mu_Delta` is a positive Borel measure on the ordered
interval `[a,c]`.  Its density is exactly

```text
rho(t)= (t-a)/(pL),       a<=t<=b,
        (c-t)/(qL),       b<=t<=c.                              (7)
```

Indeed, for `a<=t<=b`, the sublevel set `{r<=t}` is the corner triangle at
`a`, with edge ratios `(t-a)/p` and `(t-a)/L`; its coordinate area is
`(t-a)^2/(2pL)`.  For `b<=t<=c`, the complementary corner at `c` has area
`(c-t)^2/(2qL)`.  Differentiating gives (7), and
`int_a^c rho(t)dt=1/2`.

Consequently, for

```text
K(z)=int_Delta exp(zQ(r))dA,
H_ab(z)=int_a^b exp(zQ(t))dt,
H_bc(z)=int_b^c exp(zQ(t))dt,                                 (8)
```

one has the exact B-spline transform

```text
K(z)=1/L [ 1/p int_a^b (t-a)exp(zQ(t))dt
          +1/q int_b^c (c-t)exp(zQ(t))dt ].                    (9)
```

This is the complete RMK input.  Positivity and mass are preserved, but the
pushforward forgets transverse level-set geometry.  The noncancellation
below uses the additional algebraic `E`-function sidecar; it is not a
consequence of positivity.

## 3. Endpoint correction and the common first-order operator

Put

```text
h=-B/(2A),                  d_x=Q'(x)=2Ax+B,
E_x(z)=exp(zQ(x)).                                               (10)
```

On either segment `[x,y]`, integration of
`d(exp(zQ(t)))/dt` gives

```text
2Az int_x^y t exp(zQ(t))dt+Bz H_xy(z)=E_y(z)-E_x(z).            (11)
```

Substitution into (9) yields

```text
L K(z)=alpha H_ab(z)+beta H_bc(z)+U(z),                         (12)

alpha=(h-a)/p,            beta=(c-h)/q,
U(z)=1/(2Az)[(E_b-E_a)/p-(E_c-E_b)/q].                         (13)
```

Although (13) displays `1/z`, its numerator vanishes at zero.  Define the
entire corrected period directly by its segment expression

```text
M(z):=alpha H_ab(z)+beta H_bc(z)=L K(z)-U(z).                   (14)
```

Each `H_xy`, hence `M`, is an `E`-function.  Its coefficient of `z^m/m!` is
an integral over algebraic endpoints of the algebraic polynomial `Q^m`.
All conjugates grow exponentially; after clearing fixed coefficient and
endpoint denominators, the remaining common denominators divide a fixed
exponential factor times `lcm(1,...,2m+1)`.  Thus their logarithmic heights
are `O(m)`.  The differential system below supplies holonomicity.

The common segment operator is

```text
D_Q=4Az d/dz+2A+(B^2-4AC)z.                                   (15)
```

Differentiating `(2At+B)exp(zQ(t))` and integrating gives the exact identity

```text
D_Q H_xy=d_y E_y-d_x E_x.                                     (16)
```

Therefore

```text
D_Q M=S_a E_a+S_b E_b+S_c E_c,                                (17)

S_a=d_a^2/(2Ap),
S_b=-L d_b^2/(2Apq),
S_c=d_c^2/(2Aq).                                               (18)
```

For example, `S_a=-alpha d_a`; using
`alpha=1-d_b/(2Ap)` and `beta=1+d_b/(2Aq)` gives every square in (18).
The alternating sign is the second-difference signature of the linear
B-spline, while the squares retain the quadratic critical-point data.

## 4. The knot-path quotient and every endpoint collision

Equation (16) has an intrinsic path interpretation.  The ordered knots
`a-b-c` are vertices, `H_ab,H_bc` are its two oriented edges, and `D_Q`
sends an edge to its signed endpoint current `d_yE_y-d_xE_x`.  The spline
weights turn that incidence vector into (18).  This is a path, not a
tournament: introducing orientations between nonadjacent knots would add no
observable and would lose the level-set order.

Functional relations group vertices whose exponentials agree.  For a
quadratic on three distinct knots, at most one pair of endpoint values can
agree.  The three possible reflection collisions and their grouped sources
are

```text
Q(a)=Q(b):  B=-A(a+b),   S_a+S_b=-A p^2/(2q),
Q(b)=Q(c):  B=-A(b+c),   S_b+S_c=-A q^2/(2p),
Q(a)=Q(c):  B=-A(a+c),   S_a+S_c= A L^3/(2pq).                 (19)
```

Every quantity on the right is nonzero.  With no collision, some `d_x` is
nonzero and its ungrouped source in (18) is nonzero.  Thus after quotienting
the knot path by **any** equality among `E_a,E_b,E_c` (including equality
with the constant function when `Q(x)=0`), at least one exponent group has a
nonzero source.

This identifies the missing sidecar in the exponential quotient: the value
`Q(x)` determines which vertices merge, but the squared derivative `d_x^2`
in (18) prevents their currents from cancelling.

## 5. Functional independence and the two forbidden values

Take the distinct functions among

```text
M(z), E_a(z), E_b(z), E_c(z), 1.                               (20)
```

They form a homogeneous first-order rational system using (17) and
`E_x'=Q(x)E_x`; its common denominator is `T(z)=z`.  We now prove that they
are linearly independent over `Qbar(z)`.

Distinct algebraic exponentials are already independent over `Qbar(z)`.
If a relation involving `M` existed, solve it for `M` as a rational linear
combination of the distinct exponentials and substitute into (17).  For an
exponent group of value `xi`, its rational coefficient `R` would satisfy

```text
4AzR'+[2A+(B^2-4AC+4Axi)z]R=T_xi,                             (21)
```

where `T_xi` is the grouped source.  Choose the nonzero group supplied by
section 4.  If it contains a knot `x`, then

```text
B^2-4AC+4Axi=(2Ax+B)^2=d_x^2!=0.                              (22)
```

After division by `4A`, (21) is

```text
zR'+(1/2+kappa z)R=tau,       kappa!=0, tau!=0.                (23)
```

There is no rational solution.  A pole of order `n` away from zero gives an
uncancellable pole of order `n+1` in `zR'`.  At zero the leading coefficient
is `(-n+1/2)r_(-n)`, never zero for integral `n>=1`.  Hence `R` is a
polynomial; then `kappa zR` raises its degree, contradicting the nonzero
constant right side.  This proves the claimed functional independence.

Beukers, *A refined version of the Siegel--Shidlovskii theorem*, Corollary
1.4, applies to (20) at the algebraic point `z=1`, since
`1*T(1)!=0`.  Therefore their values at one are linearly independent over
`Qbar`.

But (12)--(14) give

```text
K(1)=0     => M(1)+U(1)=0,
K(1)=1/2   => M(1)+U(1)-L/2=0.                                (24)
```

Here `U(1)` is an algebraic linear combination of the endpoint exponential
values.  Each line of (24) is a nontrivial algebraic linear relation because
the coefficient of `M(1)` is one, even after all endpoint collisions.
Corollary 1.4 excludes both.  Together with the reductions following (5),
this proves (2)--(3).

For the flat-top FC(3) coefficient of THM-3116, the consequence is precise:
an algebraically specialized first projective layer that is a quadratic in
any collinear real affine coordinate cannot be the required zero period.
This is a branch exclusion, not an FC(3) proof.

## 6. Sharp hostile and failure boundary

The algebraic coefficient hypothesis cannot be dropped.  Take a real affine
coordinate with vertex values `0,1,2`.  Its area-measure pushforward is the
tent density `t/2` on `[0,1]` and `(2-t)/2` on `[1,2]`.  Hence

```text
int_Delta exp(g r)dA=(exp(g)-1)^2/(2g^2).                      (25)
```

At `g=2 pi i`, (25) is exactly zero.  This is a three-distinct-knot
cancellation witness with a positive representing measure; its coefficient
is necessarily transcendental by Lindemann--Weierstrass.  It demonstrates
both sharp boundaries:

1. RMK positivity does not prevent a Fourier--Laplace transform from
   vanishing;
2. the arithmetic `E`-function independence, not positivity, excludes the
   algebraic rank-one quadratic branch.

The hostile coordinate `r` itself is still real-valued, collinear, and
algebraic.  What lies outside the theorem is precisely the phase coefficient:
`g=2 pi i` is transcendental, so `Q(t)=gt` is not in `Qbar[t]`.  Thus (25)
challenges an extension to unrestricted complex **coefficients**, not the
real/collinear direction hypothesis proved here.

If the three complex vertex images of an affine coordinate are noncollinear,
there is no real interval order, no two-piece density (7), and no knot path
to which (18) applies.  The pushforward then lives on a genuine triangle in
the complex plane.  That is a different, still-open rank-one-as-an-algebraic-
expression geometry and is not silently included here.

## 7. Exact verification

Run

```bash
python3 04-computation/fc3_rank_one_quadratic_spline_thm3126.py
python3 -O 04-computation/fc3_rank_one_quadratic_spline_thm3126.py
```

The frozen controls verify:

* the spline mass and seven direct simplex/pushforward moments;
* fourteen coefficient recurrences for the segment ODE (16);
* the three signed-square sources (18) and the identity (22);
* every pair-collision source in (19);
* the half-integral pole obstruction in (23); and
* the exact hostile transform (25) and its zero at `2 pi i`.

The normal and optimized executions are byte-identical.  QED in the stated
real/collinear rank-one scope.

```text
source sha256 = b7157cff1e6b4278cf28e4ce0f37315c46100060714f8332eec9450c15a96d7c
output sha256 = 1b9af0130c443faa04f5f43c470e3d6b0f2ea6175dc185dc792e7a1468ff2856
```
