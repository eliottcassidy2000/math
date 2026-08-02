---
id: THM-3143
title: "FC(3) complex affine-coordinate quadratics: Cauchy--Green cycle sources cannot cancel"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let ell be an algebraic
  affine polynomial on the coordinate two-simplex and let P in Qbar[t] have
  degree at most two.  Then int_Delta exp(P(ell)) dA is nonzero, and it is
  1/2 iff P(ell) is identically zero.  For noncollinear complex vertex
  images, Cauchy--Green converts the area period into a weighted cycle of
  three edge E-functions.  Its common quadratic differential operator has
  vertex source S_j=P'(z_j)(b_(j-1)-b_j).  Locally this is the derivative
  square times the turn of the two adjacent unoriented edge directions.
  After grouping every equality P(z_j)=P(z_k), a reflection pair cancels
  only its shared edge direction and leaves two distinct outer directions;
  hence a nonzero source group survives.  Functional independence and
  Beukers Corollary 1.4 at s=1 exclude the values 0 and 1/2.  The collinear
  boundary is THM-3142 and the affine phase boundary is THM-3116.  This is
  the full degree-at-most-two affine-coordinate (univariate-aligned) branch,
  not the larger class of quadratics with rank-one Hessian and a transverse
  linear term, and not FC(3).
source: codex-2026-08-02-fc3-simplex
depends_on:
  - THM-3116
  - THM-3142
related:
  - THM-3039
external:
  - "Cauchy--Green formula."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
output: 05-knowledge/results/fc3_complex_affine_coordinate_cauchy_green_thm3143.out
---

# THM-3143 — complex affine-coordinate quadratic phases on the triangle

## 1. Statement and scope guard

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},       area(Delta)=1/2.       (1)
```

For every affine `ell in Qbar[u,v]` and every `P in Qbar[t]` with
`deg P<=2`, put

```text
K(s)=int_Delta exp(s P(ell(u,v))) du dv.                       (2)
```

Then

```text
K(1)!=0,                                                       (3)
K(1)=1/2       iff       P(ell)=0 identically.                 (4)
```

The exact factorization `q=P(ell)` is essential.  We call this the
**affine-coordinate** or **univariate-aligned** quadratic class.  A
quadratic may have Hessian rank one and still contain a transverse linear
term, in which case it need not factor through one affine coordinate and is
not covered here.  Nor does the theorem address a nonflat leading form in
the Factorial Conjecture or prove FC(3).

If the three vertex values of `ell` are collinear, algebraically
reparameterize their real affine line and apply THM-3142.  If `P` is
affine, `P(ell)` is affine and THM-3116 applies.  Constant cases reduce to
Hermite--Lindemann.  It therefore remains to prove (3)--(4) when

```text
z_j=ell(vertex j),          z_0,z_1,z_2 noncollinear,
P(t)=A t^2+B t+C,           A!=0.                              (5)
```

Relabel the vertices, if necessary, so that `(z_0,z_1,z_2)` is
counterclockwise.  This only permutes barycentric coordinates and does not
change (2).

## 2. Cauchy--Green retains the missing two-dimensional geometry

Write

```text
e_1=z_1-z_0,       e_2=z_2-z_0,
J=Im(conjugate(e_1)e_2)>0,       W=2 i J!=0.                   (6)
```

The affine map from `Delta` to the oriented triangle `T` with vertices
`z_0,z_1,z_2` has real Jacobian `J`.  Since

```text
d_bar(conjugate(t) exp(sP(t)))=exp(sP(t)),                     (7)
```

Cauchy--Green gives the exact area-to-boundary identity

```text
W K(s)=integral_(boundary T) conjugate(t) exp(sP(t)) dt.        (8)
```

This is the noncollinear replacement for THM-3142's one-dimensional spline
pushforward.  It retains all three edge directions, which will be the
sidecar preventing cancellation.

Center the quadratic at

```text
c=-B/(2A),             w=t-c,             P'(t)=2Aw.           (9)
```

Index cyclically modulo three and orient edge `j` from `z_j` to `z_(j+1)`.
Set

```text
m_j=conjugate(z_(j+1)-z_j)/(z_(j+1)-z_j),
w_j=z_j-c,
b_j=conjugate(w_j)-m_j w_j,
H_j(s)=integral_(z_j)^(z_(j+1)) exp(sP(t))dt,
E_j(s)=exp(sP(z_j)).                                           (10)
```

Along edge `j`, the conjugate-line equation is

```text
conjugate(w)=m_j w+b_j.                                       (11)
```

The term `conjugate(c)` integrates to zero around the closed cycle because
`exp(sP(t))` is entire.  Integrating the `m_j w` term by endpoints therefore
turns (8) into

```text
W K(s)=M(s)+U(s),                                              (12)
M(s)=sum_j b_j H_j(s),
U(s)=1/(2As) sum_j (m_(j-1)-m_j)E_j(s).                        (13)
```

The displayed numerator of `U` vanishes at `s=0`; (12), or the edge
integrals themselves, gives its removable continuation there.

There is an important quotient loss to record before using (13):

```text
H_0(s)+H_1(s)+H_2(s)=0.                                       (14)
```

Indeed this is the integral of an entire function around a closed triangle.
We never claim that the three edge periods are independent.  Only the one
cycle-weighted period `M` will be compared with the endpoint exponentials.

## 3. One common operator and the vertex-turn source

Put

```text
d_j=P'(z_j)=2Aw_j,
D_P=4As d/ds+2A+(B^2-4AC)s.                                  (15)
```

The polynomial identity

```text
P'(t)^2=4A P(t)+B^2-4AC                                      (16)
```

and differentiation of `P'(t)exp(sP(t))` give, on every edge,

```text
D_P H_j=d_(j+1)E_(j+1)-d_jE_j.                               (17)
```

Consequently

```text
D_P M=sum_j S_jE_j,             S_j=d_j(b_(j-1)-b_j).         (18)
```

Both edge constants in `S_j` must be evaluated at the same vertex.  From
(11), on the incoming and outgoing edges respectively,

```text
b_(j-1)=conjugate(w_j)-m_(j-1)w_j,
b_j    =conjugate(w_j)-m_jw_j.                                (19)
```

Thus, and only after this local subtraction,

```text
S_j=2A w_j^2(m_j-m_(j-1))
   =d_j^2(m_j-m_(j-1))/(2A).                                 (20)
```

The factor `m_j-m_(j-1)` is nonzero.  Equality of
`conjugate(e)/e` for two nonzero complex edge vectors says that their
directions agree modulo `pi`; adjacent sides would then be collinear, which
contradicts (5).  Formula (20) is therefore a derivative-square decorated
turning current on the three-cycle, not an assertion that the individual
`H_j` are independent.

## 4. Every exponential quotient retains a source

Functional relations identify vertices having the same endpoint
exponential.  For a quadratic, the complete collision law is

```text
P(z_j)=P(z_k), j!=k       iff       w_k=-w_j.                  (21)
```

There are two cases.

* If no endpoint values collide, at most one `w_j` is zero.  Every other
  singleton source in (20) is nonzero.
* Suppose `P(z_j)=P(z_k)`.  Then `w_k=-w_j!=0`, so the grouped source is

```text
S_j+S_k=2A w_j^2[(m_j-m_(j-1))+(m_k-m_(k-1))].                (22)
```

  In a three-cycle the shared edge direction occurs once with each sign and
  cancels.  What remains is the difference of the two outer-edge `m` values.
  Those directions cannot agree unless the triangle is degenerate.  Hence
  (22) is nonzero.

Thus, after grouping **every** equality among `E_0,E_1,E_2` and also merging
with the constant function when `P(z_j)=0`, there is a value `xi` whose
grouped source `T_xi` is nonzero.  For any endpoint in that group,

```text
B^2-4AC+4A xi=P'(z_j)^2=d_j^2!=0.                             (23)
```

The edge-direction cycle is precisely the missing coordinate in the
endpoint-value quotient: a reflection collision forgets the shared edge,
but cannot forget both outer edges.

## 5. Functional independence of the corrected cycle period

Each `H_j`, and hence `M`, is an E-function.  Indeed the coefficient of
`s^n/n!` in `H_j` is the integral of the algebraic polynomial `P(t)^n`
between algebraic endpoints.  Its conjugates have exponential growth, its
denominators have exponential growth after the factors up to `2n+1` are
cleared, and (17) supplies holonomicity.  The functions `E_j` and `1` are
elementary E-functions.

Take the distinct functions among the endpoint exponentials and `1`, and
denote them by `E_xi(s)=exp(xi s)`.  Distinct exponentials are linearly
independent over `Qbar(s)`.  We claim that adjoining `M` preserves this
independence.

Otherwise solve a rational functional relation for

```text
M=sum_xi R_xi(s)E_xi(s),          R_xi in Qbar(s).             (24)
```

Apply `D_P` and compare distinct exponentials.  For the nonzero source group
from section 4, `R=R_xi` would satisfy

```text
4AsR'+[2A+(B^2-4AC+4Axi)s]R=T_xi.                            (25)
```

By (23), division by `4A` gives

```text
sR'+(1/2+kappa s)R=tau,       kappa!=0, tau!=0.               (26)
```

No rational `R` solves (26).  At a finite nonzero pole, the derivative term
has one higher pole order and cannot cancel.  At zero, a pole of integral
order `n>=1` has leading multiplier `-n+1/2`, never zero.  Thus `R` is a
polynomial.  Its top term then produces an uncancelled term of one higher
degree through `kappa sR`; `R=0` contradicts `tau!=0`.  This proves the
functional independence.

The vector consisting of `M` and these distinct exponentials satisfies a
first-order system over `Qbar(s)` by (18), with common denominator

```text
T(s)=s.                                                        (27)
```

Since `1*T(1)!=0`, Beukers Corollary 1.4 applies at `s=1`: their values are
linearly independent over `Qbar`.

## 6. The two forbidden values

At `s=1`, equations (12)--(13) say

```text
W K(1)=M(1)+1/(2A) sum_j(m_(j-1)-m_j)E_j(1).                  (28)
```

All coefficients in (28) are algebraic and `W!=0`.

If `K(1)=0`, (28) is a nontrivial algebraic linear relation among `M(1)`
and the distinct endpoint exponentials.  If `K(1)=1/2`, subtracting
`W/2` gives such a relation after adjoining the constant value `1`.
Both contradict section 5.  This proves (3) and excludes `1/2` in the
noncollinear quadratic case.

If `P(ell)=0` identically, then (2) gives `K(1)=area(Delta)=1/2` directly.
The reductions in section 1 prove the converse and nonvanishing in every
remaining boundary case.  This completes (3)--(4).

## 7. Exact verification and failure boundary

Run

```text
python3 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
python3 -O 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
```

The frozen exact controls normalize the triangle to vertices
`0,1,x+iy`, `y>0`, and verify:

* Cauchy--Green orientation on seven holomorphic monomials;
* all conjugate-line and endpoint-correction identities;
* six coefficient instances of the closed-cycle relation (14);
* six polynomial instances of the common ODE (17);
* the locally derived source (20) at all three vertices;
* all three reflection collisions, including the surviving outer-edge
  source in (22);
* the half-integral pole obstruction in (26); and
* the exact collinear boundary `y=0`, where `W` and all three turns vanish
  together and THM-3142's spline mechanism takes over.

The normal and optimized runs are byte-identical.

```text
source sha256 = e2683b37419f6bdc32e872440db1bbbc88287e68e5c2ebc3d0ed38ad2ddbd718
output sha256 = 9a2477aa33e7d203223721be920da2933505c4ac679c317d338efe7f0e67d7ae
```
