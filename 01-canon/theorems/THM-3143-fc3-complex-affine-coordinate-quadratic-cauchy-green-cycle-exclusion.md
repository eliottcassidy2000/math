---
id: THM-3143
title: "FC(3) rank-at-most-one quadratic phases: discriminant-class edge currents cannot cancel"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every algebraic
  quadratic q on the coordinate two-simplex with Hessian rank at most one,
  int_Delta exp(q) dA is nonzero, and it is 1/2 iff q is identically zero.
  The rank-zero boundary is THM-3116.  THM-3142 and the inherited first part
  of this theorem handle quadratics aligned with one affine coordinate.  A
  rank-one phase with a transverse linear term is converted by Cauchy--Green
  (noncollinear leading coordinate) or Green (collinear leading coordinate)
  into edge E-functions.  Equal edge discriminants, rather than raw edges,
  are the correct first-order blocks.  Singleton, two-edge-path, full-cycle,
  and vertical-edge source audits show that every discriminant block survives
  every endpoint-value quotient.  The full-cycle obstruction is the fact
  that an odd three-cycle of nonzero equal-square edge increments cannot
  close.  The jointly independent discriminant-block system has sole
  denominator s, so Beukers Corollary 1.4 at s=1 excludes 0 and 1/2.  This
  closes Hessian rank at most one, not full-rank quadratics and not FC(3).
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

# THM-3143 — rank-at-most-one quadratic phases on the triangle

## 0. Promoted statement and reduction map

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},       area(Delta)=1/2,
q in Qbar[u,v],                          deg q<=2,
K_q(s)=int_Delta exp(s q(u,v))du dv.                         (0.1)
```

If the Hessian of `q` has rank at most one, then

```text
K_q(1)!=0,                                                     (0.2)
K_q(1)=1/2       iff       q=0 identically.                    (0.3)
```

The rank-zero case is the affine theorem THM-3116.  If the Hessian has rank
one, its homogeneous quadratic part is

```text
q_2=A ell^2,              A!=0,                               (0.4)
```

for an algebraic nonconstant linear form `ell`.  There are two faithful
normal forms.

* If the three vertex images of `ell` are noncollinear in `C`, put `z=ell`.
  Then `z` and `conjugate(z)` span the two-dimensional algebraic affine-linear
  space modulo constants, so uniquely

  ```text
  q=A z^2+Bz+C+lambda conjugate(z).                            (0.5)
  ```

  The case `lambda=0` is the aligned Cauchy--Green proof in sections 1--7.
  The transverse case `lambda!=0` is proved in sections 8--11.
* If the images of `ell` are collinear, divide by an algebraic direction of
  that line and write `ell=eta x`, with `x` real-valued and algebraic.  After
  choosing an independent real algebraic affine coordinate `y`,

  ```text
  q=A x^2+Bx+C+lambda y.                                      (0.6)
  ```

  The case `lambda=0` is THM-3142; the transverse case is section 12.

Thus the former affine-coordinate restriction is discharged, not silently
dropped.  Full-rank Hessians remain outside the theorem, and (0.2)--(0.3) do
not prove `FC(3)`.

## 1. Inherited aligned statement and scope

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

Sections 1--7 prove the inherited **affine-coordinate** or
**univariate-aligned** core `q=P(ell)`.  The factorization is essential for
the common-operator argument in those sections, but it is not a restriction
on the promoted theorem: sections 8--12 supply the missing transverse
linear term.  The theorem still does not address a full-rank leading form
or prove `FC(3)`.

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

The normal and optimized runs are byte-identical.  Section 13 records the
current combined source and output hashes after the transverse controls.

## 8. Noncollinear transverse normal form

Assume now that the leading coordinate has noncollinear vertex images and
that the transverse coefficient in (0.5) is nonzero:

```text
q(z,conjugate(z))=A z^2+Bz+C+lambda conjugate(z),
A lambda!=0.                                                  (29)
```

Orient its image triangle counterclockwise, write its vertices as `z_j`, and
put

```text
J=Im(conjugate(z_1-z_0)(z_2-z_0))>0,
m_j=conjugate(z_(j+1)-z_j)/(z_(j+1)-z_j),
b_j=conjugate(z_j)-m_j z_j.                                  (30)
```

On edge `j`, `conjugate(z)=m_jz+b_j`, so the restriction of `q` is

```text
Q_j(z)=Az^2+beta_j z+gamma_j,
beta_j=B+lambda m_j,       gamma_j=C+lambda b_j,
delta_j=beta_j^2-4A gamma_j.                                 (31)
```

Define the edge and endpoint E-functions

```text
H_j(s)=integral_(z_j)^(z_(j+1)) exp(sQ_j(z))dz,
E_j(s)=exp(s q(z_j,conjugate(z_j))).                          (32)
```

Since `d_bar exp(sq)=s lambda exp(sq)`, Cauchy--Green gives

```text
2 i s lambda J K_q(s)=N(s):=H_0(s)+H_1(s)+H_2(s).             (33)
```

The zero at `s=0` on the right is removable.  Unlike (14), the transverse
cycle need not vanish: the three edge restrictions are different
quadratics.

For an endpoint `z_k` incident to edge `j`, set

```text
d_k^(j)=2A z_k+beta_j.                                       (34)
```

The quadratic identity

```text
(2Az+beta_j)^2=4A Q_j(z)+delta_j                              (35)
```

gives the exact edge equation

```text
D_(delta_j)H_j=d_(j+1)^(j)E_(j+1)-d_j^(j)E_j,
D_delta=4As d/ds+2A+delta s.                                 (36)
```

The operators in (36) are generally different.  This is the point at which
the aligned proof cannot simply be copied.

## 9. The discriminant blocks and their source audit

Let `F=Qbar(s)`, and let `mathcal E` be the `F`-span of the distinct endpoint
exponentials together with `1` (with `1` identified with the endpoint class
of value zero when necessary).  For every edge-discriminant value `delta`,
put

```text
I_delta={j:delta_j=delta},       G_delta=sum_(j in I_delta)H_j. (37)
```

Equation (36) makes each `G_delta` a first-order block modulo `mathcal E`.
We first prove that no block belongs to `mathcal E`.

Suppose a block did.  After grouping equal endpoint values `xi`, write

```text
G_delta=sum_xi R_xi(s)exp(xi s),        R_xi in F.             (38)
```

Applying `D_delta` shows that a grouped source coefficient `tau_xi` obeys

```text
4AsR_xi'+[2A+(delta+4Axi)s]R_xi=tau_xi.                       (39)
```

Whenever

```text
tau_xi!=0,                 delta+4Axi!=0,                     (40)
```

equation (39) has no rational solution.  A pole away from zero gains one
order under the derivative.  At zero, a pole of order `n>=1` has nonzero
leading multiplier `2A(1-2n)`.  Thus `R_xi` is a polynomial, whose top term
then gains one degree from `(delta+4Axi)sR_xi`.  This contradicts the
nonzero constant right side.

It remains to find (40) for every possible block and endpoint-value
partition.

### 9.1 One edge

On one edge from `a` to `b`, the two endpoint derivatives differ by
`2A(b-a)!=0`.  If the endpoint values differ, at least one source is nonzero
and its square is

```text
delta+4Aq(a)=d_a^2       or       delta+4Aq(b)=d_b^2.         (41)
```

If the values agree, quadratic reflection gives `d_b=-d_a!=0`, and the
grouped source is `d_b-d_a!=0`.  A third vertex in the same exponential
class contributes no source to this singleton block.  Hence (40) holds.

### 9.2 Two adjacent edges

Write the directed path as `a -> b -> c`.  Equal discriminants give equal
squares for the incoming and outgoing derivatives at `b`.  The two edge
slopes are distinct in a nondegenerate triangle, and `lambda!=0`, so those
derivatives are distinct.  Therefore they are `d` and `-d`, with `d!=0`.
The middle source is `2d`, and its square witnesses (40).

If `q(b)` merges with `q(a)`, reflection on the first edge makes the source
at `a` equal to `d`, so the grouped source is `3d`.  The same statement holds
for a merge with `c`.  If only `q(a)=q(c)`, the middle singleton still has
source `2d`.  If all three values merge and `a!=c`, the sole grouped source
is `4d=2A(c-a)!=0`.  Thus every two-edge block survives every endpoint
partition.

### 9.3 All three edges: the odd-cycle obstruction

If all three discriminants agree, the source at vertex `j` is

```text
c_j=d_j^(j-1)-d_j^(j)=lambda(m_(j-1)-m_j).                   (42)
```

Adjacent slopes are distinct, so `c_j!=0`.  The two derivatives in (42)
have equal squares but are distinct; hence they are opposite and their
common square `delta+4Aq(z_j)` is nonzero.  With three distinct endpoint
values, every source survives.  If exactly two values merge, their grouped
source is `-c_k`, because `sum_j c_j=0`, and is again nonzero.

Only a merger of all three endpoint values could erase (42).  That merger
is incompatible with a common discriminant.  Indeed, equality of the values
at the ends of edge `j` makes its endpoint derivatives opposite, whence

```text
delta+4Axi=A^2(z_(j+1)-z_j)^2.                               (43)
```

The left side is independent of `j`.  Thus the three nonzero edge increments
are each one of the two square roots `+h,-h`.  Three signs of one nonzero
number cannot sum to zero, whereas the increments around a cycle must sum
to zero.  This odd-cycle contradiction completes the source audit and proves

```text
G_delta notin mathcal E             for every delta.          (44)
```

## 10. Joint functional independence, not total-period independence alone

On the quotient by `mathcal E`, let

```text
T=d/ds+1/(2s),             kappa_delta=delta/(4A).             (45)
```

Equations (36)--(37) say

```text
T[G_delta]=-kappa_delta[G_delta].                             (46)
```

The distinct discriminant blocks are jointly independent over `F`.  If a
relation with a minimal number of blocks existed, divide by one coefficient
and write its first coefficient as one.  Applying `T+kappa_0` removes that
block.  Minimality forces every remaining rational coefficient `R_delta` to
satisfy

```text
R_delta'/R_delta=kappa_delta-kappa_0!=0.                      (47)
```

No nonzero rational function has a nonzero constant logarithmic derivative.
Together with (44) and the usual independence of distinct exponentials,
this proves that

```text
{G_delta: delta occurs} union {exp(xi s): xi occurs} union {1} (48)
```

is linearly independent over `Qbar(s)`.

Every function in (48) is an E-function: the coefficient-and-denominator
argument of section 5 applies edgewise to every algebraic `Q_j`, and finite
sums preserve the class.  Summing (36) within a class gives a first-order
system for exactly this vector, with common denominator `s`;
the endpoint exponentials satisfy `E_xi'=xi E_xi`.  Therefore `s=1` is an
ordinary algebraic point, and Beukers Corollary 1.4 makes the values in (48)
linearly independent over `Qbar`.

There is a validity warning here.  It is not enough to prove only that the
raw total `N=sum_jH_j` is independent from the endpoint exponentials and then
invoke Beukers on that list.  When the `delta_j` differ, that list is not
closed under differentiation; an arbitrary basis completion can introduce
an apparent singularity at `s=1`.  The discriminant-class vector (48) is the
required ordinary-point repair.  The raw-`N` shortcut is forbidden.

## 11. The transverse noncollinear values

Since `N=sum_delta G_delta`, equation (33) at `s=1` gives

```text
2 i lambda J K_q(1)=sum_delta G_delta(1).                     (49)
```

If `K_q(1)=0`, (49) is a nontrivial algebraic relation among the values in
(48).  If `K_q(1)=1/2`, subtracting the nonzero algebraic constant
`i lambda J` gives another such relation after adjoining `1`.  Both are
impossible.  This proves (0.2)--(0.3) in the noncollinear transverse case.

## 12. Collinear leading coordinate and the vertical-edge boundary

In normal form (0.6), let the real algebraic `(x,y)` coordinates map `Delta`
to an oriented nondegenerate triangle of Jacobian `J>0`.  Green's formula is

```text
-s lambda J K_q(s)=integral_(boundary T) exp(sq)dx.            (50)
```

Every nonvertical edge has a finite line equation `y=m_jx+b_j`; its
restriction is exactly (31), with `z` replaced by `x`, and satisfies
(35)--(36).  A vertical edge contributes zero because `dx=0`.

If the three vertex `x`-values are distinct, all three edges are active and
sections 9--11 apply verbatim with real slopes.  In particular, the
all-value/full-discriminant collision would again force three nonzero signed
equal-square `x`-increments to sum to zero.

If two `x`-values agree, the edge between them is vertical and the other two
active edges form a directed path whose initial and final `x`-coordinates
agree.  For distinct discriminants both blocks are singletons.  For a common
discriminant, the derivatives at the shared vertex are `d,-d`, `d!=0`, and
the middle source is `2d`.  A collision with either vertical endpoint gives
the grouped source `3d`; a collision of the two vertical endpoints leaves
the middle source unchanged.  An all-three collision would require the two
signed active-edge increments both to equal `d/A`, but those increments are
opposites.  It is therefore impossible.  This covers vertical edges,
repeated leading-coordinate values, and every endpoint-value collision.

The discriminant-class system remains ordinary at one, and (50) excludes
`K_q(1)=0,1/2` exactly as in section 11.  Together with THM-3116, THM-3142,
and sections 1--7, this completes the reduction in section 0.

## 13. Extended exact verification and sharp scope

The same script now additionally checks:

* every transverse edge restriction and its discriminant ODE;
* the singleton reflection source;
* the `2d`, `3d`, and `4d` two-edge path currents;
* all endpoint partitions of the three-cycle current and the eight-sign
  odd-cycle obstruction;
* Green orientation on a triangle with a vertical edge; and
* the vertical-path collision formulas and repeated-coordinate boundary.

The verifier does not numerically sample periods and does not inject the
desired conclusion.  It checks the algebraic identities used by the proof;
the E-function value step is the cited Beukers theorem.  Full-rank Hessians
are the exact remaining quadratic boundary.

Run the two byte-comparison controls displayed in section 7.  The current
frozen hashes are:

```text
source sha256 = ab46ba7df9b66159154715b71d0099c3c875c3cd104dde2f641c27604cb43a2e
output sha256 = 273a4cc7fb40fb7c8411653cf5c894fbcc397e83a23cd7ea0235a1921b93cd2a
```
