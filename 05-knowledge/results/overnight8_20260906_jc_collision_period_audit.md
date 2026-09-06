# Independent audit of the ninth-source collision-period obstruction

**PASS: analytic proof and exact independent controls.** I read the complete
[producer report](overnight8_20260906_jc_collision_period.md)
and [producer source](../../04-computation/overnight8_20260906_jc_collision_period.py),
as well as the proved
[seventh source-module theorem](overnight7_20260906_jc_fifth_module.md)
and [THM-4046](../../01-canon/theorems/THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction.md).
No mathematical correction is required. The additional all-form converse
below is independently proved and checked. It strengthens the interpretation
of the source conditions without changing the producer's obstruction.

## 1. General filtered lemma: hypotheses and proof accepted

Work in characteristic zero with fixed source coordinates `(x,t)`, a single
smooth target germ at a common point y0, and `N>=2`. Let
`ftilde-f=t^N v+O(t^(N+1))`. All retained source branch tangents T_i must be
nonzero, and at least two tangent lines must be distinct. Both old and new
labelled section families are assumed to collide modulo `t^(N+1)`.

These are sufficient hypotheses even if the target has more than two
coordinates. At the first differing section coefficient of order j<N,
subtracting the two common-image equations gives `r_i T_i=b` for a common
vector b. Two distinct tangent lines have zero intersection; hence b=0,
and all r_i=0. Induction proves equality of the entire section prefix below
N, not merely equality of its values. At order N the equation is exactly

```
v_i+r_i T_i=b.
```

For an arbitrary regular target two-form omega, differentiating the changed
map in t first introduces `N t^(N-1)v`; its changed x derivative and the
changed coefficients of omega start at order N. Thus at retained points

```
delta J_(N-1)(x_i)=N omega_y0(T_i,v_i).
```

All lower coefficient functions agree identically in x, so all their
retained x derivatives also agree. If `sum ell_i T_i=0`, substituting the
common-motion equation annihilates the weighted change by alternation:

```
sum ell_i omega(T_i,v_i)
 =sum ell_i omega(T_i,b-r_i T_i)=0.
```

This establishes the claimed weighted period identity. It does not claim
that the last coefficient is unchanged pointwise, and it does not require
closed or decomposable forms. A single coherent target form at y0 is
essential; independently assigned branch forms would be a different object.

The tangent hypothesis has a sharp local hostile. Take two base germs
`f1(s,t)=(s,0)` and `f2(s,t)=(s,s^2)`, with old sections s=0. Perturb only
the second to `(s,s^2-t^2)` and take both new sections s=t. Both old and
new pairs collide, but both tangents are `(1,0)`. At N=2 and weights `(1,-1)`,
the weighted change in the coefficient of `ds wedge dt` is 2. The section
prefix already moves at order one. This is an explicit failure if the
distinct-line premise is dropped; the producer retains that premise.

## 2. The exceptional compiler and the source/target distinction

The inherited source map is the Russell compiler on
`q=Q(x)+t^2,w=t`, with normalized target coordinates `(C/3,E+3,w)`.
Its base tangents and normal q-directions are exactly

```
T_-=(1,-9,0), T_0=(1,4,0), T_+=(1,9,0),
U_-=(2/3,-2,0), U_0=(0,4,0), U_+=(-2/3,-2,0).
```

The independent symbolic calculation verifies the literal compiler
Jacobian `Jac_(x,q)(C,E)=6(D+1)(D^2+2D-2)`, all three vectors above, and
`det_(y,z)(T_i,U_i)=4`. With `ell=(5,-18,13)`, the complete three-coordinate
tangent relation is `sum ell_i T_i=0`, and the inherited normalization is
`Lambda=ell/18`.

The old uncorrected source sections collide through t-order nine because
the constant pencil agrees through fourth s-order and s=t^2. Therefore
the general lemma applies with N=9 to every smooth-target map perturbation
`ftilde=f0+O(t^9)` whose labelled branches still collide modulo t^10.

The entire old retained relation survives: J0 through J7 and their retained
x-jets are unchanged, while its only J8 block, Lambda(J8), is unchanged.
The nonzero constant response `S0(1)=kappa` is an inherited proved
THM-4046 fact. Thus the same old relation excludes a nonzero constant
pullback density for every regular target two-form. This audit does not
repeat THM-4046's complete 45-by-495 matrix or its quartic norm certificate.
Those are explicit proved dependencies, not finite-field conclusions here.

The ninth corrections to an actual **target pair** in THM-4046 act through
its old restriction operator. The present calculation instead changes the
**source map**. It is not a duplicate use of the target-pair obstruction.
The proof also does not assert that every formal target-map perturbation
comes from a polynomial source automorphism, or replace constant density
with a formal unit density. Those boundaries are correct in the producer.

## 3. Exact source conditions and the additional all-form converse

For an ambient source perturbation

```
(x,q,w)=(x+t^9 r,Q(x)+t^2+t^9 s,t+t^9 k)+O(t^10),
h_i=s(x_i)-Q'(x_i)r(x_i),
```

the target forcing is `r_i T_i+h_i U_i+k_i e_w`. Tangential r_i is absorbed
by a section correction. The surface collision system has six equations
and five unknowns, rank five. The full target system has nine equations
and six unknowns, rank six. Consequently solvability is exactly

```
L(h)=0,     k_-=k_0=k_+,     L=5 ev_- -18 ev_0+13 ev_+.
```

For `omega=A dy wedge dz+B dy wedge dw+C dz wedge dw`, only the base
coefficients A0,B0,C0 enter the first response. The full formula is

```
delta Lambda(J8)=2 A0 L(h)+(B0 L(k)+C0 L(e*k))/2,
e=(-9,4,9).
```

The factor9 from differentiating t^9 and the factor18 in Lambda are both
retained. For graphs k=r=0 this specializes to `delta J8(i)=36A0 h_i`
and `delta Lambda=2A0 L(h)`. It applies to the genuine polynomial graph
shears, with no assumption of descended Hamiltonian primitives.

**Converse for this source class:** vanishing of this response for every
base two-form is equivalent to collision compatibility. Indeed it requires
`L(h)=L(k)=L(e*k)=0`. The two independent rows `(ell_i)` and `(ell_i e_i)`
have common kernel exactly the constant k-vectors. This proves the converse,
including the stable coordinate, by exact linear algebra.

It is false for one arbitrarily selected form. Take h=x and k=x. The
response is `16A0+4B0+81C0`; the nonzero form with `(A0,B0,C0)=(1,-4,0)`
has zero response, although the collision splits. A second hostile takes
`k=4x^2-9x`, for which L(k)=0 while L(e*k) is nonzero. The C form slot is
necessary to detect this stable-coordinate split.

## 4. Independent exact engines and controls

The audit source imports no producer or repository code. It combines a
literal univariate series compiler, an independently differentiated
dual-number compiler, symbolic exterior algebra, and modular overdetermined
section equations. Fresh exceptional embeddings are
`(p,alpha)=(449,120),(467,169)`, different from the seventh producer's primes.

At each fresh embedding it verifies the quartic equation and solves the
old six surface equations for three branch corrections and two common-image
coefficients, successively through t-order nine. All nine stages agree in
the literal compiler, and all odd section terms vanish. The order-ten
system is inconsistent, confirming the genuine next collision debt. This
independently reconstructs the exact old contact premise of the new lemma.

Eight source controls cover zero, a constant graph normal, a nonconstant
kernel graph normal, h=x, a tangential reparametrization with common stable
motion, k=x, a stable vector in ker L but not constant, and the one-form
hidden split. The full nine-equation changed collision system agrees with
the necessary-and-sufficient conditions in every case. For compatible
controls, substitution of the solved ninth section corrections makes all
three literal target coordinates agree modulo t^10.

Separately, at each fixed retained source point, a dual-number source
compiler carries both the map value and its x derivative through t-order
ten. Differentiating its t-series and forming the three actual coordinate
wedges verifies864 density controls across the two primes. Constant and
nonconstant target coefficients `1,y,z,w,yz,w^2` are included. Every lower
density J0 through J7 is unchanged, and every J8 variation matches its
predicted base-form response. These are genuinely different arithmetic
paths from the producer's characteristic-zero local first-jet expansion.

The finite controls support the analytic argument; no characteristic-zero
or all-degree claim is inferred solely from them.

## 5. Frozen reproduction

```
python -B 04-computation/overnight8_20260906_jc_collision_period_audit.py
python -B -O 04-computation/overnight8_20260906_jc_collision_period_audit.py
```

Both runs pass **1,871 optimization-live gates** with byte-identical LF
output. The separate normal and optimized transcripts are frozen.

```
source 01600a2c46790da5716a71017202e05181add970fdab598eae76e354c0bc96c4
output 24893fd9f63d3788dab4ed677e8846a0593a423773fee2d324a8fd24a05fec18
```

All audit files remain outside the repository. The producer, shared
navigation, and Git state were not modified by this referee.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
