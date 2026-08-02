---
id: THM-3250
title: "FC(3) noncollinear pure-power phases: the two-residue vertex-turn current cannot cancel"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  Let ell be
  an algebraic affine polynomial on the coordinate two-simplex whose three
  vertex values are noncollinear in C.  For every integer d>=3 and algebraic
  A,h,C with A nonzero, the genuinely higher-degree rank-one phase
  q=A(ell-h)^d+C satisfies int_Delta exp(q)dA !=0 and !=1/2.  After setting
  z=ell-h, Cauchy--Green turns the area period into two vertex-turn primitive
  blocks G_0,G_1.  They obey `dsG_0'+G_0=-S` and
  `dsG_1'+2G_1=S`, where the grouped endpoint source is
  `S=sum c_j z_j^2 exp(Az_j^d s)` and the exact turning identity is
  `S(0)=-W!=0`.  Endpoint-power collisions therefore cannot erase every
  source.  The residues 1/d and 2/d are nonintegral and incongruent for
  d>=3, making G_0,G_1 jointly independent modulo the endpoint
  exponentials.  The E-function system is ordinary at s=1, so Beukers
  Corollary 1.4 excludes both forbidden values.  This closes every pure
  cubic and higher power in one noncollinear affine coordinate, not a
  general cubic, a collinear coordinate, nonalgebraic coefficients, or
  FC(3).
source: codex-2026-08-02-fc3-cubic-frontier
depends_on:
  - THM-3203-fc3-complex-affine-coordinate-quadratic-cauchy-green-cycle-exclusion
related:
  - THM-3039-the-FC-n-exponential-period-bridge-forced-level-is-the-simplex-volume
  - THM-3116-fc3-flat-top-simplex-singular-coefficient-and-affine-nonvanishing
  - THM-3202-fc3-real-affine-coordinate-quadratic-two-piece-spline-exclusion
external:
  - "Cauchy--Green formula."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc3_noncollinear_pure_power_turn_current_thm3250.py
output: 05-knowledge/results/fc3_noncollinear_pure_power_turn_current_thm3250.out
script_sha256: e343e9a6dcd65f6734850432908c9d00c0652a3171ebc2957cb2da4692fc4578
output_sha256: 494f98657e3cb32dd3d2ced08c79651e58a73e76ebdcb81e4e5a807ee011fa58
---

# THM-3250 -- noncollinear pure-power phases on the triangle

## 1. Statement

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},        area(Delta)=1/2.     (1)
```

Let `ell in Qbar[u,v]` be affine, and suppose its values at the three
vertices of `Delta` are noncollinear as points of `C`.  Fix

```text
d in Z, d>=3,        A,h,C in Qbar,        A!=0,              (2)
q(u,v)=A(ell(u,v)-h)^d+C.
```

Put

```text
K(s):=int_Delta exp(sq(u,v))du dv.
```

Then

```text
K(1)!=0,                                                      (3)
K(1)!=1/2.                                                    (4)
```

In particular, `d=3` closes the noncollinear affine-coordinate **pure
cubic** branch of the simplex exponential obstruction.  Translation by
`h` is allowed, but an arbitrary cubic

```text
A(ell-h)^3+B(ell-h)+C,              B!=0,                    (5)
```

is not claimed.  Thus (3)--(4) are a genuine first result beyond quadratic
degree, not all cubic simplex phases and not `FC(3)`.

## 2. Cauchy--Green keeps the turning sidecar

Put `z=ell-h`, and denote its three vertex values cyclically by
`z_0,z_1,z_2`.  Relabel the vertices, if necessary, so their image triangle
`T` is counterclockwise.  Set

```text
e_j=z_(j+1)-z_j,
J=Im(conjugate(e_0)(z_2-z_0))>0,
W=2iJ!=0,                                                       (6)
```

with indices modulo three.  The real-affine map `Delta -> T` has Jacobian
`J`.  Since

```text
d_bar(conjugate(z) exp(sAz^d))=exp(sAz^d),                    (7)
```

Cauchy--Green gives the exact identity

```text
W K(s)=e^(Cs) integral_(boundary T) conjugate(z)exp(sAz^d)dz. (8)
```

On the oriented edge `j`, define

```text
m_j=conjugate(e_j)/e_j,
b_j=conjugate(z_j)-m_jz_j,                                   (9)
```

so that `conjugate(z)=m_jz+b_j` along that edge.  The data `m_j` are the
edge directions modulo `pi`; they are the geometric sidecar that an
endpoint-value quotient would forget.

For `k=0,1`, define the path-independent entire primitives

```text
J_(k,j)(s)=integral_0^(z_j) t^k exp(sAt^d)dt.                (10)
```

The integrand is entire in `t`, so the path in (10) is immaterial.  Edge
integration followed by regrouping at vertices turns (8) into

```text
W K(s)=e^(Cs)(G_0(s)+G_1(s)),                                (11)

G_1=sum_j c_j J_(1,j),
G_0=-sum_j c_j z_j J_(0,j),
c_j=m_(j-1)-m_j.                                             (12)
```

Indeed the coefficient of `J_(1,j)` is `m_(j-1)-m_j=c_j`, while

```text
b_(j-1)-b_j=(m_j-m_(j-1))z_j=-c_jz_j.                       (13)
```

This fixes both the Cauchy--Green orientation and primitive normalization.

## 3. The two residue blocks and their common source

Differentiating `t^(k+1)exp(sAt^d)` with respect to `t` gives

```text
ds J_(k,j)' +(k+1)J_(k,j)
   =z_j^(k+1)exp(Az_j^d s),             k=0,1.               (14)
```

Consequently

```text
dsG_1'+2G_1=S,
dsG_0'+ G_0=-S,                                               (15)

S(s)=sum_j c_jz_j^2 exp(Az_j^d s).                           (16)
```

The same source occurs with opposite signs, but the two blocks have
different residues.  This is the key structural feature: Cauchy--Green does
not produce one opaque boundary period, but two adjacent Mellin residue
classes tied by one exact turning current.

At `s=0`, equations (10)--(12) give

```text
G_1(0)= (1/2)sum_j c_jz_j^2,
G_0(0)=-sum_j c_jz_j^2.                                      (17)
```

On the other hand, (11) and `K(0)=1/2` give

```text
G_0(0)+G_1(0)=W/2.
```

Therefore the vertex-turn identity is

```text
sum_j c_jz_j^2=-W!=0.                                        (18)
```

This identity is also the degree-zero Cauchy--Green formula written in the
vertex basis.  It is valid for every allowed triangle, not merely generic
endpoints.

## 4. Endpoint-power collisions cannot kill the source

Group the three endpoint exponentials in (16) by their algebraic exponents

```text
xi=Az_j^d,
tau_xi=sum_(j:Az_j^d=xi)c_jz_j^2.                            (19)
```

Equation (18) says

```text
sum_xi tau_xi=-W!=0.                                         (20)
```

Thus some grouped source survives, even if two or all three `d`th powers
coincide.  Moreover the group `xi=0` has zero source: because `A!=0`, every
vertex in that group has `z_j=0`, hence `c_jz_j^2=0`.  We may therefore fix
a source group satisfying

```text
tau_xi!=0,                         xi!=0.                     (21)
```

This handles every endpoint collision needed below, including collisions
with the extra exponential `exp(-Cs)` or with the constant function.  Those
functions are simply merged into the same distinct-exponential basis; the
nonzero source coefficient in (21) remains the coefficient on that basis
element.

## 5. Each residue block survives the exponential quotient

Let `mathcal E` be the `Qbar(s)`-span of the distinct functions among

```text
{exp(Az_j^d s):j=0,1,2} union {exp(-Cs),1}.                  (22)
```

Distinct exponentials are linearly independent over `Qbar(s)`.  We first
show that neither `G_0` nor `G_1` belongs to `mathcal E`.

Suppose `G_k`, where `k=0` or `1`, were a rational combination of the
exponentials in (22).  Apply the corresponding operator in (15) and compare
the coefficient of the nonzero source class (21).  Its rational coefficient
`R` would solve

```text
dsR'+(k+1+d xi s)R=epsilon tau_xi,                            (23)
```

where `epsilon=-1` for `k=0` and `epsilon=1` for `k=1`.
Equation (23) has no rational solution.

* At a finite nonzero pole, `dsR'` has one higher pole order than every
  other term.
* At zero, a pole of order `n>=1` has leading multiplier
  `k+1-dn`.  This cannot vanish because

  ```text
  1<=k+1<=2<d.                                                (24)
  ```

* Hence `R` is a polynomial.  Since `xi!=0`, the term `d xi sR` raises the
  degree of every nonzero polynomial, whereas the right side is a nonzero
  constant.

This contradiction proves

```text
[G_0]!=0, [G_1]!=0           in Qbar(s)-functions/mathcal E. (25)
```

## 6. The two blocks are jointly independent

Modulo `mathcal E`, put `T=ds d/ds`.  Equation (15) becomes

```text
T[G_0]=-[G_0],                 T[G_1]=-2[G_1].                (26)
```

If the two classes were rationally dependent, (25) would allow a relation

```text
[G_0]+R(s)[G_1]=0,             R in Qbar(s), R!=0.            (27)
```

Apply `T+1` to (27).  The first term disappears, leaving

```text
dsR'=R.                                                       (28)
```

A nonzero rational solution of (28) would have order `1/d` at zero.  A
rational function has integral order, and `d>=3`, so none exists.  Thus

```text
G_0,G_1, the distinct endpoint exponentials, exp(-Cs), and 1 (29)
```

are linearly independent over `Qbar(s)`, after merging any coincident
exponentials.

## 7. E-functions and the ordinary point

Every primitive in (10) has the expansion

```text
J_(k,j)(s)=sum_(n>=0)
  A^n z_j^(dn+k+1) s^n / ((dn+k+1)n!).                       (30)
```

It is an E-function.  Its coefficients are algebraic with exponential
conjugate growth.  After fixed algebraic denominators are cleared, the
common denominator through index `n` divides a fixed exponential factor
times `lcm(1,...,dn+k+1)`, which has exponential growth.  Equation (14)
gives holonomicity.  Hence `G_0,G_1` are E-functions, as are all the
exponentials in (29).

Equations (15) together with `E_eta'=eta E_eta` form a first-order rational
system whose only displayed finite denominator is `s`.  Therefore `s=1` is
an ordinary algebraic point.  By Beukers Corollary 1.4, functional linear
independence in (29) implies algebraic linear independence of their values
at `s=1`.

## 8. Excluding zero and the forced level

At `s=1`, (11) is

```text
W K(1)=e^C(G_0(1)+G_1(1)).                                  (31)
```

If `K(1)=0`, then `G_0(1)+G_1(1)=0`, contradicting section 7.  If
`K(1)=1/2`, then

```text
G_0(1)+G_1(1)-(W/2)exp(-C)=0,                                (32)
```

again contradicting section 7 because `W` is nonzero and algebraic.  This
proves (3)--(4).

The second exclusion is exactly the level needed by THM-3039's
simplex-period bridge.  It rules out this phase shape as the projective
exponential obstruction of a hypothetical homogeneous `FC(3)`
counterexample.  It does not supply the separate projective-leading-form
reduction needed to put every counterexample into this shape.

## 9. Failure boundary and connection contract

The proof boundary is sharp for this two-residue mechanism.

* At `d=2`, the `k=1` zero-pole multiplier in (23) can vanish at pole order
  one: `k+1-d=0`.  Indeed that primitive is elementary.  THM-3203 closes all
  algebraic quadratic phases by its stronger discriminant-block Euler-flux
  argument; the present proof does not replace it.
* A general depressed cubic with the linear term in (5) couples the two
  residue blocks.  The diagonal quotient law (26) is then lost.
* If the vertex values of `ell` are collinear, `W=0` and the area-to-complex-
  boundary isomorphism degenerates.  A one-dimensional spline pushforward is
  the faithful carrier, but its cubic analysis is not asserted here.
* Algebraicity is used by the E-function value theorem.  No statement is
  made for arbitrary complex coefficients.

The exact connection is

```text
source:    two-dimensional simplex period of A(ell-h)^d+C
target:    two adjacent primitive residue blocks on the image triangle
map:       Cauchy--Green, edge primitives, then vertex regrouping
preserves: the period and the forbidden values 0 and 1/2
loses:     transverse area geometry after endpoint-power grouping
sidecar:   the turning coefficients c_j and identity sum c_jz_j^2=-W
hostile:   total dth-power collision, partial collision, and z_j=0.       (33)
```

The natural carrier is a three-cycle with algebraic turning weights, not a
tournament: edge directions have no intrinsic pairwise winner relation.

## 10. Exact verification

Run

```text
python3 04-computation/fc3_noncollinear_pure_power_turn_current_thm3250.py
python3 -O 04-computation/fc3_noncollinear_pure_power_turn_current_thm3250.py
```

The normal and optimized transcripts are byte-identical.  The companion
checks:

* Cauchy--Green orientation on 24 exact monomial cells;
* the symbolic line-offset regrouping and `S(0)=-W` identity;
* 288 coefficientwise instances of (14);
* total cubic, total quartic, partial cubic, and zero-vertex endpoint
  collisions;
* 90 direct simplex-moment / primitive-boundary identities; and
* the nonintegral residues for `3<=d<=40`, with the quadratic failure planted
  as a hostile control.

The computation verifies the algebraic identities and normalization; it
does not numerically infer the Beukers transcendence step.

**End of proof.**
