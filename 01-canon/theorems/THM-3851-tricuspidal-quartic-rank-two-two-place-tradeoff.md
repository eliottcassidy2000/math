---
id: THM-3851
title: "The tricuspidal quartic trades rank-two resolvent torsion for two places at infinity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The parametrized
  tricuspidal quartic has a smooth
  bitangent whose complement has exactly two normalization places at
  infinity.  Its split-boundary quadratic-resolvent surface has constant
  units and class group (Z/3)^2, with a rank-two local cusp-address map and
  deck-anti-invariant torsion.  No line has fourfold pullback to the
  normalization, so no projective affine chart of this embedded quartic has
  one normalization place.  The reusable THM-3841 Jelonek mechanism makes
  the two-place branch hostile to a deleted-ramification plane atlas.  No
  cubic completion, Keller atlas, or Jacobian counterexample is claimed.
source: jc_zero_debt_lift / one-place S3 completion design lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn3709, 2026-08-23).  The audit
  rederived the birational normalization and all three A2 packets; checked
  the two distinct smooth bitangent contacts and both finite and infinite
  cases of the fourth-power line test; traced the weak degree-two del Pezzo,
  split-boundary, affine class-group, and unit exact sequences; and checked
  the uniform E7 quotient, local rank-two map, deck anti-action, and Kummer
  identification.  It also verified that the Jelonek conclusion in Section
  6 remains explicitly conditional on a deleted ramification completion.
  Normal and optimized runs byte-match the frozen 763-gate transcript.  A
  second independent audit identifies the exact reciprocal-cubic and
  ordered-root torus completion.  Its 48-gate companion verifies the
  Fourier and normalization changes, Sylvester discriminant, Laurent source,
  different, branch injectivity, two missing places, Vandermonde square and
  quotient diagram.  It corrects the scope: the ordinary toric S3 closure
  exists, but it is unit-rich and is not the desired plane completion.
depends_on:
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
script: 04-computation/jc2_tricuspidal_quartic_rank_two_two_place_thm3851.py
output: 05-knowledge/results/jc2_tricuspidal_quartic_rank_two_two_place_thm3851.out
script_sha256: 2023e45b5e63bb3447347c0e78ffe134224b23ea4df0114b18a453722b880c10
output_sha256: b3574cd466e4c2d635a852493d7a021f0ef3b5b6ef9666fb1e42e02f459e3ff6
semantic_sha256: 5b5e5e17390dbd80042a334cf1f53d04fd548ec80520b4fc806a900d87874d2c
toric_extension_script: 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_extension_thm3851.py
toric_extension_output: 05-knowledge/results/jc2_tricuspidal_reciprocal_cubic_toric_extension_thm3851.out
toric_extension_script_sha256: 62b297e196fd56934e57c985e4cdedd861e1508b3f901e85c29972c60abc7cfd
toric_extension_output_sha256: aaf9bf047eb068d3a1a814e2be9d862fe16e5d8e69ba73b3665389135cf7309c
toric_extension_semantic_sha256: 3d1406b37d55dbb02630abd1de69a2c0175c936a44609d5471bb555196b8ad5f
toric_independent_audit_script: 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_independent_audit_thm3851.py
toric_independent_audit_output: 05-knowledge/results/jc2_tricuspidal_reciprocal_cubic_toric_independent_audit_thm3851.out
toric_independent_audit_script_sha256: 820ac69c385ff8a8a953cd7abf0541a44a1d7d35bc8b6fe48ec8743e8b82090c
toric_independent_audit_output_sha256: 6e7b2b6aa67c994ba19cbb873a40bde93300581300fc7e339825391478bcfaf1
toric_independent_audit_semantic_sha256: 614126e3142de944089d842d7e4ac48c88242b570a899157e172b23a2f54f136
fixed_product_extension_script: 04-computation/jc2_fixed_nonzero_product_marked_root_extension_thm3851.py
fixed_product_extension_output: 05-knowledge/results/jc2_fixed_nonzero_product_marked_root_extension_thm3851.out
fixed_product_extension_script_sha256: ef681020058547340afd9e6bed8434b71e10ab4ca435d36e2f43b059213f0870
fixed_product_extension_output_sha256: b33d145b2c65cc378622388b146cc3e1301d7b509574892b5b67ce86032c48ec
fixed_product_extension_semantic_sha256: b40b22611b39e328749e52952b4e6eebe4f0543b69e4ee67d9975cf36f783780
hash_basis: raw LF bytes
---

# THM-3851 -- rank two costs the one-place branch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of
characteristic zero.  Consider the basepoint-free map

```text
nu:P1_[s:t] -> P2_[A:B:C],
[s:t] |-> [s^2t^2 : t^2(s-t)^2 : s^2(s-t)^2].                 (1)
```

Its image is the tricuspidal quartic

```text
Delta_3=(AB+AC+BC)^2-4ABC(A+B+C)=0.                           (2)
```

Let

```text
ell: A+B+C=0.                                                   (3)
```

Then `ell` is a bitangent at two distinct smooth points, and the affine
normalization of `C_3 minus ell` is `G_m`.  For the normal quadratic double
plane restricted to this affine chart,

```text
Q_ell -> P2 minus ell,                    W^2=Delta_3,          (4)
```

one has

```text
Q_ell^*=k*,                          Cl(Q_ell)=(Z/3)^2.          (5)
```

The two global three-classes span a rank-two plane inside the three local
`A2` class groups, and the quadratic deck involution acts as `-1` on both.
This is the first exact quartic in this design lane with two independent
anti-invariant Kummer directions.

The gain is inseparable from a boundary loss for this embedded curve: **no
nonzero line in `P2` pulls back under `(1)` to a fourth power of one linear
form.**  Thus no choice of projective line at infinity leaves only one
normalization place.  The bitangent `(3)` attains the minimum of two places.
This is a precise obstruction for `(1)-(2)`, not a classification of all
possible higher-degree branch curves.

## 1. Normalization and the three cusps

Substitution of `(1)` into `(2)` gives zero.  On `AB!=0`, the inverse
normalization ratio is

```text
s/t=(AB+AC-BC)/(2AB).                                         (6)
```

Indeed cross-multiplication after `(1)` is an identity.  Hence `nu` is
birational onto its image.  A general target line pulls back to a degree-four
binary form, so the image has degree four.  Since the nonzero quartic `(2)`
vanishes on the irreducible image, it is the image equation and is
irreducible.  Thus `(1)` is the projective normalization.

The three singular points are exactly the coordinate vertices.  In the
three standard target charts the singular ideals have Groebner bases

```text
A=1:  B-C, C^2;          B=1: A-C, C^2;          C=1: A-B, B^2. (7)
```

There are no other projective singularities.  The normalization addresses
are respectively

```text
s=t,                         s=0,                         t=0.  (8)
```

For example, put `u=s/t` near `u=0` and work in the chart `B=1`.  The two
local coordinates may be taken as

```text
z=C/B=u^2,
x-z=A/B-C/B=u^3(2-u)/(u-1)^2.                                (9)
```

They have leading orders two and three.  At `u=1+v` and at `u=infinity`
the corresponding exact packets are

```text
v^2,
v^3(2+v)/(1+v)^2;

v^2,
v^3(2-v)/(1-v)^2,                                             (10)
```

after choosing the analogous chart coordinates.  Therefore all three
vertices are `A2` cusps.

## 2. The bitangent has exactly two places

The line `(3)` pulls back exactly as

```text
nu^*(A+B+C)=(s^2-st+t^2)^2.                                  (11)
```

The splitting is already visible on the target line:

```text
Delta_3(A,B,-A-B)=(A^2+AB+B^2)^2.                            (11a)
```

The quadratic `s^2-st+t^2` has discriminant `-3`, so it has two distinct
roots in `P1(k)`.  It is nonzero at all three addresses `(8)`, hence both
contacts occur at smooth points.  The two target points are distinct as
well: in the chart `t=1`, a root `u` has normalized image

```text
[u-1:-u:1],                                                    (12)
```

while the other root is `1-u` and swaps the first two coordinates; equality
would force `2u-1=0`, incompatible with `u^2-u+1=0`.

Consequently

```text
normalization(C_3 minus ell)
  =P1 minus V(s^2-st+t^2) isomorphic to G_m.                  (13)
```

In particular the affine branch has exactly two normalization places at
infinity, not one projective point with two uncounted branches.

## 3. No line can give one place

Take an arbitrary line

```text
L=aA+bB+cC.                                                     (14)
```

On `t=1`, its pullback is

```text
P_L(u)=cu^4-2cu^3+(a+b+c)u^2-2bu+b.                           (15)
```

Suppose first that its whole degree-four divisor were supported at a finite
address `r`.  Then for some `kappa!=0`,

```text
P_L(u)=kappa(u-r)^4.                                           (16)
```

The coefficients of `u^4` and `u^3` give

```text
c=kappa,                         r=1/2.                         (17)
```

The constant coefficient in `(16)` then gives `b=c/16`, while the linear
coefficient gives

```text
-2b=-4c(1/2)^3=-c/2,                 hence b=c/4,              (18)
```

a contradiction.

The address at infinity is a separate boundary case, not hidden by the
dehomogenization.  Fourfold support at `[1:0]` would mean
`nu^*L=kappa t^4`, so `(15)` would be the nonzero constant `kappa`.  Its
linear coefficient forces `b=0`, whereas its constant coefficient forces
`b=kappa!=0`.  This is again impossible.

The cases exhaust `P1`.  Therefore no line has a one-point pullback.  Since
`(11)` exhibits a two-point pullback, two is the exact minimum number of
normalization places among projective affine charts of this embedded
quartic.

## 4. The split-boundary quadratic class group

Let `Qbar` be the weighted-projective double plane

```text
W^2=Delta_3                in P(1,1,1,2).                     (19)
```

It is normal: the reduced branch makes the hypersurface regular in
codimension one, and a hypersurface is `S2`.

It has precisely three `A2` rational double points over the cusps and is
smooth over the two bitangent contacts.  Its minimal resolution `S` is the
standard rational weak del Pezzo surface of degree two, with

```text
Pic(S)=ZH direct_sum ZE1 direct_sum ... direct_sum ZE7,
(H^2,E1^2,...,E7^2)=(1,-1,...,-1),
-K_S=3H-E1-...-E7.                                            (20)
```

Because `(11)` is a square, the inverse image of `ell` splits into two
rational curves

```text
B_+ + B_-=-K_S,
B_+^2=B_-^2=-1,                    B_+ B_-=2.                  (21)
```

They meet once over each bitangent contact.  All six exceptional `A2` roots
lie over affine points and are disjoint from the boundary.  Hence

```text
Cl(Q_ell)=Pic(S)/<the six A2 roots,B_+,B_->.                   (22)
```

As in THM-3844, `(20)-(22)` are the standard geometric bridge; the following
quotient computation is finite exact evidence.  Fix `B_+=E7`.  The companion
lists all 126 `E7` roots, the 72 roots in its orthogonal `E6`, all 120 `A2`
subsystems there, and every unordered triple of mutually orthogonal `A2`
subsystems.  There are exactly 40 compatible `3A2` configurations.  For
every one, the eight relation rows in `(22)` have Smith profile

```text
(1,1,1,1,1,1,3,3).                                           (23)
```

Thus `(5)` is independent of the chosen marking.

One marking exposes both generators.  Put

```text
r1=H-E1-E2-E3,                  r2=H-E4-E5-E6,
r3=E2-E3,                       r4=E1-E2,
r5=E5-E6,                       r6=E4-E5,
B_+=E7,
B_-=3H-E1-E2-E3-E4-E5-E6-2E7.                                (24)
```

The relations first give

```text
E1=E2=E3=x,                    E4=E5=E6=y,
H=3x=3y,                       E7=0,                           (25)
```

and the last boundary gives `3(2x-y)=0`.  Equivalently,

```text
3(x-y)=0,                         3(2x-y)=0,                   (26)
```

so `3x=3y=0` and `x,y` are independent.  Exact integral principalizing
relations are

```text
3x=-2r1-r2+r3+2r4+2B_++B_-,
3y=-r1-2r2+r5+2r6+2B_++B_-.                                  (27)
```

This proves `Cl(Q_ell)=(Z/3)^2` directly in the displayed marking.

The unit computation is separate.  The boundary-divisor sequence gives

```text
Q_ell^*/k*=ker(ZB_+ direct_sum ZB_- -> Pic(S)).                (28)
```

The two boundary vectors in `(24)` are independent in the torsion-free
lattice `Pic(S)`, proving `Q_ell^*=k*`.

## 5. Three local directions and two global directions

Orient the three `A2` root bases as

```text
(r1,r2),                         (r3,r4),                 (r5,r6). (29)
```

For an intersection pair `(m,n)`, use the local discriminant character
`m+2n mod 3`.  The two global classes in `(25)` restrict as

```text
x |-> (1,1,0),                    y |-> (-1,0,1).               (30)
```

They span the rank-two plane

```text
z1-z2+z3=0              inside (Z/3)^3.                       (31)
```

Thus all three cusp addresses participate, but they satisfy exactly one
global relation.  The quadratic deck involution acts on the degree-two
del Pezzo lattice by

```text
sigma(D)=(D K_S)K_S-D.                                        (32)
```

Since `-K_S=B_++B_-=0` in `(22)`, `(32)` sends both `x` and `y` to their
negatives.  Both Kummer directions are anti-invariant.

On the regular locus, normality and `(5)` give the Kummer identification

```text
H^1((Q_ell)_reg,mu_3)=Cl(Q_ell)[3]=(Z/3)^2.                   (33)
```

Hence there are genuinely two independent cyclic cubic quasi-etale cover
directions.  Anti-invariance is necessary for an `S3` descent, but an
equivariant descent cocycle and a polynomial cubic completion have not been
constructed here.

## 6. Why the two places remain fatal to the current design

THM-3841 proves a reusable deleted-divisor lemma: if a normal finite
completion deletes a prime divisor whose image is a branch curve, then any
dominant plane morphism to the complement forces that branch to be an
irreducible component of the composite's Jelonek set.  Its polynomial
uniruledness argument applies to `(13)` even more cheaply.  A dominant
polynomial curve in `C_3 minus ell` would lift through the normalization and
induce

```text
k[z,z^(-1)] -> k[q].                                          (34)
```

But `z` must map to a unit of `k[q]`, hence to a scalar, contradicting
dominance.

Consequently, **conditionally on a cubic completion deleting a prime
ramification divisor which maps dominantly to this affine branch**, the
THM-3841 mechanism forbids a dominant plane atlas.  This theorem does not
claim that `(4)` alone is such a cubic completion, nor that every rank-two
torsion branch has two places.  It isolates the exact obstruction in the
smallest quartic where the desired second torsion direction appears.

## 7. Reciprocal-cubic and toric completion

Let `omega^2+omega+1=0`.  The invertible Fourier change

```text
U=A+omega B+omega^2 C,          V=A+omega^2 B+omega C,
Z=(A+B+C)/3                                                     (35)
```

satisfies

```text
D(U,V,Z)=9Delta_3,
D=U^2V^2-4(U^3+V^3)Z+18UVZ^2-27Z^4,                           (36)
```

and sends the bitangent `ell` to `Z=0`.  On the normalization, the
Mobius coordinate

```text
q=omega^2(p+omega)/(p+omega^2)                                 (37)
```

identifies `(U/Z,V/Z)` with `(2q+q^(-2),q^2+2q^(-1))`; the two roots of
the bitangent pullback become `q=0,infinity`.

In the affine chart `Z=1`, equation `(36)` is exactly

```text
D(U,V,1)=disc_T(T^3-UT^2+VT-1).                               (38)
```

Writing the displayed cubic as `f(T)`, its root-sheet algebra is

```text
k[U,V,T]/(f)=k[U,T,T^(-1)],       V=T^(-1)-T^2+UT.             (39)
```

Indeed `T(T^2-UT+V)=1`.  The algebra is finite free of rank three over
`k[U,V]`, but `T` is a nonconstant unit.  In `(U,T)` coordinates its
Jacobian is `-f_T/T`.  On the ramification divisor the repeated root `q`
gives precisely

```text
U=2q+q^(-2),                    V=q^2+2q^(-1).                 (40)
```

This map from `G_m` is injective.  If distinct `q,r` had equal images,
cancellation would give `2q^2r^2=q+r` and `qr(q+r)=2`; hence
`(qr)^3=1` and `(q+r)^2=4qr`, forcing `q=r`.  Its derivative vanishes
only at `q^3=1`, the three cusp addresses.  Its projective closure adds
`q=0,infinity` at the distinct smooth points `[1:0:0]` and `[0:1:0]`.
Thus the reciprocal cubic independently recovers the exact two-place
normalization debt.

The ordered-root completion makes the group action explicit.  Put

```text
z=(xy)^(-1),
U=x+y+z,             V=xy+x^(-1)+y^(-1),
W=(x-y)(y-z)(z-x).                                             (41)
```

Then

```text
W^2=D(U,V,1),
det d(U,V)/d(log x,log y)=-W.                                  (42)
```

The symmetric/alternating invariant theorem therefore gives the global
affine quotients

```text
(G_m)^2/S3=A2_(U,V),
(G_m)^2/A3=Spec k[U,V,W]/(W^2-D).                              (43)
```

The nontrivial `A3` fixed locus consists only of the three triple-root
points, so over the regular quadratic resolvent this realizes one genuine
cyclic Kummer direction.  It does not realize both independent classes in
`(5)` at once.

This is the ordinary ordered-root `S3` closure; it was therefore too strong
to say that no `S3` closure exists.  What fails is the design target.  The
root sheet has unit group `k^* T^Z` and the ordered-root torus has unit group
`k^* x^Z y^Z`.  Any map from `A2` sends those character units to constants,
so neither source admits a dominant plane atlas.  The construction is
globally monogenic and unit-rich, supplies no Keller pair, and does not show
that every other Kummer line in `(5)` has the same obstruction.

The obstruction persists on the entire fixed-nonzero-product marked-root
stratum.  For `c in k^*`, put

```text
f_c(T)=T^3-UT^2+VT-c,
D_c=U^2V^2-4V^3-4cU^3+18cUV-27c^2.                            (44)
```

Then

```text
disc_T(f_c)=D_c,
k[U,V,T]/(f_c)=k[U,T,T^(-1)],
V=cT^(-1)-T^2+UT.                                             (45)
```

Its discriminant normalization is the injective map

```text
q |-> (2q+cq^(-2),q^2+2cq^(-1)),             q in G_m,       (46)
```

with cusps exactly at `q^3=c` and missing places `q=0,infinity`.  Choosing
`gamma^3=c` and scaling

```text
(T,U,V)=(gamma tau,gamma u,gamma^2 v)                          (47)
```

reduces `(44)-(46)` to `(38)-(40)` and multiplies the discriminant by
`gamma^6`.  Likewise, with `xyz=c`, the ordered-root source is again
`(G_m)^2`, its Vandermonde squares to `D_c`, and its logarithmic Jacobian is
minus that Vandermonde.  Thus varying a fixed nonzero root product creates
no new modulus and never repairs the unit or two-place debt.

As a hostile cross-control, the displayed target curve of the provisional
inverse-discriminant lane has the independently checked parametrization

```text
D(t)=-t(5t+1)(17t-3)(19t+4),
C=-D(t),                         A=-tD(t).                     (48)
```

Its normalization is `A1`, with degrees `(deg A,deg C)=(5,4)` and one place
at infinity.  No affine reparametrization of the coefficient plane can turn
the fixed-product branch into this target: `k[q,q^(-1)]` and `k[t]` have
different unit groups.  This proves only the fixed-`c!=0` boundary of the
reserved all-degree marked-root lane; variable or zero product remains open.

## 8. Exact companions and audit boundary

Reproduce with

```bash
python3 04-computation/jc2_tricuspidal_quartic_rank_two_two_place_thm3851.py
python3 -O 04-computation/jc2_tricuspidal_quartic_rank_two_two_place_thm3851.py
python3 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_extension_thm3851.py
python3 -O 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_extension_thm3851.py
python3 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_independent_audit_thm3851.py
python3 -O 04-computation/jc2_tricuspidal_reciprocal_cubic_toric_independent_audit_thm3851.py
python3 04-computation/jc2_fixed_nonzero_product_marked_root_extension_thm3851.py
python3 -O 04-computation/jc2_fixed_nonzero_product_marked_root_extension_thm3851.py
```

Both modes byte-match the frozen output.  The assertion-free companion checks
the implicit quartic, normalization inverse, complete singular support and
three `A2` packets, the two distinct smooth bitangent contacts, the finite
and infinite fourth-power line systems, all 126 `E7` roots, all 120
boundary-orthogonal `A2` systems, all 40 compatible `3A2` configurations and
their Smith profiles, the explicit order-three relations, local rank-two
address map, and deck inversion.  It reports 763 active exact gates.

The 1,376-gate extension independently repeats the lattice census while
deriving `(35)-(43)`; the separate 48-gate audit rederives the toric bridge
without that census, and the 28-gate corollary checks `(44)-(48)`.  The
weak-del-Pezzo resolution and class-group quotient
`(20)-(22)` remain human geometric inputs.  The lattice census is exhaustive
finite exact evidence after that bridge, not a computational proof of it.
All normal, optimized, and frozen replays pass.
**QED.**
