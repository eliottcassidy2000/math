---
id: THM-4141
title: "Delta-D collision-wall boundary packet and monodromy exclusion"
status: >
  PROVED RELATIVE TO THM-3827/4053/4103/4120/4122/4138/4140 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED. On the live delta-only exact-M=8 Delta_D=0
  wall, the doubled outer-edge root is one smooth index-seven tangency, not
  a secondary stratum. The source has genus eight and infinity packet
  (7,7,2,2,1), forcing degree 15 or 19. THM-4140's critical length eighteen
  and seam-independent punctured monodromy exclude both degrees. Hence this
  entire wall is empty. The two-term collision wall, M>=9, other cells,
  entry, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-cycle6-63-20260825
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
  - THM-4140-delta-d-collision-wall-affine-critical-length-eighteen
related:
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
script: 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141.py
output: 05-knowledge/results/jc23_delta_d_collision_wall_boundary_monodromy_thm4141.out
independent_audit_script: 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_delta_d_collision_wall_boundary_monodromy_thm4141_independent_audit.out
script_sha256: 38d7e6881c4628093cf807f4d0881c8cd06cb6c6cab1ce1a76924f4e7708cc1d
output_sha256: 3eb50c8d27057fc8f9b0f3e8511bd3671e4ea03d908f4cfc3e30374ff01c0e31
independent_audit_script_sha256: e14278b2ce23fa7cb88073b1d61ff72b210ccaff127bf65c4126768ce5ee52a7
independent_audit_output_sha256: b13554d6cada9a342a20370229ded0d0c81b7844b0b45bec58dfd7ab7a41a00a
semantic_sha256: 400997aaf03937bcb4adc6fd82be47305677e3671c6ab60b50c8d57586ee7b9f
independent_semantic_sha256: 6a0a29b3b106ac8288ac9f06e2ae6c3c526cb0521b3397995e76b26ddda168bd
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy calculation rebuilds the complete Newton polygon,
  proves the generic source equation primitive and irreducible, resolves the
  repeated edge in the (u,v) chart, computes the residue order, freezes the
  response degrees and both monodromy budgets, and checks three sharp
  hostiles. Normal, optimized, and hash-seeded replays byte-match.
independent_audit: >
  ACCEPT. A standard-library sparse-Laurent referee imports no primary code,
  uses the reciprocal (z,w) chart, reconstructs the polygon and local rows
  after clearing the live denominator, audits all shared-pivot commutators
  and finite-carrier capacity cases, and obtains the same packet, degrees,
  and exclusions. Normal, optimized, and two hash-seeded replays byte-match.
---

# THM-4141 -- the `Delta_D` wall is empty

**PROVED RELATIVE TO THM-3827/4053/4103/4120/4122/4138/4140 +
VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.** Work over
`C` and retain the exact reduced `(2,3)` cell and maximum-weight-eight
hypotheses of THM-4053.

## 1. Theorem and inheritance

> **Theorem.** No nonautomorphic planar Keller pair lies on the live
> delta-only collision wall
>
> ```text
> theta=0,       delta*kappa!=0,       Delta_D=phi^2-4kappa delta=0. (1)
> ```

The inheritance pass is deliberately collision-sensitive:

- the closest proved mechanism is THM-4138's fixed-sheet versus
  orbit-merger obstruction;
- the canonical hostile is THM-4134, where a repeated edge loses critical
  roots and creates a secondary wall;
- the corrected near miss is specialization of a smooth-edge ramification
  packet across a doubled root;
- the least-used sidecar is THM-4140's square top row and exact affine
  critical length `18`.

The new mechanism is that the forced cubic coefficient resolves this doubled
edge in the transverse direction. Two off-wall index-four punctures become
one smooth index-seven puncture. There is no analogue of THM-4134's secondary
node.

## 2. Exact generic source and connectedness

Use THM-4140's normalized variables

```text
t=p-s^2,
G=-s^2/(2t)+H(p,sp),

H=-3p+(8/3)p^2-(1376/135)p^3
  +K s^2p^2+Phi sp^3+Delta p^4.                         (2)
```

The live wall has the complete parameterization

```text
Delta=5696/[15(6r^2+7)],
K=Delta r^2,                 Phi=2Delta r,
r!=0,                        6r^2+7!=0.                 (3)
```

For a generic pencil value put `Q=q^-1`. Multiplying `G-q` by `Qt` gives

```text
F_Q(s,p)=(s^2-p)(1-QH)-Q s^2/2.                         (4)
```

Write `F_Q=a_0+Qa_1`. Then

```text
a_0=s^2-p,
a_1=-(s^2-p)H-s^2/2,
a_1=-s^2/2                 modulo (s^2-p).               (5)
```

Thus `gcd(a_0,a_1)=1`. Equation `(4)` is a primitive polynomial of degree one
in `Q`, so Gauss's lemma makes its total `Q`-equation irreducible after any
live coefficient specialization. This useful audit does not by itself prove
geometric connectedness after algebraic base extension.

For that stronger conclusion, use the closed-polynomial factor theorem cited
and audited in THM-3827: if a polynomial `J:C^2->C` has disconnected
geometric generic fibre in characteristic zero, then

```text
J=P(J_0)                         with deg(P)>1.           (5a)
```

Apply this to the original polynomial `J=E(A,C)`. If `(5a)` held, a root
`c` of `P'` would make the nonempty curve `J_0=c` lie in `Crit(J)`. But

```text
grad(E o F)=DF^t grad(E),                               (5b)
```

`DF` is invertible, and the target polynomial `E` has only its two nodal
critical points. Their inverse images are finite; THM-4140 computes the
total length as eighteen. This contradiction proves that the geometric
generic source fibre is connected for every live `r`. It supplies the
transitivity premise below without importing a smooth-seam claim.

## 3. Complete boundary and the doubled edge

After coincident monomials are combined, the Newton polygon of `(4)` is

```text
conv{(0,1),(2,0),(4,2),(2,4),(0,5)}.                    (6)
```

Deleting the optional interior point `(2,3)` on the special support wall
`K=epsilon` does not change `(6)`. Its edge ledger, in counterclockwise
order, is

| edge | inward normal | length | residue distance | ordinary index |
|---|---:|---:|---:|---:|
| `AB:(0,1)--(2,0)` | `(1,2)` | 1 | 1 | 1 |
| `BC:(2,0)--(4,2)` | `(-1,1)` | 2 | 2 | 2,2 |
| `CD:(4,2)--(2,4)` | `(-1,-1)` | 2 | 4 | 4,4 |
| `DE:(2,4)--(0,5)` | `(-1,-2)` | 1 | 7 | 7 |
| `EA:(0,5)--(0,1)` | `(1,0)` | 4 | 1 | 1,1,1,1 |

Here the last four roots have `s=0`, hence `(x,t)=(0,p)` with `p!=0`; they
are affine source points, not punctures at source infinity. Pick's theorem
gives

```text
2 Area=24,              boundary steps=10,
arithmetic genus=(24-10+2)/2=8.                         (7)
```

The only nonsquarefree edge on `(1)` is `CD`:

```text
K+Phi W+Delta W^2=Delta(W+r)^2.                         (8)
```

It is essential to resolve `(8)` in the complete curve. Put

```text
s=1/u,                    p=(v-r)/u,
L(u,v)=u^6 F_Q(1/u,(v-r)/u).                            (9)
```

Exact expansion gives

```text
L(0,v)=-Delta Q v^2(v-r)^2,
L_u(0,0)=-1376 Q r^3/135,
[v^2]L=-Delta Q r^2.                                   (10)
```

Both displayed coefficients are nonzero on `(3)`. The collision point is
therefore smooth, and its branch begins

```text
u=-(135Delta/(1376r))v^2+O(v^3).                        (11)
```

The Keller residue identity from THM-4103 is

```text
varphi_q^*(dU/(2V))=Q ds/(F_Q)_p.                       (12)
```

In `(9)`, `(F_Q)_p=u^-5L_v` and `ds=-u^-2du`. Substituting `(11)` into
`(12)` gives a nonzero constant times

```text
v^6 dv.                                                 (13)
```

Thus the doubled root is one puncture of ramification index seven. It
replaces the two off-wall index-four punctures without changing total
defect. All other boundary roots are simple by the inherited complete-model
gate; `(10)` proves the only missing local case. The projective source is
smooth of genus eight and its complete infinity packet is

```text
(7,7,2,2,1),                    total defect=14=2*8-2.  (14)
```

Consequently `(14)` is the entire ramification divisor of the generic map
to the elliptic target. In particular, no hidden affine branch point can
enter either monodromy budget.

## 4. Exact target response: only degrees 15 and 19

The `AB`, collided `CD`, and `DE` punctures are `C(q)`-rational and have
indices

```text
1,7,7.                                                   (15)
```

THM-4120's target-only calculation

```text
E_q(C(q))={O}                                            (16)
```

forces all three to the target origin. The `BC` edge is instead

```text
(1-Q/2)-QK W^2=0,
K W^2=q-1/2.                                             (17)
```

The valuation of `q-1/2` at `q=1/2` is odd, so `(17)` is one irreducible
quadratic closed point over `C(q)`. Its two geometric index-two punctures
respond together. No affine source point maps to projective target infinity.
Therefore the mapping degree is exactly

```text
n=15  if BC has finite image,
n=19  if BC maps to O.                                  (18)
```

This is an exhaustive response dichotomy, not a degree upper bound.

### 4.1 The finite `BC` carrier satisfies the earlier transport gates

When `n=15`, choose a constant square root of `K` and put `rho=sqrt(K)W`.
Equation `(17)` becomes

```text
q=1/2+rho^2.                                             (19)
```

This is exactly the quadratic target base change audited in THM-4138, after
the harmless target scaling used there. The proof is source-wall independent:
THM-4122 makes the normalization of a horizontal nonproperness component
`A1`; degree factorization and `(16)` force its degree over the `q`-line to
be two; and the rank-one Mordell--Weil calculation on

```text
y^2=x^3-3x+2+rho^2                                      (20)
```

has polynomial sections `+-P,+-2P,+-3P`. The first two pairs have coordinate
degree pair `(0,1)`, whereas `+-3P` has `(2,3)`. THM-4122 requires a positive
intrinsic pole pair `(2rho_C,3rho_C)`, so only the `+-3P` nodal carrier survives.
Thus the two generic finite images are the two marked points `Q_+,Q_-` of
THM-4138. This rechecks the inputs to its puncture-avoiding construction; it
does not apply that theorem's `Delta_V`-scoped conclusion by analogy.

## 5. Seam-independent fixed-sheet transport

Let `o_0,o_1` be the two target nodes and

```text
r_i=#F^(-1)(o_i) in A2.                                  (21)
```

THM-4140 proves

```text
r_0+r_1=18,                                             (22)
```

and its two universal Morse pairs show `r_0,r_1>=2`.

For each affine inverse of `o_i`, Keller etaleness identifies target and
source Morse coordinates:

```text
E-q_i=uv.                                                (23)
```

A core circle in a nearby Milnor annulus therefore has one closed
degree-one lift in each inverse neighborhood. Choose it parallel to avoid
the finite moving `BC` branch set when that set is present.

Resolve the compactified graph of the pencil map and delete the finite
additional Hurwitz discriminant from the `q`-line. After this shrinking the
morphism between the proper smooth curve families is quasifinite, hence
finite. Deleting the complete marked branch divisor `O` and, for `n=15`,
`Q_+,Q_-`, makes it finite etale. Ehresmann transport along a path to a
common reference fibre therefore carries the core and all its distinct
closed lifts. Sheet labels are changed by conjugation; no branch permutation
is postmultiplied. If a carrier meets a node in the central fibre, a nearby
parallel core still avoids its two annular intersections. This is precisely
the local and proper transport argument of THM-4138, with its hypotheses
re-established by `(14),(19),(20)`.

Let `X,Y` be the two transported vanishing permutations. Then

```text
#Fix(X)>=r_0,                    #Fix(Y)>=r_1.            (24)
```

Removing finitely many points from the connected projective source does not
disconnect it. Hence the unramified punctured cover is connected and its
monodromy is transitive.

## 6. The full-boundary degree `19` is impossible

Here every puncture maps to `O`, so the once-punctured torus group is
generated by `X,Y`, and the `O`-meridian is their commutator up to inversion.
Equations `(22),(24)` give

```text
|supp(X)|+|supp(Y)|<=2*19-18=20.                         (25)
```

Transitivity forces the two supports to cover all nineteen letters. If one
of `X,Y` is the identity, the other must be a nineteen-cycle and the
commutator is the identity, contradicting `(14)`. If both are nonidentity,
their supports must meet in exactly one pivot. Any second nontrivial cycle
of either permutation would form a separate orbit, so each permutation is
one cycle on its support. Two cycles meeting in one pivot have commutator
cycle type

```text
(3,1^16).                                                (26)
```

The actual meridian type is `(7,7,2,2,1)`, not `(26)`. This contradiction
excludes `n=19`.

## 7. The finite-carrier degree `15` is impossible

Now the origin packet is `(7,7,1)`. The two `BC` points have finite images;
their meridians are two transpositions, or one product of two disjoint
transpositions if the images collide. In either case their total permutation
index is two.

The punctured torus is generated by `X,Y` and those carrier meridians.
Equations `(22),(24)` give

```text
|supp(X)|+|supp(Y)|<=2*15-18=12.                         (27)
```

For a nonidentity permutation `sigma`,

```text
ind(sigma)=15-#Cycles(sigma)<=|supp(sigma)|-1.           (28)
```

Thus, if at least one of `X,Y` is nonidentity,

```text
ind(X)+ind(Y)<=11.                                      (29)
```

Adding the carrier index gives at most `13`, below the `14` edges needed to
merge fifteen singleton sheet orbits into one. If both are identities, the
total index is only two. Either way the generated action is not transitive,
contradicting the connected cover. Hence `n=15` is impossible.

## 8. Sharp controls, consequence, and replay

Three hostile controls mark the exact failure boundary.

1. Deleting the forced `-(1376/135)p^3` row makes `L_u(0,0)=0`; the smooth
   tangency proof genuinely depends on the complete lower model.
2. Replacing critical length `18` by `17` permits a thirteen-cycle and two
   attaching transpositions to generate all fifteen sheets with total index
   `14`.
3. Keeping length `18` but allowing carrier index three likewise reaches the
   transitivity threshold. With the actual two transpositions it misses by
   exactly one.

Both degree alternatives in `(18)` are impossible. Therefore

```text
the delta-only exact-M=8 Delta_D=0 wall is empty.         (30)
```

Together with THM-4130 and THM-4134/4138, this leaves only the two-term wall
`delta*theta!=0, delta+theta=0` inside THM-4053's exact maximum-eight
trichotomy. It does not treat that wall, `M>=9`, another reduced cell, chart
entry, `JC(2)`, or `DC(2)`.

Replay with

```text
python3 -B 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141.py
python3 -B -O 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141.py
PYTHONHASHSEED=29 python3 -B 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141.py
python3 -B 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141_independent_audit.py
python3 -B -O 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141_independent_audit.py
PYTHONHASHSEED=123 python3 -B 04-computation/jc23_delta_d_collision_wall_boundary_monodromy_thm4141_independent_audit.py
```

All streams byte-match their corresponding frozen outputs. **QED.**
