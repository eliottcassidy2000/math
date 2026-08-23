# JC(2): the one-place/C3/nonmonogenic junction after THM-3879--3885

**Status: CURRENT RESEARCH SYNTHESIS, 2026-08-23.  JC(2) remains OPEN.**
This note records the design object exposed by the proved theorems, not a
counterexample and not a theorem dependency.

## 1. The target is an intersection of contracts

The recent cubic and branch experiments should no longer be ranked by a
single feature such as rationality, a torus equation, or local cusp torsion.
A degree-three finite completion relevant to a planar Keller map must meet all
of the following contracts at once:

1. a reduced simple branch component with affine normalization `A1` and the
   correct common point at infinity;
2. a connected `S3` field with a genuine unramified `C3` layer over its
   quadratic resolvent;
3. a normal finite-flat cubic order which is globally **nonmonogenic**;
4. after deleting the ramified sheet, constant units and exactly one visible
   reduced companion over each branch component;
5. an actual dominant etale plane atlas, including the determinant/density
   cocycle forgotten by coefficient quotients.

The three closest laboratories each miss a different coordinate:

```text
THM-3874 branch: one-place, but its resolvent has no C3 packet;
THM-3844 branch: one-place and a C3 packet, but the cubic is monogenic;
THM-3879 branch: a global C3 packet, but the cubic is monogenic and the
                 normalization loses two places at infinity.
```

The next construction must not rediscover any one of these features in
isolation.  It must change the coordinate on which the corresponding no-go
acts.

## 2. Duality now has an exact place ledger

THM-3882 reframes a line on a rational dual curve as the Wronskian of a point
projection.  If `D_P` is the projection base fibre and `phi_P` the resolved
map to `P1`, then

```text
(dual normalization)^* H_P = 2D_P + Ram(phi_P).
```

This explains THM-3879's best line structurally:

```text
two node addresses * two base units + one ramification unit = (3,3).
```

Riemann--Hurwitz makes one-point support impossible for the dual of every
immersed rational curve of degree at least three.  Pluecker then removes the
whole sextic `6A2+4A1` equisingularity packet from the one-place search.  A
different contact conic or torus presentation inside that packet cannot help.

The honest exits are now typed:

- make the primal normalization nonimmersive, and account for the common
  tangent factor cancelled from the dual map;
- change the singularity packet, hence the Pluecker and Kummer ledgers; or
- construct the one-place branch directly rather than as an immersed dual.

This is a loss ledger, not a reason to abandon the `C3` packet itself.

## 3. Marked-root descent is no longer the vague part

THM-3880 and THM-3883 completely isolate affine descent once a regular
normalization root

```text
z^2=1+(2/3)AC
```

exists.  For a fixed global sign, the tests are exactly:

```text
A=0:      z has the matching sign;
z!=0:     all branch residues agree (then Hensel uniqueness descends z);
z=0:      z^3 belongs to the completed local branch ring.
```

There is no unbounded mystery family of conductor jets.  The only higher
singularity sidecar is the square/cube gap.  The sharp local controls are

```text
k[[t^2,t^3]]: z=t has z^2,z^3 in the ring       (positive),
k[[t^2,t^5]]: z=t has z^2 in but z^3 out        (hostile).
```

The surviving global difficulty is therefore upstream or downstream:
produce the normalization square root without forbidden sign gluing, and
then satisfy the projective companion/finite-order geometry.

A useful nonbirational-root hostile is

```text
z=t^3,
A=t^2+t+1,
C=(3/2)(t-1)(t^3+1).
```

The two nontrivial normalization addresses with the same `(A,C)` have
opposite nonzero `z` values, so the forced polynomial carrier fails exactly
by THM-3883.  By contrast

```text
A=t^2,
C=3t+(3/2)t^4,
z=1+t^3
```

does carry the polynomial

```text
B=(4C^2-6A^2C-9A)/18,
```

but the shear `C -> C-(3/2)A^2` exposes the already closed parabola.  This is
a positive descent control, not a new branch architecture.

## 4. The three-cusp polynomial-coefficient debt is now a filtered norm

THM-3881 replaces three arbitrary cusp-ideal coefficient polynomials and a
hidden lift ambiguity by the pair

```text
(T,f),                 T(0,0)=4f(0,0),
```

with matrix determinant `Delta` and exact sidecars

```text
C_side=L^2f,                    B_side=Pf^2-T^2.
```

The full `T=0` lane has only the base `f=0`.  THM-3884 shows that a square
survivor with nonconstant `f` must satisfy

```text
deg T >= deg f+1.
```

At equality its leading pair lies in the `(K,-a)` lift-gauge direction.
This is an associated-graded symbol, not yet a lift: peeling it changes the
mixed lift by a `Delta` multiple and creates a lower residual debt.  The
incoming sextic filtered-kernel scout has exactly the same failure anatomy:
top arms can vanish while omitted lower buckets decide liftability.

THM-3885 gives the complementary `f=0` arm packet: the `a=0` restriction is
zero or one of three constants with cube `-625/32`, and the `L=0` restriction
has a finite root-polarization grammar.  Low total degrees close; the live
problem is nonlinear interpolation between the two arms.

Thus the residual search should be organized as a Rees/strict-transform
problem, retaining the lower debt after every gauge peel.  A leading gauge
direction alone is not a recursive proof.

## 5. The next genuinely different object is the binary cubic index form

THM-3801 says that the finite-flat cubic order of a constant-unit etale open
must be nonmonogenic.  With a trace-zero basis, its intrinsic generator test
is the binary cubic index form

```text
I(X,Y)=aX^3+bX^2Y+cXY^2+dY^3.
```

The order is monogenic exactly when `I` represents a scalar unit over
`k[A,C]`.  This determinant sidecar is the order-theoretic analogue of the
incoming Russell density cocycle: quotient coefficients can agree while the
actual determinant predicate does not descend.

The easiest sufficient nonmonogenic gate is a common zero of `a,b,c,d`.
It has a useful tariff.  The binary discriminant is quartic in those four
coefficients, so it lies in the fourth power of that maximal ideal.  Any
reduced discriminant curve then has multiplicity at least four and local
delta invariant at least six.  A rational plane discriminant using this gate
must consequently have degree at least five.  In particular the one-place
quartic of THM-3844 could never have obtained nonmonogenicity from this cheap
common-zero certificate.

This gate is only sufficient.  A better construction may use an index form
whose coefficients generate the unit ideal but which nevertheless represents
no unit.  That arithmetic/geometric distinction should be retained in every
search.

## 6. Ranked construction tests

1. **Binary-order anchor.** Search small binary cubic forms for a squarefree
   irreducible one-place discriminant, an irreducible `S3` generic algebra,
   and a proved no-unit index form.  Audit normalization places before class
   groups or Kummer claims.
2. **Nonimmersed-dual niche.** Enumerate the first primal singularities whose
   cancelled tangent factor can pay the Wronskian one-place debt, then recompute
   the dual degree and the full singularity packet.  Do not assume the sextic
   torus equation survives.
3. **Residual wildcard.** Continue the THM-3884 equality seam one filtered
   layer at a time, or attack the `f=0` nonlinear arm interpolation.  Record
   the residual created by every gauge peel.
4. **Sextic normal-strip cross-check.** Treat the fifteen nonprincipal pole
   cones and higher shift jets as separate Rees charts.  The principal arm
   quotient is not the full sextic theorem.

The prize remains a complete object satisfying all five contracts in Section
1.  None is presently known, and no result in this note changes the OPEN
status of `JC(2)`.
