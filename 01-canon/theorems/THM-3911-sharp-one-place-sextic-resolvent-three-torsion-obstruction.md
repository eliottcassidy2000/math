---
id: THM-3911
title: "Sharp one-place sextic quadratic-resolvent three-torsion obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The sharp
  rational one-place sextic control from THM-3906 has finite singularity
  packet one ordinary quadruple point, two A2 cusps, and two A1 nodes.  On a
  resolution of its quadratic double plane, the split infinity boundary has
  Gram matrix [[-2,3],[3,-2]], of determinant -5.  The complete removed-
  divisor lattice is 3-saturated in the ambient integral Picard lattice, so
  the affine quadratic-resolvent class group has no 3-torsion and its units
  are scalar.  Consequently this sharp plane control cannot be the branch of
  a normal irreducible S3 cubic algebra.  Other one-place sextics, higher
  degree branches, nonnormal orders, and JC(2) remain open.
source: root / post-THM-3906 sharp sextic resolvent audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS after one proof-completeness repair.  The
  audit independently checked the hidden-node exhaustion, the double-cover
  canonical and boundary formulas, the complete exceptional/boundary
  quotient, the mod-three residue calculation, scalar units, and the
  Kummer/purity passage.  The original candidate called Pic(S) an integral
  lattice without excluding Picard torsion; the repaired proof computes
  q(S)=0 and kappa(S)=-infinity, hence S is rational and Pic(S) is genuinely
  torsion-free.  Normal and optimized companion runs byte-match the frozen
  output, raw hashes agree, and agent documentation passes.
depends_on:
  - THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
  - THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy
script: 04-computation/jc2_sharp_one_place_sextic_resolvent_three_torsion_thm3911.py
output: 05-knowledge/results/jc2_sharp_one_place_sextic_resolvent_three_torsion_thm3911.out
script_sha256: 117248d6eb4c56d0566f4c66930cbed4c9d583e248ad3b68b7a90af80053cf2c
output_sha256: 3974383d1d1c05bc9901fe3bb43dfb4a8d3c7c845c52c7e1b456ffa26612599e
semantic_sha256: 637b640a6ea58c4b80abacde95d3135fc4ac069d5a72ea31bff6fa53125d4291
hash_basis: raw LF bytes
---

# THM-3911 -- the sharp sextic has the wrong boundary arithmetic for a cubic lift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Put

```text
G=t^4+t^3-2t-5,
A=tG,                         C=(t^2+1)G.                    (1)
```

The image is the rational sextic `Gamma=V(delta)`, where

```text
delta=-A^6+6A^5-25A^4-6A^3C^2-32A^3C
      -5A^2C^3-18A^2C^2+AC^4+7AC^3+C^5+5C^4.              (2)
```

It is the sharp one-place plane control of THM-3906.  Let

```text
Q=Spec k[A,C,W]/(W^2-delta).                                (3)
```

Then

```text
Q*=k*,                         Cl(Q)[3]=0.                  (4)
```

In particular there is no normal finite-flat cubic `k[A,C]`-algebra whose
generic fibre is an irreducible `S3` cubic field and whose discriminant
divisor is `V(delta)` with generic exponent one.  Thus the sharp plane
geometry cannot be inverse-realized as the required normal cubic order.

## 1. The normalization has two hidden nodes

Eliminating `t` from `(1)` gives exactly `(2)`.  Conversely,

```text
t=-(A^3-3A^2-AC^2-6AC)/(A^2C+4A^2+AC+C^2)                 (5)
```

on a dense open.  The homogeneous normalization morphism is

```text
[T:S] |-> [TS Gh : (T^2+S^2)Gh : S^6],
Gh=T^4+T^3S-2TS^3-5S^4.                                    (6)
```

It has no base point and `(5)` makes it birational.  Hence `delta` is
irreducible and `(6)` is the normalization.

The four simple roots of `G` map to the origin.  Their tangents are distinct:
two such slopes could agree only at reciprocal parameters, while

```text
gcd(G(t),t^4G(1/t))=1.                                     (7)
```

Thus the origin is an ordinary quadruple point.  The only critical
parameters away from that fibre satisfy

```text
gcd(A'(t),C'(t))=t^2-1.                                    (8)
```

They give `A2` cusps

```text
t=1  |-> (-5,-10),                 t=-1 |-> (3,-6).         (9)
```

There are also two nodes which are easy to miss if one only counts critical
parameters.  Since `C/A=t+t^{-1}` away from `A=0`, distinct equal-image
parameters must be reciprocal.  The exact collision factor is

```text
H=t^4+t^3+3t^2+t+1,                                        (10)
```

whose four simple roots form two reciprocal pairs.  Their two targets are

```text
C=-3,                         A^2-3A+9=0.                   (11)
```

The two tangent slopes in each pair are distinct.  An independent singular-
ideal elimination has projected support

```text
C^5(C+3)(C+6)^2(C+10)^2,                                   (12)
```

and the four fibre gcds recover exactly the origin, the two points `(11)`,
and the two cusps `(9)`.  Thus no finite singularity is omitted.  Finally,
the degree-six leading form is `-A^6`, so the projective branch has one point
at infinity, `[0:1:0]`; `partial_Z(delta_h)=1` there, hence it is smooth.
The genus ledger is consequently

```text
p_a(Gamma)=10=6 (ordinary quadruple)+2 (two A2)+2 (two A1). (13)
```

## 2. The removed-divisor lattice

Blow up the ordinary quadruple point in `P2`; write `H` for the pullback of a
line and `E` for the exceptional curve.  The strict branch has class

```text
6H-4E=2(3H-2E).                                             (14)
```

Take the double cover with branch `(14)` and resolve the two `A2` and two
`A1` rational double points.  Denote the resulting smooth projective surface
by `S`.  The inverse image `D` of `E` is a double cover of `E=P1` ramified at
the four distinct tangent directions.  Therefore `D` is a smooth genus-one
curve and

```text
D^2=2E^2=-2.                                                 (15)
```

Moreover the double-cover canonical formula, unchanged by the crepant ADE
resolutions, gives

```text
K_S=pi^*(K_Bl_p(P2)+3H-2E)=-D.                              (16)
```

There is one necessary torsion check before using intersection theory as a
Picard-lattice argument.  Let `Y=Bl_p(P2)`, let `L_0=3H-2E`, and let `X` be
the normal double cover of `Y` before resolving its ADE points.  The standard
double-cover decomposition is

```text
pi_* O_X = O_Y direct-sum O_Y(-L_0).                       (16a)
```

Serre duality and `K_Y+L_0=-E` give

```text
H^1(Y,O_Y(-L_0)) dual=H^1(Y,O_Y(-E))=0.                   (16b)
```

The last equality follows from
`0 -> O_Y(-E) -> O_Y -> O_E -> 0`, because the map on constants is an
isomorphism.  Thus `q(X)=0`; the rational ADE resolutions preserve this, so
`q(S)=0`.  On the other hand `(16)` and an ample divisor show that no
positive multiple of `K_S=-D` is effective.  Hence `kappa(S)=-infinity`.
The ruled-surface classification now makes `S` rational (the base genus is
`q(S)=0`).  In particular

```text
Pic(S) is a torsion-free integral intersection lattice.                 (16c)
```

This is the bridge that prevents invisible Picard three-torsion from
evading the square calculation in Section 3.

The preimage of the line at infinity splits as two rational curves `B_+`
and `B_-`, because the restriction of the branch equation is the square
`-A^6`.  They are disjoint from `D` and the finite exceptional curves.
Adjunction with `(16)` and `(B_++B_-)^2=2H^2=2` gives

```text
B_+^2=B_-^2=-2,                       B_+ B_-=3.             (17)
```

The divisors removed when `S` is mapped to the affine normal surface `Q`
are:

```text
D;  B_+,B_-;  two A2 chains;  two A1 curves.                (18)
```

Different blocks in `(18)` are disjoint.  The affine double plane `Q` is
normal: it is a hypersurface, and the singularity packet in Section 1 makes
its singular locus finite, so `S2+R1` applies.  The standard
resolution/class-group exact sequence therefore has relation lattice `R`
with Gram blocks

```text
[-2],  [-2  3],  two copies of [-2  1],  two copies of [-2],
       [ 3 -2]                [ 1 -2]
                                                                  (19)
```

and determinant `360`, and gives

```text
Cl(Q)=Pic(S)/R.                                               (20)
```

## 3. Why the two cusp classes cannot globalize

Suppose `x in Pic(S)` has `3x in R`.  Pairing `3x` against the generators in
`(19)`, modulo three, forces the coefficients of `D`, both boundary curves,
and both `A1` curves to vanish modulo three: all those Gram blocks are
nonsingular modulo three.  In either `A2` block the only radical direction is
the difference `v_i` of its two simple roots.  After subtracting an element
of `R` from `x`, one therefore has

```text
3x=a v_1+b v_2,                  a,b in {0,+1,-1}.           (21)
```

The two blocks are orthogonal and `v_i^2=-6`.  Since `Pic(S)` is an integral
lattice, `(21)` would require

```text
x^2=-2(a^2+b^2)/3                                                (22)
```

to be an integer.  For a nonzero pair `(a,b)`, it is not.  Hence `a=b=0`;
the adjusted `x` is then killed by three, so torsion-freeness `(16c)` makes
it zero.  Therefore the original `x` lies in `R`.  The sublattice `R` is
3-saturated in `Pic(S)`, and `(20)`
proves

```text
Cl(Q)[3]=0.                                                    (23)
```

This is the structural difference from THM-3844.  There the split quartic
boundary has determinant three and can absorb a diagonal cusp three-class.
Here the sextic boundary determinant is five, so the same two local `A2`
directions have nowhere integral to glue.

## 4. Units and the cubic obstruction

A unit of `Q` pulls back to an invertible function above the affine surface,
so its divisor on `S` is supported only on `B_+` and `B_-`.  If it is
`mB_++nB_-`, intersection with the two boundary curves gives

```text
-2m+3n=0,                         3m-2n=0.                   (24)
```

The determinant is `-5`, so `m=n=0`; a rational function without a pole on
projective `S` is constant.  This proves `Q*=k*`.

Let `Q_reg` be the regular locus.  Normality and codimension two give

```text
Gamma(Q_reg,O)^*=k*,                    Pic(Q_reg)=Cl(Q).     (25)
```

The Kummer sequence and algebraic closedness of `k` now imply

```text
H^1(Q_reg,mu_3)=0.                                             (26)
```

If the normal irreducible cubic algebra excluded in the statement existed,
its Galois closure after the quadratic-resolvent base change would be a
connected cyclic cubic cover.  Generic discriminant exponent one means the
quadratic base change absorbs the transposition inertia; the cyclic layer is
unramified in codimension one.  Purity over the regular surface makes its
restriction a nontrivial `mu_3`-torsor on `Q_reg`, contradicting `(26)`.

## 5. Scope and replay

This theorem closes the particular sharp curve `(1)`, not all one-place
sextics.  It does not exclude a different singularity packet, a nonnormal
cubic order, a non-common-zero coefficient grammar, a positive-genus branch,
or degree at least seven.  It constructs no Keller map; `JC(2)` remains
**OPEN**.

Reproduce the exact algebra and integral-lattice packet with

```bash
python3 04-computation/jc2_sharp_one_place_sextic_resolvent_three_torsion_thm3911.py
python3 -O 04-computation/jc2_sharp_one_place_sextic_resolvent_three_torsion_thm3911.py
```

Both runs must byte-match the frozen output.  The companion checks the
eliminant, rational inverse, projective normalization, address collisions,
critical parameters, complete singular support, local types, smooth unique
infinity, boundary contacts, the full Gram matrix, and every mod-three
relation-lattice residue.  The double-cover resolution, divisor quotient
`(20)`, rationality/torsion bridge `(16a-c)`, and Kummer/purity passage have
been independently hostile-audited.  **QED.**
