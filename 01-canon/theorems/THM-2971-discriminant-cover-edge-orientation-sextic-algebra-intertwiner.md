---
id: THM-2971
title: "Discriminant-cover edge/orientation sextic algebra intertwiner"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For a full-S4 depressed quartic on q*Delta*J_or nonzero, the edge and
  oriented-cycle sextic algebras are generically nonisomorphic over the base
  field but become explicitly isomorphic over the discriminant quadratic
  cover.  In the exact common THM-2968 gauge the map is
  Z -> sqrt(Delta) A(Y), where E'(Y)A(Y)=-8q mod E.  Writing
  A(Y)=YB(Y^2), the cubic-base Kummer ratio is exactly
  F(U)=Delta*U*B(U)^2 mod S.  Conjugating the
  discriminant orientation negates Z.  The construction is ramified at
  Delta=0 and loses primitivity at J_or=0, so it supplies neither a Keller
  layer nor an SFC(4)/JC(2) closure.
source: codex-discriminant-cover-sextic-algebra-intertwiner-2026-07-30
audit: >
  Independent hostile audit rederived both root-level signs, the cubic-base
  Kummer identity, the actual rank-six algebra isomorphism and its
  2^30*q^5*J_or^2*sqrt(Delta) power-basis determinant, the nonconjugate
  V4_edge/C4 base-field obstruction, the identical embedded C2 stabilizer
  over A4, and every q/Delta/J_or and Keller-scope boundary.  Normal and
  optimized raw transcripts match each other; both LF-normalize to the stored
  output, and the declared LF hashes were independently recomputed.
depends_on:
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
related:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
  - THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary
  - THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary
script: 04-computation/quartic_discriminant_cover_sextic_intertwiner_thm2971.py
output: 05-knowledge/results/quartic_discriminant_cover_sextic_intertwiner_thm2971.out
script_sha256: 6ca7c7f5695876dc38dc0a9423b4c5f81d648c29a5af90cb11b287ea0169b693
output_sha256: c9c41a197031b54f12c1bda7661cdadd1a97bf5149655d718575b6f9665dba49
hash_basis: LF-normalized bytes
---

# THM-2971 -- discriminant-cover edge/orientation sextic algebra intertwiner

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

[THM-2968, edge and oriented-cycle `S4`
complements](THM-2968-quartic-edge-and-oriented-cycle-s4-complements.md)
identifies two inequivalent six-sheet actions of a full-`S4` quartic.  Their
point stabilizers are respectively an edge stabilizer `V4_edge` and an
oriented-cycle stabilizer `C4`.  The two actions have the same matching
quotient and differ by the central sign graph.  That theorem is deliberately
only an action-level bridge.

This theorem repairs its multiplication gap, but only after adjoining the
quartic discriminant orientation.

Let `K` be a characteristic-zero field and let

```text
f(X)=X^4+pX^2+qX+r                                      (1)
```

have full Galois group `S4` over `K`.  Write `Delta=Disc(f)` and

```text
S(U)=U^3+2pU^2+(p^2-4r)U-q^2,                           (2)
E(Y)=S(Y^2),                                             (3)
F(U)=16(16r-4pU-3U^2)(U^2+2pU+p^2-4r),                 (4)
O(Z)=Res_U(S(U),Z^2-F(U)),                              (5)
J_or=1024r^3+768p^2r^2-288pq^2r+27q^4.                 (6)
```

Thus `E` is the edge sextic and `O` is the oriented-cycle sextic in the
common gauge fixed by THM-2968 and [THM-2864, the two polynomial
resolvents](THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product.md).
Assume

```text
q Delta J_or != 0.                                      (7)
```

Put `D=K(v)`, where `v^2=Delta`, with the Vandermonde orientation inherited
from the ordered roots.  Since

```text
Disc(E)=64q^2 Delta^2,                                  (8)
```

`E'` is invertible in `K[Y]/(E)`.  Let `A(Y)` be the unique class of degree at
most five satisfying

```text
E'(Y) A(Y) = -8q                 mod E(Y).               (9)
```

Then

```text
Phi: D[Z]/(O)  ->  D[Y]/(E),
     Z         |-> T(Y)=v A(Y)                          (10)
```

is an isomorphism of rank-six `D`-algebras.  Equivalently,

```text
T(Y)^2 = F(Y^2)                 mod E(Y),                (11)
O(T(Y)) = 0                     mod E(Y).                (12)
```

The class `A` is odd, so there is a degree-at-most-two class `B(U)` in the
matching cubic algebra `K[U]/(S)` with `A(Y)=YB(Y^2)`.  It satisfies the
stronger Kummer identity

```text
F(U)=Delta U B(U)^2                 mod S(U).           (12a)
```

Thus the two sextics are the binary fibres `Y^2=U` and `Z^2=F(U)` over the
same ternary matching-cubic base, and their Kummer square-class ratio is
exactly `Delta`.  In the ordered power bases used in `(10)`, the coordinate
change has the sharp determinant

```text
det(Phi)=2^30 q^5 J_or^2 v.                            (12b)
```

Thus `J_or=0` is precisely the primitive-coordinate Jacobian wall.  The
nontrivial element of `Gal(D/K)` sends

```text
v -> -v,             T -> -T,                           (13)
```

so it exchanges the two inverse orientations over each matching.  Over `K`
the two generic sextic algebras are not isomorphic; over `D` they are the same
algebra in two explicit primitive coordinates.

## 2. The root formula

Let the ordered roots of `(1)` be `a_0,a_1,a_2,a_3`, with sum zero, and put

```text
v=prod_(i<j)(a_i-a_j).                                  (14)
```

For the first edge in THM-2968's common gauge set

```text
s=a_0+a_1,
d_1=a_0-a_1,              d_2=a_2-a_3,
Omega=d_1 d_2(d_1^2-d_2^2).                            (15)
```

The three roots of `S` are the squares of the three matching sums.  Direct
factorization gives

```text
v=d_1 d_2 S'(s^2),
s(d_1^2-d_2^2)=-4q,
E'(s)=2s S'(s^2).                                       (16)
```

Consequently

```text
Omega=-8q v/E'(s).                                      (17)
```

The inverse edge has coordinate `-s` and inverse orientation `-Omega`.
Relabelling by `S4` transports `(17)` to all six edge/orientation sheets in
the exact pairs fixed by THM-2968.  Hence `(9)--(10)` are not merely an
abstract primitive-element argument: they are the common-gauge sheet map.

Squaring `(17)` removes the discriminant orientation.  THM-2864's radicand
identity, or direct reduction in the cubic algebra, gives

```text
Delta A(Y)^2=F(Y^2)                 mod E(Y),            (18)
```

which proves `(11)`.  Dividing by `Y^2` and using `A(Y)=YB(Y^2)` gives the
cubic-base identity `(12a)`.  Since `(5)` is the norm polynomial of `F(U)`,
`(12)` follows.  Under `(7)`, THM-2864 proves that the six orientation values are
distinct.  After passage to a splitting field, `(10)` is therefore a
bijection on six reduced points, and hence an isomorphism of the two rank-six
etale `D`-algebras.  Equivalently, the trace-discriminant change-of-basis law
and THM-2864 give

```text
Disc(O)/Disc(E)=2^60 q^10 J_or^4 Delta=det(Phi)^2.       (19a)
```

The common gauge and the rational control in the exact companion fix the
positive sign in `(12b)`.

## 3. Why the discriminant cover is exact

The base-field obstruction is group-theoretic and sharp.  For one selected
sheet, the two `S4` stabilizers are

```text
H_edge = V4_edge,              H_or = C4.                (19)
```

They are not conjugate, so the corresponding transitive degree-six
`K`-algebras are not isomorphic.  A sheet transposition fixes two edges and no
oriented cycles, which already distinguishes their permutation characters.

The discriminant cover replaces `S4` by `A4`.  In the common THM-2968 gauge,

```text
H_edge intersect A4 = H_or intersect A4 = C2.            (20)
```

Thus the two restricted six-point `A4`-sets are literally the same coset
space.  Formula `(10)` is the polynomial realization of this equality.  The
descent cocycle is exactly `(13)`; no additional square root is hidden.

This also isolates the recurring binary/ternary asymmetry.  A three-cycle is
even and already acts identically in the edge and orientation lifts.  A sheet
transposition is odd and differs by the central three-pair flip.  Passing to
the discriminant cover removes precisely that odd `C2` move while retaining
the common `C3` move.  This is the finite `S4` realization of the two
`C2*C3` faces highlighted by THM-2968; it is not a claim that every modular
group action admits this algebraic repair.

## 4. Sharp boundaries and the Keller no-go

All three factors in `(7)` are load-bearing.

1. If `q Delta=0`, the edge discriminant `(8)` vanishes and `E'` cannot be
   inverted.
2. If `J_or=0` while `q Delta!=0`, the edge algebra may remain etale but the
   chosen orientation coordinate ceases to be primitive.  At

   ```text
   (p,q,r)=(1,4,-3),        Delta=-22000,       J_or=0,
   O(Z)=(Z^2-1280)^2(Z^2+14080).                (21)
   ```

   Thus edge separability alone does not extend `(10)` across the
   orientation-collision wall.
3. The quadratic cover is ramified along `Delta=0`.  Nothing here makes it a
   quasi-etale affine cover.

The last point blocks the tempting Keller transfer.  [THM-2769, the full-`S4`
quartic hostile](THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile.md)
simultaneously has the matching cubic, the identical quartic/cubic
discriminant, both sextic complements, the grade-three depressed/cuspidal
identities, and the intertwiner above on a dense open.  Nevertheless its monic
root cover is finite/proper, its Jelonek nonproperness locus is empty, and its
`V4` normalization ramifies.  The quotient forgets the affine source,
divisor ownership, maximal-order index, and the three `V4` Kummer residues.

Accordingly `(10)` is a genuine multiplication-preserving advance over
THM-2968, but only inside the quartic splitting algebra after an oriented
discriminant base change.  It does not identify this algebra with
[THM-2950's real six-point Artin
algebra](THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame.md),
does not supply the pair-respecting contraction missing in
[THM-2951](THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary.md),
and proves no SFC(4), degree-four point-cap Keller, `JC(2)`, `DC(2)`, GMC, or
LRC closure.

The cheapest positive affine successor is therefore divisor-wise: retain, on
one common chart, the maximal-order inertia, graph-order index, present versus
omitted sheet owner, and the three `V4` Kummer residues.  The algebraic
intertwiner alone cannot reconstruct those four sidecars.

## 5. Exact evidence

Run

```text
python 04-computation/quartic_discriminant_cover_sextic_intertwiner_thm2971.py
python -O 04-computation/quartic_discriminant_cover_sextic_intertwiner_thm2971.py
```

Both modes byte-match the LF-stored transcript

```text
05-knowledge/results/quartic_discriminant_cover_sextic_intertwiner_thm2971.out.
```

The companion checks the two `S4` stabilizers, their common `A4` intersection,
the transposition fixed-sheet hostile, the explicit inverse derivative in
`(9)`, the symbolic square and orientation-polynomial identities
`(11)--(12)`, oddness, the cubic-base Kummer identity `(12a)`, both symbolic
root identities in `(16)`, the discriminant ratio and signed determinant
`(12b)`, all six common-gauge root values for one rational quartic, and the
exact `J_or=0` boundary `(21)`.  Every truth-bearing gate is
an explicit `require`/exception; there is no Python `assert` or floating-point
decision.

LF-normalized SHA256:

```text
script  6ca7c7f5695876dc38dc0a9423b4c5f81d648c29a5af90cb11b287ea0169b693
output  c9c41a197031b54f12c1bda7661cdadd1a97bf5149655d718575b6f9665dba49
```

**QED.**
