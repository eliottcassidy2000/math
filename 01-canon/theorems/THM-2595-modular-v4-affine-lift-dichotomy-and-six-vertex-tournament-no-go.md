---
id: THM-2595
title: "Modular V4 affine-lift dichotomy and the six-vertex tournament no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For the standard quotient
  G=C2*C3 -> GL(2,F2)=S3 on V=F2^2=V4, affine lifts to
  AGL(2,F2)=S4 form eight 1-cocycles and four translation
  coboundaries, so H^1(G,V)=F2.  Both restrictions H^1(C2,V) and
  H^1(C3,V) vanish: the unique nonzero global class is a genuine
  free-factor compatibility defect.  The zero class has image S3
  and orbit partition 1+3 on V; the nonzero class has image S4 and
  is transitive.  Equivalently, each affine C2 generator fixes a
  pair, each affine C3 generator fixes one point, and the class is
  zero exactly when that point lies in that pair.  The six elements
  of the quotient S3 admit no left-translation-invariant tournament,
  because each involution reverses its edge with the identity.  This
  is pure finite group anatomy.  It neither excludes degree-four
  Keller maps nor identifies an LRC carrier.
source: codex-modular-v4-session-2026-07-27
depends_on: []
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
script: 04-computation/psl2z_v4_affine_lift_dichotomy_thm2595.py
output: 05-knowledge/results/psl2z_v4_affine_lift_dichotomy_thm2595.out
script_sha256: 82fe926bcb2814fbfcf441387ec2c1e788a25c394128072fa8c8a2ae8d19fbeb
output_sha256: 07b48abdca9d33faffaf557796a1c1c60c7815468b78bfefbe16377632a6671b
hash_basis: working-tree bytes (LF)
---

# THM-2595 -- modular V4 affine lifts see 2 and 3 only together

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem makes one part of the proposed `2/3 -- V4 -- S3 -- S4`
picture exact and sharply types what it does not say.  Its central result is

```text
locally trivial on C2 and locally trivial on C3
                 does not imply
globally translation-trivial on C2*C3.                         (1)
```

The surviving bit is precisely the dichotomy between the intransitive
point-stabilizer `S3 < S4` and the transitive full `S4` action.

## 1. The marked modular quotient

Let `V=F_2^2`, written additively, and take

```text
      [0 1]          [0 1]
S  =  [1 0],    T =  [1 1].                                  (2)
```

Then `S^2=T^3=1` and `<S,T>=GL(V)`.  Hence

```text
G = <s,t | s^2=t^3=1> = C2*C3  -->  GL(V)=S3,
s |-> S, t |-> T                                                   (3)
```

is the standard marked quotient.  The affine group is

```text
AGL(V)=V semidirect GL(V)=V4 semidirect S3 = S4,               (4)
```

where the last equality is the faithful action on the four points of `V`.

An affine lift of (3) has the form

```text
x |-> rho(g)x+c(g),
c(gh)=c(g)+rho(g)c(h).                                         (5)
```

Conjugating the lift by the translation `x |-> x+v` adds the coboundary

```text
delta_v(g)=v+rho(g)v.                                          (6)
```

Thus translation-conjugacy classes of lifts are `H^1(G,V)`.

## 2. Eight cocycles, four coboundaries, one compatibility bit

Write `a=c(s)` and `b=c(t)`.  Since the presentation (3) has no mixed
relation, the complete cocycle conditions are

```text
(1+S)a=0,                 (1+T+T^2)b=0.                        (7)
```

The first solution space is the one-dimensional line
`{(0,0),(1,1)}`.  The second is all of `V`, because
`1+T+T^2=0`.  Therefore

```text
|Z^1(G,V)| = 2*4 = 8.                                         (8)
```

No nonzero vector is fixed by all of `GL(V)`, so the coboundary map
`V -> Z^1(G,V)` is injective.  Consequently

```text
|B^1(G,V)|=4,       H^1(G,V) = F2.                             (9)
```

The mechanism is stronger than the count.  On the two factors separately,

```text
ker(1+S)=im(1+S),
ker(1+T+T^2)=V=im(1+T),                                       (10)
```

so

```text
H^1(C2,V)=0,                    H^1(C3,V)=0.                   (11)
```

Every affine lift can therefore be gauged to linear form on the binary face
and on the ternary face separately.  The nonzero class in (9) records that
the two local gauges cannot be chosen to agree.  This is the exact
free-factor **co-occurrence** phenomenon: neither prime face detects it on
its own.

## 3. Fixed loci and the S3/S4 dichotomy

For every admissible `a`, the affine involution `x |-> Sx+a` fixes exactly
two points.  For every admissible `b`, the affine order-three map
`x |-> Tx+b` fixes exactly one point.  A cocycle is a coboundary precisely
when the full affine action has a fixed point.  Hence

```text
[c]=0  <=>  Fix(S+a) intersects Fix(T+b)
       <=>  the unique T+b fixed point lies in the S+a fixed pair.       (12)
```

The exhaustive census contains four cocycles of each kind.  In the zero
class, the generated group has order six and orbit sizes `1+3`; it is a
point stabilizer in `S4`.  In the nonzero class, the generated group has
order 24 and one orbit of size four; it is all of `S4`.

There is also a conceptual proof of the second assertion.  The affine image
surjects onto `GL(V)=S3`, and its kernel inside the translation group is an
`S3`-invariant subgroup of `V`.  Since `S3` acts transitively on the three
nonzero vectors, that kernel is either zero or all of `V`.  If it is zero,
the image is an `S3` complement; every such subgroup fixes one of the four
points, so the cocycle is a coboundary.  A nonzero class must therefore have
kernel `V`, giving the full group of order `4*6=24`.

## 4. What this says about a quartic resolvent

For four labelled roots, the three partitions into two unordered pairs are
the three perfect matchings of `K4`.  The permutation action on those three
matchings is

```text
S4 --> S3,       kernel V4.                                   (13)
```

Equivalently, the matchings are the three nonzero translation directions in
`V`, and (13) is the linear-part map `AGL(V)->GL(V)`.  They are also the
three root labels of a standard cubic resolvent.  The theorem therefore
classifies affine lifts of a **marked full-S3 resolvent action** coming from
the modular presentation: up to translation gauge, it has exactly the split
`1+3` lift and the transitive `S4` lift.

This is a structural interface, not a degree-four Keller obstruction.
THM-2455 already proves that the quartic and its standard resolvent cubic
have identical discriminant.  The cubic resolvent of a graph quartic is not
thereby a three-dimensional Keller map, so the proved grade-three Keller
anatomy cannot be transferred to it without a new realization theorem.
Moreover, the marked full-`S3` hypothesis places this statement only in the
`S4` quartic-monodromy branch.  It says nothing about the `A4 -> C3` or
`D4 -> C2` resolvent branches.
The next legitimate test is more specific: construct the graph quartic of a
hypothetical degree-four point-cap map, mark its `C2/C3` local monodromies,
and ask whether the two local affine gauges can satisfy the Keller/Jelonek
sidecars simultaneously.  No such construction is supplied here.

## 5. Why the natural six-object is not a tournament

The six quotient elements form `S3`, but there is no
left-translation-invariant tournament on this vertex set.  If `u` is any of
the three involutions and the edge is oriented `1 -> u`, left translation by
`u` forces `u -> 1`; the opposite starting orientation gives the same
contradiction.  Thus a size-six tournament here requires extra, non-invariant
gauge data.  It is not intrinsic to the quotient.

On the three matching labels, an order-three element does choose a cyclic
orientation.  Every order-two element reverses it.  The typed object is
therefore the `S3` action on a three-set (or a cyclic orientation together
with its reversing gauge), not a canonical six-vertex tournament.

## 6. Tree-language guardrail

The Bass--Serre tree of `C2*C3` is `(2,3)`-biregular: the two valences come
from cosets of the finite factors.  This does **not** make every binary tree
the `C2` factor or every ternary tree the `C3` factor.  In particular, a
rooted ternary semigroup tree has three independent child moves, whereas an
element of `C3` is one invertible move of order three.  Any identification of
a concrete fraction or Pythagorean-triple tree with (3) must exhibit the
actual involution, order-three map, preserved object, and quotient losses.

## 7. Reproduction and scope

Run

```bash
python3 04-computation/psl2z_v4_affine_lift_dichotomy_thm2595.py
```

The script enumerates `GL(2,F2)`, all eight cocycles, all four
coboundaries, both restriction complexes, every affine image and orbit,
fixed-locus intersections, all 24 affine permutations, and the three
involution reversals.  Normal and `python3 -O` outputs byte-match the stored
result.

This theorem is exact finite group theory.  It does not prove the existence
of a polynomial cover with the indicated monodromy, constrain a Keller map's
coefficients, trivialize the LRC `C91` holonomy, or construct a common
physical Boolean carrier.
