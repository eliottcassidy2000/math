---
title: "Two and three are local generators; the theorem lives in the extension data"
date: 2026-07-28
status: RESEARCH SYNTHESIS FROM PROVED THM-2595/2596/2646/2680/2681/2682; NO LRC OR JACOBIAN CLOSURE
source: root-2026-07-28-modular-path-lifting-synthesis
---

# The recurring object is controlled forgetting

The prompt's `2,3,4,6` frame has a rigorous core, but it is not a universal
tree or a numerical coincidence.  Four proved mechanisms now show the same
operation at four different categorical depths:

```text
global affine lift -> factor restrictions loses a Cech/H^1 compatibility class;
physical edges -> sector adjacency    loses path lifting;
full ordered-root normalization -> all-root torus loses codimension-one valuations;
braid group    -> modular quotient    loses central height.               (1)
```

In every row, the quotient is honest for consumers which annihilate the
forgotten coordinate.  It is false for the live target until that coordinate
is restored.  This is the precise sense in which two and three are special:
they generate the smallest nontrivial free-product and matching actions, but
the proof obligation is carried by how those local pieces coexist.

## 1. The four exact failures

### Free-factor gauges do not glue

THM-2595 computes the affine lifts of

```text
C2*C3 -> GL_2(F2)=S3=Aut(V4).                              (2)
```

The restrictions of every lift to `C2` and `C3` are translation-trivial:

```text
H^1(C2,V4)=H^1(C3,V4)=0.                                  (3)
```

Nevertheless

```text
H^1(C2*C3,V4)=F2.                                         (4)
```

The nonzero class is exactly the incompatibility of the two local gauges.  It
upgrades the intransitive `S3` affine action to the transitive `S4` action.
Thus checking the binary and ternary faces separately cannot detect the
global affine carrier.

THM-2596 supplies the essential correction to the tree slogan.  The
Bass--Serre tree is `(2,3)`-biregular; the binary rooted Farey grammar appears
only after contracting the degree-two vertices.  The positive Farey children
are parabolic words in the torsion generators, and the ternary Berggren tree
is a reflected PGL prefix code, not the `C3` action.  The faithful LRC flank
state is `(U^T U,U^T w)`, including its off-diagonal Gram coordinate.

### Edge-surjectivity does not lift paths

THM-2680's physical dilation fibre products project to the two sector labels
`safe,danger`.  Every allowed edge of

```text
              [1 1]
A_Fibonacci = [1 0]                                         (5)
```

has a positive physical lift: the exact labelled counts are `17160,4488,
3696,0`.  Yet THM-2682 proves that every three-event product is empty.  If
`P` is the physical event graph and `Gamma` is (5), then

```text
E(P) -> E(Gamma) is surjective,
Paths_2(P) is empty,
Paths_2(Gamma) is nonempty.                                (6)
```

So the quotient is edge-surjective but not path-lifting.  The spectral radius
`phi` of (5) is a true statistic of the one-edge shadow and a false entropy
for the physical chronology.  A local Fibonacci grammar is not evidence for
a Bass--Serre action unless reduced words lift through at least the next
nerve level.

### A generic Kummer plane need not cross the missing divisors

THM-2681 normalizes the actual THM-1310 cubic splitting field on `c!=0` as

```text
B=C[c^+-1,r_1,r_2,r_3]/(r_1+r_2+r_3-1).                   (7)
```

On `a!=0`, every `r_i` is a unit and the ratios `r_1/r_3,r_2/r_3` give a
connected `S3`-equivariant `V4` Kummer torsor.  The tempting generic carrier
is real.  It fails globally because

```text
div(r_1/r_3)=D_1-D_3,   div(r_2/r_3)=D_2-D_3,             (8)
```

where `D_i=(r_i=0)` are irreducible divisors above `a=0`.  Their odd
valuations force codimension-one ramification.  On the full chart,
`B^*=C^*c^Z` and `Cl(B)=0`, so no standard two-dimensional Kummer plane
remains.

This separates the live quartic cases sharply.  An `A4` canonical resolvent
is cyclic Galois of degree three and therefore cannot equal THM-1310's
nonnormal cubic root field.  In the `S4` case the field type matches, but
THM-2655 requires the missing global `V4` torsor.  Same discriminant is much
weaker than same normalization plus its kernel cover.

### The modular knot quotient forgets the full twist

THM-2646 proves

```text
B3=<x,y | x^2=y^3>,       B3/<x^2>=C2*C3=PSL_2(Z).        (9)
```

The modular word and its exponent/central height are jointly faithful; the
modular word alone is not.  Braids with the same modular image can close to
the unknot, trefoil, or torus knots of unbounded genus after central twisting.
Thus the knot analogue of (6) is not a missing edge but a missing lift in a
central extension.

## 2. Transfer ledger

| source | quotient target | preserved | destroyed | necessary sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| `C2*C3` affine lift | its two factor restrictions | each local fixed gauge | gauge co-occurrence | class in `H^1(C2*C3,V4)` | intersect the two fixed loci |
| physical `D` events | sector matrix (5) | every one-edge type | composable middle clock | next Cech nerve | test one lifted length-two path |
| ordered cubic roots | generic root torus | `S3` Kummer plane | `r_i=0` valuations | divisor parity / `Cl[2]` | restore one missing divisor |
| `B3` braid | `PSL_2(Z)` word | modular conjugacy shadow | full-twist exponent | central integer height | compare closures in one fibre |
| Farey flank | endpoint defects | two scalar endpoint values | cross Gram entry | `(U^TU,U^Tw)` | compare the mediant defect |

The source-to-target map, preserved predicate, destroyed coordinate, and
test are different in every row.  The common theorem is methodological: a
finite shadow becomes a proof carrier only after its extension datum is
shown to be irrelevant or is reconstructed.

## 3. Sharp next targets

### LRC(14)

The fixed arrival-six `D` chronology is finished by THM-2682.  Decorating its
edges cannot create a missing two-simplex.  The next lawful target is a
**configuration-switching** correspondence: change the arrival chart or the
handoff itself, retain the predecessor/carry label, and ask for one positive
three-event fibre.  The reserved THM-2672 slope-seven carry atlas is relevant
only after promotion; its current signal concerns large simplices in the
carry direction, not temporal path lifting.  The cheapest experiment is a
bigraded incidence table

```text
(configuration switch, D-time edge, retained carry, endpoint phase)       (10)
```

and the decisive statistic is the first mixed two-simplex, not another
marginal edge count.

### Quartic Keller / Jacobian

The exact THM-1310 resolvent lane, under the stated target-field/base-ring
identification, is now closed for both `A4` and `S4`, but the monodromy types
remain live.  Search an unrelated cyclic/full-`S3`
resolvent normalization whose regular locus contains the standard
`F2[C3]`/`F2[S3]` Kummer plane globally.  Before imposing the Keller
equations, compute units, `Cl[2]`, and all codimension-one valuations.  A
generic all-root-unit chart is now known to be an adversarially weak test.

### Knots

Any conclusion meant to distinguish closures within a modular fibre must
either prove invariance under central twisting or retain central-height or
Markov-stable data.  A finite family
can share an optimal context by THM-2281, but that says nothing about a
positive catalyst or about stability under the full-twist fibre.  The cheap
test is to hold the modular class fixed, vary central exponent, and demand
the proposed invariant distinguish the resulting closures.

### Tournaments and graceful/partial-cube shadows

THM-2606 already shows that the six unordered pairs carry an intrinsic
incident/disjoint association scheme, while a `C6` choice spends an origin
and no invariant tournament exists.  Therefore a six-vertex tournament is
useful only when a genuine pairwise observable or oriented gauge is supplied.
The new LRC lesson is stronger: even an exact binary relation needs a
path-lifting audit before its Hamiltonian or spectral grammar is transferred.

## 4. Stopping rule

Do not infer a global theorem from any of the following alone:

```text
the right generator orders;
the right finite quotient group;
the right one-edge adjacency matrix;
the right generic Kummer characters;
the right discriminant;
the right modular knot class.                              (11)
```

Ask instead which extension datum the target consumes.  In the present
frontier those data are, respectively, the common affine gauge, the middle
clock, divisor parity, and central height.  This reframing turns the role of
the primes two and three from numerology into a precise audit of local
generators versus global compatibility.
