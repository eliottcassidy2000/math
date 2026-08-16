# The tetrahedron behind tournaments, XOR, Berggren words, and graph H1

**Status: elementary structural proof plus VERIFIED-EXACT finite atlas.**

The user's phrase “these are all essentially a tournament of size four” has
an exact positive core and a sharp boundary.  The common carrier is the
tetrahedron `K4`; the predicates live in different representations of its
automorphism group `S4`.

## Inheritance pass

- **Closest proved mechanisms.** THM-3497 identifies the Berggren matching
  quotient `S4->S3`, while the typed cohomology synthesis constructs the
  opposite-face map `Omega:C^0(K4)/1->H^1(K4)`.
- **Canonical hostile.** The matching quotient identifies the identity with
  all three double transpositions.  Any claimed reconstruction of the full
  vertex action from the three matchings fails there.
- **Corrected near miss.** XOR detects whether two `K4` edges meet; it does
  not orient either edge.  Tournament direction is an extra six-bit cochain.
- **Least-used sidecar.** The four oriented opposite triangles give the
  missing cycle coordinate and expose exactly which source information the
  12-state Berggren language discards.

## 1. Five tetrahedral representations

The same four labels produce:

```text
vertices:           P1(F3), or the four affine points of F2^2,
edges:              six weight-two vectors e_i+e_j in F2^4,
perfect matchings:  three partitions 01|23, 02|13, 03|12,
opposite faces:     four oriented triangle boundaries,
tournaments:        2^6 orientations of the six edges.
```

These sets have cardinalities `4,6,3,4,64`.  `S4` acts on all of them, but
the actions do not retain the same information.  The edge action is the
XOR/incidence carrier; the matching action has kernel `V4`; and the
opposite-face action is the sign-twisted standard representation.

This is the precise sense in which the objects are “the same size-four
thing.”  It does not identify an edge orientation with an XOR weight, a
matching with a current, or a word language with a tournament class.

## 2. Exact cut/cycle coordinates for every labelled T4

Orient the six reference edges from smaller to larger endpoint.  Let `B` be
the `6 x 4` incidence matrix, and let

```text
Omega(v)=sum_i (-1)^i v_i partial(F_i),                 (1)
```

where `F_i` is the face opposite vertex `i`.  Delete one constant column
from each map.  The exact integer determinant is

```text
det [B_0 B_1 B_2 Omega_0 Omega_1 Omega_2] = -16.        (2)
```

Thus, over every field of characteristic not two, cut and opposite-face
cycle coordinates form a direct decomposition of the six-dimensional edge
space.  The lattice index `16` records the integral denominator cost.  This
is the same split used in the typed cohomology synthesis, now interpreted on
tournament orientations.

Encode a labelled tournament by `s in {+1,-1}^6`, with `-1` meaning that the
lower-labelled endpoint beats the higher-labelled endpoint.  Define

```text
score(s)=B^t s,
face(s)=Omega^t s.                                      (3)
```

Then

```text
score_i=2 outdegree(i)-3.                               (4)
```

Each `face_i` is the oriented edge sum on the triangle opposite `i`:

```text
|face_i|=1  iff that triple is transitive,
|face_i|=3  iff that triple is a directed 3-cycle.       (5)
```

Equations (2)--(3) prove that the score vector together with the four face
circulations reconstructs every labelled tournament.  Neither half alone
does so.

The 64 labelled orientations have the four familiar unlabeled orbits:

| score sequence | orbit size | cyclic opposite faces |
|---|---:|---:|
| `(3,1,1,1)` | 8 | 1 |
| `(2,2,2,0)` | 8 | 1 |
| `(3,2,1,0)` | 24 | 0 |
| `(2,2,1,1)` | 24 | 2 |

This table also handles tournaments with missing or both-way edges once a
tie is retained as edge value zero: (3) remains a faithful linear observer
of the full edge cochain, while the object is no longer a tournament.

## 3. What the three perfect matchings forget

The action

```text
q:S4 -> Sym({01|23,02|13,03|12})=S3                    (6)
```

is onto and has kernel

```text
ker(q)={identity and the three double transpositions}=V4. (7)
```

By contrast, in odd characteristic the opposite-face action on graph `H1`
is faithful: all 24 elements of `S4` give distinct `3 x 3` matrices.  Its
character is `sign(pi)(#Fix(pi)-1)`.  The exact class ledger is

| `S4` class | size | matching class | `H1` trace | affine linear bits allowed |
|---|---:|---|---:|---|
| `1^4` | 1 | `1^3` | 3 | `0` |
| `2*1^2` | 6 | `2*1` | -1 | `1` |
| `2^2` | 3 | `1^3` | -1 | `0` |
| `3*1` | 8 | `3` | 0 | none |
| `4` | 6 | `2*1` | 1 | `1` |

The two pairs merged by the matching quotient have different `H1` traces:

```text
identity versus double transposition: 3 versus -1,
transposition versus four-cycle:      -1 versus 1.       (8)
```

Thus the opposite-face representation is an exact sidecar that restores the
full `S4` action.  It is not a scalar invariant: choosing one coordinate
still requires a face/orientation gauge.

## 4. Why A and C fail while B passes

THM-3497's projective branch permutations are

```text
A=(0 1 3), fixing 2,
B=(0 1 3 2),
C=(0 3 2), fixing 1.                                   (9)
```

So `A` and `C` are literal rotations of the face opposite their fixed
vertex.  Their matching actions are 3-cycles.  The prescribed affine linear
parts `I` and `J` have only the following translate classes:

```text
I+t: identity or double transposition,
J+t: transposition or four-cycle.                       (10)
```

No translate in (10) has a matching 3-cycle, which is the structural
three-primary obstruction for `A,C`.  `B` is a four-cycle, its matching
action is a reflection, and it lies in the `J`-affine fibre.

For a composite word, THM-3497's variable language retains only

```text
(q(rho_3(w)), epsilon(w)) in S3 x C2.                   (11)
```

It accepts matching identity with bit zero and matching reflections with bit
one.  Equation (8) shows exactly what (11) forgets: it deliberately merges
the two source conjugacy classes inside each lawful affine fibre.  That loss
is why the minimal automaton has 12 states rather than a faithful 24-state
`S4` carrier.  The fixed-drift language keeps the full source `S4` and pairs
it with `D4`, yielding its distinct 192-state automaton.

## 5. Harmonic and cohomological boundaries

The word languages become subsets of the natural numbers only after a
shortlex address map is fixed.  THM-3497's coefficients `1/3` and `17/96`
are then logarithmic densities of those index subsets.  They are not masses
of tournament orbits, nor reciprocal sums of Pythagorean triple values.

Likewise, `Omega` gives a genuine graph-`H1` class, whereas an edge gradient
remains a coboundary.  The tetrahedral atlas explains how a nonzero cycle
observer can coexist with THM-3479's zero absolute current class.  It still
does not supply a physical LRC phase, a D5 Kummer marking, or a common-
ancestry lift.

## 6. Reproduction and scope

```bash
python3 04-computation/tetrahedral_k4_tournament_berggren_representation_atlas_20260816.py
python3 -O 04-computation/tetrahedral_k4_tournament_berggren_representation_atlas_20260816.py
```

The companion checks all 24 permutations, all four affine translations for
both linear bits, and all 64 labelled tournaments.  Its complete
score/face ledger has SHA256

```text
0e3878081e5cfd5ee63a24e85554605e1f06c5dcb54d389cb8c29634e2a4a2df.
```

This is a finite representation atlas.  It proves no physical current, LRC
row exclusion, Jacobian-conjecture consequence, or all-size tournament law.
