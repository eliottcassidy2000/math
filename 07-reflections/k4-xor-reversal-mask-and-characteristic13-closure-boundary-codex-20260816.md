# The pointed P4 is a K4 reversal mask, but XOR erases its characteristic-13 seam

**Status: FINITE-EXACT BOOLEAN/TOURNAMENT SIDECAR; NO D5 OR LRC
PROMOTION.**  The user's “essentially a tournament of size four” intuition
has an exact version, but it requires one extra datum.  The bidirected/missing
relation on the four owner states is not itself a tournament.  After choosing
a transitive orientation of all six `K4` pairs, its three realized path edges
act as an XOR reversal mask and produce a genuine strongly connected
tournament on four vertices.

This recoding is useful for tournament calculations, but it is not faithful
to the relation and cannot carry the proposed D5 class.  In particular, the
Boolean closure mask is exact over `F2`, whereas its consistently oriented
all-one cochain has circulation `4!=0` over `F13`.

## Inheritance and typing

Use the physical state order

```text
A=0,       B=1,       C=3,       D=2.                (1)
```

The pointed-six sidecar proves that nonempty source intersection gives the
bidirected path

```text
A <-> B <-> C <-> D,                                 (2)
```

with generalized-tournament pair census

```text
(both-way,one-way,missing)=(3,0,3).                  (3)
```

Its alternating state/root refinement is `P7`, so both static graphs have
`H1=0`.  The missing Boolean-square edge is `D--A`; adjoining it gives `C4`
and an `F13` seam coordinate.  The canonical hostile is THM-3524's rank
obstruction: a four-row Walsh carrier does not arise merely because four
states are visible.  The least-used relevant sidecar is the Boolean XOR
operation on the six unordered pair slots.

## 1. XOR needs a complete orientation gauge

Choose the transitive tournament

```text
A -> B -> C -> D,
```

with every earlier vertex also pointing to every later vertex.  On each
unordered pair its two arc bits are either `(1,0)` or `(0,1)`.  Represent an
undirected mask edge by the symmetric bit pair `(1,1)`.  Then

```text
(1,0) XOR (1,1)=(0,1),                               (4)
```

so XOR reverses precisely the selected tournament edge.  A missing mask edge
has `(0,0)` and leaves it unchanged.  Consequently, for this fixed gauge,

```text
{six-bit undirected K4 masks} <-> {labelled T4 tournaments}              (5)
```

is a bijection of two 64-element sets.

Equation (5) is the exact sense in which every both-way/missing four-state
relation can be made into a tournament.  The transitive gauge is additional
information: forgetting it loses which pairs were originally both-way and
which were missing.

## 2. The canonical P4 mask gives the strong T4

XOR the three path edges `AB,BC,CD` into the chosen transitive tournament.
The resulting arcs are

```text
B->A, C->B, D->C, A->C, A->D, B->D.                 (6)
```

They form the strongly connected `T4` with score sequence

```text
(2,2,1,1)                                             (7)
```

in the order `(A,B,C,D)`.  Exactly two of its four triples are directed
cycles.  Thus the six pointed arcs and the six pair slots of `K4` are related
by a lawful Boolean transform, not by a literal identification.

All 12 labelled Hamilton-path masks have the same sharp parity property:
their simplicial coboundary is nonzero on exactly two of the four `K4`
triangles.  For the path in (2), in triangle order `ABC,ABD,ACD,BCD`, it is

```text
(0,1,1,0).                                            (8)
```

This two-triangle obstruction is the Boolean shadow of the missing endpoint
closure.

## 3. Adding the closure edge changes the mask class

Add `DA` to the reversal mask.  Its generalized-relation census becomes

```text
(4,0,2),                                              (9)
```

while XOR again gives a strongly connected score-`(1,1,2,2)` tournament with
two directed triangles.  The more informative change is cohomological on the
complete pair set:

```text
AB+BC+CD+DA = delta({A,C}) over F2.                  (10)
```

The eight `K4` cut masks are exactly the kernel of the four triangle-parity
tests.  The path mask is outside that kernel by (8); the closed `C4` mask lies
inside it.  Hence the closure edge repairs the Boolean two-triangle defect and
makes the reversal mask a vertex cut.

This is an exact XOR formulation of the branch closure, but it does not yet
construct a clock edge or current.

## 4. Why the XOR closure cannot be the D5 flux

Orient the closed state cycle consistently as

```text
A->B->C->D->A.                                       (11)
```

For a cochain `j`, the `F13` cohomology coordinate from the pointed-P4
sidecar is

```text
Phi(j)=j_AB+j_BC+j_CD+j_DA.                          (12)
```

The all-one oriented cochain therefore has

```text
Phi(1,1,1,1)=4 mod 13 !=0.                           (13)
```

As a Boolean XOR mask, however, the same four `0/1` entries have parity

```text
1+1+1+1=0 mod 2                                      (14)
```

and equation (10) makes that cochain exact.  The cycle remains a nonzero
`F2` *homology chain*; what vanishes is this particular all-one edge
*cochain*.  Confusing those two roles would falsely turn an XOR mask into a
flux generator.

There is no field homomorphism `F13->F2`, and naive parity reduction changes
the seam value from `4` to zero.  Therefore XOR supplies none of the three
missing arrows in the proposed D5 cospan:

1. no lawful `U_clock` realization of `DA`;
2. no map from cyclotomic response amplitudes to an `F13` word-current; and
3. no semantic identification of this four-state cycle with THM-3496's
   seven-chart cycle.

Even after those arrows are built, THM-3496's additive-flux hostile and
THM-3450's two full-germ margin hostiles remain.

## 5. Complete finite atlas

The exact six-bit enumeration gives the four labelled-tournament strata

| sorted score sequence | cyclic triples | strongly connected | count |
|---|---:|---:|---:|
| `(0,1,2,3)` | 0 | no | 24 |
| `(0,2,2,2)` | 1 | no | 8 |
| `(1,1,1,3)` | 1 | no | 8 |
| `(1,1,2,2)` | 2 | yes | 24 |

The canonical P4 and C4 masks both land in the last stratum.  Tournament
isomorphism type therefore cannot detect the closure edge; the triangle
parity/cut sidecar can.

## Connection contract

| field | exact answer |
|---|---|
| vertices | owner states `0,1,3,2` |
| intrinsic relation | bidirected `P4`, census `(3,0,3)` |
| extra datum | a complete transitive `T4` orientation gauge |
| XOR map | reverse the tournament edge on each realized undirected mask edge |
| P4 image | strong `T4`, scores `(1,1,2,2)`, two cyclic triples |
| closure effect | adds endpoint edge; turns the mask into the alternating `K4` cut |
| preserved | four vertices, selected pair set relative to the chosen gauge |
| destroyed | both-way versus missing typing without the gauge, amplitudes, roots, address, time |
| coefficient hostile | oriented seam `4 mod 13` becomes parity zero mod 2 |
| D5 status | no coefficient/current/chart realization map |
| LRC status | no clock, physical current, row exclusion, or LRC(14) conclusion |

## Reproduction

```text
python -B 04-computation/k4_xor_path_closure_tournament_atlas_20260816.py
python -B -O 04-computation/k4_xor_path_closure_tournament_atlas_20260816.py
```

Normal and optimized outputs match the stored transcript.  Script/output/
semantic LF SHA-256 are respectively
`bcfdcf20a595e2e11bb0d70221624b8d435f34284609272c514ac8474f9aea64`,
`0e003b5fedffd5382cf85a39f371ae600cf1e5ffabac2eccc8deb154916e20de`,
and `1d0dfe67dc8537728c3a71bed89c3109e9604130c4dcd24d56412f37d9fadf7d`.
