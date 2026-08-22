---
id: THM-3157
title: "Pointed resolvent C3 lift, edge-hexagon orientation, and tournament alignment"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
audit: >
  An independent hostile audit rebuilt the pointed A4 lift, edge six-cycle,
  tournament census and automorphisms, odd-deck reversal, modular/Farey
  distinction, Johnson/partial-cube typing, and THM-2681 boundary.  It caught
  and repaired two novelty/scope defects: the abstract tournament is the
  inherited THM-2597 mask-873 type, and a chosen discriminant orientation
  does orient the matching C3 even though an unmarked A4 reduction does not.
  The repaired companion adds the explicit inherited-tournament alignment;
  fresh normal and optimized runs LF-byte-match the stored 25-line transcript
  and both declared hashes.
depends_on:
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
  - THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient
related:
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2641-modular-abelianization-theta-blindness-and-637-residue-no-go
  - THM-2757-marked-reference-opposite-edge-clutch-transgression
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
  - THM-3141-quartic-v4-modular-congruence-shadow-and-gamma3-sidecar-boundary
script: 04-computation/pointed_resolvent_edge_hexagon_abelianization_thm3157.py
output: 05-knowledge/results/pointed_resolvent_edge_hexagon_abelianization_thm3157.out
hash_basis: LF-normalized bytes
script_sha256: 78b989e5eb9b0fbc2ac0ed0e8bfd65c787595faeaf675c3a66162810a6d54aeb
output_sha256: b64f85f55d849102121c65cae4544b04fbd0ef98f9a986029169dcd45855b9bc
semantic_sha256: 8209e4a7944e771b7ac2603bcad27345b11f58f49e163068ec4458415c4207f1
---

# THM-3157 -- pointed resolvent `C3` lift, edge-hexagon orientation, and tournament alignment

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 0. Inheritance and exact novelty

The finite objects in this theorem are substantially inherited:

- THM-2606 Sections 5--6 proves that an owner selects one of the four
  antipode-compatible `C6` subgraphs `C_o` on the six edge slots, but selects
  no orientation of it.
- THM-2597 proves the abstract six-vertex tournament type occurring below:
  its tile-mask-`873` tournament has score multiset `(3,3,3,2,2,2)`, `45`
  Hamiltonian paths, and automorphism group `C3`.  Its bicycle-support `C3`
  rotation is explicitly not that tournament automorphism.
- THM-2632 distinguishes the Farey permutahedral `C6` from the regular cyclic
  `C6`, and THM-2641 proves that modular abelianization does not canonically
  glue the Farey/theta and cyclic coordinates.
- THM-2753 gives the matching quotient `S4/V4=S3`; THM-2756 gives the central
  complement and its rational `3+3` eigenspace/integral-clutch split; and
  THM-2757 gives the marked-owner `2 x 3` star/opposite table.
- THM-2768 gives the modular free product, its `A4/S4` quotients, and their
  Bass--Serre cycle ranks.  THM-3141/3145 identify the pointed mod-three
  congruence shadow and the once-subdivided tetrahedral quotient.

Accordingly neither the existence of `C_o`, the abstract tournament type,
the cyclic abelianization, the matching quotient, the central complement,
nor the `3+3` split is claimed as new here.  The new payload is the
simultaneous relative selector:

```text
owner o + oriented matching generator bar R
    -> unique R_o -> h_o=cR_o^-1
    -> an orientation of the inherited C_o
    -> a canonical placement of THM-2597's tournament.      (0)
```

The canonical global hostile is THM-2681: the actual THM-1310 cubic root
field cannot be the classical quartic resolvent field under the stated
target-base identification.  Thus `(0)` is local finite organization, not a
route around that obstruction.

## 1. Statement

Let `X` be a four-element set and let

```text
E = binom(X,2),                    |E|=6,                 (1)
c(e) = X \ e.                                             (2)
```

The quotient `E/<c>` is the three-element set of perfect matchings of `X`.
The natural action gives

```text
1 -> V4 -> S4 -> S3 -> 1,
A4/V4 = C3.                                               (3)
```

Choose two pieces of relative data:

1. an owner `o in X`; and
2. one of the two generators `bar R` of the matching `C3`.

Then there is a unique `R_o in A4` which fixes `o` and induces `bar R` on
the three matchings.  On `E` define

```text
h_o = c R_o^-1.                                           (4)
```

The complement commutes with every sheet permutation, and hence

```text
h_o^2=R_o,             h_o^3=c,             h_o^6=1.      (5)
```

Moreover `h_o` is one six-cycle.  For `X={0,1,2,3}`, `o=0`, and
`R_o=(1 2 3)`, it is

```text
01 -> 12 -> 02 -> 23 -> 03 -> 13 -> 01.                  (6)
```

After orientation is forgotten, this is exactly THM-2606's inherited
owner-cycle `C_o`: its edges are the pairs `ov--vw` with `o,v,w` distinct.
Thus `(4)` does not construct a new fifth hexagon; it canonically orients the
one already selected by `o`.

There is a unique tournament `T(o,bar R)` on the six edge slots satisfying
all three relative prescriptions:

- the owner-star triangle and the complementary triangle follow `R_o`;
- every opposite pair is oriented from the owner-star edge to its
  complement; and
- the other six cross-pairs follow the oriented `h_o` hexagon.

Its score multiset is

```text
(3,3,3,2,2,2),                                           (7)
```

and its full automorphism group is the selected `C3=<R_o>`.

Abstractly, this is exactly THM-2597's mask-`873` tournament.  In the edge
order `(01,02,03,12,13,23)`, the relabeling

```text
(01,02,03,12,13,23) -> (2,4,6,5,3,1)                    (7a)
```

is an isomorphism, and it carries `R_o` to THM-2597's tournament
automorphism `(1 3 5)(2 4 6)`.  This is a new **relative alignment**, not a
new tournament isomorphism type.  In particular it does not identify
`R_o` with THM-2597's different bicycle-support `C3` rotation.

More sharply, among all tournaments invariant under `R_o`, the successive
numbers surviving the opposite-pair, two-triangle, and hexagon prescriptions
are

```text
32 -> 16 -> 4 -> 1.                                      (8)
```

There is no unpointed `A4`-invariant tournament on `E`.

Finally, if `tau in S4` is odd, then relative to a fixed initial matching
orientation,

```text
tau R_o tau^-1 = R_(tau o)^-1,
tau h_o tau^-1 = h_(tau o)^-1.                           (9)
```

Thus an odd/discriminant relabeling reverses the oriented hexagon.  With an
owner retained, the underlying unoriented hexagon transports; its orientation
does not.  Forgetting the owner as well leaves a four-member orbit, not one
canonical hexagon.

## 2. Proof of the pointed lift

The action of `S4` on the three complementary edge pairs has kernel `V4`,
giving `(3)`.  The `V4` action on `X` is regular.  Therefore the point
stabilizer `(A4)_o`, which has order three, meets `V4` trivially.  Its image
in `A4/V4=C3` is consequently an isomorphism.  A chosen nonidentity
`bar R` has exactly one lift in this stabilizer.  This proves existence and
uniqueness of `R_o` without coordinates.

The involution `c` commutes with `S4` because
`sigma(X\e)=X\sigma(e)`.  Since `c^2=1` and `R_o^3=1`, equation `(4)` gives

```text
h_o^2=R_o^-2=R_o,       h_o^3=cR_o^-3=c.                 (10)
```

Thus `h_o` has order six.  It cannot split into shorter cycles: `(10)` says
its cube is the fixed-point-free complement on all six edges, while its square
has the two three-cycles formed by the owner star and its complement.  The
only compatible cycle type is one `6`-cycle.  Directly evaluating the
owner-zero frame gives `(6)`.

## 3. Proof and census of the tournament

The six unordered pairs internal to the two `R_o`-triangles, the three
opposite pairs, and the six remaining cross-pairs are disjoint and exhaust
all

```text
3 + 3 + 3 + 6 = binom(6,2)=15                            (11)
```

tournament pairs.  The three prescriptions therefore define one relative
placement of the inherited tournament.  They are `R_o`-stable, so `<R_o>`
lies in its automorphism group.

The `R_o` action has five orbits on unordered pairs: one orbit internal to
each triangle and three cross-orbits.  Because these orbits have odd length,
each can be oriented in either direction with no self-reversal.  This gives
`2^5=32` invariant tournaments.  Fixing the opposite orbit leaves `16`;
fixing both triangle orbits leaves `4`; fixing the two remaining cross-orbits
by the `h_o` cycle leaves `1`.  This proves `(8)`.

The owner-star vertices have score three and the other three have score two,
so every automorphism preserves both triangles.  Its restriction to the
owner triangle is a rotation.  After composing with a power of `R_o`, suppose
that it fixes the owner triangle pointwise.  Each complementary vertex has a
different singleton win-signature against that triangle, so all three are
fixed as well.  Hence

```text
Aut(T(o,bar R))=<R_o>=C3.                                (12)
```

For the unpointed obstruction, each complementary edge pair is swapped by
two elements of the `V4` kernel.  No orientation of that pair can be invariant
under such a swap.  Consequently no tournament on `E` is `A4`-invariant.

## 4. Why the hexagon is not quartic monodromy

The complement `c` acts on `E` with cycle type `2^3`.  The edge action of an
involution in `S4` has cycle type `1^2 2^2`; the other `S4` edge-action cycle
types are `1^6`, `3^2`, and `2 4`.  Thus

```text
c notin S4 acting on E.                                  (13)
```

It is central, so adjoining it gives the order-`48` octahedral edge group

```text
<S4,c> = S4 x C2.                                        (14)
```

The order-six motion `h_o` lives in `(14)`, not in the quartic sheet
monodromy.  This is the first failed implication in any attempted Keller
application: a canonical six-edge motion does **not** follow from an `S4`
quartic fibre alone.  It requires both the owner and the external complement
operation.

Conjugation by an odd sheet permutation reverses `A4/V4=C3`, yielding `(9)`.
If the odd permutation fixes `o`, it sends `h_o` literally to `h_o^-1`.
Accordingly the unoriented hexagon is the sharp descent survivor after the
orientation is forgotten, while the directed tournament is not.

## 5. The modular `C6` shadow is abelianization, not the `A4` quotient

Write the modular group in free-product form

```text
Gamma=PSL2(Z)=<s,r | s^2=r^3=1>=C2*C3.                  (15)
```

THM-2632/2641 already identify

```text
Gamma^ab=C2 x C3=C6.                                    (16)
```

Once the orientation of `h_o` is fixed, the inherited quotient is realized
on the six edge slots by

```text
s |-> h_o^3=c,           r |-> h_o^2=R_o,
s r^-1 |-> h_o.                                          (17)
```

Because `h_o` has order six, `(17)` is precisely the full abelianization
action.  Reversing the chosen `C3` generator reverses its cyclic gauge.  This
is the cyclic `C6` side of THM-2632, now realized on the owner-selected
`C_o`; it is not a new proof that `(16)` exists.

The nonabelian mod-three map of THM-2768/3141 is different:

```text
s |-> S in V4,           r |-> R_o,
<S,R_o>=A4.                                              (18)
```

In `(17)`, the binary image is the external complement `c`, which is not in
the sheet `S4`; in `(18)`, it is a genuine `V4` sheet permutation.  The
first images commute and generate `C6`; the second do not commute and
generate `A4`.  There is no compatible factorization `A4 -> C6`: the
abelianization of `A4` is `C3` and kills its `V4` binary generator.

The information loss is concrete.  The distinct reduced words

```text
s r^-1,                         r^-1 s                    (19)
```

have the same image `h_o` in `(17)`.  With the standard integral matrices
of THM-2641 they are, however,

```text
s r^-1 = [1  0],                r^-1 s = [1 -1],
          [1  1]                          [0  1]          (20)
```

and hence are distinct even projectively in `PSL2(Z)`.  The abelianization
kernel is the commutator subgroup (a torsion-free index-six free group of
rank two); all commutator order, Bass--Serre path, and Farey ancestry inside
it is forgotten.  Equations `(19)--(20)` are the cheapest hostile to treating
the oriented edge hexagon as a faithful modular tree action.

The Farey six-state graph of THM-2632 is abstractly another `C6`, and there
are `12` graph isomorphisms from it to a fixed `C_o`.  But its modular left
action is regular `S3` and has no one-step order-six element.  The `h_o`
rotation is instead the regular cyclic abelianization action.  Owner plus
oriented `bar R` supplies one particular cyclic realization on `E`; it does
not retroactively make the Farey gluing canonical in the unpointed data of
THM-2641.

## 6. Johnson octahedron versus the partial-cube incidence frame

Identify the six edge slots with the two-subsets in the Johnson graph

```text
J(4,2):  e~f iff |e intersect f|=1.                      (21)
```

This is the octahedral graph: it has six vertices, twelve edges, degree four,
and eight triangular faces.  Complement is its antipodal involution; for
each vertex `e`, `c(e)` is its unique distinct nonneighbor.  Every consecutive
pair in `(6)` meets in one point, so `h_o` is a Hamilton cycle of `J(4,2)`.
The two `R_o` triangles use the other six Johnson edges, while the three
opposite tournament pairs are precisely the three nonedges of `J(4,2)`.

This must not be confused with THM-3145's partial cube.  The once-subdivided
`K4` is the singleton/two-subset inclusion graph

```text
{ {i}:i in X } disjoint_union binom(X,2),                 (22)
```

embedded in `Q4` by characteristic vectors.  Inclusion changes one bit, and
graph distance is exactly Hamming distance; its four theta classes each have
three edges.  By contrast an `h_o` step joins two weight-two vectors and has
Hamming distance two.  Hence `(6)` is a Hamilton cycle in the Johnson
octahedron, **not** an edge cycle in the partial-cube incidence graph.
Indeed `J(4,2)` contains triangles and is not itself a partial cube.

For representation bookkeeping only, THM-2756's central complement gives

```text
Q[E]=Q[E]_+ direct_sum Q[E]_-,              dimensions 3+3, (23)
```

spanned by opposite-pair sums and differences.  This is a `Z/2`-graded
`S4` representation, and `h_o` preserves both grades because it commutes
with complement.  No multiplication, parity-changing action, or superbracket
has been supplied.  Thus `(23)` is not a Lie superalgebra, Lie supergroup, or
physical supersymmetry claim.

## 7. Geometric contract and stopping boundary

This theorem is finite fibre combinatorics.  A geometric application has the
following exact contract.

```text
source:     a separable quartic four-sheet fibre and its matching cubic;
target:     the six pair slots and their three complementary classes;
map:        owner-fixed lift of a chosen matching-C3 generator, followed by c;
preserved:  matching classes, opposition, selected owner, selected C3 direction;
destroyed:  owner and direction after unpointed/discriminant descent;
needed:     common base-ring identification, physical owner, oriented cubic
            generator, and a common nonzero current/atom;
test:       24 sheet relabelings, 32 tournaments, and odd-deck reversal.
```

The sheet-orientation torsor and the orientation torsor of the three matching
classes are canonically identified because
`sign_S4=sign_S3 o mu`.  Thus a **chosen** square root/trivialization of the
quartic discriminant does orient the matching `C3` and selects `bar R`.
What does not select a direction is only the bare statement that the
discriminant is a square, or an unmarked reduction to `A4`; it forgets which
trivialization was chosen.  Neither version chooses a sheet owner, and no
version supplies the required common physical current.

Most importantly, THM-2681 excludes the actual THM-1310 cubic root field from
being the classical matching-resolvent field of a quartic `S4` Keller map
under the stated global target-coordinate/base-ring identification; its `A4`
field type is excluded there as well.  Therefore the depressed cubic law,
discriminant identity, Jelonek lead, and cuspidal law cannot simply be grafted
onto `(4)`--`(12)`.  A usable Keller bridge would first need a different
resolvent, a weaker explicitly preserved shadow, or a scoped base change,
together with the owner/orientation/current data above.

No quartic Keller exclusion, `G1`, LRC gate, degree-four point-cap theorem, or
physical modular-group action follows here.  The positive result is exact but
local and conditional: once a common quartic owner and an oriented matching
generator are genuinely present, they canonically organize the six pair slots
as one oriented hexagon and one relative tournament.

## 8. Exact evidence

The self-contained companion is

```text
04-computation/pointed_resolvent_edge_hexagon_abelianization_thm3157.py
```

and its stored transcript is

```text
05-knowledge/results/pointed_resolvent_edge_hexagon_abelianization_thm3157.out
```

It uses explicit raises rather than Python `assert`.  Normal and optimized
runs byte-match after LF normalization.  It enumerates all `24` sheet
permutations, the `S4/V4=S3` and `A4/V4=C3` actions, all eight pointed/oriented
frames, all `2^15` tournaments at every frame, all `720` six-vertex
automorphisms, the `32/16/4/1` census, the absence of an `A4`-invariant
tournament, the explicit THM-2597 mask-`873` isomorphism and automorphism
alignment, the `48`-element octahedral extension, and all `96`
sheet-permutation/owner reversal cases.  It additionally checks `(17)` and
the noncommuting `A4` contrast, the two exact matrices in `(20)`, all four
inherited THM-2606 owner cycles, the `12` abstract Farey/owner-cycle graph
isomorphisms, the Johnson and triangular-face censuses, the full `10`-vertex
`Q4` distance table and four theta classes, every `h_o` Hamming-two step, and
the complement `3+3` eigenspace split.

Final LF-normalized evidence hashes are

```text
script 78b989e5eb9b0fbc2ac0ed0e8bfd65c787595faeaf675c3a66162810a6d54aeb
output b64f85f55d849102121c65cae4544b04fbd0ef98f9a986029169dcd45855b9bc
```

The reproduction commands are

```text
python 04-computation/pointed_resolvent_edge_hexagon_abelianization_thm3157.py
python -O 04-computation/pointed_resolvent_edge_hexagon_abelianization_thm3157.py
```
