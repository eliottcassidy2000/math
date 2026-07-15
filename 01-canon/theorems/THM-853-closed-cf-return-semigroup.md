---
id: THM-853
title: The closed continued-fraction lift/face returns form a 46,126-map identity-adjoined contraction monoid with six coordinate horizons
status: PROVED EXACT COORDINATE MONOID/CAYLEY FLOW + FINITE-EXACT Q FUTURE-CONGRUENCE CLASSIFICATION
source: codex-2026-07-15-S13g
depends_on: [THM-828, THM-838, THM-840, THM-842, THM-846]
related: [THM-810, THM-813, THM-843, THM-851, THM-852, THM-854, HYP-6880]
verification:
  - 04-computation/n9_closed_cf_return_semigroup_codex_S13g.py
  - 05-knowledge/results/n9_closed_cf_return_semigroup_codex_S13g.out
  - 05-knowledge/results/n9_closed_cf_return_semigroup_codex_S13g.json
---

# THM-853 — the closed CF return monoid

Let `A,B,C:X_9->X_9` denote THM-846's three closed returns

```text
F_R=R_10 Phi_(9->10).                                   (1)
```

Read words from left to right, so `AB=B composed A`.  The identity-adjoined
transformation semigroup—that is, the monoid of the empty word and all
nonempty A/B/C words—has exactly

```text
46,126 coordinate-copy maps.                              (2)
```

The identity cannot be represented by a nonempty word because every generator
has rank below 28.  Hence the strict nonempty-word semigroup has 46,125 maps.

It is a large but sharply contracting recursive object.  Its right Cayley
graph has six rank-one terminal sinks, supported on six central staircase
coordinates.  On THM-828's 58 canonical rows, B is a Q-congruence after
*every* monoid state, while A and C each fail at exactly 76 states.  The
present Q equality partition takes only 144 values, but equality under every
future A/B/C continuation has 20,419 classes.  Thus a compact static quotient
can be two orders of magnitude too small as a recursive state.

## Scope: return words, not increasing-size CF words

Every letter in (1) lifts from size nine to ten and immediately returns to
size nine.  A product in this theorem repeats that fixed return.  It is a
genuine endomorphism monoid of the n=9 staircase, but it is **not** the
increasing-size continued-fraction cocycle

```text
X_9 -> X_10 -> X_11 -> ... .                              (3)
```

Interpreting (2) as arbitrary consecutive CF digits would require the ambient
size, Christoffel phase, and substitution token action.  Those data are not
present here.  This guardrail is part of the theorem.

## 1. Exact transformation-monoid census

Breadth-first enumeration from the identity, in generator order A,B,C, gives
shortest representatives of length at most 22.  The coordinate ranks are

```text
rank   maps
 28       1
 27       1
 21       4
 20       4
 15      16
 14      14
 11      12
 10      54
  9      32
  8      36
  7      78
  6     180
  5     318
  4   1,038
  3  25,146
  2  19,186
  1       6.                                               (4)
```

In particular, no rank 22 through 26 map exists, and 44,332 of 46,126 maps
have rank two or three.  There are 3,149 idempotents:

```text
rank 28:1, 27:1, 4:2, 3:130, 2:3009, 1:6.                 (5)
```

The simplest contraction laws are

```text
rank(A^k)=rank(C^k)=21,15,11,8,5,4 for k=1,...,6,
A^6=A^7,             C^6=C^7,
B^2=B,
rank((AC)^k)=rank((CA)^k)=15,6,3 for k=1,2,3,
(AC)^3=(AC)^4,       (CA)^3=(CA)^4.                       (6)
```

The two-letter rank table is

```text
        next A   next B   next C
A          15       21       15
B          20       27       20
C          15       21       15.                          (7)
```

The monoid contains only 131 distinct visible-coordinate supports, even
though it contains 46,126 different copy maps.  Multiplicity and placement of
the copied outputs therefore carry most of the transformation data; support
alone is not a recursive address.

Every map commutes with bitwise complement.  Reflection conjugation exchanges
A and C and fixes B.  It has only twelve fixed maps, so there are 23,069
reflection orbits.  This is an exact left/right symmetry of the coordinate
monoid.  It will coexist below with a strong B versus A/C dynamical
asymmetry.

The sorted map certificate has SHA-256

```text
1df8e46ee522452ebf91841ae3d420abfa75f4fd056a3233b1453c64632aadc7. (8)
```

## 2. Right Cayley flow and the six horizons

The right Cayley graph has 46,126 vertices and 138,378 labelled edges.  Rank
never increases.  Its 31,075 strongly connected components have size census

```text
31,058 x 1,       6 x 216,       11 x 1,252.              (9)
```

Each 216-component consists of rank-three maps and can reach three terminal
sinks.  Each 1,252-component consists of rank-two maps and can reach two.
The six terminal components are singleton rank-one idempotents.  Their
surviving source coordinates and shortest words are

```text
source tile   shortest word
(5,3)         AABAAACC
(6,3)         ABACCAACC
(6,4)         AAACCCCC
(7,3)         BACCAACC
(7,4)         CBACCAACC
(7,5)         CCBCAACC.                                  (10)
```

Reflection exchanges `(5,3)<->(7,5)` and `(6,3)<->(7,4)`, while fixing
`(6,4)` and `(7,3)`.  The two fixed coordinates are exactly the two entries
of THM-846's B seam

```text
x_(6,4) xor x_(7,3).                                     (11)
```

Thus the one-step seam is the central fixed pair of the entire asymptotic
return geometry, not an isolated coordinate accident.

For each map, count the rank-one sinks reachable by further words.  The exact
histogram is

```text
reachable sinks   maps
       6             1
       5             1
       4             4
       3         1,746
       2        31,342
       1        13,032.                                  (12)
```

The identity is the unique six-horizon map, B the unique five-horizon map,
and the four four-horizon maps have shortest words `A,C,AB,CB`.  Hence every
nonidentity return has already discarded at least one possible asymptotic
coordinate.  Rank records present information; the horizon set records which
information can survive any future continuation.

The generator imbalance is exact.  A and C each have 394 self-loops in the
Cayley graph, while B has 32,253.  Reflection enforces equality of the A/C
counts but says nothing that would make B comparable to them.  This is the
closed-return version of “blue symmetry with black directional imbalance”:
left/right endpoint roles mirror, while the central gap role is dynamically
sticky.

## 3. Complete Q audit on the 58-row bank

Across every monoid map and the 58 canonical representatives there are
6,138 output masks, or 8,940 after adjoining all-free-coordinate complements.
The number of distinct literal Q cells of one map ranges from one to 58.  Only
ten maps are Q-injective on the bank; their shortest words are

```text
I, A, B, C, AB, BA, BC, CB, BAB, BCB.                    (13)
```

Only `I` and `B` are simultaneously Q-injective and constant across all 58
old collision pairs.  There are 1,752 maps with full old-pair descent, but
all the other 1,750 have already collapsed some representative cells.  This
separates three predicates:

```text
literal coordinate rank,
injectivity on present Q cells,
descent through the old collision relation.              (14)
```

The complete two-step boundary is

| word | coordinate rank | Q fibre profile | old-pair descent |
|---|---:|---|---:|
| AA, CC | 15 | `24x1 + 9x2 + 4x4` | 0 |
| AB, CB | 21 | `58x1` | 0 |
| AC, CA | 15 | `38x1 + 10x2` | 4 |
| BA, BC | 20 | `58x1` | 0 |
| BB=B | 27 | `58x1` | 58 |

The A/C equality is forced by reflection.  Order matters: AB has rank 21
while BA has rank 20, although both remain Q-injective.  B is the unique
nonidentity map in the enumerated monoid that is both Q-injective and has full
old-pair descent.

## 4. Operation kernels and future equivalence

Apply THM-840 to the current Q equality partition at every Cayley vertex.
For a generator R, ask whether

```text
ker(Q after w) subset ker(Q after wR).                    (15)
```

The exact split-source-cell histograms over all 46,126 states are

```text
R=A: 0 splits on 46,050 maps; 1 on 60; 2 on 10; 5 on 6,
R=B: 0 splits on all 46,126 maps,
R=C: 0 splits on 46,050 maps; 1 on 60; 2 on 10; 5 on 6.   (16)
```

Thus at every monoid state, B induces a well-defined map from that state's
current row quotient to its B-successor row quotient.  A and C each fail on
precisely 76 states, in mirrored fashion.  This does not by itself make the
set of 144 present partitions a B-Markov quotient, because distinct monoid
states with the same current partition need not share one labelled successor.
It is nevertheless far stronger than THM-846's one-step B purity and gives an
exact algebraic sense in which transitivity flows through the middle face.
It remains a finite-bank statement: Q values outside the 58 rows are not
classified here.

Now forget the numerical Q labels and retain only the equality partition of
the 58 rows.  There are just 144 present partitions.  Refine two monoid
states whenever an A, B, or C successor has already been distinguished.
The exact Moore-refinement class counts are

```text
144, 502, 1,741, 4,367, 10,983, 17,598, 20,397,
20,419, 20,419.                                           (17)
```

Consequently there are 20,419 future-Q classes: two monoid maps are equivalent
precisely when their present equality partitions, and the equality partitions
after every further return word, agree.  Seven effective lookahead rounds
suffice.  The 144-state static partition summary is therefore not remotely a
Markov state for the full A/B/C action, even though each state-specific B
kernel inclusion in (16) holds.

This is a finite Myhill--Nerode construction on proof obligations.  It does
not assert that 20,419 numerical labels are intrinsically necessary for some
other predicate; it is the minimal deterministic state for this named
partition-output action.

## 5. The defect coordinate and the node boundary

Only 462 maps preserve rank four on THM-832's defect basis.  Exactly those 462
also keep all eleven occupied sectors distinct.  The complete defect-rank
census is

```text
rank 4:462, rank 3:236, rank 2:21,374,
rank 1:22,308, rank 0:1,746.                              (18)
```

Thus B's one-step pointwise fixation of the defect bank does not make defect
rank stable under arbitrary endpoint returns.  The six rank-one coordinate
horizons, the Q future classes, and the four-bit defect image are three
different summaries of the same monoid action.

This computation deliberately stops before classifying all 8,940 output
masks into size-nine tournament isomorphism classes.  THM-834/846 provide
exact local node maps, and THM-851 is independently targeting node-coloured
constant-composition refinement.  Extending the subset-DP canonicalizer to
the 8,940-mask return closure is the next node-level computation.  No bare
node claim is inferred from Q equality.

## 6. Tournament Analysis and preservation boundary

Tournament Analysis uses seven information carriers as vertices: coordinate
rank, visible support, reachable sink set, defect rank, Q-cell count, Q
equality partition, and old-pair descent.  Pairwise separation of monoid
maps is the observable.  Switching from raw retention to retention per
logarithmic cell cost flips eight of 21 edges.  Both gauges are transitive,
with no directed triangle, singleton SCCs, and one tie Hamiltonian path.
This ranks candidate state summaries; the carrier tournament is not the
return monoid or an LRC proof object.

The theorem preserves exact coordinate-copy composition, complement,
reflection conjugation, finite-bank Q equality, and the stated defect images.
It defers merged-node identity on the full return closure and destroys every
LRC gap, owner, wall, and metric-loneliness field.  It proves no implication
for LRC(14). ∎
