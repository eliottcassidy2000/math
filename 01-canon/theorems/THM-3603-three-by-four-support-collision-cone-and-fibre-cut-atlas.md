---
id: THM-3603
title: "Three-by-four additive-support collision cones and fibre-cut atlas"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For positive three-point and four-point supports, the complete
  connected additive-collision locus consists of 83 rational cones and 149
  oriented sum-fibre words.  Simultaneous reversal has three fixed words and
  73 two-cycles and is retained, not quotiented.  Every connected word of
  sumset size eight or nine loses projection connectivity after deletion of
  one collision fibre; the unique size-six word and exactly four size-seven
  words have cut debt two.  The only coefficient deletion used is the inherited
  deletion of whole disconnected zero-bracket components; no deletion follows
  merely from a one-fibre cut.  No new Darboux nonentry or JC(2) conclusion is
  claimed.
source: kps-s188 + agent Anscombe / open exponent-two 3x4 support sector, 2026-08-21
audit: >
  PASS.  An independent SymPy reconstruction recovered 786 flats, maximal
  connected dimension two, 83 cones, 149 words, the full cone/chamber and
  reversal counts, and every fibre-cut, cut-debt, cut-size, and exposure count.
  An independent gap-at-most-five scan realizes exactly all 149 words and all
  83 closure masks.  Ordinary and optimized companions are byte-identical to
  the stored transcript after LF normalization.  The audit caught and repaired
  an expository scope contradiction: whole disconnected components are lawfully
  deleted, but a one-fibre projection cut is never treated as coefficient
  deletion.
depends_on:
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
related:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_three_by_four_additive_support_atlas_thm3603.py
output: 05-knowledge/results/jc2_three_by_four_additive_support_atlas_thm3603.out
script_sha256: 85c8a5e8200e0999b0e671f2419f3b9a52f03c862370cfa5a6afc7aaf8715c74
output_sha256: 531e7edf2b8e2af53f6b74266d36a1bd9fcaa3ddd8a3ac374665513299fbfc44
flat_semantic_sha256: 2c5e4c1f448cfb96ff6f4f240d615339f26a853ea7406e7801cf3f581e5cbeb9
cone_semantic_sha256: a29b4a617a677ba6b6e1778e5a1e94eb47de8f9d4cb8a0df8fd25f1017c9e8f0
word_semantic_sha256: 1c32e4a56baaaecc56014899bf70f8f4c051ec0a31f377c2400b7c67c544eca8
bounded_semantic_sha256: 0f8ace8fa981601cdc0a51e9d03e10db232b0c5ee23a99a6b999bbbd7d8dfc6d
hash_basis: LF-normalized bytes
---

# THM-3603 -- three-by-four additive-support collision cones and fibre-cut atlas

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This classifies the additive skeleton of the first still-open
seven-piece exponent-two Darboux sector.  It deliberately stops before the
coefficient and arm-regularity equations.

## 1. Supports, collision arrangement, and orientation

Translate independently and write

```text
A={0,X,X+Y},
B={0,U,U+V,U+V+W},                 X,Y,U,V,W>0.          (1)
```

The three positive differences of `A` are

```text
X, X+Y, Y,
```

and the six positive differences of `B` are

```text
U, U+V, U+V+W, V, V+W, W.                              (2)
```

A collision between two distinct cells of the `3 x 4` address rectangle is
equivalent to equality of one difference from each list.  Hence all collision
types lie in the central arrangement of the 18 rational hyperplanes

```text
A_ij=B_kl.                                               (3)
```

The three-point support is named first.  Simultaneous reversal acts by

```text
(X,Y;U,V,W) -> (Y,X;W,V,U).                              (4)
```

It is retained as an involution.  THM-3592's failed first audit shows why
this matters downstream: reversal of abstract weights need not preserve the
signed regularity floors on a Danielewski surface.

## 2. Exact global cone and word atlas

For each sum `s`, let

```text
F_s={(i,j):a_i+b_j=s}.                                   (5)
```

Join two `A` indices when some collision fibre contains both, and define the
analogous projection graph on the four `B` indices.  Call a word connected
when both projection graphs are connected.

The companion enumerates the row-space lattice of `(3)` over `Q`.  There are
exactly

```text
786 rational flats.                                     (6)
```

For every flat it computes the exact closure mask, rejects disconnected
projection graphs, and intersects its nullspace with the strict positive
orthant.  Every surviving connected flat has dimension at most two.  In
dimension one the positive ray is immediate.  In dimension two, the
remaining collision hyperplanes and the five coordinate walls form a central
line arrangement; one exact primitive point between every adjacent pair of
rays samples every oriented chamber.  Thus this is a global chamber
enumeration, not stabilization inferred from a bounded search.

The result is

```text
83 positive connected rational cones,
149 oriented sum-fibre words.                            (7)
```

Their distribution by sumset size `m=|A+B|`, flat rank, and dimension is

| `m` | rank | dimension | cones | oriented words |
|---:|---:|---:|---:|---:|
| 6 | 4 | 1 | 1 | 1 |
| 7 | 4 | 1 | 7 | 7 |
| 8 | 4 | 1 | 11 | 11 |
| 8 | 3 | 2 | 12 | 26 |
| 9 | 3 | 2 | 52 | 104 |

Equivalently, 25 cones have one chamber-word, 50 have two, and eight have
three.  No connected positive word has sumset size at least ten.

Reversal `(4)` has exactly three fixed words, one each at sizes six, seven,
and eight.  The other 146 words form 73 two-cycles, for 76 reversal orbits in
all.

## 3. Collision identity

Put

```text
C(A,B)=sum_s binom(|F_s|,2),
T(A,B)=#{s:|F_s|=3}.                                    (8)
```

No fibre has more than three cells.  Comparing the twelve addresses with
the number of occupied fibres gives the exact identity

```text
|A+B|=12-C(A,B)+T(A,B).                                 (9)
```

Indeed a double fibre loses one address and contributes one to `C`, while a
triple loses two addresses and contributes three, corrected by `+T`.

The unique size-six word is the common-step ray

```text
(X,Y;U,V,W)=(1,1;1,1,1),
00 | 01=10 | 02=11=20 | 03=12=21 | 13=22 | 23.         (10)
```

The seven size-seven words are exactly the rays represented by

```text
(1,1;1,1,2), (1,2;1,1,1), (1,2;1,1,2),
(2,1;2,1,1), (1,1;1,2,1), (2,1;1,1,1), (1,1;2,1,1).   (11)
```

The transcript records every closure mask, cone, sample, fibre word, profile,
and reversal partner.

## 4. Fibre-cut theorem

Delete a collection of collision fibres from `(5)` and retain all other
collision edges in the two projection graphs.  Define `d_A`, `d_B`, and `d`
as the least number of deleted collision fibres needed to disconnect the
`A` graph, the `B` graph, or either graph, respectively.  Then

```text
d=min(d_A,d_B).                                          (12)
```

The exact word census is:

| sumset size | words | `d=1` | `d=2` | exceptional cut debt |
|---:|---:|---:|---:|---|
| 6 | 1 | 0 | 1 | `(d_A,d_B)=(3,2)` |
| 7 | 7 | 3 | 4 | two have `(2,2)`, two have `(3,2)` |
| 8 | 37 | 37 | 0 | none |
| 9 | 104 | 104 | 0 | none |

Therefore:

> **Every connected oriented `3 x 4` word with `|A+B|>=8` has a
> one-collision-fibre projection cut.  The only cut-debt-two words are the
> common-step size-six word and four of the seven size-seven rays.**

Among the size-eight words, a one-fibre cut can be a double fibre in 27 words
and a triple fibre in ten.  Among the size-nine words, 86 admit a double cut
and 18 admit both double and triple cuts.  A stricter rectangle-exposure test
holds for 58 of the 104 size-nine words and for none at size eight.

This language is intentionally exact: deleting a *collision fibre* only
removes its collision edges from the combinatorial projection graph.  It does
not delete an input weight, a coefficient polynomial, or a bracket equation.
Rectangle exposure is likewise only a candidate filter.  Neither operation
is yet a lawful Darboux reduction.

## 5. Independent bounded control

As a hostile check, the companion separately scans every primitive positive
five-gap vector with all gaps at most 16.  It finds

```text
1,012,441 primitive vectors,
4,617 connected vectors,
connected counts by |A+B|: 1, 7, 1265, 3344,
distinct connected words: 1, 7, 37, 104.                (13)
```

The bounded word set equals the global 149-word chamber atlas exactly, and
every global word has a primitive witness whose largest gap is at most five.
This confirms the arrangement computation by a disjoint representation.  It
does not prove exhaustiveness; exhaustiveness comes from the rational-flat and
chamber argument in Section 2.

## 6. Darboux boundary and the next gate

THM-3592 and THM-3583 leave `(3,4)` and its transpose among the minimal
seven-piece exponent-two sectors.  Their component-deletion argument gives a
precise corollary here.  If either projection graph is disconnected, retain
the unique component containing the scalar row and delete every other whole
component, whose bracket is zero on its disjoint set of rows.  A scalar fibre
contains at least two distinct indices on each side.  Hence a disconnected
three-vertex `A` graph leaves exactly two `A` weights in its scalar component
and reduces to the `2 x <=4` theorem THM-3583.  A disconnected four-vertex `B`
graph leaves either two `B` weights, again covered by THM-3583 after exchange,
or three, covered by the `3 x 3` theorem THM-3592.  Consequently every
putative exponent-two `3 x 4` Darboux pair has one of the 149 connected
oriented words in this atlas.

This is the only coefficient deletion used in the theorem.  In particular,
the one-fibre cuts of Section 4 are not components and are not deleted.

This removes the unbounded additive-support search from that cell.  What
remains is finite but algebraic: for each word, impose the bracket
cancellation equations, the scalar-arm address condition, and the all-arm
Hermite--Pade divisibility of THM-3600.

The cut theorem suggests the right order of attack.  At `m>=8`, first isolate
a cut fibre and eliminate its bracket equation.  The projection graph on one
side then splits, exposing a candidate Euler factor or a multiarm
Hermite--Pade degree invoice.  At `m<=7`, attack the eight explicit rays
directly.  A proof still has to show that coefficient cancellation respects
this combinatorial cut; that implication is open and is not asserted here.

## 7. Verification

Reproduce with

```bash
python3 04-computation/jc2_three_by_four_additive_support_atlas_thm3603.py
python3 -O 04-computation/jc2_three_by_four_additive_support_atlas_thm3603.py
```

Both runs are byte-identical after LF normalization and match the stored
transcript.  The four semantic digests separately pin the flat lattice, cone
atlas, oriented-word ledger, and bounded hostile census.

This theorem proves no coefficient-system emptiness, no polynomial Darboux
nonentry on a Danielewski surface, no planar Keller pair, and no result on the
planar Jacobian conjecture.

**QED.**
