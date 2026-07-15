---
id: THM-832
title: The n=9 false-palindrome defect is a rank-four constant Cech core on every nonempty face intersection
status: PROVED LINEAR CORE + FINITE-EXACT WITNESS-FIBRE AUDIT
source: codex-2026-07-15-S13
depends_on: [THM-801, THM-828]
related: [THM-553, THM-810, HYP-3234, HYP-6880]
verification:
  - 04-computation/n9_rank_four_cech_core_codex_S13.py
  - 05-knowledge/results/n9_rank_four_cech_core_codex_S13.out
---

# THM-832 — the rank-four defect is the common Cech core

Let `V` be the binary span of THM-828's four survivor-difference generators

```text
b0=0192486, b1=08c2c0c, b2=11b4600, b3=4483414.          (1)
```

Let `A,B,C` be the three size-eight B3 faces of the size-nine staircase.  For
every nonempty intersection `I` of these faces, restriction

```text
r_I: V -> F_2^I                                             (2)
```

is injective.  In particular every face, every pairwise overlap, and even the
ten-cell triple overlap sees all four defect coordinates.  After identifying
each image with `V`, the resulting Cech complex is the constant rank-four
system on a filled triangle:

```text
dim H^0=4,       H^1=H^2=0.                              (3)
```

Thus the four-dimensionality is neither created by one face nor hidden in a
boundary corner.  It is the redundantly stored common core of the three-face
recursion.

## Exact face and Venn ranks

The face sizes are 21, the pairwise intersections have size 15, and the
triple intersection has size 10.  Direct binary elimination gives

| chart `I` | `A` | `B` | `C` | `AB` | `AC` | `BC` | `ABC` |
|---|---:|---:|---:|---:|---:|---:|---:|
| cells | 21 | 21 | 21 | 15 | 15 | 15 | 10 |
| `rank r_I(V)` | 4 | 4 | 4 | 4 | 4 | 4 | 4 |

The exact seven Venn regions separate the location of the information:

| exact membership | `A` | `B` | `C` | `AB` | `AC` | `BC` | `ABC` |
|---|---:|---:|---:|---:|---:|---:|---:|
| cells | 1 | 1 | 1 | 5 | 5 | 5 | 10 |
| restricted rank | 0 | 0 | 0 | 3 | 2 | 3 | 4 |

Every defect vanishes on the three face-exclusive corner cells.  The
ten-cell common intersection by itself reconstructs it.  Numerically the old
full recursion is therefore literal at the rank level:

```text
A+B+C-AB-AC-BC+ABC = 4+4+4-4-4-4+4 = 4.               (4)
```

Rank is not an inclusion-exclusion valuation in general; (4) is a special
consequence of all seven intersection maps being injective on this `V`.

## An explicit decoder from every chart

Write a sector as

```text
D=c0 b0+c1 b1+c2 b2+c3 b3,
z=(c3 c2 c1 c0)_2.                                     (5)
```

For each chart, read the four listed tile bits of `D`, in order, as a nibble
`y=(y0,y1,y2,y3)`.  If the decoder word has hexadecimal nibbles
`r0 r1 r2 r3`, then

```text
ci = parity(ri bitand y).                               (6)
```

The following are exact left inverses of (2):

| chart | four pivot tiles | decoder |
|---|---|---|
| `A` | `(8,1),(7,1),(6,1),(7,2)` | `1487` |
| `B` | `(8,1),(7,1),(6,1),(7,2)` | `1487` |
| `AB` | `(8,1),(7,1),(6,1),(7,2)` | `1487` |
| `C` | `(9,2),(7,2),(6,2),(5,2)` | `182f` |
| `AC` | `(7,2),(6,2),(5,2),(4,2)` | `f418` |
| `BC` | `(9,2),(7,2),(6,2),(5,2)` | `182f` |
| `ABC` | `(7,2),(6,2),(5,2),(6,3)` | `941e` |

Equations (5)--(6) give a canonical function from any one restricted
difference word to its global defect-sector index.  The verifier exhausts
all sixteen words in `V` against every decoder.

## Exact Cech calculation

Trivialize the seven images using (1).  Then

```text
C0=V^3,        C1=V^3,        C2=V,
delta0(a,b,c)=(a+b,a+c,b+c),
delta1(p,q,r)=p+q+r,                                      (7)
```

where subtraction equals addition over `F_2`.  Exact elimination gives

```text
rank(delta0)=8,       rank(delta1)=4.                    (8)
```

The composition vanishes, so (3) follows.  There is no linear phase
holonomy: a compatible local defect has one global four-bit coordinate.
THM-828's four punctures must therefore be imposed by the nonlinear
metagraph-face and literal-histogram predicates, exactly as its gate audit
found, not by failure to glue the linear difference carrier.

## Mapping sectors to their literal fibres

For a THM-828 collision pair `{u,v}`, the forward map is

```text
{u,v} -> D=u xor v -> z(D),                              (9)
```

where `z` can be decoded from any chart by (6).  Conversely, for a sector
`z`, its exact literal fibre is obtained by filtering the stored witness table
for the corresponding `D`.  The eleven occupied sector rows are

| `D` | `z` | pair mass | weight | exact-region weights `(AB,AC,BC,ABC)` |
|---|---:|---:|---:|---|
| `0192486` | `0001` | 6 | 8 | `(2,0,2,4)` |
| `08c2c0c` | `0010` | 2 | 8 | `(2,0,2,4)` |
| `095088a` | `0011` | 2 | 8 | `(2,0,2,4)` |
| `11b4600` | `0100` | 2 | 8 | `(0,2,0,6)` |
| `18e4e8a` | `0111` | 4 | 12 | `(2,2,2,6)` |
| `1976a0c` | `0110` | 4 | 12 | `(2,2,2,6)` |
| `4483414` | `1000` | 4 | 8 | `(2,2,2,2)` |
| `4511092` | `1001` | 2 | 8 | `(2,2,2,2)` |
| `4c41818` | `1010` | 26 | 8 | `(2,2,2,2)` |
| `5c67a9e` | `1111` | 2 | 16 | `(4,4,4,4)` |
| `5df5e18` | `1110` | 4 | 16 | `(2,4,2,8)` |

The other four nonzero words have coordinates

```text
1026286:0101, 4dd3c9e:1011, 54a5692:1101, 5537214:1100 (10)
```

and empty final fibres.  This realizes an objective cube order and a pair of
exact functions between collision pairs and sector fibres.  It does **not**
identify a sector with one tournament isomorphism-class node: THM-828 proves
that all 58 ordered six-node face signatures are distinct.  A sector groups
several address cells by a shared reflection-difference direction.

## The nonlinear no-go

The exact Venn support profile does not determine fibre mass.  Profile
`(2,0,2,4)` occurs with masses two and six; profile `(2,2,2,2)` occurs with
masses two, four, and 26.  In particular neither the four coordinates nor
their regional Hamming weights replace the literal `H_8` kernel rows and raw
`S2` base realization.

This precisely separates the three recursive levels:

1. the rank-four difference coordinate is a constant, holonomy-free Cech
   core;
2. the nonlinear `H_8` relation decides which local bases can carry it;
3. literal `S2` decides which compatible bases realize equal histograms.

All 58 realized pairs are pure black, but the linear space `V` itself has no
blue/black edge type.  Colour belongs to the witness fibre, not to a bare
difference word.

## Tournament Analysis and preservation boundary

Use the eleven occupied sectors as the alternate vertex set.  The structural
and empirical gauges of THM-828 order the observable
`(first skew layer,weight,mass)` differently, with hexadecimal `D` as the tie
Hamiltonian path.  Both tournaments are transitive; 22 edges flip.  This
orders the sector fibres without mislabelling cube-coordinate adjacency as a
metagraph flip.

The carrier preserves the difference word, its four coordinates, every face
restriction, and the exact witness fibre supplied by (9).  A bare sector
destroys the base tiling, endpoint-node identities, `H_8` kernel truth, raw-S2
realizability, curvature, and the LRC metric predicate.  The challenged
assumption is again that recursive vertices must be runners or tournament
nodes: here the useful vertices are defect sectors, but only when decorated
by literal fibres.

The next computation is now smaller and sharper.  Apply the centered-CF
coordinate-copy action to these four generators and their 58 literal
reflection-orbit witnesses; test separately whether the linear core is
invariant and whether its nonlinear fibres descend. ∎
