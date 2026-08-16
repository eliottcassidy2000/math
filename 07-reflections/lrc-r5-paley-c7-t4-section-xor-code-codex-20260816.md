# A folded heptagon gives an affine cube of tournaments, not one tournament

**Status: finite-exact structural sidecar to the folded-`C7` transporter
audit.**  This note sharpens the distinction between the `K4` carrier, its
six edges, a chosen tournament orientation, and a genuine current or ancestry
coordinate.  It adds no LRC(14) claim.

## Source, fold, and lost coordinate

Start with the Paley tournament on `F_7`, oriented by

```text
r -> s iff s-r in {1,2,4}.
```

Folding by negation gives four blocks

```text
O0={0}, O1={1,6}, O2={2,5}, O3={3,4}
```

and arc-count matrix

```text
0 1 1 1
1 1 2 2
1 2 1 2
1 2 2 1.
```

The off-diagonal support is `K4`, but the matrix is symmetric.  Negation has
erased orientation.  Thus the exact quotient is a weighted looped symmetric
digraph, not a tournament.

An orientation can be restored only by a section: choose one representative
from each nonzero pair.  Encode the three choices by `(x,y,z) in F_2^3`, where
zero selects `(1,2,3)` and one selects its negative.  In edge order

```text
01,02,03,12,13,23,
```

the induced tournament has the affine six-bit code

```text
(x,y,1 XOR z,x,z,y).                                  (1)
```

This is the precise sense in which the four-vertex tournament and its six
edges are the same object here: the eight lawful orientations occupy only an
affine three-cube inside `F_2^6`, not all `2^6` edge orientations.  Global
negation adds `(1,1,1)` to the section and reverses all six arcs.

## The XOR invariant

The eight sections contain six strongly connected tournaments, one 3-cycle
with a source, and its globally reversed 3-cycle with a sink.  No lawful
section is transitive.  If

```text
E=z XOR xy XOR xz XOR yz,                              (2)
```

then `E=1` exactly at `(0,0,1)` and `(1,1,0)`, the source/sink pair.  For the
skew adjacency matrix `S` and the number `c3` of cyclic triples,

```text
det(S)=1+8E,                 c3=2-E.                   (3)
```

Thus all eight skew matrices are invertible; the exceptional pair has
determinant nine and the six strong sections determinant one.  The section
type is governed by a quadratic Boolean polynomial, not a single linear
Walsh character.

## Why the dimension-three coincidence is not H1

Subtract the base orientation from (1).  The three gauge generators have
edge supports

```text
{01,12},       {02,23},       {03,13}.                 (4)
```

They are two-edge paths, not cycles.  Over `F_2`, their boundaries are the
three independent reduced vertex vectors

```text
0+2,           0+3,           0+1.
```

Consequently the boundary map is an isomorphism from the section-difference
space `G_section` to the reduced zero-chain space.  Since the `K4` cycle
space also has dimension three and intersects `G_section` trivially,

```text
C1(K4;F_2)=Z1(K4;F_2) direct-sum G_section.            (5)
```

Equation (5) is the typing answer: the section bits form a boundary
complement, not hidden absolute `H1` flux.  Their dimension agrees with the
cycle rank only because both halves of the six-edge chain space have
dimension three.

The three section paths in (4), the three Walsh perfect matchings

```text
01|23,          02|13,          03|12,
```

and the `O1 -> O2 -> O3 -> O1` multiplication-by-two cycle are three
different cardinality-three grammars.  None canonically identifies the other
two without an additional map.

## Transporter boundary

The folded-transporter audit proves that every selected THM-2594 source bank
has rank at most three while the U_full Walsh target has rank four.  A chosen
tournament matrix `S` cannot repair this: even though `S` is invertible,

```text
rank(S X C) <= rank(X)=3
```

for every common right operator `C`.  The section orientation is a left gauge
mixer, not a fourth common-base response channel.  The quadratic scalar `E`
classifies that gauge; it does not create ancestry, Boolean coupling, a
physical current, absolute `H1`, a bispectrum, or a scalar exclusion.

## Exact artifact

```text
04-computation/lrc_r5_paley_c7_t4_section_xor_code_audit_20260816.py
05-knowledge/results/lrc_r5_paley_c7_t4_section_xor_code_audit_20260816.out
```

The integer/Boolean-only companion exhausts all eight sections and all
sixty-four edge chains.  Normal and optimized transcripts are byte-identical;
its semantic digest is

```text
85640c02220b246781bab2e447a6a2f9fec7ca912d27c56377fbcd943bf041bc.
```
