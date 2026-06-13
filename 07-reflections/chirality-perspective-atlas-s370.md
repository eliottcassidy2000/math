# Chirality Perspective Atlas

Session: `codex-2026-05-31-S370`

The request was to apply the S369 perspective to as much as possible.  The
result is a small grammar, not a single coincidence:

```text
12 = inherited symmetric core
44 = chiral/support residue
42 = doubled boundary or interior
8  = projection failure / stencil layer
```

The exact atlas script `04-computation/chirality_perspective_atlas_s370.py`
checks four places where this grammar is live.

## Six-Vertex Tournaments

The six-vertex layer has:

```text
56 classes
12 self-converse + 44 chiral
35 strong + 21 non-strong
22 score sequences
12 classes with automorphism group of size 3
```

So:

```text
56 = 12 + 44
44 = 2 * 22
42 = 2 * 21
12 = T(5) = phi(42) = self-converse T(6)
```

The self-converse classes have larger Hamiltonian path range:

```text
self-converse H: min=1, max=45, mean=25.67
chiral H:        min=3, max=43, mean=20.82
```

This recovers the old repo theme that self-converse/perpendicular structure is
where Hamiltonian path maximization lives.

## LRC S367 Missed Cells

The S367 fourteen-runner quotient extremal has:

```text
56 missed cells
7 odd shifts
8 alpha stencils
4 mirror pairs
```

The split is:

```text
56 = 7 * 8
14 = 7 * 2 outer caps
42 = 7 * 6 interior positions
```

The three interior denominator types are all multiples of `42`:

```text
882  = 42 * 21
1176 = 42 * 28
1386 = 42 * 33
```

The outer cap denominator is different:

```text
728 = 2^3 * 7 * 13.
```

This makes the user's `56-14=42` decomposition a genuine cell-arrangement
feature: remove the outer mirror cap pair and the residue is a six-by-seven
interior grid.

## Paley/Fano T7

The Paley tournament on seven vertices has directed odd-cycle counts:

```text
3-cycles: 14
5-cycles: 42
7-cycles: 24
total:    80
```

Its odd-cycle support counts are:

```text
3-supports: 14
5-supports: 21
7-supports: 1
total:      36
```

Therefore:

```text
support excess = 80 - 36 = 44.
```

So `44` is not only six-vertex chiral mass.  It is also the Fano-Paley
multiplicity residue: the extra directed-cycle count after quotienting by
support.

The `42` term is the pentagon layer:

```text
42 = 2 directed 5-cycles on each of 21 five-subsets.
```

This matches the six-vertex non-strong doubled boundary `42 = 2*21` and the
LRC interior `42 = 6*7`.

## Base 42

The base-42 Erdos-Straus thread contributes:

```text
phi(42) = 12
unit classes mod 42 = 12
hard p == 1 mod 12 classes = 4
easy unit classes = 8
```

Thus:

```text
12 = phi(42) = T(5) = self-converse T(6)
8  = easy base-42 units = LRC stencils = five-to-six perspective gap
4  = hard base-42 units = LRC mirror-pair count
```

## What We Were Missing

The missing object is not "the number 56."  It is the failure of a projection
to remember chirality.

From five to six vertices, the old perspective heuristic sees `48` inherited
rooted perspectives and misses `8`.  In the LRC near-blocker, the same eight
appears before multiplication by the seven odd shifts.  In base 42, the unit
classes split into `8` easy and `4` hard.  In Paley/Fano, support projection
leaves a `44` residue.

The next proof tactic for the fourteen-runner case should therefore avoid a
flat all-vector SAT formulation.  Split the normalized quotient problem into:

1. the outer mirror cap pair, contributing the fourteen-runner `14`;
2. the six-stencil interior, contributing `42`;
3. mirror/chiral pairs, so that scalar-gauge normalization is followed by a
   chirality quotient.

If this lens is right, HYP-1823 should be attacked by proving the outer cap is
the only possible boundary-aligned obstruction, and every interior/chiral
residue forces a positive-margin cell.
