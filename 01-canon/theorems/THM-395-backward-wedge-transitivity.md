---
id: THM-395
name: backward-wedge-transitivity
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S550
depends_on:
  - HYP-2040
---

# THM-395: Transitivity is exactly zero backward-wedge debt

## Statement

Let `T` be a finite tournament.  For every oriented edge `X -> Y`, define the
backward-wedge set

```text
W_T(X,Y) = { Z distinct from X,Y : Z -> X and Y -> Z }.
```

Then the following are equivalent:

1. `T` is transitive: whenever `X -> Y` and `Y -> Z`, then `X -> Z`.
2. Every oriented edge has empty backward wedge:

```text
W_T(X,Y) = empty for every edge X -> Y.
```

Moreover, if `c_3(T)` is the number of directed 3-cycles in `T`, then

```text
sum_{X -> Y} |W_T(X,Y)| = 3 * c_3(T).
```

Equivalently, each accepted transitive edge `X -> Y` carries the hidden local
clearance

```text
not (Z -> X and Y -> Z)
```

for every third vertex `Z`.

## Proof

Assume `T` is transitive and fix an edge `X -> Y`.  If some third vertex `Z`
satisfied `Z -> X` and `Y -> Z`, then transitivity applied to

```text
Z -> X -> Y
```

would force `Z -> Y`, contradicting `Y -> Z`.  Hence `W_T(X,Y)` is empty.

Conversely, assume every backward-wedge set is empty.  If `X -> Y` and
`Y -> Z`, then tournament completeness gives exactly one of `X -> Z` and
`Z -> X`.  If `Z -> X`, then `Z` lies in `W_T(X,Y)`, contradicting emptiness.
Therefore `X -> Z`, so `T` is transitive.

For the count, a membership

```text
Z in W_T(X,Y)
```

is exactly the directed triangle

```text
X -> Y -> Z -> X.
```

Each directed triangle contributes one such membership for each of its three
oriented edges, and a transitive triple contributes none.  Therefore the total
backward-wedge mass is `3*c_3(T)`.

## Cascade Corollary

For any cascade whose accepted relations are represented by tournament edges,
the scalar clearance product can be refined by edge-local anti-wedge factors.
For an accepted edge `X -> Y`, set

```text
C_{X,Y} = product_{Z distinct from X,Y} 1[(X -> Z) or (Z -> Y)].
```

Then

```text
C_{X,Y}=1  iff  W_T(X,Y)=empty.
```

A transitive cascade has every such factor equal to `1`; a zero factor is
exactly a directed-triangle debt attached to the accepted edge.  Thus a product
of conditional clearances does not only multiply scalar survival factors.  If
the proof route uses a transitive or hierarchical lift, each accepted clearance
also exports no-backward-wedge obligations to the remaining objects.

## LRC Use

This theorem is deliberately neutral about what the tournament vertices are.
For LRC, the vertices need not be runners.  HYP-2040 records the useful
candidate lifts: clearance obligations, endpoint owners, wall-crossing events,
p-adic zero branches, cover arcs, Fourier/Gabor modes, and proof obligations.
HYP-2041 identifies the same anti-wedge condition with no-return or 3-term
resonance debt.  HYP-2042 records the complementary warning: this 3-cycle layer
is necessary, but a full LRC proof must still control the whole order-`n`
conditional-clearance ladder.

In such a lift, a conditional LRC factor

```text
P_k = P(next obligation clears | previous obligations cleared)
```

should be paired with the anti-wedge ledger from this theorem.  AP/regular rows
are expected to align the wedge debt near the final runner or compact wall,
while non-AP rows should distribute it across earlier conditional clearances.
