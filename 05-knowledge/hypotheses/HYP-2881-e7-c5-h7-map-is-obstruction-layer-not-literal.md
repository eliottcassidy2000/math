---
id: HYP-2881
title: E7 metagraph C5 holes live one quotient level above the directed C5 support that realizes the H=7 K3 obstruction
status: SUPPORTED by exact S100 audit, integrated with S37. Directed C5 support = H=7 K3 support in cycle space; E7 metagraph C5-hole = five-class quotient cycle, not that single support object.
source: codex-2026-06-22-S100
tags: [tournaments, even-graphs, e7, odd-holes, h7, ocf, strong-components, euler-totient, multiplicative-functions, tournament-analysis]
related:
  - HYP-2878
  - HYP-2879
  - HYP-+2879
  - HYP-+2880
  - HYP-2880
  - HYP-2877
  - HYP-2872
  - HYP-2629
  - HYP-2660
  - THM-200
  - THM-201
---

# HYP-2881: E7 metagraph C5 holes are quotient-level analogues, not the directed H=7 C5 support

The S31e lead asked whether E7's `C_5` holes might map, through the
even-graph/tournament-cycle bijection, onto the H=7 `Omega=K_3` obstruction.
After S37, the answer has two levels:

```text
directed C5 support = H=7 K3 support in the cycle space;
E7 metagraph C5 hole != that support object.
```

The S100 audit proves the second line.

The fixed-path gauge is:

```text
0 -> 1 -> ... -> 6,
free bit (i,j), j>i+1, equals 1 iff j -> i.
```

A reversed free arc creates the directed cycle

```text
i -> i+1 -> ... -> j -> i,
```

whose undirected support is exactly the path-fundamental cycle used in the E7
metagraph.  This is the right labelled cycle-space bijection for the test.

Exact audit `e7_c5_h7_obstruction_map_codex_s100.py` verifies:

```text
E7 classes: 54
E7 metagraph edges: 951
Chordless C5 holes: 1496
C5 holes touch classes: 48/54
```

The H=7 obstruction data in the same fixed-path cube is:

```text
forbidden_h7_point alpha=(3,0): 0 masks, 0 classes
near_h7_P3:                         0 masks, 0 classes in this gauge
k3_forces_pentagon:                54 masks, 5 classes
```

Alignment with E7 `C_5` holes:

```text
k3_forces_pentagon classes hit 835/1496 C5 holes
flagged incidences: 1180/7480 = 15.78%
no C5 hole is entirely inside the flagged classes
```

The directed pentagon support itself maps to a single E7 class:

```text
undirected 5-cycle support class: class 3
edge_count=5
C5-hole incidence=209
```

So the metagraph-hole equality cannot be correct.  An E7 `C_5` hole is a
five-class cycle in the quotient metagraph.  The directed pentagon forced in
THM-200 is a single even-graph support class under the same cycle-space map.
This is compatible with S37: the directed C5 support itself is the XOR of three
pairwise vertex-conflicting triangles, hence the literal H=7/K3 cycle-space
support.

## What survives

The thematic bridge remains useful, but it has to be phrased with quotient
levels explicit:

```text
scalar forgetting -> cycle-space closure -> first odd-cycle obstruction.
```

For H=7, scalar OCF data wants the missing point `alpha=(3,0)`, equivalently
`Omega=K_3`.  Tournament incidence closure blocks it: three mutually-conflicting
triangles force a directed pentagon support, hence extra odd-cycle mass.

For E7, the even-graph quotient first becomes non-chordal at the same apex
prime and produces metagraph pentagons.  These are not the directed pentagon
support; they are quotient cycles through support classes.  Both mark the same
first odd-cycle obstruction layer after scalar or quotient data forgets
incidence.

This also explains why HYP-2872 remains the guardrail: even-graph minors do not
preserve `H`.  The E7 metagraph is an address quotient for obstruction
geometry, not the `H`-complete object.

## Multiplicative-function reading

Since `H` is multiplicative over strong components, the correct closure for
forbidden values is a strong-atom Euler product, not an even-graph minor order.
The older Euler-copy/totient work gives the right warning.

In the LRC residue setting,

```text
sum_{d|n} phi(d) = n
```

means `phi(d)` counts exact-denominator packets, and the main term in
`N(S,D)` is scaled by `phi(D)`.

In the tournament setting, strong components give the analogous multiplicative
packet law:

```text
H(T) = product_i H(C_i).
```

But scalar arithmetic factorization is not the atom factorization.  S36 verifies
`49` and `75` are single irreducible strong atoms, even though their integer
values are composite.  Conversely, `21=3*7` is absent because the strong atom
`7` is absent.  Thus the useful multiplicative function is not ordinary
integer primality; it is the Dirichlet-convolution/Euler-product ledger of
strong-H atoms.

The next arithmetic object should be an atom-count function `a(h)` for strong
irreducible `H` values, with product decompositions handled by multiplicative
partitions of `h`.  This is the H-spectrum analogue of the exact-period
totient packet law: keep the packet labels before scalarizing.

## Tournament Analysis / Assumption Challenge

Candidate vertices considered:

```text
E7 C5 metagraph holes
directed C5 support class
Omega=K3 forbidden point
near-H7 P3 replacement
k3-forces-pentagon classes
strong-component atoms
Euler/totient exact-period packets
raw H scalar
```

Pairwise observable: whether the quotient preserves the predicate "blocks
`H=7` by forcing extra odd-cycle mass."

Resulting order:

```text
strong-component atoms
> OCF alpha/incidence packets
> k3-forces-pentagon classes
> directed C5 support class
> E7 C5 metagraph holes
> Euler/totient packet analogy
> raw H scalar.
```

Challenged assumption: an E7 metagraph pentagon hole and the H=7 directed
pentagon support are the same object under the cycle-space bijection.  The
audit says no: they live at different quotient levels.  The preserved insight
is that the directed support is the H=7/K3 cycle-space object, and the E7
metagraph hole is an odd-hole incidence object one layer above it.

## Artifacts

- `04-computation/e7_c5_h7_obstruction_map_codex_s100.py`
- `05-knowledge/results/e7_c5_h7_obstruction_map_codex_s100.out`
- `07-reflections/e7-c5-h7-obstruction-layer-and-euler-totient-atoms-codex-s100.md`
