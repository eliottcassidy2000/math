---
source: codex-2026-06-03-S575
status: reflection + parity-ladder audit
tags: [tournaments, simplex, polygon, dihedral, LRC, nplus2, round, boundary]
---

# Tournament Simplex/Polygon Ladder

The user's sentence is right, with one important precision:

```text
tournaments are simplex meshes when all edge/face choices matter;
tournaments are polygon outsides when circular order determines the arcs.
```

The regular polygon is the switch.

For odd `m`, the regular `m`-gon has no antipodal pairs.  Its half-turn
orientation is a genuine tournament, the rotational `R_m`.  Rotations are
automorphisms.  Reflections reverse all arcs, so they are anti-automorphisms.
Together they form the dihedral automorphism-plus-converse object.

For even `m`, the regular `m`-gon has antipodal pairs.  The outside circle is
not enough.  There are `m/2` tie diagonals, and a full tournament needs mesh
choices for those diagonals.

## LRC Parity

LRC at parameter `n` uses `m=n-1` runners.

So:

```text
even LRC n -> odd m -> clean regular polygon tournament -> dihedral/converse
odd LRC n  -> even m -> regular polygon wall -> tie-resolution mesh
```

That is the every-other-`n` occurrence of the dihedral line.  The recursion
`n -> n+2` keeps you on the same line:

```text
clean line: D_ext grows by 4
wall line:  tie pairs grow by 1 and mesh choices double
```

This also explains why the same tournament object can feel like an outside in
one argument and a mesh in another.  Open LRC times see the round necklace body:
the outside circular order is enough.  Boundary/tight times see tie and
reflection data: the simplex mesh comes back.

## Connection Back

HYP-1998: open runner tournaments are round necklaces.

HYP-2086: the hard regime is Burnside fixed/self-converse boundary.

HYP-2087: the time-word dual lives in a dihedral/cosine fixed sector.

HYP-2089: the strong lens sees the tight object as a regular rotational
encirclement.

The incoming round/converse quotient note says open-time reversal sends a
round runner tournament to its converse, so even the open polygon body wants a
reflection quotient before attaching to the merged graph.  That is the same
outside symmetry seen from the class-counting side.

The compact HYP-2090 simplex/polygon parity note states the dichotomy.  HYP-2091
adds the explicit `n+2` ledger: clean polygon behavior on even LRC `n`,
mesh/tie behavior on odd LRC `n`, and the growth rules on each ladder.

**Artifacts:** `04-computation/tournament_simplex_polygon_recursion_s575.py`,
`05-knowledge/results/tournament_simplex_polygon_recursion_s575.out`, HYP-2091.
