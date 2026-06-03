---
id: HYP-2091
status: SYNTHESIS - parity ladder audited; proof use open
source: codex-2026-06-03-S575
related:
  - HYP-1998
  - HYP-2090
  - HYP-2089
  - HYP-2087
  - HYP-2086
  - HYP-1881
---

# HYP-2091: tournaments alternate between polygon-outside and simplex-mesh behavior along the n+2 ladder

A tournament has two compatible readings.

```text
simplex/mesh:   orient every edge of the complete graph, the 1-skeleton of an
                (m-1)-simplex; boundary ambiguities require interior edge data.

polygon/outside: place m points on a circle and orient by the half-turn order;
                 round tournaments are circular necklaces.
```

For the regular `m`-gon these readings separate by parity.

If `m` is odd, the polygon outside determines the rotational tournament `R_m`
with connection set `{1,...,(m-1)/2}`.  Its automorphism group contains the
rotations `C_m`, and its reflections are anti-automorphisms, so the
automorphism-plus-converse symmetry is dihedral of size `2m`.

If `m` is even, the regular polygon is not a tournament: it has `m/2`
antipodal tie diagonals.  To turn the outside into a tournament, the simplex
mesh must choose orientations for those ties, giving `2^(m/2)` labelled
tie-resolution choices.

Since LRC `n` uses runner tournament size `m=n-1`, clean dihedral polygon
behavior occurs exactly at even LRC `n`, every other value.  The recursion
`n -> n+2` preserves the ladder:

```text
even n: clean polygon ladder, m odd, D_ext grows by 4 each step
odd n:  wall/mesh ladder, m even, tie pairs grow by 1 and choices double
```

## Evidence

`tournament_simplex_polygon_recursion_s575.py` audits `n=4..18`.

Clean polygon ladder:

```text
n = 4,6,8,10,12,14,16,18
m = 3,5,7,9,11,13,15,17
D_ext = 6,10,14,18,22,26,30,34
tie choices = 1 throughout
```

Wall/mesh ladder:

```text
n = 5,7,9,11,13,15,17
m = 4,6,8,10,12,14,16
tie pairs = 2,3,4,5,6,7,8
tie choices = 4,8,16,32,64,128,256
```

The same audit records simplex dimension and edge growth:

```text
n -> n+2 means m -> m+2, simplex dimension -> +2,
and simplex edge count increases by 2m+1.
```

The round necklace body becomes tiny inside the full simplex mesh.  At the
LRC `n=14` row (`m=13`), `round/all = 10^-11.19`; at `n=18` (`m=17`),
`round/all = 10^-22.80`.

## Interpretation

This refines HYP-1998, HYP-2086, HYP-2087, and HYP-2089.

HYP-1998 says the open-time runner tournament lives in the round necklace body.
HYP-2086/2087 say the hard boundary lives on the reflection/self-converse
fixed side.  HYP-2089 says the strong-lens tight object is the regular
rotational encirclement.

The incoming round-body converse-merged quotient note
`HYP-2089-lrc-round-body-converse-merged-quotient.md` fits this same picture:
open-time reversal sends `T_runner(t)` to `T_runner(t)^op`, so the outside
polygon body naturally enters through a converse quotient.

The incoming `HYP-2090-lrc-simplex-polygon-parity.md` gives the compact parity
dichotomy: the dihedral face appears on odd runner sizes and hence even LRC
`n`.  HYP-2091 extends that base statement with the explicit `n -> n+2`
ledger:

```text
the regular rotational encirclement exists as a clean tournament exactly on
the even-LRC n ladder.
```

On the odd-LRC ladder, the regular outside is a wall; proof data must include
the mesh/tie-resolution layer.  This is why the same tournament can sometimes
behave like only its outside circular order, and sometimes like the full
simplex of edge/face obligations.

## Tournament Analysis

Vertices were LRC `n`-ladder entries.

Pairwise observable:

```text
(clean polygon flag, extended dihedral size, tie choices, round/all fraction,
regular 3-cycle mesh count).
```

Switch/gauge:

```text
cleaner polygon quotient wins; within a parity ladder, larger extended
dihedral size wins, then fewer tie choices, then smaller round/all fraction.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Assumption Challenge

Possible vertices considered:

```text
tournaments, simplex faces, polygon vertices, antipodal tie diagonals,
boundary resolutions, dihedral orbits, SCC blocks, and LRC proof obligations.
```

Chosen quotient:

```text
LRC n-ladder entries with runner size m=n-1.
```

Predicate preserved:

```text
regular polygon parity, clean outside-vs-wall status, dihedral
automorphism/anti-automorphism size, tie-resolution mesh choices, and the
round necklace fraction inside all tournament classes.
```

Information destroyed:

```text
specific runner ownership, wall endpoint provenance, exact non-regular
iso-class identity, and speed-set arithmetic.
```

Challenged assumption:

```text
dihedral symmetry is an ordinary tournament automorphism at every n.
```

It is not.  For odd `m`, reflections are anti-automorphisms of `R_m`; for even
`m`, the regular polygon is a wall rather than a tournament.

## Files

- `04-computation/tournament_simplex_polygon_recursion_s575.py`
- `05-knowledge/results/tournament_simplex_polygon_recursion_s575.out`
- `07-reflections/tournament-simplex-polygon-nplus2-ladder-s575.md`
