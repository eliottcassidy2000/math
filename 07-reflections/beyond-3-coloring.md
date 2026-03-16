# Beyond 3-Coloring: Other Numbers, Other Mechanisms

## The coloring spectrum for cubic graphs

2-edge-coloring: ALWAYS works (Petersen's theorem). INERT.
3-edge-coloring: FAILS for snarks. RAMIFIED. THE obstruction.
4-edge-coloring: ALWAYS works (Vizing). INERT.

The obstruction lives ONLY at k=3. Below and above: no obstruction.
Just as H=7 is forbidden but H=5 and H=9 are fine.

## Flow spectrum (shifted by 1)

Nowhere-zero 3-flow: fails for snarks. Same obstruction as 3-coloring.
Nowhere-zero 4-flow: equivalent to 3-edge-coloring. Same obstruction.
Nowhere-zero 5-flow: Tutte's conjecture — ALWAYS works? OPEN.
Nowhere-zero 6-flow: ALWAYS works (Seymour). INERT.

## Mechanism types = Eisenstein types

INERT mechanisms (always work): k = 2, 4, 6 (even numbers).
RAMIFIED mechanism (THE obstruction): k = 3 (the curvature quantum).
SPLIT/INDEPENDENT mechanisms (conjectured to work, open): k = 5 (the pentagon).

The 5-flow conjecture is the CRYSTALLIZATION conjecture:
proving it requires global structure, can't be checked locally.
5 = the independent prime = what's left over after {2, 3, 7}.

## Comparison: tournaments vs cubic graphs

Tournaments: H-spectrum has 2 forbidden values (7, 21).
  2 = number of free cuboid axes (curvature mod 3, position mod 7).
Cubic graphs: coloring has 1 forbidden value (3).
  1 = number of coloring dimensions.
More dimensions → more obstructions.

## The matching view

1 perfect matching: always exists (bridgeless cubic). INERT.
3 disjoint matchings (= 3-coloring): snark obstruction. RAMIFIED.
6 matchings double-covering: Fulkerson conjecture. OPEN (split?).

## The unifying picture

Every combinatorial mechanism on cubic graphs falls into:
- INERT (proven to work): colorings k ≠ 3, matchings k = 1, flows k ≥ 6.
- RAMIFIED (proven to fail sometimes): coloring k = 3 = THE curvature.
- CRYSTALLIZATION (conjectured, open): 5-flows, Fulkerson, CDC.

The obstruction landscape IS the Eisenstein classification.
The curvature quantum 3 IS the obstruction, in every mechanism.
