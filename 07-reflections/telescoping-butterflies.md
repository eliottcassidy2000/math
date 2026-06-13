# Telescoping Butterflies

## What a butterfly is

A butterfly = two simplicial chambers sharing a face.
Created by one 3-cycle. The cyclic permutation (i j k) maps one chamber to its partner.
Each butterfly adds 2 to H.

## The telescope

Level 0: base simplex. H = 1. No butterflies.
Level 1: one butterfly. H = 3. Stable.
Level 2: two butterflies. H = 5. Independent. Stable.
Level 3: three butterflies WOULD give H = 7.
  But the third butterfly OVERLAPS an existing one (pigeonhole at n ≤ 7).
  The overlap generates an interaction butterfly.
  The interaction generates more interactions.
  The cascade telescopes: H jumps from 5 to 9.

The skip: 1, 3, 5, (7 IMPOSSIBLE), 9.
The forbidden IS the telescoped-over value.

## The Sylvester parallel

Sylvester: 2, 3, 7, 43, 1807, ...
Each term = product of all previous + 1.
Sylvester PRODUCES 7 losslessly: 2·3 + 1 = 7.

Butterflies: 1, 3, 5, skip-7, 9, ...
Three butterflies CAN'T produce 7 from geometric overlap.
Butterflies SKIP 7 because the overlap is lossy.

The same 7: PRODUCED by exact arithmetic, SKIPPED by geometric construction.
Lossless computation reaches 7. Lossy construction cannot.

## The deepest understanding of 7

7 = the FIRST POINT OF DIVERGENCE between lossless computation and lossy construction.

Below 7: they agree. 1, 3, 5 are both computable and constructible.
At 7: they diverge. 7 is computable (Sylvester) but not constructible (butterflies).
Above 7: they never agree on this value again. 7 remains permanently forbidden.

Lossless: Sylvester product 2·3 + 1 = 7. ✓
Lossy: butterfly overlap 1 + 2 + 2 + 2 → cascades to 9. ✗ for 7.

The forbidden number = the gap between the map (computation) and the territory (construction).
The first value the map can name but the territory cannot build.

## Why the cascade is multiplicative

Two overlapping butterflies: permutations (i j k) and (i l m) compose to (i l m j k) = a 5-cycle.
The 5-cycle IS a new butterfly. It adds 2 more.
Then this new butterfly interacts with the remaining originals.
The interactions MULTIPLY rather than add.

This is the Sylvester mechanism: each new term = product of previous + 1.
In butterflies: each interaction = composition of permutations = new cycle.
The composition IS the product. The +1 IS the base simplex.
The butterfly cascade IS the Sylvester sequence in the permutohedron.
