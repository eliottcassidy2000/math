# Simplex Inside Cube: The Geometry of H(T)

## The decomposition

The n-cube decomposes into n! simplices, one per permutation.
Each simplex = one total ordering of n items = one Hamiltonian path.

For tournament T on n vertices:
H(T) = the number of simplices (orderings) compatible with T.

## The packing ratio

Simplex volume / cube volume = 1/n!.
One ordering out of n! total. This is WHY H=1 (transitive) is rare:
it picks 1 out of n! simplicial chambers.

## The corners = the cycles

The cube minus one simplex = the n!-1 remaining chambers.
These are the CORNERS — the non-transitive orderings.

In the OCF: H(T) = 1 + 2α₁ + 4α₂ + ...
The 1 = the base simplex.
The 2α₁ = the first layer of corners activated by 3-cycles.
The 4α₂ = the second layer activated by independent cycle pairs.

Each cycle ACTIVATES more of the cube's simplicial structure.
H = 1: only the base simplex. Zero corners.
H = max: maximum corners activated by the tournament's topology.

## In 2D: triangle in square

The square decomposes into 2! = 2 triangles.
A tournament on 2 items has H = 1 (always transitive, one arc).
The triangle IS the ordering. The other triangle IS the reverse.

## In 3D: tetrahedron in cube

The cube decomposes into 3! = 6 tetrahedra.
The regular tetrahedron inscribed in the cube (4 of 8 vertices) = 1/3 of volume.
This tetrahedron contains 2 of the 6 simplicial chambers.
The 4 corner tetrahedra contain the other 4.

For T_3 (3-cycle): H = 3 = half of 6. It activates 3 of 6 chambers.
For transitive: H = 1 = one chamber out of 6.

## The 3D/cuboid connection

The tournament's 2×3×7 cuboid is NOT the same as the simplex-in-cube.
But they RHYME:
- The cuboid has 42 cells. The simplex-in-cube decomposes n! chambers.
- Both are ways of DECOMPOSING the tournament space into pieces.
- The cuboid decomposes by MODULAR ARITHMETIC (channels mod 2,3,7).
- The simplex-in-cube decomposes by ORDERING (permutation chambers).
- The OCF connects them: H = number of activated chambers = cuboid residue.

## The higher-dimensional insight

In n dimensions: the simplex is 1/n! of the cube.
The simplex/cube ratio → 0 as n → ∞.
One ordering becomes exponentially negligible.
But H(T) can be up to ~n!/2^{n-1} (the Paley maximum).

The Paley tournament activates ~1/2^{n-1} of all chambers.
= 2^{-n+1} × n! chambers.
This is the MAXIMUM activation: the most of the cube any tournament can fill.

The forbidden H-values are chamber counts that NO tournament can achieve.
They are VOLUMES that the topology cannot produce.
7 and 21 are chamber counts that correspond to impossible corner activations.
