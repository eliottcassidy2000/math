# The Cuboid and the Forbidden: z=0 Floor Analysis

## The z=0 floor

The forbidden values 7 and 21 both have z = 0 (divisible by 7) in the 2×3×7 cuboid.

Odd points on the z=0 floor mod 42: 7=(1,1,0), 21=(1,0,0), 35=(1,2,0).
7 and 21 are forbidden. 35 is achievable.

## The curvature channel distinguishes them

7: curvature y = 1 (unit curvature). FORBIDDEN.
21: curvature y = 0 (zero curvature). FORBIDDEN.
35: curvature y = 2 (two curvature units). ACHIEVABLE.

The forbidden values have curvature ≤ 1. The achievable has curvature 2.
Insufficient curvature on the z=0 floor → forbidden.

## The channel mismatch

For H = 1 + 2·α₁ + 4·α₂ + ...:
H mod 3 depends on α₁ mod 3 (and higher terms).

At H = 7: H mod 3 = 1. But the required α₁ = 3 has α₁ mod 3 = 0. MISMATCH.
At H = 21: H mod 3 = 0. The required α₁ = 10 has α₁ mod 3 = 1. MISMATCH.
At H = 35: H mod 3 = 2. The required α₁ = 17 has α₁ mod 3 = 2. MATCH.

The forbidden values are where the curvature channel (mod 3) of H
disagrees with the curvature channel (mod 3) of α₁.
The two channels SEE DIFFERENT THINGS. The disagreement = the obstruction.

## Implication

Forbidden-ness requires:
1. z = 0 (divisible by 7 — on the forbidden floor)
2. Curvature mismatch (H mod 3 ≠ α₁ mod 3 for all valid decompositions)
3. The OCF cascade (the topological obstruction that forces the mismatch)

The cuboid gives conditions 1 and 2. The topology gives condition 3.
All three together: only {7, 21} survive as permanently forbidden.
