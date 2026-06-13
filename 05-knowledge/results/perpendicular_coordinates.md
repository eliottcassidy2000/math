# Perpendicular Coordinate Analysis Results
## opus-2026-04-03-S27

### The Two Perpendicular Families
- x+y = c (parallel to hypotenuse): controls SC/NS, blue/black
- x-y = c (perpendicular to hypotenuse): controls skip/wiggly structure
- Grid reflection R: SWAPS x+y layers, PRESERVES x-y layers

### THEOREM: Cross-Pair Swap Counting
At distance d from hypotenuse with k reflection pairs:
  Grid-symmetric (GS) choices = 2^k
  x+y-profile-symmetric (XPY) choices = C(2k, k)
  Cross-pair swap extras = C(2k,k) - 2^k

| k pairs | GS = 2^k | XPY = C(2k,k) | Extras |
|---------|----------|---------------|--------|
| 0 | 1 | 1 | 0 |
| 1 | 2 | 2 | 0 |
| 2 | 4 | 6 | 2 |
| 3 | 8 | 20 | 12 |
| 4 | 16 | 70 | 54 |
| 5 | 32 | 252 | 220 |

### Predicted xpy_sym / gs ratio
| n | GS | XPY-sym | Ratio |
|---|-----|---------|-------|
| 3 | 2 | 2 | 1.000 |
| 4 | 4 | 4 | 1.000 |
| 5 | 16 | 16 | 1.000 |
| 6 | 64 | 96 | 1.500 |
| 7 | 512 | 1152 | 2.250 |
| 8 | 4096 | 23040 | 5.625 |
| 9 | 65536 | 921600 | 14.062 |

### The XPY-sym tilings
These are tilings with symmetric layer COUNTS but not necessarily pairwise
bit matching. They represent "weakly symmetric" configurations that maintain
the statistical balance across the hypotenuse without the strict local
matching that grid-symmetry requires.

At n<=5: XPY = GS (no multi-pair layers exist, so count = bit matching).
At n>=6: XPY > GS, with the gap growing rapidly. The cross-pair swaps
represent a new intermediate symmetry class.

### Tiles per distance layer
```
n\d    0    1    2    3    4    5    6
 3     1
 4     1    2
 5     2    2    2
 6     2    4    2    2
 7     3    4    4    2    2
 8     3    6    4    4    2    2
 9     4    6    6    4    4    2    2
```

### Skip counts per skip value (tiles per x-y value)
```
n\s    1    2    3    4    5    6
 3     1
 4     2    1
 5     3    2    1
 6     4    3    2    1
 7     5    4    3    2    1
 8     6    5    4    3    2    1
```
Formula: skip s has n-s-2 tiles.

### H vs Coordinate Profiles (n=5)

Correlation of H with coordinate projections:
  corr(H, Hamming weight) = 0.569
  corr(H, hypotenuse bits) = 0.394
  corr(H, periphery bits) = 0.418
  corr(H, skip-1 bits) = 0.282

H is non-monotone in the (hyp, per) decomposition:
  (hyp=0, per=0) → H=1 (transitive)
  (hyp=0, per=4) → H=15 (maximum)
  (hyp=2, per=0) → H=15 (also maximum!)
  (hyp=2, per=4) → H=9 (anti-transitive)

Maximum H achieved by EITHER all bits on hypotenuse OR all bits off hypotenuse.
The middle (hyp=1) has intermediate H. This suggests H is determined by how
"concentrated" the bits are along one perpendicular axis vs the other.

Skip profile determines H more tightly than distance profile.
Some skip profiles give unique H: (0,0,0)→1, (0,2,1)→15.
Others allow 2-3 values: (2,1,0)→{5,9,13}.

### Vertex-Layer Anatomy (n=5)

Vertex 1 (endpoint): touches dists [0,1,2], skips [1,2,3] — ALL layers
Vertex 3 (middle): touches dists [2,2], skips [1,1] — ONLY periphery
Vertex 5 (endpoint): touches dists [0,1,2], skips [1,2,3] — ALL layers

Endpoints see the full triangle. The middle vertex only sees the periphery.
This is why vertex-star flips at the center are more neutral: they only
affect periphery bits, leaving the hypotenuse backbone intact.

### H-Landscape on the (hyp, per) Plane

Maximum H is achieved at BOTH extremes of the hypotenuse axis:

n=5: max H=15 at (hyp=0,per=4) AND (hyp=2,per=0) — opposite corners
n=6: max H=45 at (hyp=0,per=8), (hyp=2,per=2..6), and many intermediate points

The minimum H=1 (transitive) is always at (0,0) — all tiles in default direction.

The anti-transitive tiling (all bits=1) has:
  n=5: (hyp=2,per=4) → H=9 (NOT maximum)
  n=6: (hyp=2,per=8) → H=17 (NOT maximum)

Key insight: the anti-transitive tournament is NOT the regular tournament!
The regular tournament (max H) is achieved by concentrating flips in a
structured pattern, not by flipping everything.

The H-landscape has a SADDLE POINT structure: high H at both
(hyp-full, per-empty) and (hyp-empty, per-full), lower H at the corners
(all-0 and all-1).

At n=6, max H=45 is spread across many (hyp,per) combinations,
suggesting the regular tournament has multiple tiling representations
with different coordinate decompositions.
