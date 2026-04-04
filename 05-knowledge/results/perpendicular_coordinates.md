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
