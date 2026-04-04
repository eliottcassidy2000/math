# Edge Overlap and Structural Analysis (n=3..7)
## opus-2026-04-03-S27

### Pure-blue counts (VERIFIED, explorer conventions)
| n | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| Pure-blue | 2 | 1 | 3 | 2 | **4** |
| Mixed | 0 | 1 | 5 | 10 | **84** |
| Pure-black | 0 | 2 | 4 | 44 | **368** |

### Transpose == True-complement: PROVED through n=7
The grid reflection on tiling space induces T→T^op on iso classes.
Both operations produce identical class pairings at every n tested.
SC = Transpose-self at every n.

### Edge overlap (wiggly d=1 vs complement d=m)
| n | Pairs | Wiggly | Comp | Both | Wig-only | Comp-only | Neither |
|---|-------|--------|------|------|----------|-----------|---------|
| 3 | 1 | 1 | 1 | 1 | 0 | 0 | 0 |
| 4 | 6 | 5 | 3 | 3 | 2 | 0 | 1 |
| 5 | 66 | 30 | 24 | 15 | 15 | 9 | 27 |
| 6 | 1540 | 290 | 362 | 99 | 191 | 263 | 987 |
| 7 | 103740 | 4086 | 12312 | 680 | 3406 | 11632 | 88022 |

### Density trends
| n | Wig/Total | Comp/Total | Both/Total | Neither/Total |
|---|-----------|------------|------------|---------------|
| 3 | 1.000 | 1.000 | 1.000 | 0.000 |
| 4 | 0.833 | 0.500 | 0.500 | 0.167 |
| 5 | 0.455 | 0.364 | 0.227 | 0.409 |
| 6 | 0.188 | 0.235 | 0.064 | 0.641 |
| 7 | 0.039 | 0.119 | 0.007 | **0.848** |

Key: the metagraph becomes EXTREMELY sparse. At n=7, 85% of class pairs
have NO direct connection by either operation.

Complement edges overtake wiggly at n≥6 (362 vs 290 at n=6, 12312 vs 4086 at n=7).
The complement operation reaches ~3x more class pairs than wiggly at n=7.

### Mixed complement edges at n=7
50 mixed edges out of 12312 (0.4%). All have pattern 1-2 blue + 2-4 black lines.
One edge has 2blue+2black (class 82↔277). All others are 1blue + kblack.
