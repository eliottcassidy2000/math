# Speedup Ideas from Master Identities
## opus-2026-04-03-S27

### Speedup comparison (operations to enumerate all classes)

| n | Full 2^m × n! | GS-only 2^gs × n! | Speedup | Half-triangle |
|---|---------------|-------------------|---------|---------------|
| 7 | 165M | 2.6M | 64× | 2× additional |
| 8 | 84.6B | 165M | 512× | 2× additional |
| 9 | 97.4T | 23.8B | 4096× | 2× additional |
| 10 | 249P | 3.8T | 65536× | 2× additional |

### IDEA 1: GS-Only Enumeration (for SC classes)
Enumerate only 2^gs_dim grid-symmetric tilings to find ALL SC classes.
Exponential speedup: 2^(m - gs_dim) = 2^{(m - floor((n-1)/2))/2}.
Already used for n=7 (512 tilings → 88 SC classes in 15s)
and n=8 (4096 tilings → 176 SC classes in ~12min).
At n=9: 65536 GS tilings (feasible), vs 268M full (not feasible).

### IDEA 2: Half-Triangle (THM-280 reflection pruning)
Each non-GS tiling pairs with its reflection under R.
Both tilings give isomorphic tournaments (THM-280).
So: enumerate only one from each pair → halves the work.
Combined with GS-only: enumerate source-half configurations.

### IDEA 3: |Aut| is always odd and small
|Aut| divides the odd part of n! = n!/2^v_2(n!).
At n=6: odd part = 45, possible |Aut| = {1,3,5,9,15,45}.
Observed |Aut| = {1,3,5,9} only. Most classes have |Aut|=1.
Consequence: for |Aut|=1 classes, tiling_count = H (no division needed).
73% of classes at n=6 need no automorphism computation at all.

### IDEA 4: Score-based early pruning
Score sequences are fast to compute (O(n) per tournament).
Different score sequences → different iso classes (guaranteed).
Same score sequence → MIGHT be same class (need full canonicalization).
At n=6: 19 distinct score sequences for 56 classes.
Average: ~3 classes per score, much cheaper to disambiguate.

### IDEA 5: ME4 gives exact averages without enumeration
E[H | uniform labeled] = n!/2^(n-1). No enumeration needed.
E[t3 | uniform labeled] = C(n,3)/4. No enumeration needed.
These give exact first moments of tournament invariants.
Second moments (variance) would need more work.

### IDEA 6: Nauty/Traces for canonicalization
Replace brute-force n! canonicalization with nauty/Traces.
For directed graphs, nauty runs in quasi-polynomial time.
Python wrapper: `pynauty`. Would turn each canonicalization
from O(n!) to O(n^c log n) for small c.
At n=8: from 40320 to ~few hundred operations per tournament.
Combined with GS-only: 4096 × ~1000 ≈ 4M ops (seconds, not hours).
