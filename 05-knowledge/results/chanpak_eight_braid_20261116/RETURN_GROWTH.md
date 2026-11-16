# Uniform lower growth in the repeatable return language

The exact mod-32 controller uses 168 eight-letter return blocks. The lower-growth inequality for **all** constant-rank literal stretches is open, but this particular infinite rank-four language admits a direct uniform lower estimate without a cone LP:

| Ordered blocks in the isolated induced suffix | Number of choices | Minimum actual H |
|---:|---:|---:|
| 1 | 168 | 55 |
| 2 | 28,224 | 27,143 |
| 3 | 4,741,632 | 13,460,883 |

This is an exact finite matrix calculation for each row; choices can repeat. The isolated suffix starts from the empty profile, so any backward span directed beyond its first vertex is correctly omitted. This is not the H obtained by extending the seed.

For literal tournaments A and B occupying two disjoint vertex intervals of a larger ribbon, H of the union tournament is at least H(A)H(B): use the independent old odd-cycle packs on the two supports, by [Grinberg–Stanley](https://arxiv.org/abs/2307.05569). No forward orientation across the boundary is required for the inequality. Partition an arbitrarily long return schedule into groups of three blocks. With m return blocks appended to any physical prefix T this gives the all-length estimate
```
H(final) ≥ H(T) · 13,460,883^floor(m/3).
```
In particular the asymptotic per-vertex rate of this lower bound is 13,460,883^(1/24) ≈ 1.98173. The two-block online controller has the simultaneous uniform guarantee H(after request)≥27,143 H(before request). The proof works for every choice of blocks in the bank, not merely a favorable periodic path, and it is compatible with staying at exact rank four.

This strengthens the ribbon growth picture in one highly restricted language. The same lower bound cannot just be transferred to an arbitrary constant-rank stretch; the known alias and slow-pattern issues still apply. This rapid multiplicative growth also makes the online controller very different from a selector of narrowly spaced absolute counts.

Replay: after generating `return_automaton.json` with `lattice.py`, run `python check_return_growth.py`. It builds one literal transfer matrix per block and compares complete one-, two-, and three-block induced counts; integer-overflow bounds are checked before each matrix approach.
