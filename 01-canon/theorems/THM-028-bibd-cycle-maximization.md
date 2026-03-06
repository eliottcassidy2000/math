# THM-027: BIBD Arrangement Maximizes Directed Cycles, Not Disjoint Pairs

**Status:** PROVED (exhaustive n=7, algebraic formula)
**Date:** 2026-03-06
**Author:** kind-pasteur-2026-03-06-S18h

## Statement

Among all regular tournaments on p=7 vertices:

1. **All have the same 3-cycle count**: c_3 = p(p-1)(p+1)/24 = 14 vertex sets, each supporting exactly 1 directed 3-cycle. (Moon's formula; depends only on score sequence.)

2. **All have the same 5-cycle and 7-cycle vertex sets**: 21 five-vertex subsets support directed Hamiltonian cycles, and 1 seven-vertex set (the full tournament). But the NUMBER OF DIRECTED CYCLES per set varies.

3. **The BIBD arrangement (Paley T_7, uniform lambda=2) MINIMIZES alpha_2** (disjoint 3-cycle pairs): alpha_2 = 7, vs 10 or 14 for other arrangements. This is because D = C(b,2) - p*C(r,2) + sum C(lambda_e, 2) and uniform lambda minimizes the convex sum.

4. **The BIBD arrangement MAXIMIZES total directed odd cycles**: 80 = 14 + 42 + 24, vs 65 = 14 + 36 + 15 or 59 = 14 + 28 + 17 for other classes.

5. **H-maximization is driven by total directed cycles**: H = 1 + 2*alpha_1 + 4*alpha_2 where alpha_1 = total directed cycles. The 2*alpha_1 term dominates 4*alpha_2.

## Proof

### Part 1: BIBD formula for alpha_2

For a tournament T with c_3 triples forming a 2-(p,3,lambda) BIBD:
- b = c_3, r = per-vertex count, lambda = per-pair count
- Number of blocks disjoint from a given block B:
  d(B) = b - 1 - (3r - 3*lambda) = b - 1 - 3r + 3*lambda
- Total disjoint pairs: alpha_2 = b * d(B) / 2

For general (non-BIBD) arrangement with same b, same d(v)=r (regular):
D = C(b,2) - p*C(r,2) + sum_e C(lambda_e, 2)

Since C(x,2) is convex and sum lambda_e = 3b is fixed:
sum C(lambda_e, 2) >= C(p,2) * C(lambda, 2)  with equality iff all lambda_e equal.

So the BIBD MINIMIZES sum C(lambda_e, 2), hence MINIMIZES D = alpha_2.

### Part 2: Exhaustive verification at n=7

All 2640 regular n=7 tournaments fall into exactly 3 classes:

| lambda distribution | alpha_2 | dc3 | dc5 | dc7 | total_dir | H |
|---|---|---|---|---|---|---|
| [2]*21 (BIBD) | 7 | 14 | 42 | 24 | 80 | 189 |
| [1]*7+[2]*7+[3]*7 | 10 | 14 | 36 | 15 | 65 | 171 |
| mixed | 14 | 14 | 28 | 17 | 59 | 175 |

H = 1 + 2*total_dir + 4*alpha_2 verified for all 2640.

### Part 3: Why BIBD maximizes directed 5-cycles

Each 5-vertex subset of T_7 supports 0 or more directed Hamiltonian cycles.
For the BIBD group: EVERY 5-vertex subset has exactly 2 directed Ham cycles (each is the regular tournament T_5, which has exactly 2 directed Ham cycles).
For non-BIBD groups: some 5-vertex subsets have fewer than 2 directed Ham cycles.

The BIBD's uniform lambda forces each 5-vertex subtournament to be as "regular" as possible, maximizing its directed cycle count.

## Implications

1. **H-maximization is NOT about maximizing disjoint cycle pairs** (alpha_2). The BIBD actually minimizes alpha_2.

2. **H-maximization IS about maximizing total directed odd cycles** (alpha_1). The BIBD's uniform distribution of 3-cycles creates the most regular subtournaments, maximizing directed 5-cycles and 7-cycles.

3. **The linear term dominates**: In H = 1 + 2*alpha_1 + 4*alpha_2, the 2*alpha_1 contribution (160 for Paley) dwarfs 4*alpha_2 (28 for Paley).

4. **Structural explanation for Paley maximizer**: The quadratic residue structure creates a 2-(p,3,(p+1)/4) BIBD of 3-cycles, which forces maximum regularity of subtournaments, maximizing directed cycle counts at all odd lengths.

## Verification Scripts

- `04-computation/bibd_h_maximization.py` — exhaustive alpha_2 vs H at n=7
- `04-computation/bibd_disjoint_proof.py` — BIBD minimizes alpha_2 among regular
- `04-computation/paley_directed_cycles.py` — directed cycle counting by DP
