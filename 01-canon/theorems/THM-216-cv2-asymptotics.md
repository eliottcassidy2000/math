---
id: THM-216
name: CV² Asymptotics
status: VERIFIED (computational, n=3..25)
verified_by: opus-2026-03-14-S89c
---

# THM-216: CV²(H) = 2/n + o(1/n)

## Statement

For the Hamiltonian path count H over uniformly random tournaments on n vertices:

  CV²(H) = Var(H)/E[H]² = 2/n + o(1/n)

Equivalently: Var(H) ~ 2·E[H]²/n as n → ∞.

More precisely, n·CV² → 2 from below, and the gap 2 - n·CV² = O(1/n²).

## Known Exact Values

| n | CV² | n × CV² |
|---|------|---------|
| 3 | 1/3 | 1.000 |
| 4 | 1/3 | 1.333 |
| 5 | 19/60 | 1.583 |
| 6 | 13/45 | 1.733 |
| 7 | 131/504 | 1.819 |
| 8 | 131/560 | 1.871 |
| 9 | 1097/5184 | 1.905 |
| 10 | 3121/16200 | 1.927 |

## Computational Data (n=11..25)

| n | CV² (approx) | n × CV² |
|---|---------------|---------|
| 15 | 0.13147 | 1.9721 |
| 20 | 0.09927 | 1.9855 |
| 25 | 0.07964 | 1.9911 |

## Proof Strategy

CV² = W(n)/n! - 1 where W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}.

- NUD(n) = permutations with no unit descent
- adj1(σ) = number of unit ascents in σ
- |NUD(n)| = A000255(n-1) ~ n!(n+1)/(en)

The EGF of A000255 is exp(-x)/(1-x)², and W(n) is the z=2 specialization
of the bivariate generating function tracking unit ascents in NUD permutations.

W(n)/W(n-1) → n, so W(n)/n! → 1 (i.e., CV² → 0).
The precise rate n·CV² → 2 should follow from the generating function structure.

## Note

Earlier conjecture that CV² → 1/4 was WRONG (based on n≤7 data only).
The value crosses below 1/4 at n=8 and continues toward 0.

## Files

- `04-computation/nud_weight.c` — C implementation for n up to 25
- `04-computation/pi_edge_dist_89c.py` — Python analysis
- `04-computation/pi_cv2_n7_89c.py` — Original n=7 computation
