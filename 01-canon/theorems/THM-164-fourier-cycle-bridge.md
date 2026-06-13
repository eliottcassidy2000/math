# THM-164: Fourier-Cycle Bridge

**Status:** VERIFIED (n=5 exhaustive, n=6 partial, n=7 regular)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

The Walsh-Hadamard Fourier decomposition of H(T) on tournament space {-1,+1}^m has a direct correspondence with the cycle count decomposition:

### Degree-2k Fourier term encodes (2k+1)-cycle counts

| Fourier degree | Cycle length | Score-determined? | Energy fraction (n=5) |
|:---:|:---:|:---:|:---:|
| 0 | none (constant E[H]) | trivially | 75.95% |
| 2 | 3-cycles (c3) | YES | 22.78% |
| 4 | 5-cycles (c5) | NO | 1.27% |
| 6 | 7-cycles (c7) | NO | 0% at n<=6 |

### At n=5, within the only varying score class (1,2,2,2,3):

H = 9 + 2 * c5_dir, exactly.

All tournaments in this class have c3_dir = 4 (score-determined).
c5_dir = {1, 2, 3} distinguishes H = {11, 13, 15}.

H_4 = 2 * c5_dir - 3 (exact linear relationship).

## The Vitali Structure

The Fourier decomposition reveals a hierarchy of "measurability":

1. **Level 0 (measurable)**: The constant E[H] = n!/2^{n-1}. Universal.
2. **Level 2 (score-measurable)**: H_2 = c_2 * n * (m - 2*Var_s). Determined by score sequence. Captures 3-cycle contribution. This is the "Lebesgue-measurable" part.
3. **Level 4 (cycle-measurable)**: H_4 encodes 5-cycle count beyond score determination. This is the first "non-score" contribution. It's the first "Vitali-like" component — you need more structure than scores to determine it.
4. **Level 2k**: (2k+1)-cycle counts. Each successive level captures finer cycle structure.

## Key Relationships

- c3 is ALWAYS score-determined: c3 depends only on the degree sequence via the formula c3 = C(n,3) - sum_v C(s_v, 2) + ... (Rao's formula). For regular tournaments: c3 = n(n-1)(n+1)/24.

- c5 is NOT score-determined: Within a fixed score class, c5 can vary. The variation is captured by H_4.

- AA^T variance (common-neighbor non-uniformity) correlates with c5_dir: lower AA^T variance -> more 5-cycles -> higher H_4 -> higher H.

## Connection to Overlap Weight

The overlap weight W(C_i, C_j) = |V(C_i) ∩ V(C_j)| operates in the same space:

- At n=5: all directed cycles pairwise conflict (n too small for disjoint cycles)
- H = 1 + 2 * (c3_dir + c5_dir + c7_dir + ...) since alpha_2 = 0
- H_2 captures the c3 contribution to 2*alpha_1
- H_4 captures the c5 contribution to 2*alpha_1

At larger n (n >= 6):
- Disjoint cycle pairs become possible
- H_4 encodes BOTH the c5 contribution AND the disjoint pair structure
- corr(H_4, disj_3_pairs) = 0.548 at n=6 (moderate, not perfect)

## Effective Dimension

H has effective dimension ~1.8 DOF across all tested n (3-6):
- DOF 1: Score variance (Var_s)
- DOF 2: Cycle structure beyond scores (captured by H_4)

This means the tournament space {-1,+1}^m (dimension up to 15 at n=6) collapses to essentially a 2D manifold when viewed through H.

## Verification

- vitali_deep_structure.py, degree4_overlap_bridge.py, H_fourier_deep.py
- All results in 05-knowledge/results/
