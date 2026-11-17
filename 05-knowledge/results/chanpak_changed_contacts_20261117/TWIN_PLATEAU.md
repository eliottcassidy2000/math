# The twin-seam lower witness cannot climb monotonically to the upper count

Let T_m be the literal [slope-129 negative witness](TWIN_SEAM.md) with exterior masks (286,276,8,402,420,463,31,327,68), marked q=8, and B=TT_m attached as D−q→B→q, m≥1. Its path count is N_m=903·2^m−7. The paired construction gives the count N_m+14 on the same number of vertices.

**Proved:** there is no nondecreasing sequence of single arc reversals from this T_m to *any* tournament on those vertices with count N_m+14. This allows neutral moves, changes to all contacts, reversals inside B, and arbitrary detours. It is a property of these starting witnesses, not a theorem about all tournaments with count N_m. The only plateaus reachable without overshooting are transitive relabelings, and at m=3, one genuine extra exterior move.

## Internal and exterior moves

An internal reversal in a transitive block either swaps adjacent vertices, which only relabels the block, or creates a directed triangle. For the latter, the odd-cycle packing expansion gives ΔP₂≥2ΔH>0 from the original transitive block: every nonempty pack contributes 4^c·2^s to ΔP₂ and 2^c to ΔH, with c≥1. Because the exterior reads 903P₂−7H even for an arbitrary modified B, the increase is at least (1806−7)ΔH≥3598, too large.

For the 36 exterior moves use the exact [closed-band calculation](TWIN_SEAM.md#a-physical-four-arc-path-but-not-a-monotone-one) at r=P₂/H=2^m. The only move with count in [N_m,N_m+14] is {1,6} when m=3, and it leaves the count at N_m.

## All-order certificate for moves that change a contact

Fix a vertex x of B and reverse its arc to an exterior vertex i. Write a,b for the numbers of B-vertices before and after x in transitive order. An entire run in B can exit only to q or, using x, to i. Thus there are at most three B-runs. Within any run, the pieces before and after x follow the transitive order, so each of the two transitive subblocks appears at most once in that run. A collapsed path therefore reads at most three ordered nonempty slots of each subblock. For a positive subblock length s, the respective allocation counts are

    F₁(s)=1,  F₂(s)=2^s−2,  F₃(s)=3^s−3·2^s+3.

On the nonempty/nonempty stratum the count is a bilinear form in the two triples (F₁,F₂,F₃). The one-empty strata use a single triple and the both-empty stratum a constant. Subtracting the old count still has this form (2^s=F₂(s)+2F₁(s) for s>0). The evaluations at s=1,2,3 uniquely determine the coefficients, since the basis matrix is triangular with diagonal 1,2,6. Thus the literal subset-DP grid and interpolation in the replay are an exact algebraic certificate, not a test at scattered heights.

For every i≠q all difference coefficients are nonnegative in this basis, with constant coefficient at least 168 on each of the four strata. For i=q the certificate simplifies to

| positions | exact change |
|---|---:|
| a,b>0 | 199·2^a−903·2^b |
| a>0,b=0 | −903·2^a |
| a=0,b>0 | 199·2^b−984 |
| a=b=0 | −984 |

None equals zero or lies in [1,14]. In the first row, a−b≤2 makes the difference negative, while a−b≥3 gives a positive difference of at least 689·2^b. In the third row b≤2 is negative and b≥3 gives at least 608. The replay also verifies literal controls outside the interpolation grid.

## Closing the one nontrivial plateau

At m=3 apply the neutral exterior reversal {1,6}. This is now a literal order-twelve tournament with count 7217. Direct subset DP checks all 66 single-arc neighbors. None increases by 1 through 14. The only equal-count moves are undoing {1,6}, or swapping adjacent vertices of the still-transitive block. Up to block relabeling, a nondecreasing path staying below or at the target could therefore never leave the starting plateau. This proves the claimed barrier, including both prime-twin heights m=3 and m=5.

## The exact price of a shortest fixed-endpoint route

Use the [published marked relabeling and four edits](TWIN_SEAM.md#a-physical-four-arc-path-but-not-a-monotone-one), numbered {0,4}, {0,5}, {4,7}, {5,7}. For arbitrary B write P=P₂(B), h=H(B). A shortest route must use each of these four arcs once. If the excursion of a route is its largest distance outside the closed endpoint band [903P−7h,903P+7h], the minimum is

    94P+12h     (order 1,3,4,2, below the lower endpoint).

If a route must stay at or above the lower endpoint, the least overshoot is

    191P−8h     (order 2,3,1,4).

Indeed the four possible penultimate states have affine coefficients (W,K) given by (703,7), (809,−19), (1132,17), and (1094,−1). The first two are below the lower endpoint, the other two above the upper one. Their distances and P≥2h give the two lower bounds. The two displayed orders attain them by checking their intermediate reads (the replay lists the complete four-edit cube). This quantifies the difference between a small endpoint gap and the cost of a physical route. It does not measure arbitrary longer routes or all count-equivalent witnesses.

Replay: `node replay_twin_plateau.js`.
