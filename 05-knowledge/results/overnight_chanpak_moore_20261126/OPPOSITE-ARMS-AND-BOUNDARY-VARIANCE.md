# Opposite-arm forks and integer boundary variance

2026-11-26. Follow-up to the [overnight checkpoint](README.md). The two results preserve the physical edit coordinates: a tournament's two marked spokes must belong to the same tournament, and a Moore cut must account separately for the sites that are forced to carry its changed edges. Neither open problem is resolved.

## 1. Opposite-arm forks: four heights through address 2²⁵

In F1, let H(b,c)=bc+F(b−1)F(c−1)−[b odd][c odd] with the addresses and tournament construction of the checkpoint. Search for central addresses with one unit-spoke reversal in each binary block, with opposite signs. Marked weights are powers of two, a weight ℓ acting on address b by b→1+((b−1) xor ℓ). Both orientations are on the same labelled vertex set (pad at the high end if needed).

The [replay](opposite_arm.js) exhausts **all central pairs through 33,554,432 on both arms**. Up to exchanging the two sides, it yields precisely these marked rows. Each row is also independently checked by subset DP after building the actual orientations:

| b | c | first weight | second weight | first shift | central count | order |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | −2 | 3 | 4 |
| 31 | 4 | 2 | 1 | −2 | 141 | 9 |
| 127 | 4 | 4 | 1 | −2 | 573 | 11 |
| 182 | 29 | 8 | 2 | +2 | 6,681 | 15 |
| 91 | 58 | 4 | 4 | +2 | 6,681 | 15 |

For fixed c and ℓ, let g=F(c−1). A zero-to-one flip on the first arm changes the count by cℓ+γ+g(Z+ε−β), with γ=[c odd] if ℓ=1 and zero otherwise. Here Z is the reverse complement of the log₂ℓ-bit suffix, ε is Z mod 2 (zero when ℓ=1), and β=q+[q odd]. Test both signs t=±2. Necessarily g divides cℓ+γ−t, after which each Z gives at most two admissible q (β or β−1) and a unique decoded suffix. Test both central endpoints, negate the sign at the upper bit endpoint, keep b≥c with b≤N, and directly evaluate the second-arm flips. The loop stops when cℓ>g(⌊N/(2ℓ)⌋+1)+2, a necessary bound from β≤⌊N/(2ℓ)⌋+1 that stays violated as ℓ increases. A second-arm bit weight m≥N is newly leading and increases the scalar by at least bm, so it cannot yield a unit exit. The scan made 1,247,465 first-unit readings after applying b≥c. The central bounds do not require a flipped endpoint to be within the same address box.

**An exact transfer explains the duplicate height 6,681.** If b is even and its marked weight ℓ≥2, replace (b,c;ℓ,m) by (b/2,2c;ℓ/2,2m). The product is unchanged and F(2v+1)=F(v) preserves both unary readings; the parity correction is zero before and after, also on each of the two flipped neighbors. Thus all three totals, and not just the center height, are preserved. The two 6,681 rows are exactly one such transfer apart. Combining this with exchange of arms, one can normalize a transferable fork until its donor address is odd or its donor spoke is the lowest bit. The transfer is an identity between separate physical witnesses, not an assertion that the tournaments are isomorphic.

### All orders when one address is at most four

The only central heights in this lane are 3,141,573. Here is the all-order argument, with h=F(b−1).

* At c=1 or c=2 the weight-one change has magnitude b+[b odd]. Any higher bit either is newly leading or has magnitude too large (at c=2 all higher bits are newly leading). Thus only b=1,2 can contribute, giving just the fork of height 3.
* At c=3, changing to 1 has magnitude 2(b+h)>2, and higher bits add at least 4b. The change 3→4 is b+[b odd]−2h. For b=2u the only unit case is +2 at b=4: h≤u, and h=u−1 requires u even and h=F(u/2−1)≤u/2. Direct substitution gives no opposite exit at (4,3). For odd b=2u+1, a unit exit requires −2 with F(u)=1 and u odd, hence b=2^k−1, k≥2 (F=1 precisely for all-one words or the empty word). On the b arm all non-lowest shifts are multiples of 3. Its lowest shift at those addresses is 4−3·2^(k−1), never +2.
* At c=4 the change 4→2 is −2b, which is unit only at b=1 where there is no opposite exit; higher changes add at least 4b. The change 4→3 is −b+2h−[b odd]. The even case gives only b=4 with −2, and in (4,4) the marked unit shifts cannot have opposite signs. In the odd case the shift is 2(F(u)+[u odd]−1), so +2 forces u odd, F(u)=1, b=2^k−1. On the b side, the lowest-bit change is 4−2^(k−1), never −2. At a set bit i≥1 of b−1 below its leading bit, a −2 change requires β=4·2^i+2^(i−1)+[i=1]−2. The left side β is a power of two (or zero); for i≥3 the right side lies strictly between consecutive powers 4·2^i and 8·2^i. The i=1,2 equations give k=5,7, respectively. At the leading bit β=0, and a newly leading spoke cannot decrease the count. This leaves precisely (31,4) and (127,4).

## 2. Moore separation with frozen exits: at least 1,523 removals and additions

This strengthens the [earlier 1,408 bound](../moore57_quadratic_tensors_20261125/PRIMITIVE-HOLONOMY-AND-EDGE-SEPARATION.md) for two distinct hypothetical full Moore charts agreeing on at least three entire stars. Retain its notation: U is the support of the adjacency difference, u=|U|, r_x is the removed degree, a=min r_x, M is the number of removals, b_x is the number of old neighbors outside U. Its identity for inactive vertices implies b_y≤r_x on removed edges xy, and on added edges too. All active leaves have b_x≥4 (a head and three frozen exits) and 4≤r_x≤53.

**Boundary charge.** Keep the four compulsory old neighbors when using the earlier moment inequality:

    a Σ b_x + 4(2M−au) ≤ Σ r_x b_x ≤ Σ r_x² ≤ 2(53+a)M−53au.

Thus if C is the boundary cut inside the 3024-vertex induced exterior (so Σ b_x=4u+C), then 2(49+a)M≥a(53u+C). This improves the comparison even before using a spectral correction.

**Eight supported sites.** The earlier generic support bound still applies for every a. At a=4 take a minimum-loss row x with four removed neighbors N and four added neighbors P. All eight have exactly four inactive neighbors. For each member of N, its four inactive witnesses of the difference identity meet four distinct members of P, since two vertices have at most one common neighbor. Thus every N–P pair has an inactive common neighbor, has no old edge, and has no active common neighbor. The old union lower bound 415 covers only their active neighborhoods. All four members of N lie outside that union, since they have no old N–N or N–P edges. If e is the number of old P–P edges, their overlap loss is at most 6−e: adjacent pairs share no neighbor. Each of the i isolated vertices in the old graph on P also lies outside the union. So u≥419+e+i≥421, using e+i≥2 on four vertices.

**Integer spectral moments.** Let L be the exterior adjacency, n=3024, w=1_U−u/n, α=46u/n, and z=(7I−L)w. The centered spectrum lies in [−8,7] by extending by zero to the full graph, so

    δ := C − 46u(1−u/n) = wᵀ(7I−L)w ≥ ‖z‖²/15.

The sums of z on U and its complement are δ and −δ, since z is centered. On the 2a removed and added neighbors of a minimum-loss row, each has at least 57−a old neighbors inside U. Thus z≤−K there, where K=50−a−α. If K>0, set m=2a, pay mK² in the norm and use a sum lower bound δ+mK on the remaining u−m inside sites; otherwise take m=0 and keep all inside sites. The outside sites, after reversing their signs, have sum δ. This second use of the spectral estimate keeps both group sums, not just the m large coordinates.

More precisely, every coordinate of n z is congruent to 46u modulo n. For k coordinates of a single congruence class g modulo n, with sum at least a nonnegative S, bound their squared sum by taking k times the piecewise-linear interpolation of consecutive lattice squares at mean S/k. If a closest lattice point to zero lies to the right of that mean, use its square instead. This rule is valid because the interpolant is convex and nondecreasing to the right of a minimizing point. Apply it to the unmarked inside and sign-reversed outside groups, with S=n(δ+mK) and nδ respectively. The bound for the outside group uses the *negative* residue class. If u=n there is no exterior boundary.

Start with the spectral lower bound on C rounded up to parity u (the exterior has degree 53). Substituting this lower bound for C in the displayed sum constraints is valid because the numerical bound on squared sums is monotone in S. Round again to parity u and repeat twelve times or until stationary; any finite number of repetitions remains valid. The [BigInt replay](frozen_integer_support.js) covers every integral support size allowed by the old union bound and the eight-site refinement, and every a=4,…,53. It also drops candidates exceeding the trivial cut capacity. The minimum certified M is **1,523** at the lower-bound inputs a=4,u=421,C≥18,021. The next four minimum-loss degree cases give 2,190, 2,939, 3,729 and 4,540. These are conditional lower bounds, not optimality claims: even the row minimizing the numerical constraints need not be realizable. The same result bounds additions since degree is preserved. It does not say that two partial charts must differ by that many edges.

Replay both packets with Node: 'node opposite_arm.js' and 'node frozen_integer_support.js'.
