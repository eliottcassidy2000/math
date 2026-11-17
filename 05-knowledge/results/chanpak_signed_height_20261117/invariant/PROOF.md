# Invariant-changing exact +2: two-adic decoder and a second physical branch

The migration and reader are the literal ones in the preceding branching packet: a weak padded word (k,u) with leading bit 0 sends (h,R) by M=(F(u),u+1; 2^k−rev_k(u),2^k), preserves the unique full ending at the root and preserves S=0. A leading-one last word makes the rooted reader strong, and its converse used as the only return side gives R P₂+h H. No arithmetic changes below are interpreted as freely chosen nonphysical matrix entries.

## Independent small slope-five grammar

Take a transitive triple rooted at its last vertex, with (h,R,S)=(1,4,0). On 5h−R=1, the following two physical letters preserve the invariant and branch:

    A=(8,50): h -> 292h−51
    B=(8,84),(4,2), chronological: h -> 5976h−1023.

The second letter reverses the chronological order of the older slope-three two-step block. Finish by (6,50); then (h_f,R_f)=(292h−51,365h−64). The old same one-arc source profiles (187,3736) and (197,3728) give output difference 10h_f−8R_f=2, so the lower count is 1418244h−248641. The first rooted finish is checked by direct DP: (h_f,R_f)=(241,301) and exactly one full path ends at the root, giving lower count 1169603. The construction is strong and the source edge switch is internal, independent of which reader word was used.

This alternate grammar has an exact inverse on all words. Put x=291h−51, initially 240 with valuation 4. A sends x to 292x and raises valuation by two; B sends x to 5976x+7032 and resets valuation to 3 because every x is divisible by 8. An even valuation means a pure A word, with length (v₂(x)−4)/2. An odd valuation means strip (v₂(x)−3)/2 trailing A's and undo B, then repeat. Each undo is exact affine division or rejection, and every successful path decreases h to the seed. The lower target first selects h by exact division modulo 1418244, so this is an actual source decoder on the stated strand. Branch growth is genuine (the two reset valuations distinguish them), not a repeated image of the slope-three language. It remains a sparse fixed-progression bank.

## Fast inverse for the whole changed-token graph

Let level i have q=2^i, a=q+1, with state a h−R=1 and h odd. All levels i≥2 have their own physical seed h=1 given by a transitive (i+1)-vertex tournament rooted at its last vertex. For an allowed chronological edge j→i, put Q=2^j, b=Q+1, and

    2≤i≤j+1, t=1+Q/2−q/4, k=2it, D=2^k, T=(D−1)/a.
    u=T−1,  h' = m h−T,  m=(bD+C)/a,  C=D/2+1+T.

The main token-transfer proof establishes (a,−1)M=(b,−1) at this width and address. One can read C directly: in base q the address has digits q−1 in the even blocks except the last low digit q−2, and zero in each odd block. Hence its k-bit reversal is qT−D/2. Using D=aT+1 also gives

    a²m = [a(b+1/2)+1]D+q.

Every edge grows h by at least a factor of two on h≥1, since m−T≥2 and m h−T≥(m−T)h. The following congruences turn the full graph's inverse into a deterministic valuation calculation, not a search over previous levels.

**Recover the last level and climb bit.** For every nonseed h at level i, v₂(h−1) is 2i for a non-climb, or 2i−1 for a climb. Proof by induction. On a climb j=i−1 the formula simplifies to h'=1+q[(Q+1)h−1]. By induction, the predecessor is either h=1 or h−1 divisible by 2^(2j−1)>2^j; the inequality holds for j≥2. Thus the bracket has valuation j and the result has valuation i+j=2i−1. For j≥i the exact formula gives

    h'−1 == q(h−a)/a² (mod 2^(k−1)),

where odd denominators are inverted. By induction h−1 is divisible by 2^(2j−1)>2^i or is zero. Thus v₂(h−a)=i. Here k−1>2i, so the leading valuation is exactly 2i. Recover level i as ceil(v₂(h−1)/2); its parity tells whether to undo a climb from i−1 or a non-climb.

**Recover non-climb predecessor.** Compute the exact integer Y=(a²h'−a)/q. It satisfies

    Y = h + (D/(2q))*([a(2b+1)+2]h−2a).

The bracket is odd, and k−i−1>2j for j≥i≥2 (already strict at j=i, and the left side then grows exponentially in j). If the predecessor h>1, then e=v₂(Y−1)=v₂(h−1) is 2j or 2j−1, and j=ceil(e/2) is the only candidate. If the predecessor is the seed h=1, then e=k−i−1; recover the only candidate k=e+i+1, then test k=2it and solve 2^(j−1)=t−1+2^(i−2). Both branches use one exact affine inversion and their appropriate seed/nonseed checks.

They cannot both be accepted. For an actual seed branch at j, its e>2j would suggest a larger j' on the fake nonseed branch. But at fixed destination the minimum possible h' for a non-climb, namely m−T, strictly increases with j: explicitly

    m−T = ([b+1/2−q/a]D + [1+q/a])/a,

and both b,D increase with j. Thus the larger level's edge cannot produce the smaller seed output. Conversely any ambiguity would already include such an accepted seed, so this excludes all ambiguity. To avoid constructing a huge rejected candidate D, check its k against the bit length of h' first; on a non-climb m−T>D/2. The two arithmetic candidates and bit sizes are polynomial in the actual input length. Exact validation proceeds recursively until h=1; every accepted inverse can be run forward along the uniquely recovered path. The final reader family in the main proof computes its h from the prescribed output before this inverse is applied.

**An unrestricted binary subgrammar at every level.** Use A=j→j and B=(j→j+1→j) at fixed j. These are distinct productive moves; the full inverse already proves no word collisions. For a direct single-state check put x=(m_A−1)h−T_A. With z=j, q=2^z, the initial valuation is 2z; A raises it by z. Let W be the short climb and N the return. The climb has m_W=2q a, T_W=2q−1. The wide A and N have m == q/a² and T == −1/a modulo 2^(k_A−1). Hence the constant shift of B=N∘W in x is

    (m_N m_W−1)T_A − (m_A−1)(m_N T_W+T_N)
      == −3q³/a⁴ (mod 2^(k_A−1)).

The modulus exponent is larger than 3z. This resets x to valuation 3z because the composite multiplier has valuation 2z+1. Thus if v≥3z strip (v−3z)/z trailing A's, then distinguish A(seed) from a B output by exact equality: min B>A(seed), since N(1)≥2q T_N>a T_A≥A(1). At v=2z accept only the actual seed. This is a uniform two-letter exact branch as well as the independent smaller slope-five branch above.

Run `python replay.py` for the literal source and rooted finish, 255 slope-five words with its own inverse, 35 general edges with their physical migration matrices, 635 uniform binary-loop words and 300 arbitrary paths with exact inverse and nearby ghost controls. Normal and optimized transcripts agree.
