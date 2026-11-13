# A cyclic ternary address tree made out of actual Hamiltonian counts

The [19-adic selector](ADDITIVE_PORTS.md) gave a binary address tree, without an intrinsic ternary orientation. Here three physically interchangeable ports realize both the ancestry and the cyclic orientation in the **actual counts**, using the character of the leading 11-adic digit. It is a literal counterpart of the synthetic [Session 312 chart](../session312_cluster_pythagorean_symmetry_20261113/CLUSTER_CHART.md).

Use the ordinal chain and apex response law from [ADDITIVE_PORTS](ADDITIVE_PORTS.md). The following three local choices all share $(h,q,s,s^*)=(33,720,0,127)$, as computed by subset DP. Outgoing masks and cuts use zero-based vertices.

| Label | Outgoing rows | Apex out-cut mask | $(e,e^*,b)$ |
|---|---|---|---|
| A | (52,33,2,7,14,28,63) | 40 | (24,509,109) |
| B | (56,17,3,6,12,30,63) | 14 | (12,305,71) |
| C | (52,33,2,7,14,28,63) | 17 | (26,515,111) |

All choices are independent in the total count, because the start responses and both unrooted responses agree. Freeze the first site. If the latest differing site of two count words is $i\ge2$, their difference has valuation $n-i$. To prove it, for each differing site use the overlapping-pair formula with $j<i$. The lowest valuation contribution at site $i$ is the unique $j=1$ term,
$$
\Delta e_i^*\,127\,720^{i-2}\,33^{n-i}.
$$
Every other term at this site has higher valuation: $j>1$ adds factors of 33, the local and terminal terms have factor $33^{n-1}$, and the adjacent contribution vanishes since $s=0$. Earlier sites likewise have higher valuation. The three $e^*$ residues are $3,8,9$ modulo 11, so the leading coefficient cannot vanish. The factors $720$ and $33/11$ reduce to squares modulo 11; 127 reduces to the nonsquare 6. Therefore, orienting a pair by the quadratic character of the leading unit of the signed count difference gives a cyclic triangle on $A,B,C$: the differences $8-3,9-8,3-9$ are squares, and the common nonsquare factor reverses all three arcs. At each first differing digit it is exactly the same directed triangle. This gives a ternary cyclic tower of height $n-1$ on $3^{n-1}$ distinct integer counts, before adding any spacers.

## Simultaneous real order and adic ancestry

For a positive ordered chart, take $n=2d$ and activate only the even sites $2,4,\ldots,2d$, freezing all the rest. Order digits by $B<A<C$. All three local coordinates $e,e^*,b$ increase in that order. Put $\rho=720/33$ and $G_i=\sum_{t=0}^{i-2}\rho^t$. Ignoring a common positive factor $33^{n-2}$, the site-$i$ contribution is $127 e_i^*G_i+33 b_i$, plus $33 e_i$ at the terminal site. Its smallest successive gap is at least $127\cdot6G_i$. An earlier active site $j$ is not terminal and its full span is at most $(127\cdot210+33\cdot40)G_j=27990G_j$. Since $G_{j+2}>\rho^2G_j$ and $\rho^2>475$, the sum of all earlier spans is less than $27990G_i/(\rho^2-1)<127\cdot6G_i$. Thus ordinary real count order is the right-to-left lex order on the active sites.

In that order the positive Plücker coordinates are literally $H_v-H_u$ (columns $(1,H_u),(1,H_v)$). Their valuations are twice the length of the common address prefix, and their leading units still orient the cyclic split. Prefix classes are consecutive in the real order, so the minimum formula for four ordered points gives exact tropical exchange with no excess valuation. This lets us recover the same tree in three ways: by reading the chosen port sequence, by looking at actual count differences, or by forming minors of the actual count coordinates.

The outputs can be taken strong without changing these counts. Vertex 6 of each local block is a source within that block, and is outside all three out-cuts. Delete vertex 6 of the first block, a source of the whole tournament. The remaining first block is strong, as is the six-vertex interior of each other block (a reachability check on the displayed rows). The apex reaches the first interior, all later vertices are reachable through the ordinal arcs, and every vertex can reach the apex. Deleting a global source does not change the count. This is an orientation on a family **of counts**, not a claim that a single tournament realizing a count has the tower as its automorphism group.

The result gives a literal three-way physical address compiler, including the signed digit, but still not a dense interval: at a fixed length only three leading choices at each stage are available. Nothing here covers all odd targets.

## A neutral two-vertex pad makes every site usable

Append one sink to each block, placing it in the complement of the apex out-cut (so it beats the apex). The response transformation is
$$
(h,q,s,s^*,e,e^*,b)\longmapsto(h,2q,s,2s^*,h,e^*+q,b+s^*).
$$
For two-slot covers, split according to which slot contains the sink. It contributes at the end of that slot, may be a slot by itself, and is an allowed end but never an allowed start. This gives the displayed transformations, including the cases of an empty other slot. Two successive sinks therefore give $(h,q,s,s^*)=(33,2880,0,508)$, with $e=33$ constant, $e^*$ increased by 2160 and $b$ increased by 381 for all three choices. The rotation signatures are unchanged: 2880 reduces to the square 9 modulo 11 and 508 to the nonsquare 2. The same argument now gives simultaneous real and 11-adic prefix order with **every** site $2,\ldots,n$ active. Indeed $G_{i+1}>(2880/33)G_i$ and
$$
\frac{508\cdot210+33\cdot40}{2880/33-1}<508\cdot6 .
$$
So the full earlier span is less than a single current digit gap. The 11-adic valuation is the common prefix length itself, no longer doubled. One sink already suffices for real order, but reverses the cyclic convention on alternating levels; two make the same convention uniform.

The length-bypass calculation also applies here. For a uniform cut and $n\ge k+2$, the count modulo $11^k$ is
$$
s^*e^*q^{n-2}(1-33/q)^{-2}.
$$
The unit $q=2880$ has order five modulo 11 and $q^5\equiv23\bmod121$, so it generates the square subgroup at every precision by binomial lifting. Choices B and A have $s^*e^*$ respectively nonsquare and square. Thus two fixed uniform choices, varying only the length, realize every unit modulo each $11^k$; ordinal appendages of the 33-path block supply nonunits. On the other hand, every depth-$n$ count is of order $2880^n$, by the positive endpoint term and the same geometric upper bound as in [LENGTH_BYPASS](LENGTH_BYPASS.md). The combined family with appendages has exactly order $X^{\log_{2880}3}$ distinct counts up to $X$: one depth gives the lower bound, and summing the at-most-$3^n$ words over depths and appendages gives the matching upper bound. In the upper bound a fixed choice at the first site or finitely many uniform cuts changes only constants. Once again complete modular coverage leaves a thin absolute spectrum.

## A literal Zeckendorf rank selector, without a conflict-path assumption

Label the real-order palette $B<A<C$ by $0,1,2$. With $d=n-1$ active sites, restrict to words in which no two consecutive digits are nonzero. There are exactly $J_{d+2}$ admissible words and hence exactly that many distinct physical counts. Because ordinary count order is address order, the [typed Jacobsthal rank rule](FIBONACCI_CONFLICT.md) is now a rank-to-tournament algorithm: choose the $0,1,2$ branch with subinterval lengths $J_{d+1},J_d,J_d$, subtract its rank offset, and recurse. This constructs the tournament of any specified increasing rank without scanning the family. The 11-adic difference chart still recovers the restricted trie. The rank is not the count; this avoids the impossible full conflict path while retaining the same arithmetic recursion as a genuine selector of ordered physical objects.

[Standalone replay](ternary_replay.py) checks the padded port profiles by subset DP, the lex order and leading characters through four active sites, and the length formula through modulus $11^3$.
