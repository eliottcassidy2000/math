# Symmetry consumption at a nonuniform aperture

Let $H(B)$ count Hamilton paths, $F_j(B)$ ordered covers by $j$ nonempty paths, and $Q_B(2)=2H(B)+F_2(B)$ ordered covers by two possibly empty slots.

## 1. Marked-stabilizer contact budget

Suppose a tournament consists of disjoint nonempty blocks $B_i$ with uniform arcs between each pair of blocks, and $k$ additional apex vertices with arbitrary incidence. For a prime $p$, put $g_i(s)=\min_{|S|\le s} v_p(|\mathrm{Aut}(B_i)_{(S)}|)$, where the subscript denotes pointwise stabilizer. Then
$$
v_p(H(T))\ge \min_{\sum s_i\le 2k}\sum_i g_i(s_i).
$$
Partition the paths by their full block/apex position word and by the names of block vertices immediately adjacent to apexes. If a vertex lies between two apexes, mark it once. In each part, the product of the pointwise stabilizers of these names acts on paths by relabelling within their existing positions. All used apex contacts are fixed, and between-block arcs are uniform. The action is free because the path uses every vertex. Its orbit size divides the summand, and each part has at most $2k$ marks in total. A chosen subgroup can replace the full automorphism group throughout.

In the ternary rotation tree of depth $d$, the iterated rotation group has 3-exponent $E_d=(3^d-1)/2$. Fixing specified leaves kills one rotation factor per internal ancestor in their union: this follows recursively since a rotation at an ancestor cannot fix a specified leaf below it. Consequently the least exponent after $s$ marks is
$$
g_d(s)=E_d-\sum_{j=0}^{d-1}\min(3^j,s).
$$
The sum is an upper bound on ancestors hit, attained simultaneously by distributing marks as evenly as possible at each branch. One mark may cost $d$ factors; multiple marks share some ancestors.

## 2. A single apex consumes four factors of three, sharply

Take $B=C_3[C_3,C_3,C_3]$ on nine vertices. Its rotation subgroup has order $3^4$ and is vertex-transitive. The replay gives
$$
h=H(B)=3159,qquad q=Q_B(2)=48924,qquad S=(q-h)/9=5085.
$$
Their 3-valuations are respectively 5, 4, 2. The quantity $S$ counts ordered two-slot covers whose second slot starts at a specified vertex, with the first optional. Sum over the starts to get $q-h$, then use transitivity. The dual optional-second-slot end response is also $S$, and each specified Hamilton start/end count is $h/9$.

For $r\ge2$, take the ordinal join of $r$ copies of $B$ and add apex $x$. It beats just a chosen vertex $a$ in the first block, loses to every middle vertex, and beats all of the last block except a chosen vertex $b$. This tournament is strong: $x$ enters the first block; that block and the intervening blocks reach the last, whose vertices can all reach $b$ and hence $x$. Its exact count is
$$
H_r=S h^{r-1}+h^{r-2}(h-h/9)(h-[r=2]h/9)+S(S+h)q^{r-2}.
$$
Proof by the path slots on the sides of $x$. Endpoint and same-block seams together give $2Sh^{r-1}$: the first block can start at $a$ with optional pre-slot, and the last can end at $b$ with optional post-slot; marked vertices already used cannot also appear in the other slot. A forward seam must join adjacent blocks, and the only possible successor block is the last. It contributes the middle term above. A backward seam must start the post-slot at $a$. Ending the pre-slot at $b$ gives $S^2q^{r-2}$; ending in strict middle block $i$ gives $S(q-h)q^{i-2}h^{r-i}$ since that block's pre-slot is nonempty. Summing telescopes by
$$
(q-h)\sum_{i=2}^{r-1}q^{i-2}h^{r-i}=h(q^{r-2}-h^{r-2}).
$$
The last term of $H_r$ has exact valuation $4r-4$, since $S+h$ has valuation 2. The first has valuation $5r-3$, and the middle $5r-2$ if $r>2$, or 6 at $r=2$. Thus $v_3(H_r)=4r-4$ exactly, attaining the marked-stabilizer floor. This demonstrates genuine arithmetic loss of four factors, not just a reduction of the group order. A shallow $C_3$ instead has $Q(2)-H=9$, so that block's analogous marked response loses no factor.

## 3. Strong homogeneous amplification floor

If $Q$ is strong on at least three vertices, $v$ a vertex of $Q$, and $U$ a nonempty tournament of order $n$, then
$$
H(Q[v\gets U])\ge H(Q)H(U)+F_2(U)\ge(n+2)H(U).
$$
Quotient Hamilton paths give the one-run term. A directed triangle through $v$ exists: if no out-neighbor of $v$ beats an in-neighbor, the in-neighbors dominate both $v$ and its out-neighbors, contradicting strongness. Start the walk $v,a,b,v$. Insert every unused quotient vertex: prepend if it beats the first, append if the last beats it, otherwise a consecutive pair gives the required two-arc insertion by a sign transition. This produces a two-run itinerary for $U$, filled by any nonempty ordered two-cover. Also $F_2(U)\ge(n-1)H(U)$ by splitting each Hamilton path at any edge, injectively. A directed triangle gives $H(Q)\ge3$ by extending its three Hamilton paths by insertions.

More generally, $k$ directed excursions from $v$ back to $v$ with pairwise disjoint interiors concatenate to a $k+1$-run itinerary after inserting all unused vertices. This contributes at least $\binom{n-1}{k}H(U)$ by splitting at $k$ edges; disjoint run counts yield disjoint sets of paths. In particular no fixed-factor strong multiplication (including 3 or 9) can indefinitely keep the entire old tournament as one homogeneous block. The quotient may change with the order: the bound is uniform.

## 4. Arithmetic aperture capacity even without automorphisms

If the interblock order is ordinal and there are $k$ apexes, define
$$
d_i(k)=\gcd_{1\le j\le\min(k+1,|B_i|)} F_j(B_i).
$$
The product of $d_i(k)$ over untouched blocks divides each marked-itinerary summand of $H(T)$. Every apex-free segment travels monotonically through the blocks, so an untouched block has at most $k+1$ maximal runs. For a fixed itinerary and the names of used contacts, it fills its $j$ runs independently by any ordered nonempty cover, giving factor $F_j(B_i)$. Thus for any prime $p$, the sum of $v_p(d_i(k))$ minus its largest $\min(2k,r)$ entries is a lower bound for $v_p(H(T))$. For any marked block one can instead use its actual pointed count or its stabilizer. Do not multiply two divisibility claims on a single factor without separate justification.

For a sharp rigid example take $B$ with outgoing masks $(14,28,24,16,1)$, using bit $j$ for an arc to $j$. It is strong via $(0,1,2,3,4,0)$. Its degree classes have size at most two and no directed pair can be swapped, hence it is rigid. Here $h=9$ and $q=84$ and $d_i(1)=3$. At vertex 1, the start count is 4, the end count is 1, and the optional-slot start/end responses are respectively $S=25$ and $E=10$. Put $a=b=1$ in the respective first/last blocks of the same one-apex construction. The seam summation above, allowing different $S,E$, gives
$$
H_r=E h^{r-1}+h^{r-2}(h-s_b)(h-[r=2]e_a)+S(E+h)q^{r-2}
   =10\cdot9^{r-1}+5(9-[r=2])9^{r-2}+475\cdot84^{r-2}.
$$
The first and second terms have higher valuation than $r-2$ when $r\ge3$; the last has exactly $r-2$. At $r=2$ the count is 605, also of valuation zero. Thus $v_3(H_r)=r-2$ exactly. The replay independently evaluates the literal tournaments for $r=2,3$. In any one-apex ordinal grammar with these rigid blocks, $r\ge3$ forces a factor 3 and the insertion lower bound $H\ge9^r$, excluding both prime and $7p$ targets for prime $p$ within that grammar.

These are limitations and sharp exact families for particular physical interfaces, not an unrestricted prime exclusion or proof of Chan–Pak.
