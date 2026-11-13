# Selective shared-apex fusion and the even clock obstruction

## Exact formula

Let nonempty tournaments $B_1,\dots,B_r$ be in ordinal order (all interblock arcs point right). Add an apex $v$ with **arbitrary** arcs into each block. For block $i$, let $P_i$ count its Hamilton paths, $S_i$ count those whose first vertex is beaten by $v$, and $E_i$ those whose last vertex beats $v$. Let $F_i$ count ordered covers by two nonempty directed paths. Of those covers, $A_i$ counts those with first path ending at a vertex beating $v$, $J_i$ those with second path starting at a vertex beaten by $v$, and $D_i$ those satisfying both. Write lowercase for division of each coordinate by $P_i$ and $q_i=2+f_i$. Then the exact normalized Hamilton path count is

$$
\frac{H}{\prod_i P_i}=s_1+e_r+\sum_i d_i+\sum_{i=1}^{r-1}e_i s_{i+1}
+\sum_{i<j}(s_i+j_i)(e_j+a_j)\prod_{i<k<j}q_k.
$$

Here $i$ in the last sum is the block immediately after the apex, and $j$ the block immediately before it. If $i=j$, that block has a constrained ordered 2-cover. If $i>j$, there is no overlap and nonemptiness forces $i=j+1$. If $i<j$, block $i$ can lie wholly to the right of the apex or be split, giving $S_i+J_i$; block $j$ dually gives $E_j+A_j$. Every intermediate block lies wholly on either side or splits arbitrarily, giving $2P_k+F_k$. All remaining blocks are forced to their side of the apex. Paths beginning or ending at the apex give the first two terms. Each choice is independent since all cross-block edges point right. This proves the formula by a reversible decomposition without assumptions about strengths or uniformity.

With $r\ge2$ every Hamilton cycle must pass the apex; cutting there gives a path in the old spine whose first and last blocks are distinct. Thus the cycle count (unrooted up to rotation) is $(\prod P_i)s_1 e_r$. For $r=1$ instead use the number of block paths satisfying both endpoint conditions, not the product of their marginals.

## Exact coarse graining for uniform incidence

When each block is entirely beaten by the apex (letter 1) or beats it (0), the formula simplifies to

$$
\frac{H}{\prod_iP_i}=1+\sum_{i<j:\,w_i=1,w_j=0}
\left(\mathbf1_{j=i+1}+(q_i-1)(q_j-1)\prod_{i<k<j}q_k\right).
$$

Indeed (d_i=0), a pair contributes only for a 1-to-0 inversion, with endpoint factors from the optional split at each end, and $s_1+e_r+\#01=1+\#10$. With singleton blocks $q_i=2$, recovering the S302/S304 weighted inversion rule.

An ordinal merge of blocks of the same letter preserves the whole tournament, multiplying their $P$ and $q$. For two labelled path slots each block may use either slot alone or both. Choices concatenate slotwise in ordinal order, proving $2P_{A\to B}+F_{A\to B}=(2P_A+F_A)(2P_B+F_B)$; the $2P_{A\to B}$ term counts the two ways to leave one slot empty throughout.

For an arbitrary nonempty block, $F$ is even (swap the two paths freely) and $P$ odd by Rédei. Thus every actual gap clock $q=(2P+F)/P$ is an even rational with odd denominator, with strictly positive 2-adic valuation. Any normalized monomial built by inserting nonempty clock blocks has that valuation strictly increased. In particular one cannot realize a normalized pure ×3 (or any odd unit) simply by adding gap blocks, even if their internals and order change at each step. This claim concerns a monomial transfer with unchanged coefficients, not the full sum or exact cancellation between multiple changing terms.

In the simple strong pattern $C_3[A,B,1]$, corresponding to the word 10 with two blocks, the normalized count is $2+(q_A-1)(q_B-1)$. Thus it is exactly 3 only if both blocks are singletons: any block of order at least two has a two-cover by cutting a Hamilton path. Unlike an ordinal triangle factor, the shared apex admits revisiting its neighbors.

## Source trail and scope

Session 302/304 supplies the singleton weighted inversion and the fixed apex clock. The new use here is a direct seven-statistic composition across arbitrary ordinal blocks with independent selective apex ports. For comparison, THM-4208 retains endpoint response jets to compose arbitrary contexts; THM-4754 gives an exact arc-level Collatz embedding in a 3-adic tournament, but not control of path-count values; THM-4716 changes clock carries via adjacent swaps. The analogue here is to change a path-cover clock or its port without changing the rest of the seam. No assertion is made that this identity or the seven raw statistics have not appeared elsewhere in the literature or repository. No full readiness closure or conjecture proof follows.

## A strict order law for swapping two adjacent apex letters

Fix the blocks and their order, and change only two adjacent uniform apex letters at positions p,p+1 from 01 to 10. Write x=q_p and y=q_(p+1). Let L be the sum over earlier 1 blocks i<p of (q_i−1) times the product of clocks strictly between i and p. Let R be the dual sum over later 0 blocks j>p+1, with the product strictly between p+1 and j. Let l be 1 if the immediate preceding letter is 1, and r be 1 if the immediate following letter is 0 (missing neighbors give 0). The exact increase in the normalized Hamilton count is

$$
[xy-2x+1]L+[xy-2y+1]R+1+(x-1)(y-1)-l-r .
$$

The earlier 1 blocks now pair with p+1 rather than p: their multiplier changes from x−1 to x(y−1). The later 0 blocks now pair with p rather than p+1: their multiplier changes from y−1 to y(x−1). These are all differences except the newly created central inversion (which is adjacent) and the lost adjacent bonuses on either side. Since x,y≥2, both coefficients are at least 1, L≥l and R≥r. Therefore the normalized increase is at least 1+(x−1)(y−1)≥2. It is exactly 2 if and only if the two blocks are singletons and all nonzero earlier/later terms consist solely of their immediate singleton neighbors, if present.

This lifts the elementary adjacent exchange on parity words to a real apex-arc move, but with the opposite lesson from the fixed 3-smooth Collatz ladder (THM-4716): the increment is a pair of endpoint-dependent sums, not one universally known monomial. In particular there is **no nonempty count-preserving sequence of these increasing swaps**. A true neutral compensation must either switch direction or use an interface change outside this family. When both endpoint letters are unchanged, the Hamilton-cycle count stays unchanged, so the move cleanly separates the cycle and path observables.

## Physically changing a block without changing its one-apex response

If an ordinal block is uniform to the apex, its entire normalized response is specified by its letter and $q$, even when all *other* blocks have arbitrary selective incidences. Indeed on letter 1 the seven normalized coordinates are $(1,s,e,f,a,j,d)=(1,1,0,q-2,0,q-2,0)$, with the reflected statement on 0. Hence changing its internals at fixed $q$ multiplies the total by the ratio of the old and new local $P$ (new divided by old), in every such context. Any number of such switches compose independently. This is a change of actual tournament, not merely replacing a scalar. For two candidate blocks, necessity of equal $q$ for equality of every normalized response is detected with one opposite uniform singleton in a two-block context.

There are order-seven blocks with $(P,F,q)=(29,522,20),(33,594,20),(37,666,20)$. Using the three in the same uniform slot produces the exact arithmetic progression $29K,33K,37K$ in every shared-apex context, including selectively attached other blocks. Here $K$ is an integer: the count ratios are exact and factoring only the switched block’s $P$ leaves integer coefficients, since $q=20$ is integral. Unlike a switch between orders, this edit can be made on the same named vertex set without changing any exterior arc. A second same-clock bank has $q=8$ and counts $1,3,13,15$; the respective block codes are $(n,c)=(3,0),(4,2),(5,41),(5,76)$. So fixed-clock rewrites can insert or remove a factor 3 in the first two of those choices, or replace it with a factor involving 5 or 13. They do not make an arbitrary exterior count integral after division by 3: divisibility is supplied by the actual local block.

The first bank also gives a sharp *context change* caution. Its third ordered-cover counts are respectively $3402,3714,4026$, with different ratios to $P$. In the cyclic blow-up $C_3[B,TT_2,TT_2]$, the $B$ block may occupy three runs. Enumerating the three possible initial letters of the cyclic itinerary gives contributions $11P+17F+4F_3$ (where $F_3$ counts ordered covers of $B$ into three nonempty directed paths). Thus equality of the normalized one-apex response does not survive this enlarged context. The relevant obstruction is the number of available separated runs, not just the total number of vertices.

Encoding for the order-seven examples: the bits of $c$ are read least significant first on the lexicographic list of vertex pairs $(i,j)$, $i<j$; bit 1 means $i\to j$. The codes are $245600,168491,1104621$ in the order of the display. The companion replay counts the first three ordered covers by disjoint subsets and checks each identity.

### A neutral transport of a real tripling factor

Take two uniform slots (they need not have the same letter), each with
$q=8$. The two realizations $U=TT_3$ and $V=C_3\to TT_1$ have identical
normalized response and $P_U=1,P_V=3$. Using $U,V$ in the first and second
slot, respectively, or $V,U$ there, gives the same total $H$ and the same
total order, with exactly the same other blocks and all their apex ports.
This transports an actual cyclic module from one slot to the other, rather
than just relabelling a scalar 3 in a formula. The triangle is peelable in
either position: outside the triangle, travel from its out-neighbors back
to its in-neighbors requires the single apex, by ordinal order. There can
be at most two visits. Replacing that triangle by $TT_2$ returns the
slot to $U$ and divides $H$ by 3, leaving a markable $TT_2$ for the
inverse move. In particular, a readiness condition that requires a
triangle at a specified one of the two slots can be repaired at exact
target, provided that the original factor is available in the other slot.

This is not a source of additional factors 3. With $m$ such slots,
the swaps preserve the number of factor-3 slots and ordinal positions
outside the slots; they only relocate the existing budget. They need not
preserve all old arcs or the individual sizes of the two slots. A stronger
return would have to replenish the budget, or connect this transport to a
separate dyadic move which consumes one factor and prepares another.
