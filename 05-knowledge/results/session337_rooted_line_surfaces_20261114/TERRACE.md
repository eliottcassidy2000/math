# The entire first strong triangle terrace is a Boolean cube

Use the repository's literal line metagraph: a tiling names a Hamiltonian path and a line flips all pairs **not on that path**. A blue line has a reflection-fixed tiling as its presentation. The unique blue edge out of the transitive node and its endpoint are inherited facts: at $n\ge3$ the endpoint is the rigid self-converse full-tiling tournament, with tribonacci path count. See [THM-787](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-787-transitivity-flow-on-the-metagraph.md) and [THM-4731](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md). The next result is on the **whole** constant triangle level of that endpoint, not just its line neighbours.

## Classification and the right adjacency

For every $n\ge4$, a strong tournament has at least $n-2$ cyclic triangles. The extremal score multiset is unique:
$$
       1,\ 1,\ 2,\ \ldots,\ n-3,\ n-2,\ n-2.
$$
The first blue mate has exactly this multiset. To prove the extremality, sort the outdegrees increasingly. The sum on the first $k<n$ vertices is at least $\binom k2+1$ because the subtournament has $\binom k2$ internal arcs and strongness requires at least one outgoing arc from it. The displayed score sequence attains equality for every $k$. Hence it majorizes every strong score sequence; strict convexity of squared scores gives a unique maximal sorted list. The triangle identity $C_3=\binom n3-\sum_v\binom{s_v}{2}$ gives the claimed minimum.

Alternatively sort the vertices in **descending** score order, and label them $0,\ldots,n-1$. A class on this terrace has exactly one backward arc across every cut between successive sites, by the same prefix-sum equality read from the other end. The intervals of backward arcs therefore tile the successive cuts. Thus the only backward arcs are those of a chain
$$
  s_k\to s_{k-1}\to\cdots\to s_0,\qquad
  0=s_0<s_1<\cdots<s_k=n-1.
$$
Conversely reversing just those chain edges in the transitive order gives the required scores and a strong tournament. For a single reversed edge there is a Hamiltonian cycle. For two or more the chain anchors lie on a directed cycle and all remaining sites sit between two successive anchors, so they enter that same strong component.

The only equal-score choices are swaps of the two extreme pairs. Swapping positions 0 and 1 toggles whether site 1 is an anchor, and the last-pair swap likewise toggles site $n-2$. To see this, an edge ending at old site 1 transfers its backward endpoint to new site 0 when site 1 was selected; when it was absent it transfers in the other direction and the end-pair edge becomes backward. No other arc changes its relative position. This gives a unique normalized representative per class with neither site 1 nor $n-2$ selected. The remaining $n-4$ internal bits are invariant under graph isomorphism. There are therefore exactly $2^{n-4}$ classes on the first strong triangle terrace; converse reverses the bits, and there are $2^{\lceil (n-4)/2\rceil}$ self-converse classes.

Inserting an additional anchor $b$ between consecutive anchors $a,c$ is exactly a directed triangle reversal of $a\to b\to c\to a$. It replaces the reversed arc $c\to a$ by two reversed arcs and leaves every score fixed. Thus the terrace has a **separate, score-preserving triangle-move cube**. It is not a claim that those edges are blue or black lines. In fact the old blue parity law forbids level blue lines at odd orders. Rooting this cube at either the empty or full internal subset gives a precise branching picture on the terrace without claiming the actual line graph is a tree.

## How much of the terrace is reached by the first two actual lines?

For $n\ge5$ the radius-two ball of the **actual** line graph (both colors) meets this entire terrace in exactly two classes: the tribonacci mate and the single long-backward-arc class of count $1+2^{n-2}$. Indeed the [second-shell formula](branch/ROOTED_BLUE_SHELL.md) for a composition with endpoints $a,b$ forces its total interior size to be $[a=1]+[b=1]$ in order to stay at triangle height $n-2$. If both endpoints are at least two, there are exactly two blocks, and their score-matched tournament is isomorphic to the one-long-arc order below (its two end-pair positions are the only ambiguous marks). If exactly one endpoint is a singleton, the required single interior singleton is adjacent to it, forbidden. If both are singletons, the only legal interior has size two and the order is four. The two remaining order-$n$ classes are distinct for $n\ge5$ by their normalized anchor sets: all internal sites versus none, even when their path counts coincide at $n=5$.

For the two-block isomorphism, its remaining ascending seam arcs in the descending model are the first-to-last edge and the edge across the adjacent boundary. The latter switches back by exchanging those two adjacent vertices in the background **order**, without affecting any other pair's direction relative to that order; the long reversed edge remains. This is a change of rank order, not an additional physical move.

## Exact projection onto interval packings

Put $r_i=s_{i+1}-s_i-1$ for the nonanchor gaps, $0\le i<k$. A directed cycle spanning the anchor interval $[i,j]$ must use all backward chain arcs there, because they are the unique backward crossings of the cuts. The rest of the cycle is a forward route from $s_i$ to $s_j$ through a subset of the intervening nonanchors in their unique ascending order. For one chain edge a direct return arc is unavailable, but odd cycles already require at least one selected site in that case. If $R=\sum_{\ell=i}^{j-1}r_\ell$, the number of directed odd cycles supported on this interval is
$$
 c_{ij}=\begin{cases}
 2^{R-1} & R>0,\\
 1 & R=0,\ j-i\text{ even},\\
 0 & R=0,\ j-i\text{ odd}.
 \end{cases}
$$
Two cycles are vertex-disjoint exactly when their inclusive anchor intervals are disjoint. The interval weights, not just the triangle height or the degree sequence, are the information a surface projection must keep.

By the odd-cycle packing formula of [Grinberg and Stanley](https://arxiv.org/abs/2307.05569), the partition polynomial of these intervals at $z=2$ equals the ordinary Hamiltonian count. Count on the first $j+1$ anchors by paths in a weighted interval DAG:
$$
 P_{j+1}=P_j+z\sum_{i<j}c_{ij}P_i,\qquad P_0=1.
$$
The unit term skips the last anchor and every other term appends the last occupied interval. Equivalently choose a subset of chain gaps. Each maximal selected run contributes the product of $2^{r_i}$ over its gaps if it contains a positive gap. A run consisting only of zero gaps instead contributes 2 for even length and 0 for odd length. This gives a four-state transfer (closed, open containing a positive gap, open pure-zero odd, open pure-zero even) in the replay.

## Direct unbounded constructions, now strong

1. **All gaps positive.** There is no restriction on the selected subset of gaps, so
   $$H=\prod_{i=0}^{k-1}(1+2^{r_i}).$$
   This makes an arbitrary such product a **strong** tournament on $1+k+\sum r_i$ vertices, saving $k-1$ vertices compared with an ordinal sum of the single-back-arc factors. It stays at the very first strong triangle level. Varying the order of gaps can change the class without changing the target. For given product factors, no search over tournaments is needed.
2. **A single interior zero between two positive gaps $a,b$.** All other positive gaps still factor out. The contribution from this window is
   $$ 2(1+2^a)(1+2^b)-1. $$
   It comes from counting the zero unselected freely and counting it selected unless both neighbors are skipped. This permits prime values as well as composite values. For example $a=1,b=3$ gives 53. It has exact inverse on this sector: given $N$ test $(N+1)/2$ for a product of two such factors. In either product case with two factors, $a+b=\lfloor\log_2 M\rfloor$ for product $M\ne9$; $M=9$ is $(1,1)$. Then $D=M-1-2^{a+b}=2^a+2^b$ recovers the two exponents from its binary expansion (one carry when equal). This is a prescribed-target rule, not just a count formula.

## Why a connected branching picture does not fill the spectrum

Even the score-preserving cube permits both signs of count change along anchor insertion. For $n\ge6$, put $X=2^{n-5}$. The insertions $\varnothing\to\{3\}\to\{2,3\}$ (interior anchors in zero-based notation) give
$$
                 1+8X\quad\longrightarrow\quad5+5X
                 \quad\longrightarrow\quad5+6X .
$$
The first step decreases and the second increases, by amounts unbounded in $n$. Both $(C_3,\text{scores})$ stay fixed throughout. And even all orders of this terrace cannot realize 11. Indeed $H\ge1+2C_3=1+2(n-2)$, so a hypothetical terrace witness of 11 would have $n\le7$; limiting to seven orders' elementary interval recurrence excludes it, whereas the two-hub word 110 has count 11. Thus even a full description of this first terrace and two exact unbounded inverses does not decide arbitrary target existence.

Run **python replay_terrace.py** for comparison against literal path subset DP, the isomorphism counts in small orders, and the exact product and shifted-product decoders.
