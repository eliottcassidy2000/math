# An exact mixed-power surface two moves out from the transitive root

Use the actual tiling-line move $K_P$, which retains the arcs of a marked directed Hamiltonian path and reverses all other pairs. Start in the transitive tournament on $0,\ldots,n-1$, and let $P$ be its unique increasing path. By [THM-4731(6)](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md), this is the unique blue edge out of the root; its partner $T_n^*$ keeps increasing adjacent arcs and has all nonadjacent arcs decreasing. For $n\ge4$ the reflection $i\mapsto n-1-i$ is its unique anti-automorphism.

Split the increasing path into three consecutive blocks of sizes $p,q,r\ge2$. Read them **in reverse block order, increasing within each block**, to make a path $Q$ in $T_n^*$. Indeed, a between-block seam is decreasing by at least 2. The two moves sum modulo 2 to $E(P)\triangle E(Q)$. Consequently the endpoint is the original transitive tournament with exactly four arcs reversed: the two adjacent block cuts and the two extreme arcs from each neighboring pair of blocks. This is a literal metagraph walk of length two (a distance upper bound).

## The middle endpoints are the hubs

Take $x$ and $y$ to be the first and last elements of the middle block. We have $x\to y$, and they cover all four reversed arcs. Removing them leaves an ascending spine. Put $a=p-2$ and $c=r-2$; the normalized word, with bit convention $2[x\to v]+[y\to v]$, is
$$
                         1\,0^a\,2^q\,3^c\,1.
$$
For instance, the first and last vertices of the left block see bits 1 and 2 because their reversed arcs come from $y$ and $x$, respectively. The right block's first and last vertices see 2 and 1. The intermediate $q-2$ sites of the middle block also have digit 2, so the three pieces combine as shown. This word is strong: digit 1 gives the return $y\to1\to x$ and $x\to y$; the first spine site is reached from $y$ and the last reaches $x$.

On the second line $Q$ is blue exactly when $p=r$ (equivalently $a=c$). Indeed, reversal of positions together with the unique anti-automorphism of $T_n^*$ reverses the sequence of increasing blocks; the blocks are the maximal increasing-by-one segments of the path, hence recovered uniquely. Thus the two-variable chart has a genuine blue diagonal, not merely two arbitrary roles called left and right.

## Exact ordinary count and inverse

Put $A=2^a,C=2^c$. The whole count, including every unmarked path, is
$$
 H=1+2^{q+1}(A+C)+2\cdot3^qAC-[a=0]C-[c=0]A.
 \tag{1}
$$
To prove the all-length formula, hold the occupied slots of each equal-digit block fixed, together with the hub order and assignments outside that block. The internal placements are surjections onto one, two or three slots, with counts $1$, $2^s-2$, $3^s-3\cdot2^s+3$ at positive length $s$. Hence the dependence is separately in the span of $1,2^s,3^s$. Use three consecutive positive lengths per parameter on the interior and per positive parameter on the faces $a=0$ or $c=0$; verify each grid by independently counting paths in the tournament with four literal reversed arcs. The replay proves (1) on these four strata. It also checks the two explicit metagraph moves and the color of the second line.

The **boundary** $a=0,c\ge1,q\ge3$ is exactly the count family in [S336's cubic/overlap inverse](https://github.com/eliottcassidy2000/math/blob/8d2d70178f9d17b191b93dcaffa6c1a8e25dbcb6/05-knowledge/results/session336_overlap_inverse_and_intercyclic_20261114/cubic/BALANCE_CUBIC_OVERLAP.md): set $\beta=q-2,\alpha=c$ in that count formula. This is an equality of the parametrized counts, not a claim that the S336 tournaments and the four-reversal tournaments are isomorphic. It provides a second construction at literal metagraph distance at most two for the whole range of that inverse.

For the new interior $a,c\ge1$, fix the tournament order $n$ and middle width $q$. Then $AC=P_0=2^{n-q-4}$ is determined. For a prescribed target $N$, test whether
$$
 S=(N-1-2\cdot3^qP_0)/2^{q+1}
$$
is integral, and whether the roots of $t^2-St+P_0=0$ are **both** integer powers of 2 at least 2. Exactly then a rooted walk realizes the prescribed count at this $n,q$ in the interior; the two roots give $2^a,2^c$. The repeated-root case is exactly the blue diagonal. Without fixed order, enumerating the two short widths $q,a$ and solving the linear equation for $C$ is enough: each is bounded logarithmically by the positive summands in (1). The boundary is covered by its earlier valuation inverse, while the corner $a=c=0$ has explicit count $-1+2^{q+2}+2\cdot3^q$. In particular the union of this whole family has only $O((\log X)^3)$ values at most $X$, so this exact inverse cannot guarantee coverage of every admissible target.

## What the projection forgets

Lock the four sites in increasing order: the initial 1, first 2, last 2, final 1. The small extremes provide outgoing contacts from the two hubs, and the large extremes admit the incoming matching. The interior gaps have $\alpha=a+c$ doubly eligible and $\beta=q-2$ triply eligible sites. Thus at fixed $n,q$, the overlap monomial $2^{a+c}3^{q-2}$ does not change at all. Nevertheless the whole count sees the split $a,c$ through $A+C$. Up to reflection it has exactly $\lfloor(a+c)/2\rfloor$ distinct values in the interior: $2^a+2^{s-a}$ strictly decreases with $a$ for $1\le a\le\lfloor s/2\rfloor$, $s=a+c$. At $n=11,q=3$ the split $(a,c)=(1,3)$ gives 1025 over a black line and $(2,2)$ gives 993 over a blue line, with the same locked weight 48. This is an exact rooted example of a scalar projection losing information even after fixing the order, the overlap, and the number of moves.

Run **python replay_root_surface.py**.

## The entire two-hub cone after two moves is sparse

The preceding construction used three blocks, but the essential limitation holds for **every** second move. Every directed Hamiltonian path $Q$ of $T_n^*$ is obtained by reversing the order of a composition of the increasing chain into consecutive increasing blocks, with no two successive parts of size 1. To see this, once a path reaches vertex 0 it can only continue through $1,2,\ldots$ until no next site is available, ending there. The earlier path is on a contiguous upper interval, so use induction. Conversely any such reversed composition is legal: the seam jumps downward by the sum of its neighboring sizes minus 1, at least 2. The maximal increasing runs recover the composition uniquely.

Following these two moves, each seam of $Q$ has the same short-plus-long cut as above, and there are no reversals between nonadjacent blocks. A seam between sizes $u,v\ge2$ contains a directed triangle on the left block's first and last vertices and the right block's last. If one size is 1 and the other at least 3, use the far endpoint and an interior vertex on the long side, and the singleton, to get a directed triangle. The only seam types not guaranteeing a triangle are $(1,2),(2,1)$; at either type all arcs between those two blocks reverse. Therefore **any three consecutive blocks contain a triangle**: if both seams are of the exceptional types, the three blocks cyclically dominate. Disjoint triples of blocks give disjoint directed triangles. If the endpoint is a two-hub tournament, both hubs jointly hit all these triangles, so the number of blocks cannot exceed 8.

For any block of size $s\ge3$ in a composition with at least two blocks, choose any neighbor. Retain the size-$s$ block and the opposite far endpoint across their seam, omitting the size-$s$ block's near endpoint if the neighboring block is a singleton. This induced board of at least $s$ vertices is transitive except for the single reversed arc from last to first. On $t$ vertices that board has $1+2^{t-2}$ Hamiltonian paths: the increasing one plus one for each allocation of the intermediate vertices before or after the long reversed arc. Insertion of the deleted vertices injects those paths into paths of the endpoint, so $H\ge1+2^{s-2}$. Hence at \(H\le X\), every block size is \(O(\log X)\), independently of the ambient order and any singleton padding.

Together these facts give
$$
 \#\{H\le X:\text{some two-hub endpoint at distance at most 2 from its transitive root}\}
                         =O((\log X)^8).
$$
The root itself only gives 1, and the one-move partner has a directed triangle on every consecutive triple, so if it is two-hub its order is at most 8. This is an obstruction to a fixed-depth **rooted search through two moves**, not to the full two-hub class or to successive depths that grow with the target.

Moreover it misses many **already realized** counts even in the much smaller strong one-hub class. For a binary apex word $u$ of length $k$, let $A(u)$ be its binary integer and $\epsilon(u)$ its last bit. The [one-apex inversion formula](https://github.com/eliottcassidy2000/math/blob/1b91742ce/05-knowledge/results/session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md#3-exact-binary-product-approximation) says an inversion at distance $d$ contributes $2^{d-1}$ to $H-1$, plus 1 when $d=1$. In particular
$$
 H(u0^t)=H(u)+(2^t-1)A(u)+\epsilon(u)
$$
for $t\ge1$: each new last zero sees the binary integer of its old prefix, plus 1 only when that prefix ends in 1. Allowing every possible pair in the underlying sum gives the elementary upper bound
$$
 H(u)\le 1+(k-1)+\sum_{d=1}^{k-1}(k-d)2^{d-1}=2^k-1.
$$
Take $u$ through all $k$-bit words starting in 1 and set $t=k+1$. These are strong one-hub boards (first bit 1, last bit 0). Their $2^{k-1}$ whole counts are distinct: the coefficient of $A(u)$ is $2^{k+1}-1$, strictly larger than the maximum variation of $H(u)+\epsilon(u)$. Every count lies below $2^{2k+1}$. Thus the number of already-attained strong one-hub counts at most $X$ is $\Omega(\sqrt X)$, growing faster than the polylogarithmic radius-two bound, and the rooted cone misses $\Omega(\sqrt X)$ such counts. The replay checks the padding identity against literal orientations at small orders.
