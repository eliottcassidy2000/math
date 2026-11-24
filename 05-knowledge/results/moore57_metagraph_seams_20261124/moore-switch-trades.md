# Exact Moore-switch obstructions

For a degree-$k$ Moore graph of diameter two, $A^2+A=(k-1)I+J$. If $A$ and $A'=A+\Delta$ are two such adjacency matrices on the same labels, their degree is equal, so $\Delta\mathbf1=0$, and
$$A\Delta+\Delta A+\Delta^2+\Delta=0.$$

## One bipartite complement never gives a genuine trade for degree at least four

Suppose all changes complement the cross edges between disjoint nonempty sets $X,Y$, leaving every other edge unchanged. Equal degrees require that the old bipartite adjacency $B=A_{X,Y}$ has half the entries in every row and column. Put $x=|X|,y=|Y|$, both even. A pair of $Y$-neighbors cannot occur at two different $X$-vertices by girth, and symmetrically, so
$$x\binom{y/2}{2}\le \binom y2,\qquad y\binom{x/2}{2}\le\binom x2.$$
If $x,y\ge4$, the only possibilities are $4\times4$, $4\times6$, and $6\times4$.

For an exterior vertex $z$, the Moore difference identity in row $z$ says $A_{z,X}K=0$ and $A_{z,Y}K^\top=0$, for $K=J-2B$. In the $4\times6$ case the six columns of $B$ are all six pairs from $X$. Their sign columns span $\mathbf1^\perp$, forcing an exterior vertex to meet all or none of $X$. It cannot meet all, since every pair of $X$ already has its unique common neighbor in $Y$. Likewise $X$ must be independent. Its vertices therefore have full degree three. Transpose for $6\times4$.

In the $4\times4$ case $B$ is an eight-cycle (two four-cycles are forbidden). Around that cycle, the equations $A_{z,X}K=0$ say that each consecutive pair of $X$-coordinates has sum half the total, so the coordinates are constant on each of the two opposite pairs. An $X$-vertex can only have an internal neighbor at its opposite site, and any exterior neighbor also neighbors that site. At most one exterior neighbor exists by no four-cycles, and the internal and exterior possibilities cannot both occur by no triangles. Thus the full degree is at most three.

When either side is two, say $Y=\{u,v\}$, its cross neighborhoods in $X$ are complementary; the exterior equation forces $A_{zu}=A_{zv}$ for every exterior $z$, since the two nonzero change rows are negatives. The switch is simply transposition of $u,v$. Thus for $k\ge4$ a single bipartite complement never changes the isomorphism class. A single-cell [Godsil–McKay switch](https://aeb.win.tue.nl/preprints/isomGM.pdf) is a special case, with $Y$ the set of exterior vertices meeting exactly half the cell. Degree-three feasibility is sharp: in the Petersen graph on two-subsets of a five-set (disjointness adjacency), take $X$ the four subsets containing 0 and $Y$ the other six. The accompanying finite calibration also finds fifteen $4\times4$ exact Petersen swaps, although no switch can change its isomorphism class.

## Degree 57: at least 101 changed vertices for any exact trade

Let $S$ be the set of changed-edge endpoints. The nonzero subspace $U=\operatorname{im}\Delta$ is supported on $S$ and contained in $\mathbf1^\perp$. For $v\in\ker\Delta$, the difference identity gives $\Delta Av=0$; by symmetry both the kernel and its orthogonal complement $U$ are $A$-invariant, and $A'$-invariant as well. Both have eigenvalues among $7,-8$ on $U$ and the same trace there, since $\Delta$ vanishes on its kernel and has zero diagonal globally. If $A$ had no eigenvalue $7$ on $U$, trace equality would force both restrictions to be $-8I$, contradicting $\Delta\ne0$. Pick a centered $7$-eigenvector from $U$.

On either sign set $P$ or $Q$, discarding the negative contribution from the opposite sign yields $A_P|v_P|\ge7|v_P|$. Each therefore has spectral radius at least seven. The row sums of $A_P^2$ are at most $|P|-1$ in a girth-five graph: the diagonal contribution is the degree and the off-diagonal contributions count distinct vertices at distance two. Hence each sign has at least 50 sites.

If $|S|=100$, both signs have exactly 50 sites. For $u=|v_P|>0$ we have $A_P^2u\ge49u$. Summing coordinates and using symmetry and row sums at most 49 forces equality and every row sum to be 49. Moreover $A_P(A_Pu-7u)+7(A_Pu-7u)=0$ with both terms nonnegative, so $A_Pu=7u$. Likewise for $Q$. Each induced sign graph therefore has diameter at most two, and the original eigenvector equation gives no cross edges. An outside vertex touches at most one site per sign. Every cross pair $x,y$ has its unique common neighbor outside $S$, so the block $A_{V\setminus S,S}$ contains every row $e_x+e_y$ for cross pairs $x,y$. The exterior matrix equation then gives $\Delta_x+\Delta_y=0$ for every cross pair, where $\Delta_x$ is the row of $\Delta$ at $x$. Thus all $P$-rows are identical and all $Q$-rows are their negatives. Since $S=P\sqcup Q$ and $\Delta$ has zero diagonal, the common row vanishes. This contradicts the trade and proves $|S|\ge101$.

The underlying lower bound for a general degree is $2(\lceil r^2\rceil+1)$ for the positive centered eigenvalue $r$, with the same strict saturation argument if $r^2$ is integral. A transposition of two degree-57 vertices touches exactly 114 vertices, so no labelled lower bound can exceed 114.
