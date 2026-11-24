# Pointed-spread tests for a Moore graph with a Hoffman-tight coclique

The premise is a *construction choice*, not a condition known for every putative Moore57: try building one with a coclique attaining Hoffman's bound, size 400. The associated Steiner $2$-$(400,8,1)$ design and candidate PG$(3,7)$ have been proposed before; see [Problems in Algebraic Combinatorics](https://doi.org/10.37236/1224) and [A note on Moore graphs](https://ajc.maths.uq.edu.au/pdf/77/ajc_v77_p373.pdf). The following are stronger precise construction and screening rules.

Let $N$ be the point–block incidence matrix of any such design. There are 2850 blocks, with 57 blocks through each point. Let $D$ be the adjacency on disjoint blocks, so $D=J+7I-N^\top N$ and its degree is $2401$. A symmetric zero-diagonal binary matrix $H$ gives a completion by point–block incidence exactly when
$$H^2+H=49I+D.$$
Indeed the diagonal forces degree 49, and an $H$-neighbor is disjoint from its center. Two such neighbors must also be disjoint, as otherwise they already have a common point and $D$ is zero on their pair. Thus these 49 blocks partition the 392 points outside the center. This gives $NH=J-N$, so the other two blocks of the Moore identity follow from the design. Conversely they are necessary by the Moore identity.

Consequently select one spread $S_B$ of 50 blocks containing each block $B$, with reciprocity $C\in S_B\iff B\in S_C$. Require only that every pair of disjoint blocks cooccur in at least one selected spread, counting pointed selections with multiplicity. This covering condition suffices: reciprocity puts each block into exactly 50 selected spreads. For fixed $B$, its $50\cdot49=2450$ paired incidences already require one each for its 2401 disjoint peers plus one extra each for its 49 reciprocal neighbors, whose two own spreads both contain the pair. These exhaust all pair occurrences. Thus a disjoint pair cooccurs exactly once unless reciprocal (then twice), intersecting pairs never cooccur, and $H^2+H$ is exactly as above. The 2850 selected spreads are necessarily distinct, since two closed stars cannot share 50 vertices.

This gives a Boolean/ILP test over a chosen bank of spreads: $x_{B,S}$ for $S\ni B$, one selected per center; use one linear reciprocity equality per disjoint pair and one linear covering lower bound per disjoint pair (no upper constraints necessary). A completely unpointed bank of 2850 distinct spreads can be tested using multiplicities: require each block in 50 spreads, each disjoint pair in one or two, and require that each spread contains one block whose pairs with all others in the spread have multiplicity two. Repeat pairs have degree49 and every such center then identifies its unique spread, recovering the pointing.

## Spectral screen

For a symmetric reciprocal spread selection put $H_{BC}=1$ when $C\in S_B\setminus\{B\}$. On the constants $H$ has eigenvalue49, and on the centered image of $N^\top$ it has eigenvalue $-1$ by $NH=J-N$. Thus on $\ker N$ the traces of $H$ and $H^2$ are fixed by degree and zero diagonal. Their sum is $56\dim\ker N$. Therefore if the extreme eigenvalues of $H$ on $\ker N$ lie in $[-8,7]$, all its eigenvalues there must be $-8$ or $7$ (sum the nonnegative terms $56-h^2-h$). This is necessary and sufficient. Equivalently choose an orthogonal projector $E$ of rank 1330 contained in $\ker N$, so that
$$H=-8I+\tfrac18N^\top N+15E$$
is a zero-diagonal binary matrix. The diagonal of $E$ is $7/15$, entries on intersecting blocks are $-1/120$, and entries on disjoint blocks are $0$ or $1/15$. For these formulae use $NN^\top=56I+J$.

## Resolution-bank obstruction and calibration

A resolution partitions the block deck into spreads. A single fixed resolution could give only a union of 50-cliques, not a Moore completion. More generally, if a bank of $s$ resolutions contains all selected pointed spreads, assigning each center to a resolution containing its spread gives a proper $s$-coloring of the disjointness graph $D$: disjoint centers' selected spreads intersect, so they cannot belong to the same resolution. For PG$(3,7)$, the [Chowdhury–Godsil–Royle line-coloring theorem](https://arxiv.org/abs/math/0507319) gives $\chi(D)=7^2+7=56$, hence at least **56 separate resolutions** are needed. This is not a sufficiency result.

The [standalone checker](spread-check.py) constructs a graph from the 15 points and 35 lines of PG$(3,2)$ by the pairwise compatibility CSP. There are 56 spreads and 240 resolutions. The 35 chosen star spreads have four internal edges and maximum load four. The other 21 all have zero internal edges and maximum load two; no resolution uses only these 21. Of the 240 resolutions, 30 contain seven chosen star spreads and 210 contain four. This is not an enumeration for PG$(3,7)$.

A general contact-count obstruction applies to any such design. For any spread $R$ of size $h+1$ (where $h$ is the degree of $H$), set $\ell_z=|N_H(z)\cap R|$. Then $\sum_z\ell_z=h(h+1)$ and $\sum_z\binom{\ell_z}{2}=\binom{h+1}{2}-e_H(R)$, because a disjoint pair is either adjacent or has exactly one common line-neighbor. If $n_j$ counts loads equal to $j$, subtract to get
$$
n_1=2e_H(R)+\sum_{j\ge3}j(j-2)n_j.
$$
Thus a resolution of independent spreads cannot have every load at most two when $h=49$ is odd: then there are no load-one entries, and every vertex would have even degree. Even with loads capped at three, an independent resolution must have at least 713 load-three entries across its frames: each frame has $n_1=3n_3$, while each line has an odd load in at least one frame. This is conditional on the existence of the stipulated resolution.
