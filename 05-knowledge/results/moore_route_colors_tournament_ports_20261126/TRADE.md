# Sharp Moore trade support

**Theorem.** Let $A,A'$ be distinct adjacency matrices of diameter-two Moore graphs of degree $d\ge2$ on the same labeled vertex set. Put $S=\{u:A_{u,*}\ne A'_{u,*}\}$. Then $|S|\ge2d$. Swapping two vertex labels realizes equality. Moreover, for $d\ge4$ not divisible by six, **every** equality case is a single label transposition. In particular degree 57 has a 114-vertex minimum support, at which no new isomorphism class is possible. All consequences for degree 57 are conditional on completed graphs.

## Outside colors from the matrix equation

Write $\Delta=A'-A$, $s=|S|$. It is symmetric with zero diagonal and row sums zero and supported on $S\times S$. Subtracting the quadratic equations gives
\[
 A\Delta+\Delta A+\Delta^2+\Delta=0.
\]
Taking rows outside $S$ yields $A_{\bar S,S}\Delta_{S,S}=0$. Thus for every exterior center, the sum of the edit rows of its neighbors in $S$ is zero. In particular a center cannot see exactly one vertex of $S$. Call an exterior center **clean** when it sees exactly two, and those two vertices a clean pair. Clean pairs form a simple graph on $S$, since there is at most one common neighbor. Along a clean edge the two edit rows are negatives of one another.

For $u\in S$ write $a_u=|N_A(u)\cap S|$ and $p_u=\sum_{w\in N_A(u)\cap S}(a_w-1)$. The latter counts precisely the other vertices of $S$ reached by a two-edge route with center in $S$. For each exterior neighbor $z$ of $u$, let $r_z=|N_A(z)\cap S|\ge2$. Radius-two uniqueness gives
\[
 s-1=a_u+p_u+\sum_{z\notin S,\ z\sim u}(r_z-1).
\]
Set $\ell=s-1-d$, so $\ell=p_u+\sum_z(r_z-2)\ge0$. If $c_u$ is the clean degree of $u$, its number of nonclean exterior centers is at most $\ell-p_u$, so
\[
 c_u\ge d-a_u-\ell+p_u,\qquad
 \sum_{u\in S}c_u\ge s(d-\ell-1)+\sum_{u\in S}(a_u-1)^2.
\]
The sum uses $\sum p_u=\sum a_u(a_u-1)$.

## One loss, or many

Let $t_v$ count the neighbors lost at $v$; it also counts gains. If any $t_v=1$, write $u$ for its lost neighbor and $w$ for its gained neighbor. The exterior equation in column $v$ says $N_A(u)\setminus S=N_A(w)\setminus S$. These identical sets have size at most one.

- If $u,w$ are adjacent they have no common neighbor. Their full neighborhoods lie inside $S$ and are disjoint, giving $s\ge2d$.
- If they are nonadjacent with a common neighbor outside $S$, their two sets of $d-1$ inside neighbors are disjoint. Neither contains $u$ or $w$, so again $s\ge2d$.
- If they are nonadjacent with the common neighbor inside, their union inside has $2d-1$ vertices and omits both $u,w$, giving $s\ge2d+1$.

If instead every $t_v\ge2$, the clean graph is a matching. Indeed a path of two clean edges forces its two endpoints to have identical edit rows. They would share at least two *lost* neighbors in $A$, impossible. The clean-degree bound now gives $\ell\ge d-2$, hence $s\ge2d-1$. At $s=2d-1$ it requires $\sum c_u\ge s$, impossible for a matching on an odd number of vertices. This proves $s\ge2d$ in both cases.

Transposing two labels changes the rows at the two swapped vertices and exactly $2(d-1)$ other rows: the exclusive members of their neighborhoods after the swapped vertices are excluded. This uses the zero common neighbors for an edge, and the one common neighbor for a nonedge.

## Equality types

Suppose $s=2d$ and $t_v=1$ occurs. In either possible case above, let $L,R$ be the exclusive neighbors of $u,w$ after excluding $u,w$ themselves. These are disjoint size-$d-1$ sets, and $S=L\sqcup R\sqcup\{u,w\}$.

- In the adjacent case no cross edge $L$–$R$ is possible (it would give a quadrangle). Any cross pair has its unique common neighbor outside $S$; the center is clean because it can meet at most one neighbor of each of $u,w$ and neither of $u,w$.
- In the nonadjacent case cross edges form a matching (two edges from a site into the other side would give a quadrangle). Nonadjacent cross pairs again have no possible interior common neighbor. The exterior common neighbor of $u,w$ cannot touch either side, by girth; all other exterior centers cannot see $u,w$, and can see at most one vertex on each side. Hence every nonadjacent cross pair is clean.

Thus the clean graph on $L\sqcup R$ contains $K_{d-1,d-1}$ minus at most a matching and is connected for $d\ge4$. Its paired row signs propagate the edit row at $v$ to all of $L,R$ (the vertex $v$ belongs to $L$). They all lose their old parent and gain the other. All other edits are thereby determined except possibly $uw$, which is unchanged by zero row sums. This is exactly the transposition $(u\,w)$.

Now suppose $s=2d$ and every $t_v\ge2$. Then $a_v\ge2$, and with $\ell=d-1$ the clean-degree inequality forces $\sum c_v=s$ and $a_v=2$ everywhere. Thus the clean graph is a perfect matching, the entire internal $A[S]$ is removed, and both $A[S]$ and $A'[S]$ are 2-regular. Every exterior nonclean center meeting $S$ sees exactly three of its vertices: at each $v$, $p_v=2$, so $d-3$ nonclean exterior centers exhaust the slack $d-3$ and must all have $r_z-2=1$.

Contract each clean pair just for indexing. Opposite edit rows and symmetry imply that each nonzero block between two pairs is a checkerboard of signs; each zero diagonal block is forced by the zero diagonal. Thus on the $d$ pairs there is a simple quotient support graph $Q$ of degree two. Along each $Q$ edge, $A[S]$ uses one matching between the two pairs while $A'[S]$ uses the complementary matching. A quotient triangle is impossible: switching all three matchings reverses its two-sheet holonomy, and one of $A[S],A'[S]$ would contain triangles.

Each nonclean center picks three distinct pairs (it cannot contain both vertices of a pair already sharing a clean center). Modulo two its zero edit-row sum implies that at each quotient vertex the two cyclic neighbor flags for that three-pair selection agree. A nonempty selection on an odd cycle must be the whole cycle; on an even cycle it is one parity or the whole cycle. A selection of size three is therefore either a full triangle or one parity of a six-cycle. The first is excluded. For $d>3$ every vertex of $S$ has a nonclean center, so all cycles of $Q$ have length six. Consequently $6\mid d$, completing the rigidity claim when $6$ does not divide $d$.

This is a separation theorem for two finished graphs on a common labeling. It does not bound the number of edits needed for an unfinished graph or disprove the existence of a first completion.
