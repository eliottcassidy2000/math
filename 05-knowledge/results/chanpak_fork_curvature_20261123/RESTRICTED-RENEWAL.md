# Exact two-reading renewal on two marked wrapper towers

Write $Q_t(A)$ for the number of ordered $t$ path slots of a tournament $A$, allowing empty slots, with every vertex assigned exactly once. The two fork readings are $D=Q_1$ and $Q=Q_2$. A uniform graft with reentry cannot update those readings for completely arbitrary fillers. The following is an unbounded restricted family on which there is an exact inverse, and hence a physical renewal rule on the readings.

## A two-parameter restricted family, with an optional seed bit

Let $T_m$ be the transitive tournament on $m$ vertices, and define $W(A)$ by $A\to b\to c\to A$ for two new singleton anchors. Its exterior matching size relative to $A$ is one. If $F(t)=Q_t(A)$, its slot polynomial is
$$
\mathcal L F(t)=tF(t+1)+t^2F(t)-tF(t-1).
$$
To see this, if $b,c$ are in different slots their contribution is $t(t-1)F(t)$; the slot containing $b$ ends in $b$, and the one containing $c$ begins in $c$. If they are in the same slot in order $b,c$, they are consecutive, and the interior can occupy two portions of that slot independently, giving $tF(t+1)$. If their common order is $c,b$, there must be a nonempty path of $A$ between them, giving $t(F(t)-F(t-1))$. There are no other possibilities. Thus the expressions agree at all positive integers and are a polynomial identity. In particular $D(W(A))=D(A)+Q_2(A)$ for nonempty $A$.

Take $m=1$ or $2$, $k\ge1$ and $s\ge0$. The filler is the ordinal join
$$
M(m,k,s)=T_s\to W^k(T_m).
$$
Its strong nucleus is exactly $W^k(T_m)$; the remaining $s$ components are singleton sources in ordinal order. Put $F_{m,k}(t)=\mathcal L^k(t^m)$. Its readings are
$$
D=F_{m,k}(1),\qquad Q=2^sF_{m,k}(2).
$$
A physical **pad** operation adds one more ordinal source, i.e. $s\mapsto s+1$. A physical **wrap** operation adds the two anchors around the marked strong nucleus only, while the old singleton sources continue to beat both new anchors, i.e. $k\mapsto k+1$. Both are repeatable in arbitrary order. Wrapping is uniform relative to the nucleus but not to the whole padded filler when $s>0$; this distinction is necessary. The factor-seven host can keep its fixed exterior contacts and its original two marked exit arcs under either operation, and the resulting actual fork is again evaluated at the new pair.

**Exact inverse.** If $m=1$, the first wrap forms a directed triangle that remains a homogeneous module at every subsequent stage. Its order-three rotation acts freely on Hamilton paths, so $3\mid D$. If $m=2$, the entire nucleus always has even order and $Q_2\equiv D\pmod3$ by the odd-cycle packing evaluation at $t=2$. Each wrap therefore sends $D$ to $2D\pmod3$ and keeps it nonzero modulo three. These cases identify $m$ uniquely. For either fixed seed $F_{m,k}(1)$ strictly increases, by at least a factor three at each wrap, because $Q_2\ge2D$. Enumerate $k$ until its height first equals or exceeds the input $D$; reject on strict excess. Finally require $Q/F_{m,k}(2)$ to be an exact nonnegative power of two, and read its exponent $s$. This recovers or rejects the full state using only the given two integers. With positive tests it constructs the filler and computes both exact successor pairs. The search has $O(\log D)$ wrap stages; the polynomial degree is $2k+m$ and its coefficients can be updated by binomial expansion. This is an unbounded restricted decoder, not a universal two-reading map for arbitrary fillers, and it does not infer arbitrary absolute fork heights.

## A sharp four-vertex matching threshold example

For a concrete minimum-size exterior of matching depth two, take $B=\{b_0,b_1\}$ beaten by the filler and $C=\{c_0,c_1\}$ beating it. The only arcs of $B\to C$ are $b_0\to c_0,b_1\to c_1$; reverse the other two, and orient $b_0\to b_1$ and $c_0\to c_1$. Four is the minimum possible exterior order for two disjoint bridges. Write $P_j$ for ordered nonempty $j$-path covers of the filler. The exact output reads are
$$
D'=13P_1+10P_2+2P_3,\qquad
Q'=104P_1+116P_2+48P_3+6P_4.
$$
The maximum occupied indices are three and four by the matching-depth bound. The coefficients follow by testing transitive fillers of sizes $1,2,3$ and $1,2,3,4$ respectively and using triangular interpolation $P_j(T_m)=j!S(m,j)$; the companion code uses independent subset path DP. The factor two on $P_3$ achieves the absolute minimum: a three-run path uses both matched arcs consecutively as two bridges, and the two bridges can be arranged in either order with no other bridge choices.

For the two nine-vertex equal-$D,Q$ hidden-state fillings, the output pairs are $(746195,17968288)$ and $(746675,17982688)$, so already the first reading differs by exactly 480. This gives a size-minimal literal test of the transition from matching depth one to two.

Run `python3 matching-restricted.py`. It additionally checks every small tower with both seeds, four wrap depths and three pad widths against exact subset DP.
