# Exact port elimination: a missing edge, a usable detour, and a changed pairing

This develops the endpoint-aware direction of the earlier [square-sum port work](../session32_square_monster_20261023/ports/PORT_TRANSFER.md). A port-to-port walk may be individually possible yet unusable alongside another. The right transfer retains its **interior support**, not just its endpoint count or existence.

## Squarefree elimination

For a directed graph or a bidirected version of an undirected graph, divide the vertices into kept ports $P$ and interior $I$. Work over the finite ring
$$
R_I=\mathbb Z[x_i:i\in I]/(x_i^2:i\in I).
$$
With adjacency matrix $A$ and $X_I=\mathrm{diag}(x_i)$, the effective arc matrix is
$$
E_A=A_{PP}+A_{PI}X_I(I-A_{II}X_I)^{-1}A_{IP}.
$$
The inverse is a finite series because every entry of $A_{II}X_I$ has positive degree. The coefficient of $x_J$ in $(E_A)_{uv}$ counts paths from port $u$ to port $v$ whose interior vertices are exactly $J\subset I$, with no other ports used. For distinct endpoints, this follows immediately from the walk expansion: revisiting any interior vertex makes the term zero. For a diagonal entry, the same expansion instead records returns to the one fixed port.

There is also a literal one-site recurrence. At a stage when an interior vertex $i$ is eliminated, replace every remaining effective arc by
$$
A_{uv}\longmapsto A_{uv}+A_{ui}x_i A_{iv}.
$$
The Schur pivot would give an extra inverse factor $(1-A_{ii}x_i)^{-1}$, but all its positive-degree terms vanish when multiplied by the displayed $x_i$. Thus this is exact even when effective loops have appeared. The resulting port matrix is independent of elimination order. Products across a specified itinerary of distinct ports count only disjoint interior portions, since every repeated variable vanishes; extracting the full interior monomial makes the itinerary span. One can sum over all port itineraries to obtain all Hamilton paths with prescribed endpoint ports. Or, to count all Hamilton paths at once, use a dummy source and sink as the only two kept ports, with the source pointing to every original vertex and every original vertex pointing to the sink.

A scalar reachability transfer forgets this ring. In particular substitution $x_i=1$ is not a ring map from $R_I$: it would allow paths to reuse a vertex that had already been allocated.

## Two signed pairings and a literal cycle exchange

Let $s_0,s_1,t_0,t_1$ be distinct kept ports. Form the $2$ by $2$ submatrix of $E_A$ with rows $s_0,s_1$ and columns $t_0,t_1$. The coefficient of $x_I=\prod_{i\in I}x_i$ in its permanent counts all two-path spanning covers with those port sets, and the determinant coefficient is straight pairing count minus crossed pairing count. This needs no acyclicity assumption, because products already kill every shared interior vertex.

If the four arcs $t_j\to s_i$ are available, one can require a straight or crossed external matching of those arcs. For either fixed matching, a spanning two-path cover with the opposite matching closes to a single Hamilton cycle; one with the same matching closes to two cycles. Indeed the permutation on the two starts obtained by following an internal path and then a closing seam is a transposition or identity respectively. These are counts under the stated port itinerary and required seam edges, not a claim about all cycles without those conditions. Thus the determinant difference is exactly the difference between the two choices of seam in the one-cycle channel, with sign depending on which seam choice is subtracted.

### Square-sum witness

Use the square-sum graph with starts $(6,8)$, ends $(13,15)$ and interior $\{1,3\}$. Its effective matrix is
$$
\begin{pmatrix}x_3&x_1x_3\\x_1x_3&x_1\end{pmatrix}.
$$
The full-support determinant coefficient is $+1$: the only spanning pair is $6-3-13$ and $8-1-15$. Add the interior vertex 10. Only the top-right entry changes, gaining $x_{10}$ from $6-10-15$. The full-support determinant coefficient is now $-1$: the only spanning pair is $6-10-15$ and $8-1-3-13$. All four individual port reachabilities were positive both times. In the same square-sum rule with interior $\{1,3\}$, the ports 4 and 5, and 2 and 7, have effective entry 1 by their direct edges; 2 and 4 have entry zero, whereas 8 to 6 has entry $x_1x_3$. These are different obstructions: no route at all, one actual detour, and several locally available routes competing for the same interior.

### Tournament witness

Use seven vertices with zero-based outgoing masks
$$
(82,48,107,67,108,9,34)
$$
and starts $(0,1)$, ends $(2,3)$. All four closing arcs from $\{2,3\}$ to $\{0,1\}$ exist. With interior $\{4,5\}$, in the induced six-vertex tournament the only full cover is $(0,4,2)$ and $(1,5,3)$. In the full seven-vertex tournament with interior $\{4,5,6\}$, the only full cover is $(0,6,5,3)$ and $(1,4,2)$. All four individual reachabilities are positive both times. With straight closure the first cover gives two directed triangles, the second one directed Hamilton cycle. With crossed closure the first gives a Hamilton cycle and the second gives cycles of lengths 4 and 3. No port count or reachability test can treat the two induced graphs as interchangeable for this transfer.

## Why the cycle parity needs a second channel

In the squarefree ring on all vertices, write $D_A=\det(I-AX)$, where $X$ is the diagonal vertex-variable matrix. The formal trace expansion gives
$$
-\log D_A=\sum_{C}x_{V(C)}w(C),
$$
over directed simple cycles counted up to cyclic rotation. Repeated vertices kill a walk, and the factor $1/k$ cancels the $k$ choices of a root on each length-$k$ cycle. This can be derived over the rationals; the stated cycle coefficients are integral. Taking the difference between $A$ and $-A$ then exponentiating yields
$$
\frac{\det(I+AX)}{\det(I-AX)}
=\prod_{\text{odd simple cycles }C}(1+2x_{V(C)}w(C)).
$$
Here each monomial in a single cycle squares to zero; cross products with intersecting supports vanish. In particular for a tournament, the inherited odd-cycle expansion for Hamilton paths gives the second exact expression
$$
H(T)=[x_{V(T)}]\left(\prod_{v\in V(T)}(1+x_v)\right)
\frac{\det(I+AX)}{\det(I-AX)}.
$$
The singleton factors supply unused vertices. This last equality with Hamilton paths is special to tournaments; the determinant ratio is an odd-cycle packing polynomial for arbitrary digraphs, not their Hamilton-path enumerator.

Under Schur elimination the determinant factors into the eliminated interior determinant and the retained port determinant:
$$
\det(I-AX)=\det(I-A_{II}X_I)\det(I-E_AX_P).
$$
The other sign channel uses $E_{-A}=-E_A(-x_I)$, since an internal path visiting $k$ interior vertices uses $k+1$ arcs. One can implement both channels using the one-site recurrence, accumulating each pivot factor $1-A_{ii}x_i$. Thus path reachability, complete support, and parity of route length are three distinct data. The identity gives an exact local elimination law for the signed cycle side of the tournament count, not a scalar continuation rule or a proof of the whole count spectrum.

The companion replay checks both literal witnesses by subset enumeration, the order-independent port recurrence, and the determinant-ratio Hamilton identity on every labelled tournament through order five.

## More than two ports: a cycle polynomial for alternating seams

For $k$ distinct starts and $k$ distinct ends, suppose all end-to-start closing arcs are available, and count only the stipulated closures of a disjoint internal path family. Let $P$ and $D$ denote the full-support coefficients of the permanent and determinant of its port matrix. For a given internal matching $\pi$, any outside matching $\sigma$ gives as many closed cycles as the permutation $\sigma\pi$ has cycles. As $\sigma$ varies freely, $\sigma\pi$ runs over $S_k$. Therefore the polynomial over all such seam matchings, recording resulting cycle count in $y$, is
$$
P\,y(y+1)\cdots(y+k-1).
$$
The polynomial weighted by the sign of the *exterior* matching is
$$
D\,y(y-1)\cdots(y-k+1).
$$
Indeed $\operatorname{sgn}(\sigma)=\operatorname{sgn}(\sigma\pi)\operatorname{sgn}(\pi)$, and $\operatorname{sgn}(\tau)=(-1)^{k-c(\tau)}$ for a permutation with $c(\tau)$ cycles. The rising factorial follows directly by inserting the new largest element in one of the old cycle slots or as a singleton. The falling factorial is its signed version. In particular the half-sum gives exact counts when the outside seam lies in the alternating group (literal $A_5$ for five starts). The coefficient of $y$ in the unsigned case gives $(k-1)!P$, and in the signed case $(-1)^{k-1}(k-1)!D$. This requires the full set of closing arcs as well as disjoint internal supports; it does not extend to forbidden exterior matchings by merely substituting their number.

One can retain parity of each resulting cycle, not merely the number of cycles. Fix an internal $k$-path cover with $p$ paths of odd vertex count and $q=k-p$ even. An exterior permutation cycle has odd total length exactly when it contains an odd number of odd sites. Write $Z_{p,q}(u,v)$ for the sum over all exterior permutations with weight $u$ per even cycle and $v$ per odd cycle. As above, composing with the fixed internal matching simply reindexes the exterior permutations. The two-sort exponential generating function is
$$
\sum_{p,q\ge0} Z_{p,q}(u,v)\frac{x^p z^q}{p!q!}
=(1-z-x)^{-(u+v)/2}(1-z+x)^{-(u-v)/2}.
$$
For a cycle of $\ell$ sites, using the formal colors $x,z$, its even/odd projection is the half-sum or half-difference of $(z+x)^\ell$ and $(z-x)^\ell$. A cycle contributes with factor $1/\ell$ in the labelled exponential generating function; summing those logarithms and exponentiating for the unordered set of cycles gives the identity.

In the odd-cycle-forest specialization $u=0,v=2$ the series reduces to
$$
\frac{1-z+x}{1-z-x}=1+\frac{2x}{1-z-x}.
$$
For $k>0$ its exact value is consequently $Z_{p,q}(0,2)=2p(k-1)!$, including zero if $p=0$. Thus allowing all seams converts the weighted odd-only closure problem into a **linear readout of the number of odd internal paths**; the weighted sum of full-support closures across all internal covers is $2(k-1)!$ times the sum of their odd-path counts. This is a conditional composition law on the specified seam family, not a way of discarding parity or the all-seams hypothesis.

The parity of the allowed exterior matching is also forced for every such odd-only closure. The number of output cycles has parity $p$, since each contains an odd number of odd sites. Thus $\operatorname{sgn}(\sigma\pi)=(-1)^{k-p}=(-1)^q$ and $\operatorname{sgn}(\sigma)=(-1)^q\operatorname{sgn}(\pi)$. In the alternating-seam channel, retain exactly the internal covers with $\operatorname{sgn}(\pi)=(-1)^q$, then apply the same weight $2p(k-1)!$; the others contribute zero. This is stronger than a naive split in half between even and odd seam matchings.
