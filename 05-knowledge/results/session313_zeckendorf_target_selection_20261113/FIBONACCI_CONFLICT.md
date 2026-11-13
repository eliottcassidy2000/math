# Where a Fibonacci conflict path stops being a tournament

The May [Zeckendorf/simplex/cuboid note](../../../07-reflections/eureka-zeckendorf-simplex-cuboid.md) and [THM-1425](../../../01-canon/theorems/THM-1425-the-ocf-is-the-weight-two-shifted-pascal-identity.md) identify the right symbolic recurrence: counting independent sets of a path gives Fibonacci numbers, and giving each chosen vertex two types gives Jacobsthal numbers. The assertion about odd-cycle conflict paths in THM-1425 is conditional. Here is an exact physical boundary of that condition, strengthening the old [full-path exclusion at four](../../../01-canon/theorems/THM-202-p4-omega-impossibility.md).

Let $\Omega(T)$ have one vertex per directed odd cycle, edges for intersecting vertex sets. In every nonsingleton strong component $S$ of $T$ there is an odd cycle of length $|S|$ or $|S|-1$, by [Moon's cycle theorem](https://doi.org/10.4153/CMB-1966-038-7). That cycle meets every odd cycle in $S$. Thus each nonsingleton SCC supplies one connected component of $\Omega$ with a universal vertex; arcs between SCCs cannot participate in a cycle.

If all directed triangles in a tournament are pairwise disjoint, each triangle is a module. An exterior vertex with mixed orientation to it would give a second triangle sharing vertices: around a directed triangle find a transition from an out-neighbor of the exterior vertex to an in-neighbor. Contract those disjoint modules. Any triangle in the quotient would lift to another directed triangle, so the quotient is transitive. In particular a strong tournament on at least four vertices must have two intersecting triangles.

Now in a strong component of size at least five the universal odd cycle has length at least five, and the two intersecting triangles are two other vertices of $\Omega$. They give a conflict triangle. For order three, $\Omega=K_1$. For order four there are exactly two directed triangles: all scores are 1 or 2 with two of each (strongness excludes scores 0 and 3), and the transitive triangles count is the sum of $\binom{\text{outdegree}}2$, namely 2. The two directed triangles intersect, so $\Omega=K_2$.

**Consequences.** A tournament's full odd-cycle conflict graph is triangle-free iff all SCCs have size at most four. Its only possible connected components in this case are $K_1,K_2$, and every disjoint union of those graphs occurs by ordinal joining directed triangles and strong four-vertex tournaments. In particular a full conflict path $P_m$ exists only for $m=1,2$ (with the empty case $m=0$ transitive). The actual counts in the triangle-free subclass are exactly $3^a5^b$ with $a,b\ge0$, by the odd-cycle expansion or by concatenating Hamiltonian paths in the SCCs. This does **not** exclude $J_{m+2}$ as a path count through another conflict graph. The smallest such warning is 11: the one-apex word $1100$ attains 11 with a different conflict graph (four triangles and a directed pentagon, all through the apex). The graph $P_3$ and that graph have the same independence-polynomial value at 2 but are not the same graph.

## Two typed digits recover a gapless rank language, not a tournament

Take $J_i=(2^i-(-1)^i)/3$. In a length-$n$ binary path allow nonadjacent selected sites, and give each selected site type A or B. Read the highest site first in order $0,A,B$. The counts satisfy $N_0=N_{-1}=1$, $N_n=N_{n-1}+2N_{n-2}$, so $N_n=J_{n+2}$. The first branch consumes one site, the other two consume two (the following site is forbidden). Their rank offsets are
$$
0,\quad J_{n+1},\quad J_{n+1}+J_n=2^n .
$$
Hence every integer in $[0,J_{n+2}-1]$ has a unique recursive typed decoder. This is the exact analogue of choosing a Zeckendorf branch by its interval length, with the fugacity changed to two. The untyped language has growth $\varphi^n$, not $2^n$. The types here label choices in a packing language, not the historical diagonal stripe or a choice of actual tournament with that count.

A different arithmetic bypass allows every finite binary subset of the weights $J_1,J_2,\ldots$ (including both units). Map a nonnegative integer $u$ to the sum $W(u)$ using its binary digits. Incrementing $u$ with $t=v_2(u+1)$ clears the first $t$ bits and sets bit $t+1$, so
$$
W(u+1)-W(u)=J_{t+1}-\sum_{i=1}^tJ_i
=\begin{cases}1&t\text{ even},\\0&t\text{ odd}.\end{cases}
$$
Thus $W$ is unbounded and skips no integer; a plateau is isolated since $t>0$ forces $u$ odd. There is a monotone inverse by bisection with at most two binary witnesses. The proof is arithmetic, and an additive tournament compiler for those weights is still missing.

The conflict obstruction rules out a tempting shortcut: one cannot realize this whole gapless language just by making the **full** cycle-conflict graph an arbitrarily long path. A physical compiler would need a retained path of selected interactions plus explicit control of the extra odd cycles, or a different type of transfer entirely.
