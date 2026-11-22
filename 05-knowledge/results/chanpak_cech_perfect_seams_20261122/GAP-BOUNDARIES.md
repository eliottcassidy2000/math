# The two gaps, their shifted coordinates, and the arcs that cross them

## Why the gaps survive all orders

**CITED + PROVED alternative small-core reduction.** The existing [THM-1370 proof](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) and the [original bof–Royle discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament) already exclude 7 and 21. The new small-core explanation below uses the same strong-component floor, but replaces the full small tournament census by odd cycle packings and three marked four-vertex shapes.

Strong components concatenate in ordinal order and their counts multiply. Thus 7 would have to occur strongly, while 21 either occurs strongly or needs a factor 7. [Busch's exact strong minima](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) for orders 3 to 7 are 3,5,9,15,25 and then increase. The only strong forms of orders 3 and 4 have counts 3 and 5. This excludes 7 and leaves orders 5 and 6 for a strong 21. Write $c_j$ for the number of directed $j$-cycles counted up to rotation, and $D$ for the number of unordered disjoint pairs of triangles. At order at most six the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) reads $H=1+2(c_3+c_5)+4D$. Put $e=(H-1)/2$.

* On five vertices $c_3=10-\sum_v\binom{d(v)}2\le5$. For any fixed vertex of degree $d$, a directed 5-cycle specifies its successor (an outneighbor) and predecessor (an inneighbor). The last two vertices have only one possible order. Thus $c_5\le d(4-d)\le4$. So $e\le9$.
* In a strong tournament of order six whose degrees are 2 or 3 there are three of each. Then $c_3=8$, and every vertex deletion is strong since it has neither source nor sink (at order five two nonsingleton strong components are impossible). Thus $c_5\ge6$ by Camion. If there is degree 1 but not degree 4, the score is $(1,2,3,3,3,3)$. Then $c_3=7$, and all deletions except possibly that of the degree-1 vertex's sole outneighbor are strong, so $c_5\ge5$. The degree-4-only case is dual.

In the remaining case let $x$ have outdegree 1, $y$ indegree 1, and $Q=T-\{x,y\}$. If $x\to y$, then $y\to Q\to x$. With $c=c_3(Q)$ this gives $(c_3,c_5,D)=(4+c,4+2c,c)$. For example the 5-cycles are paths through core triples closed via $x,y$. Hence $e=8+5c\ne10$, with $c\in\{0,1,2\}$.

Otherwise $y\to x$. Mark $a$ the sole outneighbor of $x$, $b$ the sole inneighbor of $y$. With core outdegree $d$, spanning path start/end counts $S_a,F_b$ in (Q), and $R_{ab}$ counting $a\to z\to b$ there, the exact reads are
$$
 c_3=c+d(a)+3-d(b)+[a=b],\qquad
 c_5=S_a+F_b+R_{ab}+2[a\to b]+3[a=b].
$$
The three terms for $c_3-c$ are triangles with $x$ alone, $y$ alone, or both. Cycles with just $x$ or $y$ use a core spanning path. Those with both $x,y$ either contain $y\to x$ (giving the two-step core path), or use the two marks and a third core vertex; the latter give two choices of positions when $a\to b$ and three choices of core vertex when $a=b$. For $D$, one term is $[a=b, Q-a\text{ is cyclic}]$; add the number of ordered choices $u,v$ with ${a,b,u,v}=Q$ and $a\to u, v\to b$, selecting separate triangles through $x,y$. This reduces the exclusion to a $4\times4$ marked table; core dual does not change the set of outcomes.

| (Q), up to dual | $e$ for $c_3=4$ | $e$ for $c_3=5$ | $e$ for $c_3=6$ |
|---|---|---|---|
| transitive | 7,8 | 11 | 16 |
| source then triangle | 7,8 | 9,13 | 15 |
| strong with arcs 01,02,12,23,30,31 | 7,8 | 9,11 | 11,12,14 |

There are no $c_3=3$ or $c_3\ge7$ entries in these tables, and entries $c_3\le2$ have $e\le4$. In particular no value is 10. This yields the two distinct mechanisms: the strong floor exhausts the small admissible orders, the marked core misses 21 internally, and multiplication cannot repair either hole.

## All orders of single-arc crossings over each exceptional gap

**PROVED from a finite-exact primitive packet.** The number of Hamilton paths shared by both orientations is denoted $u$ in the generating functions by its exponent. For the ordered lower-to-upper bridges, counted up to simultaneous isomorphism, the exact two-variable series are
$$
  \frac{2z^5(u^3+u^4)}{(1-z)^2}\quad(5\to9),\qquad
  \frac{2z^6(u^{11}+u^{14}+u^{17})}{(1-z)^2}\quad(19\to23).
$$
The [primitive seam splitting](../chanpak_support_birth_grammar_20261122/PRIMITIVE-SEAMS.md) works for any different counts. The primitive pair has a strong endpoint; the common exterior here is transitive by coprimality. Busch gives order at most 5 for the first primitive pair and at most 6 for the second. The upper endpoint cannot occur at smaller orders: on order 4 the maximum is 5, and the 5-vertex bound above gives $H\le19$. So only the two indicated orders need checking.

The path-bitset replay [check_bridges.py](check_bridges.py) enumerates all labeled orders 5 and 6 and checks the literal unchanged paths, orienting lexicographic edge bit 1 as $i\to j$ for $i<j$. It finds 240 labeled pairs per overlap 3 or 4, and 1440 per overlap 11, 14 or 17. Simultaneous automorphisms act freely on the paths of each host (a fixed ordered path fixes the vertices), so their order divides both coprime counts and is one. Divide by (5!) and (6!) respectively. Adding a common transitive prefix/suffix preserves the overlap exactly and gives the two geometric factors. Across the larger hole, the same lower mask 32747 supplies all three overlaps: flip $(3,4),(1,3),(0,4)$ for 11,14,17, with losses/gains $(8,12),(5,9),(2,6)$. So even at fixed endpoints the route across the hole is not determined by scalar heights.

## The coordinates and their arithmetic bridge

Under $q=(H-1)/2$ the two gaps are ${3,10}$; under $m=(H+1)/2$ they are ${4,11}$. Their distance is seven, and reflection through $q=13/2$ becomes reflection through $m=15/2$ by translation. [Session332](../session332_hidden_digits_and_arithmetic_diamonds_20261114/ODD_COORDINATE_PREDECESSORS.md) shows ${4,11}$ is literally the greedy inverse $3x-1$ cycle (exponents 3,0), and that its induced map on the $H$-line sends attained counts 15,31,63 to the forbidden 21. So it cannot preserve the tournament spectrum. The [base-packing comparison](../packing_bases_repunits_20261108/README.md) and [CRT/reflection comparison](../session297_parity_cutoff_20261113/BRIDGES.md) are additional arithmetic shadows, not explanations by transfer.

There is now a physical reason for these numbers to matter to an exponent induction: they are actual boundary values in the one-slot specialization $p-1=1$, which forbids exact uniform charts of that type regardless of how late their intended exponent range starts. See [the one-port proof](PORT-MODES.md). The independent norm test explains why divisible-by-7 upper heights such as 497 can still exist, and why keeping two independent growing factors will not provide their uniform construction. See [the dyadic irreducibility proof](NORM-AND-PORTS.md).
