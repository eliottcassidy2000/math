# Twenty-four odd primes, two boundary contacts, and the six-cycle prism

2026-11-26. **PROVED**: prime boundary repair mechanism, the Moore hexagon shell identities, the order-54 prism sign obstruction. **FINITE-EXACT**: the prime windows, Hoffman–Singleton control, and the pinned degree-57 chart. **CONDITIONAL**: all statements about a completed degree-57 graph. No existence or nonexistence claim.

Inherits the [square-sum arrival and endpoint kernel](../session32_square_monster_20261023/square-sum/README.md), [Monster level filtration](../monster_levels_square_sum_7_21_20261022.md), [first-25 equal-fifth rails](../first25_spread_target_returns_20261116/QUINTILES.md), [first-25 prime-value matching test](../prime25_square_rails_moore_20261126/README.md), and [coupled Latin matrices](../moore57_coupled_latin_matrices_20261126/README.md). The live distinction is between a permutation defect, a spanning-path join, a square reflection, and a missing common-neighbor slot. They are not interchangeable.

## 1. Two ways to get 25, one way to repair 26

Let $P$ be the first 24 odd primes, from 3 through 97. Use **prime values** as labels, not ranks. Square-sum edges join distinct vertices $p,q$ with $p+q$ a square. An odd–odd edge goes between 1 and 3 mod 4; all neighbors of 2 are 3 mod 4, so put 2 on the 1 side. This bipartition works for all primes, not just the window. The sides at orders 24 and 26 have sizes $11+13$ and $13+13$.

| Window | Square edges | Isolates | Minimum nonsquare entries in a derangement | Optimal derangements | Minimum nonsquare edges in a spanning path |
|---|--:|---|--:|--:|--:|
| $P$ | 20 | 37, 43, 67, 79 | 4 | 9 | 6 |
| $P+2$ | 24 | 37, 43, 67 | 3 | 2 | 5 |
| $P+101$ | 21 | 37, 67, 79 | 3 | 2 | 6 |
| $P+2+101$ | 25 | 37, 67 | 2 | 1 | 5 |

Each row in the first table has a unique square-supported permutation on its nonisolates. At 24 this is already forced by peeling leaves (or by the exact permanent replay). At 25 and 26 the new pivots are forced by the leaves 79 and 101: $2\leftrightarrow79$ and $43\leftrightarrow101$. In particular, the other three contacts from 2, to 7, 23 and 47, add options to paths but cannot appear in an optimal permutation. At optimum the isolates must permute **among themselves** without fixed points. At 26 the sole nonsquare transposition is $37\leftrightarrow67$ (sum 104). This proves the two boundary repairs commute on the optimal permutation, and that the new forced pairs did not unlock any old choices.

The path objective is different: the minimum joins are the minimum number of paths in a square-only spanning linear forest minus one. Deleting nonsquare edges cuts a path into such a forest, and concatenating optimal forest paths attains the bound. In the 24 case the 16-vertex component has three leaves and needs at least two paths; in the 26 case the 22-vertex component has five leaves and needs at least three. The receipt gives optimal covers for both and the intermediate windows. Adding 101 merely attaches one isolated component to the new leaf; it does not lower the path-cover count. This is unlike the archived consecutive-integer transition at 25: the Pell reflection frame on $[1,25]$ and its first Hamiltonian-path onset after 23 are not preserved by rank-to-prime relabeling. In fact the first square-sum 4-cycle through 2 arrives at 167 in cyclic order 2–7–29–167–2; the first all-odd one arrives at 193, in cyclic order 3–13–131–193–3. The minimality tests check common-neighbor intersections, not cycle samples.

### A Golay guard with actual repair coordinates

The four isolated prime coordinates $D_0=\{37,43,67,79\}$ play a concrete role when they label four coordinates of the extended binary Golay code. From the [octad Steiner property of Conway–Sloane's Golay construction](https://www2.kenyon.edu/Depts/Math/Aydin/Teach/Sp09/328/Lexicodes.pdf), exactly five octads contain $D_0$; their extra tetrads partition the other twenty coordinates, since each extra point specifies a unique fifth point. After the two repairs, the wanted toggle is $\{37,67\}$ and the guards are $\{43,79\}$. There are exactly

    77 - 21 - 21 + 5 = 40

octads toggling the pair while avoiding both guards. These are the minimum-weight nontrivial codewords for this constraint (weight 8, not 2). If one literally wants a length-26 code from these 24 labels, the explicit boundary map $c\mapsto(c,c_{79},c_{43})$ attaches the copied coordinates 2 and 101 along their forced partners. The same 40 minimum masks leave all four guarded coordinates unchanged. This uses the code as a selector, **not** as a prime-graph automorphism or a proof that the resulting tournament moves preserve a path count.

## 2. A hexagon gives a canonical 24 in Hoffman–Singleton

Let $G$ be a degree-$k$ diameter-two Moore graph with $k\ge3$, and $C=(c_0,\ldots,c_5)$ any hexagon. Every hexagon is induced by girth five. Let $V_i$ be the vertices contacting only $c_i$, $D$ the three common neighbors of opposite pairs, and $U$ the vertices with no edge to $C$. No outside vertex can contact a distance-one or distance-two pair of cycle vertices, by girth. Thus

```
|D|=3,  |V_i|=k-3,  |V|=6k-18,  |U|=k²-6k+10.
V = disjoint union of all V_i.
```

Write $W=V\sqcup D$, its induced adjacency $A_W$, the cycle adjacency $A_C$, the $6\times|W|$ matrix $S$ of contacts from $C$ to $W$, and the matrix $B$ of contacts from $U$ to $W$. The Moore identity $A^2+A=(k-1)I+J$ gives

```
A_W Sᵀ = J-Sᵀ-Sᵀ A_C,       BSᵀ = J.
```

For a member of $V_i$ the first identity says it has exactly one neighbor in each of the three $V_j$ where $j\notin\{i-1,i,i+1\}$, and none in $D$. Indeed a vertex of $D$ has two opposite parents on the cycle, and the two parent columns and four adjacent columns exhaust the right side of its equation, leaving no neighbor in $W$. Thus $V$ is a $(k-3)$-sheet cover of the triangular prism (the complement of $C_6$); $D$ is isolated in $W$. The prism has nine base edges, so tree gauge leaves just four sheet permutations. Its two triangular and three square loop permutations must be derangements by girth.

For **Hoffman–Singleton** this shell is a literal 24, partitioned as six fibers of four; adjoining one or two opposite mediators gives 25 or 26 vertices with the new sites isolated in $W$. This is a structural 24 next to the prime 24, not an identification. In all 5,250 hexagons of the [pentagon/pentagram construction](https://www.math.auckland.ac.nz/~hafner/his.pdf), each triangular loop is one 4-cycle on the sheets, and each square loop is two 2-cycles. In **degree 57** there are 54 sheets. The product of the signs of the three square loop permutations equals the product of the signs of the two triangle permutations: in the product the vertical edges occur twice, and the six horizontal triangle edges once. If all square loop permutations were 27 transpositions, the two triangular signs would have to be **opposite**. Therefore a model that preserves equivalent triangular loops (e.g. makes them conjugate or inverse) must give at least one square loop a cycle of length three or more. Also 54 is not a multiple of four, so no triangle loop is a disjoint union entirely of the Hoffman–Singleton 4-cycles. These are conditional necessary tests for copied prism patterns, not contradictions for degree 57.

### The 55-by-55 boards and their canonical holes

Each column of $S$ is either a single parent or an opposite pair. The identity $BS^T=J$ gives **six actual partitions of $U$** by the cliques of pairs with a shared external neighbor. A partition contains $k-3=54$ blocks of size $k-4=53$ and one distinguished block of size $k-2=55$ in degree 57. Two different blocks (different external vertices) intersect at most once.

- Opposite parents share their distinguished block. On its complement, the remaining $54\times54$ intersection table has exactly one hole per row and column. These holes are the direct matching $V_i\leftrightarrow V_{i+3}$.
- Nonopposite parents have different distinguished blocks. The table is $55\times55$ and binary, with its distinguished row and column entirely ones. Delete them: every ordinary row and column of the $54\times54$ table has exactly two holes. For neighboring parents on the hexagon the holes are the two 2-step routes in the prism; for parents distance two they are the direct edge and the 2-step route through the third triangle class.

Proof of the hole rules: vertices of different $V_i$ cannot share a parent. If they are adjacent, or their unique common neighbor lies in $V$, the two blocks are disjoint; otherwise the common neighbor lies in $U$ and the blocks intersect once. The prism routes therefore specify the complements of **all fifteen** intersection boards, not just their row totals. This is an alternative construction frontier to the prior 56-by-56 edge-deletion grid: 3 one-hole boards and 12 two-hole boards, joined by just four gauge-dependent primitive permutations for the shell. Completing these boards jointly on actual contacts remains open.

## 3. What the square-prime obstruction measures

For $u\in U$, let $\delta(u)$ be its number of neighbors in $D$. From the Moore identity, summing the cycle-contact counts of all its neighbors gives six, so $d_U(u)=k-6+\delta(u)$. Each vertex in $D$ has $k-2$ neighbors in $U$, and each pair of vertices in $D$ has a unique shared neighbor in $U$. The three pairwise shared neighbors are either all equal (**concurrent**) or all distinct (**separated**). In the first case the hexagon plus $D$ and the common apex is an induced Petersen graph. In the second case no Petersen graph contains that hexagon: its opposite pairs already have uniquely designated mediators. In fact there is then another induced hexagon on $D$ and the three distinct pairwise shared neighbors. Distinct mediators share no edge; any chord among the three new neighbors would be a triangle through their shared mediator.

In degree 57 the degree multiset on $U$ must be $51^{2754},52^{162},54^1$ when concurrent, and $51^{2755},52^{159},53^3$ when separated. This comes from counting the three 55-element neighbor sets of the mediators and their single pairwise overlaps.

In Hoffman–Singleton, the unique Petersen extension per hexagon is known (see [Hafner's Petersen-subgraph construction](https://www.math.auckland.ac.nz/~hafner/his.pdf)); the exact echo here is the graph induced by its antineighborhood $U$: in all 5,250 cases it is **two disjoint hexagons and a four-leaf star**. The 11-prime induced subgraph consisting of

    hexagon: 3–61–83–17–47–97–3
    fork: edges 5–31, 5–59, 5–11, 11–89

with no contacts between the components therefore cannot embed **induced** in Hoffman–Singleton. Ignoring nonedges loses that obstruction: the receipt gives a 26-prime edge-preserving embedding there which adds 24 extra edges.

On the actual [fourteen-site degree-57 fixture](../moore57_quadratic_tensors_20261125/planted_stars.py), 291,536 root hexagons are present: **1,770 concurrent** (295 distinct rooted Petersen graphs), **27,168 separated**, and 262,598 not decided by the two routes tested. To compute this, a rooted hexagon using center branch $a$ and other branches $b,c$ at sheet $x$ can close into Petersen exactly when the two opposite triangle transports return to the same $a$-sheet. In a completely matched triple these are $F=P_{ba}\circ P_{cb}\circ P_{ac}$ and its inverse, so the local test is $F^2$x$=x$, and girth excludes $F$x$=x$. The three full reference branches give 168 concurrent outcomes, despite 27,168 already separated elsewhere. Thus retaining only the full reference stars hides real non-Hoffman–Singleton behavior.

The replay explicitly forces the **induced** 11-prime obstacle in the pinned chart. It maps the displayed cycle to (0, 2, 115, 58, 1654, 29), and (5, 31, 59, 11, 89) to (60, 1542, 2019, 3196, 1124). The opposite mediators are (1,118,1651), and the witnesses 61 and 111 for the first mediator with the other two are distinct, so this hexagon cannot become Petersen. All required nonedges follow from saturated neighborhoods or a pre-existing path of length two or three, which would make a forbidden triangle or quadrangle if closed. This rules out a **uniform** Hoffman–Singleton-style extension in this fixture, not a Hoffman–Singleton subgraph somewhere in a hypothetical completion.

## 4. An exact rank test for the 2,917-vertex antineighborhood

Let $X$ be its adjacency. From the Moore identity,

```
X²+X = (k-1)I+J-BBᵀ,         XB = J-B-BA_W.
```

At degree 57 the only column relations of $B$ are the five from $BS^T=J$. Here is a rank proof. The identity
$B^{\mathsf T}B=(k-1)I+J-A_W-A_W^2-S^{\mathsf T}S$ is positive on $\ker S$ with lower bound $k-13$, since $A_W$ is cubic on $V$ and null on $D$. The Gram matrix preserves the span of the two type-constant vectors $\mathbf1_V,\mathbf1_D$: each term preserves it, using $S\mathbf1_V=(k-3)\mathbf1$ and $S\mathbf1_D=\mathbf1$. Its restriction is nonsingular there since the images $6\mathbf1-2\delta$ and $\delta$ are independent. In degree 57, $\delta(u)$ is sometimes positive and sometimes zero, by its sum 165 across 2,917 vertices, where $\delta$ denotes the column vector of all $\delta(u)$. A vector in $\ker B$ is therefore orthogonal to both constants. Its orthogonal projection $S^T u$ onto the row space of $S$ has $\sum u=0$ since $S^T\mathbf1$ is type-constant and $SS^T\mathbf1=(k-1)\mathbf1$. Then $BS^T u=0$, and subtracting this projection puts the remaining null vector in $\ker S$, where positivity kills it. This gives rank $327-5=322$.

The orthogonal complement of $\operatorname{im}B$ has dimension 2,595 and contains only eigenvalues 7 and −8 of $X$, since $\mathbf1\in\operatorname{im}B$. On the column quotient, $A_W$ has trace 3: it acts on $\ker B$ as $-I-A_C$ on the sum-zero cycle subspace, whose trace is −3, and $\operatorname{tr}A_W=0$. The rank-one term in $XB$ has quotient trace $k-2$, obtained by writing it as $B S^T e_i\mathbf1_W^T$. Thus the trace of $X$ on $\operatorname{im}B$ is $55-322-3=-270$, and the 2,595-dimensional complement has spectrum $7^{1402},(-8)^{1193}$.

There is also a metric consequence. Two distinct vertices of $U$ are distance three in $U$ exactly when their unique common neighbor is outside $U$. Such a pair has $k-1$ three-edge paths in $G$, and the first steps and the last steps are separately distinct; at most six paths are lost per end outside $U$. Therefore at least 44 three-edge paths remain for $k=57$. Distance-three pairs are precisely the edge-disjoint union of the cliques supported by columns of $B$, so there are $324\binom{53}{2}+3\binom{55}{2}=450,927$ of them. This supplies both a diameter-three certificate and five exact dependencies for a reduced search that does not forget the external common-neighbor sites.

Replay (standard Python, about one second), from this directory:

```sh
python3 prime_contacts.py ../moore57_quadratic_tensors_20261125
```
