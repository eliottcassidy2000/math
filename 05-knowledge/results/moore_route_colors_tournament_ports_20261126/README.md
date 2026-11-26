# What colors retain: Moore route words, tournament path ports, and A000568

2026-11-26. Continues the [six-partition fit](../moore57_singer_repunit_six_partitions_20261126/README.md) and the [independent tournament arc responses](../independent_responses_and_two_holes_20261126/README.md). **PROVED**: sharp Moore edit support and its rigidity at degree 57, port-color Latin criterion, equitable-color quotient constraints, and the two exact enumeration transforms. **FINITE-EXACT**: quotient profile enumeration, counts through tournament order 7, and the doubly regular and Hoffman–Singleton controls. Chan–Pak exact coverage and degree-57 existence are not claimed.

## 1. The shared obstacle is losing the location of an interaction

The Moore equation $A^2+A-(d-1)I=J$ sums over middle vertices. A nonedge has exactly one middle vertex, so it has a canonical **route color**, its middle vertex. For a chosen partition of centers let $D_c$ be its diagonal indicators. The matrices $AD_cA$ count two-edge routes of each color; they are off-diagonal 0/1 and have disjoint off-diagonal support. Add one color for direct edges. Equivalently, the edges of the complete graph are partitioned into the original edges and the cliques on the $d$ neighbors of every Moore vertex. This is an edge coloring of unordered *pairs*, not a proper edge coloring of the Moore graph.

A published matrix-multiplication result does not by itself address the inverse problem here. [OpenAI's November 24 paper and Lean scope](https://github.com/openai/math/tree/main/preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026) report an exact exponent $9/4$ for formal square multiplication, with the linked formalization covering division-free circuits. The [independent upper-bound exposition](https://arxiv.org/abs/2611.09716) treats the upper half. None of the bounds below assumes that result. On a proposed Moore graph, marking all endpoints of nonbacktracking length-two routes checks the equation directly in $O(nd^2)=O(n^2)$, since $n=d^2+1$. This sparse and structured verifier is faster asymptotically than general dense multiplication at exponent $9/4$. Filling an unknown adjacency matrix is not an application of that verifier.

For tournaments the comparable data must remember *visited sets*, not just middle vertices. The walk-to-path lift is to compute with $x_v^2=0$: for adjacency $M$ and $D=\operatorname{diag}(x_v)$ the top squarefree coefficient of $\boldsymbol{x}^T(MD)^{n-1}\mathbf1$ counts Hamiltonian paths. Arithmetic contraction is possible, but erasing these tags turns simple paths back into walks. This also explains the need for a common-set index in a marked arc response: for ordered endpoints $u,v$, the number of paths through a newly admitted $u\to v$ is the sum over splits $S\ni u,v\notin S$ of the number ending at $u$ on $S$ times the number beginning at $v$ on its complement. Reverse response subtracts the opposite split count. Two disjoint arc edits can still appear in the same Hamiltonian path, so their mixed response need not vanish. For the Moore matrix defect, by contrast, edits on disjoint vertex sets have zero mixed term in $\Delta^2$. This is why the old four-orientation independent-response tests have to be measured on tournaments, not inferred from disjoint endpoints.

### A sharp local repair theorem for Moore graphs

[Full proof](TRADE.md): two distinct completed degree-$d$ Moore graphs on one labeling must differ at **at least $2d$ vertices** (a vertex is counted when its neighborhood changes). A label transposition attains equality. For $d\ge4$ and $6\nmid d$, equality can only be a label transposition. In particular a degree-57 search cannot move between genuinely nonisomorphic completions by changing the neighborhoods of 114 or fewer vertices.

The proof colors outside guards by how many edited vertices they see. A one-contact guard is impossible, a two-contact guard forces opposite edit rows, and a larger guard consumes a unique-route budget. The equality case with no one-neighbor edit would force all components of an auxiliary 2-regular quotient to be six-cycles. This is a lower bound on *transitions between completions*, not a bar on finishing a partial board.

### A hostile tournament control

For a doubly regular tournament of order $4m+3$, degree $2m+1$ and adjacency $M$, $MM^T=(m+1)I+mJ$ and $M+M^T=J-I$, hence $M^2+M=(m+1)(J-I)$. This is a close directed relative of the Moore quadratic law.

Using [McKay's published doubly regular tournament data](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html) and exact subset recurrences:

- At order **19** the two inputs obey the same matrix law, have the same quadratic adjacency algebra and the same total walk counts, and have **1,172,695,746,915** and **1,172,732,961,771** Hamiltonian paths. Difference: **37,214,856**. The common law is a *parameter law*, not equality of their adjacency matrices. It determines traces of adjacency powers since the equations reduce powers to the span of $I,M,J$, whose traces are fixed.
- At order **15** the two inputs have the same total **198,335,025**, the same multiset of vertex-deletion totals, and the same multiset of all 105 unmarked single-arc responses. Yet their per-terminal path-count multisets differ. Thus even a scalar total plus two plausible local refinements can lose the boundary information needed to join gadgets.

## 2. Three different notions of Moore coloring

**Ports on edges.** A proper $d$-edge coloring of a $d$-regular graph on an even order is a one-factorization. Write its fixed-point-free involutive permutations $P_1,\ldots,P_d$. The order-$d^2+1$ Moore condition is then equivalent to a Latin-square test: the rows indexed by the identity, the $P_i$, and the $P_iP_j$ with $i\ne j$ are permutations on the vertices, and in each column their values visit every vertex once. This is an exact factorization
\[
 I+\sum_iP_i+\sum_{i\ne j}P_iP_j=J
\]
in permutation-matrix notation. Conversely the disjoint-support identity of these word matrices gives degree $d$, no length-two collisions and diameter two. An explicit seven-color certificate gives all fifty rows on Hoffman–Singleton. The Petersen replay checks all six perfect matchings and each leaves two 5-cycles, so three ports are impossible there. This is an optional, stronger encoding, not a claim that a degree-57 graph must be one-factorable.

**Colors from a root.** The $d^2-d$ vertices at distance two split into $d$ branches of size $d-1$. For each pair of branches their edges form a perfect matching (the unique two-step route from a leaf to the opposite branch parent). Group the branch pairs using the round-robin coloring of $K_d$ for odd $d$, assigning the same color to every edge in a branch-pair matching. This properly colors the leaves, but each leaf in branch $i$ misses the *same* color, namely $i$. In any extension to 57 colors on the whole degree-57 graph, at least 55 of the 56 radial edges in each branch would need a different color, which must be freed by changing an old leaf edge at that leaf. Each changed leaf edge serves at most two leaves. Thus this template requires at least **1,568** leaf-edge recolorings, or **1,539** even with a 58th color. It is a warning against using a root-symmetric coloring as if it already provided global ports.

**Equitable vertex colors.** Let $s_i$ be the class sizes and $B_{ij}$ the number of neighbors in color $j$ of every vertex of color $i$. For a genuinely equitable partition, $B^2+B-(d-1)I=\mathbf1 s^T$ and $s_iB_{ij}=s_jB_{ji}$. For a proper equitable coloring the diagonal of $B$ is zero. Degree 57 has restricted eigenvalues 7 and −8. If a proper equitable quotient has $k$ colors and 7 has multiplicity $a$, then $15a=8k-65$. Hence **$k\equiv10\pmod{15}$**: an ordinary nine-coloring is not ruled out, but an equitable nine-coloring is. (The ordinary spectral lower bound is nine.)

For $k=10$, write the symmetric quotient as $-8I+65pp^T+15vv^T$ with $p_i=\sqrt{s_i/3250}$. Its diagonal gives $v_i^2=(400-s_i)/750$. For all $s_i<400$, the squarefree parts of $s_i(400-s_i)$ must agree: write $s_i(400-s_i)=Dq_i^2$ and choose signs for the integers $q_i$ so that $\sum q_i=0$. To see the common squarefree part, the off-diagonal numbers
\[
 B_{ij}=\frac{s_is_j+Dq_iq_j}{50s_i}
\]
must be integral, so the unsigned product under the square root for any two nonzero coordinates must be a square. Include $s_i=400$ as $q_i=0$. Enumerating multisets of ten signed classes under the equation, pair integrality in both directions, $1\le s_i\le400$, and $\sum s_i=3250$ gives only four candidate quotients up to global sign and reindexing:

| Color class sizes | Number of integer quotients |
|---|---:|
| $250^4,375^6$ | 2 |
| $325^{10}$ | 1 |
| $1,57,399^8$ | 1 |

All four quotient matrices satisfy the complete quadratic equation; none is asserted to lift. The last is exactly the natural *root + first layer + eight independent leaf colors* ansatz. An eight-coloring of the leaf graph would have all classes of size 399 by the least-eigenvalue bound ($3192\cdot8/(56+8)=399$); spectral equality gives eight neighbors of each other leaf color at every leaf. Adjoining the root to a leaf class saturates the full graph bound of 400, giving eight neighbors outside that enlarged class. Each first-layer site sees the root already, so it has seven neighbors in every leaf color. Conversely the candidate with its singleton class has precisely these layers and implies an eight-coloring of the leaf graph.

## 3. A000568 vs. the Hamilton spectrum: what to retain in enumeration

[Davis's cycle-index count](https://link.springer.com/article/10.1007/s10801-022-01197-0) groups label permutations by cycle type. A fixed tournament cannot have an even vertex cycle because its antipodal pair comes back reversed. For a type with all odd parts $\ell_i$, the number of unordered pair orbits is $\sum_i(\ell_i-1)/2+\sum_{i<j}\gcd(\ell_i,\ell_j)$, and every such orbit has two orientation choices. This reduces A000568 to a sum over partitions. The weighted Burnside sum retaining $z^{H(T)}$ still exists, but counting two choices per orbit no longer works: the choices influence a single global path count. The same problem recurs when fifteen two-column Moore boards are each correct but must arise from common six-tuples.

For self-converse counts use the twisted Burnside average. In an unordered pair orbit of length $\ell$, record whether the ordered representative returns swapped. Under global converse followed by the permutation, it must return swapped exactly when $\ell$ is odd. Each compatible orbit again has two orientation choices. Averaging ordinary and twisted counts gives converse-merged counts, retaining the distinction between one fixed class and a pair of genuine classes.

An exact multiplicative refinement avoids saving all labeled tournaments. Let $Z(x,s)=\sum_{n,[T]}x^nH(T)^{-s}$, and $S(x,s)$ be its nonempty strongly connected part. Every tournament is a unique ordered join of strongly connected components, and every Hamiltonian path visits them in order, so
\[
 Z(x,s)=(1-S(x,s))^{-1}.
\]
If $C(x,s)$ counts only self-converse strong tournaments, the self-converse series is
\[
 Z_{\rm sc}(x,s)=(1+C(x,s))/(1-S(x^2,2s)).
\]
Indeed opposite SCC positions must be converse classes, giving squared Hamilton counts per opposing pair, and there is zero or one self-converse central strong component. Converse-merged series is their half-sum. Treat products of weights formally as $[a][b]=[ab]$ if preferred; this is exact by count, not just analytic moments.

The marking replay uses one fixed Hamiltonian path. Two tilings code the same tournament precisely when renormalizing one along one of its other Hamiltonian paths produces the other. The fiber size is $H/|\operatorname{Aut}|$ (the vertex permutation fixing a path positionwise is identity). This yields a direct enumeration check against Burnside and the strong/converse transforms. At order seven: **456** classes, **353** strong classes, **88** self-converse classes, **272** converse-merged classes, and **77** distinct Hamilton counts. The heaviest Hamilton fiber is 91, with fifteen classes. At order ten the lightweight Burnside counts already give 9,733,056 classes and 8,784 self-converse classes without computing their Hamilton spectrum.

Thus the useful transfer is not that either problem counts the same objects. Both need a refinement that keeps *addresses of collisions*: edit-row signs and port words on the Moore side, endpoint and used-set data for tournament seams. The trade bound says a locally exact Moore repair can still be globally forbidden. The doubly regular controls say why a uniform quadratic law cannot be used as a surrogate for a Chan–Pak selector.

## Replay

Run **python3 census.py**, **python3 drt.py**, **python3 ports.py**, and **python3 q10.py**. For the 19-vertex comparison compile **g++ -O3 hist19.cpp -o hist19** and run **./hist19 < drt19.txt**. The text records the primary data source; the recurrences independently verify the matrix law and path totals.
