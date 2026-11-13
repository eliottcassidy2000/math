# How the placement of marks controls every intermediate symmetry order

This is a follow-up to [the recursive neutral-flip atlas](WREATH_ATLAS.md), §5. Fix the monochromatic depth-$d$ tower formed by cyclically joining three copies at each internal node, starting from a rigid leaf tournament. A marked leaf is fixed individually, not merely as part of a set of marked leaves; one can equivalently name a vertex inside that rigid leaf. The marks do not alter the tournament or its full path-cover polynomial.

## Extremal budgets, with nested optimal placements

Let $E_d=(3^d-1)/2$. If $k$ of the $3^d$ leaves are marked and $a_j$ is the number of depth-$j$ nodes with a marked descendant, the pointwise stabilizer has order $3^{E_d-\sum_{j=0}^{d-1}a_j}$. Indeed the three-way rotation at an occupied node is blocked, and every entirely unmarked subtree retains its full rotation group. This identifies the group itself as the direct product of the rotation groups of the maximal entirely unmarked subtrees.

For every $0\le k\le3^d$, the smallest and largest possible exponents are
$
 e_{\min}(d,k)=E_d-\sum_{j=0}^{d-1}\min(k,3^j),\qquad
 e_{\max}(d,k)=E_d-\sum_{j=0}^{d-1}\left\lceil k/3^{d-j}\right\rceil.
$
A mark covers at most one node on each level; conversely an occupied depth-$j$ node can accommodate at most $3^{d-j}$ marks. Both bounds are attained simultaneously at every level, and in nested order as the budget increases. In ternary digit-reversal order (least significant digit first from the root), the first $k$ addresses visit exactly $\min(k,3^j)$ prefixes of length $j$. In ordinary lexicographic order, they fill contiguous subtrees of capacity $3^{d-j}$, visiting exactly $\lceil k/3^{d-j}\rceil$ prefixes of length $j$. These are the dispersion and concentration extremes of the same address system.

## There are no gaps between the two exponents

For $k>0$, an integer profile $(a_0,\ldots,a_d)$ comes from a marking if and only if (a_0=1, a_d=k), and $a_j\le a_{j+1}\le3a_j$ for each $j$. Necessity follows by counting occupied children. For sufficiency, give each occupied node one occupied child, then distribute $a_{j+1}-a_j$ extra children with capacity two per parent. Do this recursively and mark exactly the occupied leaves.

Write $b_j=\min(k,3^j)$, the componentwise maximal profile. Starting with the componentwise minimal profile $a_j=\lceil k/3^{d-j}\rceil$, increase one internal coordinate at a time until $a=b$. To see a move always exists, look at the deficient indices $j$ with $a_j<b_j$. An increase at $j$ can only be blocked to the left by $a_j=3a_{j-1}$ or to the right by $a_j=a_{j+1}$. In either case the blocking neighbor must itself be deficient, by the corresponding inequality for $b$. If every deficient coordinate were blocked, pointing each to one blocker would give a directed cycle on a path, necessarily a two-cycle. A two-cycle would require $a_{j+1}=a_j=3a_j$, impossible for a positive profile. Thus an increment is always possible. Each increment changes the stabilizer exponent by exactly one. The empty marking has only the full-group exponent. The profile path proves existence; it does not assert that a single leaf exchange realizes each step for arbitrary presentations.

At depth three, three marks in different top branches leave order $3^6$; three in one bottom triple leave $3^{10}$; all three intervening powers also occur, without changing either the tournament or its path-cover polynomial. At fixed budget, the full profile controls the order and the shapes of the unmarked pieces determine the product decomposition. This is symmetry inside a fixed Hamiltonian-count fiber, not a new count realization or a Chan–Pak coverage argument.

A finite control used exact recursive sets of attainable exponents for every budget through depth four, and direct marking masks through depth two. The all-depth result above is independent of these checks.
