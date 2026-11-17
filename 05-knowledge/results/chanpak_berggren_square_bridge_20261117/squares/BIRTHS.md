# A literal first-birth tree inside the square-sum graph

In the square-sum graph on labels 1,…,N, retain just edges between **coprime perfect-square labels**. Such an edge is exactly {a²,b²} for a primitive Pythagorean triple a²+b²=c². It first occurs at N=max(a²,b²). In particular the 3–4–5 edge {9,16} first occurs at **16**, whereas the 3–4–5 reflection path first spans the interval through **15**. These are different birth objects. The children of the edge root are born at 144, 441, and 225.

Both legs strictly increase under each classical Berggren move: for A the added terms are 2(c−b) and 2(a+c−b), for B all additions are positive, and for C use the reflected A argument. Thus the birth cutoff selects an ancestor-closed tree. A nonempty prefix of k nodes has exactly 2k+1 frontier children because its k−1 internal edges occupy that many of its 3k child slots. A priority queue keyed by max(a²,b²) therefore produces the born edges in increasing order without enumerating other square-sum edges. This parent relation is a graph on these **edge objects**, not ordinary adjacency of vertices in the square-sum graph.

Keeping only the sum loses the birth size: 63²+16²=33²+56²=65², but the literal edges are born at 3969 and 3136. This is the analogue of retaining the marked pair rather than a scalar Hamilton total.

## Number of such births

Use the unique primitive Euclid parameters m>n>0, coprime and of opposite parity. Scale both parameters by N^(1/4). In polar coordinates (r,θ), 0<θ<π/4, the scaled region satisfies

    r² max(cos(2θ),sin(2θ)) ≤ 1.

Its area is (1/2)∫ min(sec(2θ),csc(2θ)) dθ = (1/2)log(1+√2). The two opposite-parity lattice classes have total density 1/2, and excluding odd common divisors multiplies this by ∏_(odd prime p)(1−p^(−2))=8/π². Using Möbius inversion over odd divisors, the bounded piecewise smooth boundary gives the sum of lattice errors O(N^(1/4) log N); the omitted density tail beyond a constant times N^(1/4) is smaller. The number of distinct edges is therefore

    [2 log(1+√2)/π²] √N + O(N^(1/4) log N).

The branch and Euclid enumerations agree through N=10^8 (1,788 edges).

## Why the old scalar lift remains thin on these exact addresses

The [previous square-tree bank](https://github.com/eliottcassidy2000/math/blob/87d2bbcd29280356a85292d39f8612982bda70f1/05-knowledge/results/chanpak_graceful_ternary_20261115/SQUARE_TREE.md) reads H_j(x,y)≥432xy and proved O(X^(3/4)) support over arbitrary square-sum addresses. Restricting to these direct primitive edges sets (x,y)=(a²,b²), so ab≤Y=√(X/432). In Euclid variables, ab=2mn(m−n)(m+n) is bounded below by a constant times m³ min(n,m−n): separate n≤m/2 from n≥m/2. The allowed pair count is then

    O(Σ_m min(m,Y/m³)) = O(Y^(1/2)) = O(X^(1/4)),

splitting at m=Y^(1/4). The two label orientations and finitely many window choices change only the constant. Hence a perfect address tree in the graph does not make that tournament bank an absolute selector. The older [norm-one absorber bridge](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/square_sum_berggren_absorber_bijection_20261019.md) instead maps the same words to a physical path substitution; it is not the same edge projection.

Run `python births.py` for both exact finite enumerations, the three child births, and the same-sum collision. The asymptotic proof above does not rely on the finite run.
