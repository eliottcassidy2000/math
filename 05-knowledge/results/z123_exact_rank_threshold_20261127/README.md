# Rank six is the exact threshold for failure of matroid Z-log-concavity

2026-11-27. This follows up [guoriyue/z123](https://github.com/guoriyue/z123)
and our [Gabidulin–uniform rank-six construction](../z123_rank_six_uniform_spike_20261127/README.md).
The new ingredient is a rank-three incidence Gram bound, applied separately to
each corank-three flat. It rules out *every* lower-rank log-concavity
counterexample. A second construction gives a smaller witness at rank six
and separates real-rootedness from log-concavity, even under 3-connectivity.

## 1. Hyperplane-pair bound (all ranks)

Write $z_i=[t^i]Z_M(t)$. For every loopless matroid of rank $r\ge4$,
$$
z_2\le\binom{z_1}{2},
$$
with equality exactly when the lattice of flats is Boolean.

Indeed, write $W_i$ for numbers of rank-$i$ flats and $f_{ij}$ for
incidences among flats of those ranks. The definition
$Z_M(t)=\sum_F t^{\operatorname{rk} F}P_{M/F}(t)$, the degree bound on
Kazhdan–Lusztig polynomials and palindromicity give
$$
 z_1=W_{r-1},\qquad
 z_2=W_{r-2}+f_{r-3,r-1}-f_{r-3,r-2}.
$$
To see the second identity, read the coefficient at degree $r-2$.
Only corank two (constant term) and corank three (linear term) flats
contribute. In a rank-three contraction the linear KL coefficient is
its number of lines minus its number of points: compare degrees one and
two of its palindromic Z-polynomial. Parallel elements in a contraction
cause no change: its points are rank-one flats.

For a fixed corank-three flat $X$, let $h_X,p_X$ be the numbers of
lines and points in the simplification of $M/X$. Let $e_X$ count
unordered *disjoint* pairs of lines there. Then $h_X-p_X\le e_X$:
select one endpoint of each disjoint pair and discard the selected lines.
The remaining lines intersect pairwise. Two distinct surviving lines
meet in exactly one point, and each line contains at least two points.
Their point-incidence Gram matrix is
$J+\operatorname{diag}(|L|-1)$, positive definite. Thus at most $p_X$
lines survive, while at most $e_X$ were discarded.

Each of these $e_X$ pairs is a hyperplane pair of $M$ with meet
exactly $X$. Each corank-two flat $Y$ also consumes a hyperplane
pair: its rank-two quotient has at least two points, and any two
hyperplanes through $Y$ meet exactly in $Y$. These budgets are
disjoint, because a pair has a unique meet. Therefore
$$
 z_2=W_{r-2}+\sum_X(h_X-p_X)
 \le W_{r-2}+\sum_X e_X\le\binom{W_{r-1}}2 .
$$
This is the same sort of local-to-global bookkeeping used in the
[Moore Gram and last-star notes](https://github.com/eliottcassidy2000/math/commit/f65d0a4cee45698e76836c2adda807933d2aa1a7): exact meet labels, not just
the number of contacts, prevent reuse of the same pair.

If equality holds, every corank-two flat is in exactly two hyperplanes,
since three would use three pairs where only one was budgeted. Simplify
and fix a basis $B$. If a point $e\notin B$ exists, its fundamental
circuit contains at least two basis points $b,c$. In the rank-two
quotient by $\operatorname{cl}(B\setminus\{b,c\})$, the points
represented by $b,c,e$ are all distinct: the fundamental circuit
implies $e\notin\operatorname{cl}(B\setminus\{b\})$ and likewise for
$c$. This contradicts the two-hyperplane property. Conversely, the
Boolean lattice gives $Z=(1+t)^r$.

**Consequence.** For ranks four and five, this proves the first
nontrivial log-concavity inequality $z_1^2\ge z_2$; the others follow
from palindromicity and [Braden–Huh–Matherne–Proudfoot–Wang's
unimodality theorem](https://arxiv.org/abs/2010.06088), Theorem 1.2(2).
Ranks at most three follow directly from palindromicity and unimodality.
Deleting loops does not change the relevant lattice or polynomial.
Thus the minimum rank for failure is six. In fact the ranks in which
failure occurs are exactly $r\ge6$: adjoin any fixed number $k=r-6$
of coloops **after** the truncation in Section 2. Z-polynomials multiply
under direct sums. The three coefficients have leading terms
$(7q^8/2,73q^8,8q^9)$ before adjoining coloops; afterward $z_1$
has the same leading term, $z_2=O_k(q^8)$, and $z_3\sim 8q^9$.
Hence $z_2^2-z_1z_3\sim-28q^{17}$ at each such fixed rank.
This does not settle the minimum rank for failure of real-rootedness.

## 2. Smaller witnesses, including a 3-connected family

For odd prime powers $q$, use the simple rank-four affine packing
matroid $A_q$ from [Bérczi–Coron–Eur–Li, Definition 4.10 and
Corollary 4.11](https://arxiv.org/html/2610.17135#S4.SS1), and the
Fano plane $F_7=\mathrm{PG}(2,2)$. The two families are
$$
M_{q,0}=\mathrm{Tr}(A_q\oplus F_7),\qquad
M_{q,1}=\mathrm{Tr}^2(A_q\oplus F_7\oplus U_{1,1}).
$$
Their sizes are $q^4+q^2+q+8+j$ for $j=0,1$.
The second family is simple and 3-connected: truncation from rank
eight to rank six keeps both parts of a partition at their original
ranks if neither spans, so their ranks sum to at least eight; if one
spans, the other has rank at least two when it contains at least two
elements. For singleton partitions the complement has original rank
at least seven and still spans after truncation. Thus neither 1- nor
2-separations occur. In contrast the summand partition of $M_{q,0}$
has connectivity value one.

| Family | First $q$ with nonreal zeros | Size | First $q$ without log-concavity | Size |
| --- | ---: | ---: | ---: | ---: |
| $M_{q,0}$ | 29 | 708,159 | 173 | 895,775,151 |
| $M_{q,1}$ | 31 | 924,522 | 181 | 1,073,316,072 |

Here “first” refers to allowed odd prime powers in each family.
Each family is strictly log-concave at its first nonreal-root witness,
with exactly two negative real roots and four nonreal roots from its
first threshold onward. For instance at $q=29,j=0$ the polynomial is
palindromic with $z_1=1882165290897$, $z_2=38793405082335$,
and $z_3=206241535134128$. Its index-two gap is
$1116747618911365135096619409>0$. Thus the two failures are not
equivalent even inside rank six. The root claim concerns
[Proudfoot–Xu–Young's ordinary Z-polynomial](https://arxiv.org/html/1706.05575#S5),
not the Chow, augmented Chow, or dual Chow polynomial.

## 3. Exact certificate and control

Set $z=(1,a,b,c,b,a,1)$. In any rank-six matroid,
$$
c=W_3+f_{25}-f_{23}+f_{14}-f_{13}+f_{135}-f_{134}-f_{125}+f_{123}.
$$
This is read at degree three in the defining sum. The corank-four
linear KL coefficient is $W_3-W_1$ for the contraction; the
corank-five quadratic KL coefficient is $W_3-W_2+f_{24}-f_{23}-f_{14}+f_{12}$,
by comparing degrees two and three in its palindromic Z-polynomial.

Truncation preserves strict flags not reaching the target rank.
In a direct sum enumerate weakly increasing coordinate rank strings
summing to each target rank string, independently choose the component
flags, and drop repeated ranks and the forced endpoints. The Fano
block has seven points, seven lines, and 21 point-line incidences.
The replay verifies the resulting explicit linear forms in the eight
affine flag coordinates against rank splitting on a full coordinate
basis. It uses the affine flags of Corollary 4.11, splitting off the
$q^2+q+1$ auxiliary points. This needs no enumeration of the very
large matroids.

For the reciprocal cubic
$$
 F(x)=x^3+a x^2+(b-3)x+c-2a,qquad Z(t)=t^3F(t+t^{-1}),
$$
the Newton expansions of its discriminant and $F(-2)$ have strictly
negative and positive coefficients, respectively, at the integer bases
29 and 30 for the two families. (The Newton basis is
$\binom{q-q_0}{j}$, with coefficients computed by forward differences.)
On the earlier integer ranges the discriminant and the coefficients
of $F(x-2)$ are all positive. This places all three cubic roots
below $-2$ before the switch, and its unique real root below $-2$
afterward. Each real cubic root below $-2$ produces two negative
real roots of Z. Nonreal cubic roots cannot come from real $t$.
Similarly $b^2-ac$ has negative Newton coefficients at integer
bases 171 and 180 and positive values at every smaller integer
$q\ge2$; $a^2-b>0$ and $c-b>0$ have positive Newton
coefficients at base 2 in both families.

Finally, $\mathrm{Tr}^2(A_q\oplus F_7)$, with no extra coloop,
has rank five and Z-coefficients $(1,a,b,b,a,1)$. Its gamma form is
$(1+t)^5(1+(a-5)u+(b-3a+5)u^2)$, where $u=t/(1+t)^2$.
Both nonconstant gamma coefficients and their quadratic discriminant
$(a+1)^2+4-4b$ have strictly positive Newton coefficients at 2.
The two gamma roots are negative, giving four negative reciprocal Z
roots in addition to $-1$. Thus this specific reduction in rank
removes the root failure; the universal rank-five result above concerns
log-concavity only.

Run `python3 replay.py` for the integer certificates and
`replay.out` for the expected output.
