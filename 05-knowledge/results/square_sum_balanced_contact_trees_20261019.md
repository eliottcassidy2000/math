# Polynomial-bit growing balanced AP contact trees from a free two-map code

**Status: PROVED, pending independent audit.** Exact symbolic identities and finite controls are supplied separately. This is a growing-tree result; its shape includes one subdivision per prescribed tree edge and one stem edge from the prescribed seed. It is not a polynomial-bit construction of an arbitrarily long contact chain.

## 1. Statement and graph shape

Let `B` be a finite rooted binary tree with `n>=1` vertices and height `D` (maximum root-to-vertex edge count). Let `T_0` be any admissible primitive Pythagorean petal of hypotenuse `C`, and prescribe a port `epsilon in {-1,+1}`. Write `L=1+ceil(log_2 C)`.

**Theorem 1.** There is a deterministic exact algorithm producing `2n-1` further distinct admissible primitive petals, all different from `T_0`, with the following contact subgraph:

1. subdivide each edge of `B` exactly once;
2. join its root by one additional stem edge to `T_0`;
3. orient every edge away from `T_0`, using port `epsilon` at its parent and `-epsilon` at its child.

Every specified edge is a rational square contact. Thus every root-to-leaf path is port-compatible. After a common integer square scaling, the `2n` petals and their common center give `4n+1` distinct positive integer labels, realizing all these contact edges and all AP-petal triangles.

For an absolute effective constant `A`, the final labels have

\[
 O\!\left(n^5 e^{AD}L^2[L+D+\log(n+1)]\right)           \tag{1}
\]

bits, and construction takes `exp(O(D)) poly(n,L)` bit operations. In particular, for any fixed constant `c`, trees with `D<=c log(n+1)` have polynomial output-bit complexity and polynomial running time in `n,L`.

The extra edges of the ambient square-sum graph are unrestricted. Neither a graceful labeling nor a Hamiltonian cycle through all vertices of a branching contact tree is asserted.

## 2. The two sections and an exact local code

We inherit the fixed-count cusp sections from [square_sum_fixed_count_contacts_20261019.md](square_sum_fixed_count_contacts_20261019.md). On the plus curve, with its positive-sheet anchor `Q` and `P=(0,k_+(t))`, the odd section

\[
 R_d=dQ+(d-1)P
\]

supplies a parameter `F_d(t)` and the square identity

\[
 1+f(t)-f(F_d(t))\in\mathbb Q(t)^2.
\]

Use only `d=3,5`. The first is the degree-11 map already audited in the plus atlas; the second has degree `33`. Exact rational-function group-law arithmetic gives

\[
\begin{aligned}
 F_3(1+h)&=1+3h+2h^3+O(h^4),\\
 F_5(1+h)&=1+5h+O(h^4).                                \tag{2}
\end{aligned}
\]

The companion checker derives both maps from the group law, reduces their numerators and denominators, and verifies these four Taylor coefficients by integer-polynomial division. In particular the quadratic coefficients vanish for both maps, and the cubic coefficient of `F_5` is exactly zero. This is an exact identity calculation, not a numerical fit.

Set

\[
 A=F_3\circ F_5,\qquad B=F_5\circ F_3.
\]

Both have degree `363`. In the coordinate `h=t-1`, their jets are

\[
 A:15h+250h^3+O(h^4),\qquad
 B:15h+10h^3+O(h^4).                                  \tag{3}
\]

For a word `w` of length `d` in these two block maps, listed in application order, its derivative is `15^d`, its quadratic coefficient is zero, and its cubic coefficient is

\[
 c_w=15^{d-1}\sum_{i=1}^d c_i225^{i-1},\qquad
 c_i\in\{10,250\}.                                    \tag{4}
\]

Two words of equal length have different cubic coefficients: subtract the fixed baseline `10 sum 225^j`; the remaining base-225 digits are `0` or `240`, and dividing their difference by `240` reduces uniqueness to digits `0` or `1`. Different word lengths already have different derivatives. Hence `A,B` generate a free semigroup of rational maps.

### Intermediate nodes do not collide either

Each `A` block has intermediate parameter `F_5(t)` and each `B` block has intermediate parameter `F_3(t)`. These intermediate functions after a word `w` have derivatives `5*15^d` and `3*15^d`, respectively. These numbers never coincide with one another across distinct types or levels, or with any full-word derivative `15^e`. For the same type and level, a difference of word cubic coefficients remains nonzero after composition with `F_3` or `F_5`.

Thus any two distinct full or intermediate prefix functions first differ at Taylor order `1` or `3`, with a nonzero **integer** coefficient. This stronger statement, rather than mere formal freeness, will give a uniform real interval on which every desired vertex is distinct.

For the minus port use

\[
 F_d^-(t)=1/F_d(1/t).
\]

The inherited reciprocal section duality gives its minus contact identity and its regular cusp behavior. Formal conjugation by `t->1/t` preserves equality, the first order where two maps differ, and the integrality of the relevant jets. If a plus prefix has derivative `a`, zero quadratic coefficient, and cubic coefficient `c`, its minus conjugate has coefficients

\[
 a,\quad a^2-a,\quad a^3-2a^2+a+c.
\]

Consequently the same distinction at order `1` or `3` holds for the minus prefix family. Its real small positive cusp chart is established independently below; no admissibility claim is transported from `t<1` without checking it.

## 3. One effective analytic chart for all four fixed maps

There are effective absolute rational `0<r_0<=1/64` and integer `K>=6` such that each of `F_3,F_5,F_3^-,F_5^-`:

- is analytic on the closed complex disk `|t-1|<=r_0`;
- satisfies `|F(t)-1|<=K|t-1|` there;
- on the positive real interval satisfies `t<F(t)<5/4` and `F'(t)>0`.

Here is an exact, finite construction of these constants. Express each map as `N(h)/D(h)` with integer coefficients and `D(0)>0`, using the removable cusp chart. Their value and derivative at zero are `1` and `3` or `5`. Apply the elementary coefficient cutoff

\[
 r\le\min\left(1,\frac{p(0)}{2\sum_{i>=1}|p_i|}\right)  \tag{5}
\]

to each denominator and to the polynomials

\[
 (N-(1+h)D)/h,\quad (5/4)D-N,\quad N'D-ND'.
\]

Their constant terms are positive. Zero tails impose no restriction. Take the minimum of these cutoffs and `1/64`. For a denominator, (5) proves the stronger complex inequality `|D(h)|>=D(0)/2`.

Put `U=(N-D)/h` and choose

\[
 K\ge\max_F\frac{2\|U\|_1}{D(0)},\qquad K>=6.          \tag{6}
\]

Then `|F(1+h)-1|<=K|h|` on that disk. The real sign polynomials prove monotonicity, increasing parameter, and positivity. The inherited square identities prove every required local contact. All constants depend only on four fixed rational functions, not on `n,D,C`.

## 4. A polynomially small collision-free interval for a balanced tree

Every prefix needed for the once-subdivided binary tree has at most `2D` elementary maps. Define

\[
 \rho_D=r_0 K^{-2D},\qquad
 \delta_D=\min\left(\rho_D/2,\frac{\rho_D^4}{8r_0}\right).
                                                               \tag{7}
\]

These are exact rationals, with `log(1/delta_D)=O(D+1)`.

Every such prefix function `G` is analytic for `|h|<=rho_D`, and its displacement `G(1+h)-1` has absolute value at most `r_0` there. This follows by iterating (6): all intermediate arguments remain inside the disk of radius `r_0`.

Take two distinct prefix functions. Their difference is bounded by `2r_0` on the disk and first has a nonzero coefficient `a_j h^j` at order `j=1` or `3`, with `|a_j|>=1`. For `0<h<=rho_D/2`, Cauchy's coefficient bound gives

\[
 \left|G_1(1+h)-G_2(1+h)-a_jh^j\right|
 \le \frac{4r_0 h^{j+1}}{\rho_D^{j+1}}.                 \tag{8}
\]

When `h<=delta_D`, the right side is at most `h^j/2`, because `rho_D<=1` and `j<=3`. Thus the difference cannot vanish. This proves simultaneous vertex distinctness for **every** positive root displacement in `(0,delta_D]`; no root-separation estimate exponential in the expanded polynomial degree is needed.

All prefixes are strictly increasing real functions on this interval, because every elementary map is increasing on every intermediate real argument. Their output parameters lie in `(1,1+r_0)`, where spread is strictly increasing.

## 5. Attach an arbitrary prescribed seed without a tiny-seed penalty

Let `q=2n-1` be the number of full and intermediate vertices in the subdivided tree. Apply the boundary-start Lemma 2 of [square_sum_growing_chain_bit_bound_20261019.md](square_sum_growing_chain_bit_bound_20261019.md) to the prescribed seed and port, requesting

\[
 r=q+1=2n,\qquad\eta=\delta_D.
\]

It returns `2n` distinct low-parameter root candidates with spread in `(0,delta_D)`, each joined to `T_0` at the prescribed port. Since `f(1+h)>=h` in the low cone, all their root displacements are in the collision-free interval.

For each one of the `q` prefix functions, strict monotonicity means at most one candidate root makes that vertex's spread equal to `f_{T_0}`. Hence at most `q` candidates are forbidden. At least one of the `q+1` candidates avoids the seed at **every** descendant vertex. This avoids shrinking the root interval below `f_{T_0}`, which could have cost exponential dependence on the seed bit length.

Keep one candidate passing these exact tests. The free-prefix argument makes all new petals mutually distinct, and the candidate exclusion makes them different from `T_0`. The two-step block maps supply exactly the once-subdivided binary-tree contact edges; the boundary certificate supplies the stem edge. Common-center scaling gives the asserted integer graph.

The choice of “first child” and “second child” at each binary vertex may be arbitrary but fixed. A missing child simply omits its block and all descendants. No prescribed labeling of the abstract tree's nodes is needed.

## 6. Quantitative bound and reproducibility

By the boundary lemma, the maximum translation index is

\[
 Q=O(n^2 L/\delta_D)=O(n^2 e^{O(D)}L).
\]

The root's rational bit length is consequently

\[
 O\!\left(n^4 e^{O(D)}L^2[L+D+\log(n+1)]\right).         \tag{9}
\]

Every full prefix has degree `363^d` for `d<=D`, and every intermediate prefix has degree at most `33*363^(D-1)<=363^D` when `D>=1`. Repeated evaluation of the fixed homogeneous forms bounds each output height by `O(363^D)` times the root height plus an absolute constant. Primitive normalization only changes this by an absolute factor. The common hypotenuse lcm has bit length at most the sum of the `2n` primitive lengths. Integer labels are below four times its square, proving (1).

The number of tested orbit points, rational operations, and operand lengths is `exp(O(D)) poly(n,L)`. Schoolbook arithmetic therefore gives the claimed exact running time. For balanced input trees, `D=O(log n)` makes every factor polynomial in `n,L`.

The companion script derives the maps and jets, checks their contact identity at exact rational parameters, computes valid constants `r_0,K`, and tests the base-225 prefix distinction at bounded depths by a second direct jet-composition implementation. Those finite tests supplement the all-word proof (4), and do not replace it.

## 7. What changed and what remains

Source: bounded-height rooted binary tree shapes plus an arbitrary primitive seed. Target: a sparse AP-paired integer square-sum graph containing the specified subdivided contact tree and seed stem. Map: a free two-block rational-map code, uniform analytic separation, and one boundary-root selection. Preserved: distinct vertices, directed port contacts, common-center AP triangles, and the original seed after integer square scaling. Lost: exact unsubdivided tree shape, consecutive-interval coverage, small numerical labels, and graceful edge differences.

The three ingredients fix separate problems: different cubic jets prevent generic branch merging; Cauchy bounds prevent specialization collisions at a polynomially small root scale; boundary-start recurrence makes that scale reachable from an arbitrary seed in polynomial bits when the depth is logarithmic.

A long chain has depth linear in its number of vertices, so this theorem's estimate is exponential in that size. The companion chain theorem improves that case quantitatively but leaves the polynomial-bit long-chain problem open. The exact common-center capacity obstruction rules out polynomial numerical labels even before contacts are imposed. No equivalence to Collatz or to the graceful-tree conjecture is claimed.


## Reproduction from the repository root

The exact control requires Python 3 and SymPy. Run `python3 05-knowledge/experiments/square_sum_balanced_contact_trees_20261019.py`. The matching receipt is `05-knowledge/results/square_sum_balanced_contact_trees_20261019.json`; the script prints a fresh receipt.
