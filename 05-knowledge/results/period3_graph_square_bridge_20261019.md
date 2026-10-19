# Period-three capacity, graph automorphisms, and a minimal square-sum orbit

2026-10-19. **PROVED, elementary:** the real order-three realization criterion, the affine rationality obstruction, and an infinite family of square-sum embeddings of the rational three-cycle. **FINITE-EXACT:** the minimum square-sum AP triangle and the graph automorphism inventories. **INHERITED:** the quadratic period-three examples, uniqueness of the AP three-cycle, and the bipartite sum/difference gauge. No global novelty, Hamiltonicity, graceful-tree, planarity, or LRC claim is made.

Companion: `04-computation/experiments/period3_graph_square_bridge_20261019.py`; transcript: `05-knowledge/results/period3_graph_square_bridge_20261019.out`.

## 1. Inheritance and the question that actually transfers

Closest proved mechanisms:

- [THM-4139, rational three-cycle/order-six lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md), corrected/scoped by [THM-4146, horizontal divisor/fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md): the rational orbit of `x²−29/16` is the unique AP-supported three-cycle of a centered monic quadratic.
- [Quadratic observable-return note](quadratic_clock_return_20261012a.md): at `c=−7/4` the period-three cubic is squared in the dynatomic polynomial, but still has only three distinct roots. Its multiplier is one. At `c=−29/16` the rational cycle multiplier is `35/8`; the secant product is one in both cases and is not the multiplier.
- [THM-2761, edge-sum discriminants and graceful sign gauge](../../01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md): a bipartite sign change exchanges edge sums and oriented differences, but does not provide a graceful labeling.
- [THM-4701, finite unit-circle/Petersen graphs](../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md): Petersen as `KG(5,2)`. The geometry there is not the scalar dynamics below.

The canonical hostile is counting a repeated dynatomic root as a second vertex. The corrected near miss is treating a square-sum path, a quadratic orbit, and a graph automorphism as the same object. The least-used sidecar here is the **cycle type of the proposed vertex permutation**. The concept board is: period-three roots / permutation cycles / square-sum triangles / affine normalization / graceful edge collisions.

The concrete question is: given a finite graph `G` and a specified automorphism `σ` of exact order three, can distinct real numbers `ℓ(v)` satisfy

\[
\ell(\sigma v)=f_c(\ell(v)),\qquad f_c(x)=x^2+c? \tag{1}
\]

This preserves the graph's adjacency through the already specified automorphism. It does not require adjacency to be defined by distances or square sums of these labels.

## 2. A complete real capacity criterion

Write `F(σ)` for the number of fixed vertices and `T(σ)` for the number of three-element permutation orbits.

**Proposition 1.** For an exact-order-three automorphism, (1) has an injective real labeling precisely as follows:

- If `c<−7/4`, then `F(σ)≤2` and `T(σ)≤2`.
- If `c=−7/4`, then `F(σ)≤2` and `T(σ)≤1`.
- If `c>−7/4`, no such labeling exists.

Here `T(σ)≥1` is automatic from exact order three.

**Proof.** Put

\[
P_s(X)=X^3-sX^2+(c-s-1)X-cs-c-1,
\qquad c=-s^2-s-2.
\]

Direct polynomial identities give

\[
f_c^3(X)-X=(f_c(X)-X)P_s(X)P_{-1-s}(X), \tag{2}
\]

\[
\operatorname{disc}(P_s)=(4s^2+6s+9)^2>0,\qquad
\operatorname{Res}(P_s,f_c-X)=4s^2+6s+9>0, \tag{3}
\]

and `P_s(f_c(X))` is zero modulo `P_s`. Therefore, for real `s`, its three roots are distinct and real, none is fixed, and the map cyclically permutes them. The trace parameter solves

\[
s^2+s+c+2=0.
\]

For `c<−7/4` there are two distinct real traces and hence two distinct three-cycles: shared roots would force the same full three-cycle, whose trace cannot have two values. At equality the two cubics coincide. The two fixed roots are distinct in either case.

For completeness, the absence of real three-cycles above the boundary follows directly from interpolation rather than a reality assumption on the factors. Normalize any oriented triple to `0→1→t→0`, where `t` is real and neither zero nor one. Let its quadratic interpolant be `q(X)=AX²+BX+C`. The normalized parameter `AC+B/2−B²/4` obeys

\[
c+\frac74=-\frac{(t^3-2t^2-t+1)^2}{4t^2(t-1)^2}\le0. \tag{4}
\]

The interpolant is genuinely quadratic because `A=−(t²−t+1)/(t(t−1))≠0`. Affine conjugacy preserves the normalized parameter. Necessity now follows by assigning different permutation orbits to different dynamical cycles. Conversely, choose the needed distinct cycles and fixed roots, then label each permutation orbit in its cyclic order. Since `σ` was already a graph automorphism, this construction proves (1). QED.

The elementary cardinality bound `|V|≤deg(f_c³−X)=8` holds even over the complex numbers. At the parabolic parameter the number of distinct roots falls to five. This is the first check, before using graph structure.

## 3. K5, K3,3, and Petersen: what the graph side contributes

The exact-order-three automorphism types are:

| Graph | Permutation type | Number of such automorphisms |
|---|---|---:|
| `K5` | `1+1+3` | 20 |
| `K3,3` | `1+1+1+3` | 4 |
| `K3,3` | `3+3` | 4 |
| Petersen | `1+3+3+3` | 20 |

For `K5` this is immediate from permutations of five vertices. An odd-order automorphism of `K3,3` preserves its two intrinsic bipartition classes; an order-three permutation of each class is either identity or a three-cycle, giving the two rows.

For Petersen, regard vertices as the two-element subsets of `{1,…,5}`, adjacent when disjoint. The only independent four-sets are the five stars consisting of all pairs through one base element. Indeed, a pairwise-intersecting family with no common element has size at most three: after `{a,b},{a,c},{b,c}`, no fourth distinct pair intersects all three. Consequently every automorphism permutes those five stars faithfully. Base-set permutations already give 120 automorphisms, so the automorphism group is exactly the induced `S5`. A base three-cycle fixes the pair on the two remaining base elements, cycles the three pairs inside its support, and gives one three-cycle of pairs through each outside element. Its type is `1+3+3+3`.

Thus:

- `K5` supports (1) for every `c≤−7/4`, using two fixed points and one three-cycle.
- `K3,3` supports (1) exactly when `c<−7/4`, necessarily using the `3+3` type. It fails at equality.
- Petersen never supports an order-three action of a quadratic on distinct scalar labels. **The simpler 10>8 degree obstruction already proves this**; the automorphism type additionally records the three-cycle demand explicitly. This is not an obstruction unique to Petersen and is not about its non-Hamiltonicity or nonplanarity.

At the two requested parameters, write

\[
P(X)=X^3+\tfrac12X^2-\tfrac94X-\tfrac18,
\]

\[
Q(X)=X^3+\tfrac14X^2-\tfrac{41}{16}X+\tfrac{23}{64},\qquad
R(X)=(X+\tfrac74)(X-\tfrac54)(X+\tfrac14).
\]

Then `Φ3(X,−7/4)=P(X)²` and `Φ3(X,−29/16)=Q(X)R(X)`. Both `P` and `Q` are irreducible over `Q` and have three real roots. The first case therefore realizes `K5` over the reals but not `K3,3`; the second realizes both over the reals.

**Rational firewall.** Neither graph realization at these parameters has all rational labels. At `−7/4`, there are no rational roots of `f³−X`; at `−29/16`, its only rational roots are the three roots of `R`. The fixed-point discriminants are respectively `8` and `33/4`, both nonsquares. In particular the `K3,3` construction at `−29/16` needs the other, irreducible cubic cycle. It is not two rational three-cycles.

## 4. A genuine square-sum embedding of the rational orbit

Let `Q_N` be the loopless graph on `1,…,N` in which distinct vertices are adjacent exactly when their sum is a square. The inherited rational orbit is

\[
-\tfrac74\longmapsto\tfrac54\longmapsto-\tfrac14\longmapsto-\tfrac74.
\]

The affine map

\[
L(X)=1920X+3842 \tag{5}
\]

sends it to

\[
482\longmapsto6242\longmapsto3362\longmapsto482,
\]

and the underlying triangle really lies in `Q6242`:

\[
482+3362=62^2,\qquad482+6242=82^2,\qquad3362+6242=98^2. \tag{6}
\]

The transported map is explicitly

\[
(L\circ f\circ L^{-1})(Y)=\frac{(Y-3842)^2}{1920}+362. \tag{6a}
\]

Its critical point is 3842. In contrast, `X²+c` is strictly increasing on positive labels and cannot itself cycle three of them. The map (5) preserves the orbit only together with the transported quadratic; its scale and translation are essential data. Also, this three-vertex subgraph says nothing about a Hamiltonian path through all 6242 vertices.

**Finite-exact minimum.** This is the unique positive integer square-sum triangle in arithmetic progression whose largest vertex is at most 6242. Consequently 6242 is the smallest possible largest vertex of an affine integer square-sum realization of this orbit.

Here is the complete finite certificate. For `a<b<c`, `a+c=2b`, write

\[
a+b=p^2,\quad a+c=q^2,\quad b+c=r^2.
\]

Then

\[
p^2+r^2=2q^2,\qquad
(a,b,c)=\tfrac12(2p^2-q^2,\ q^2,\ 3q^2-2p^2). \tag{7}
\]

If `c≤6242`, all pair sums are at most `2·6242−1`, so all three positive square roots are at most 111. Enumerate `1≤p<q≤111`, test `r²=2q²−p²`, then the integrality, positivity, strict ordering and upper-bound tests in (7). The only survivor is `(p,q,r)=(62,82,98)`. A separately coded enumeration starts from the even square `2b=q²` and reproduces the same unique answer. The bounded search plus the exhibited witness establishes the global minimum, not a conjectural extrapolation.

**Infinite family, proved.** For every integer `k≥4`, put

\[
p=2k^2-1,\quad q=2k^2+2k+1,\quad r=2k^2+4k+1,
\]

\[
(a,b,c)=\bigl(2(2p^2-q^2),\ 2q^2,\ 2(3q^2-2p^2)\bigr). \tag{8}
\]

The identity `p²+r²=2q²` makes this an AP with pair sums `(2p)²,(2q)²,(2r)²`. It is positive: writing `k=j+4`,

\[
\frac a2=4j^4+56j^3+272j^2+508j+241>0.
\]

The inequalities `p<q<r` imply `a<b<c`. At `k=4` this gives (6). Since every nonconstant AP triple is an affine image of the rational orbit, (8) provides infinitely many genuine square-sum embeddings.

**Parabolic obstruction.** No real affine transformation takes all three points of the `−7/4` orbit to rational numbers, hence none takes them to a square-sum triangle with integer labels. If an affine conjugacy gave three distinct rational labels, the quadratic interpolating their cyclic action would have rational coefficients `A,B,C`. Normalizing by `Y↦AY+B/2` would give a rational three-cycle of `X²−7/4`, contradicting the irreducibility of `P`. The normalized parameter is unchanged by the original affine conjugacy. Equivalently, equality in (4) requires a root of the irreducible cubic `t³−2t²−t+1`, so the relative shape is irrational.

## 5. What survives toward square-sum paths and graceful trees

For a square-sum path `v0,…,vm`, put `y_i=(−1)^i v_i`. Then

\[
|y_{i+1}-y_i|=v_{i+1}+v_i.
\]

This is exactly the inherited bipartite sign gauge. It preserves edge magnitudes; it does not make them distinct, consecutive, or within the graceful-label range. For the unique square-sum Hamiltonian path at 15,

`8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`,

the fourteen edge magnitudes use only `{9,16,25}`. Thus collisions persist after the gauge. The path graph itself is graceful, but this transported labeling is not. A second hostile is the AP three-cycle: either two-edge path through its middle has equal gaps, so affine AP geometry alone cannot produce a graceful labeling of that path.

The useful common quantity is the quadratic secant

\[
f_c(u)-f_c(v)=(u-v)(u+v).
\]

Its product around a nontrivial periodic orbit is one after telescoping difference ratios. A square-sum constraint prescribes `u+v` additively, whereas the dynamical condition also requires the image differences to match. Positive original labels make every edge sum greater than one, so they cannot themselves form a nontrivial quadratic cycle. The affine conjugation in (5), with its changed quadratic coefficient and linear term, is an essential part of the successful bridge.

## 6. Reproduction, provenance, and remaining boundaries

Run from repository root:

```bash
python3 04-computation/experiments/period3_graph_square_bridge_20261019.py
python3 -O 04-computation/experiments/period3_graph_square_bridge_20261019.py
```

The two outputs byte-match. The script uses exact rational/symbolic arithmetic, explicit failure exceptions, the complete square-root search with its independent middle-label check, direct graph-adjacency automorphism backtracking, and a separate comparison with all 120 induced `S5` actions on `KG(5,2)`. It checks all factorization/discriminant identities and the all-`k` family symbolically. Finite graph counts are not being used to infer an unbounded theorem.

Primary literature checked on 2026-10-19: Patrick Morton and Serban Raianu, [Arithmetic properties of 3-cycles of quadratic maps over Q, arXiv:2105.07435v5](https://arxiv.org/abs/2105.07435v5). Their arithmetic minimality theorem for `−29/16` and numerator graph are established prior work; the present note does not claim them. The elementary period-three algebra and graph constructions above are proved here, with the already-canonical overlap identified explicitly. A narrow search did not establish prior publication of the 6242 minimum; no priority claim is warranted.

**Not established:** a new square-sum Hamiltonian-path theorem, a graceful-tree existence theorem, a new Petersen obstruction beyond the explicitly stated scalar action, or any connection preserving planarity/unit-distance geometry. The open useful next question is whether other specified rational dynamical cycles admit low-height square-sum embeddings; the AP construction settles only this one affine-shape class.
