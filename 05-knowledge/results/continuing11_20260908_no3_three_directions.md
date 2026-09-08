# A third-direction repair charge and the limit of flow duality

**Status: PROVED analytical statements; FINITE-EXACT bounded controls; independently audited.** This is a scoped advance in the deletion cost of repairing a two-regular grid board. It is not a bound on the extremal no-three-in-line number, a proof that all saturated boards fail, or a running-time bound for an adaptive repair algorithm.

The inherited proved mechanism is the unit-capacity bipartite flow and uniform fixed-row permutation-avoidance argument in `05-knowledge/results/continuing10_20260907_no3_repair.md`, independently checked in `05-knowledge/results/continuing10_20260907_no3_repair_audit.md`. The retained sidecar is which *original* cells touch an originally overfull line. Losing that information causes the naive third-direction charge to double-count deletions. The least-used representation is the exact min-cut dual, which reveals both a stronger safe charge and a small obstruction to extending integral flow to three directions.

The live concept board is: original-cell retention; capacity-two line constraints; dual line selection; isolated third-direction triples; uniformly sparse forbidden permutation boards. The connection changes the first two from a two-family problem to a three-family problem while preserving an independently chargeable subset of original cells. It also identifies exactly where the integral bipartite representation stops being valid.

## 1. Objects and exact source-to-target map

Write cells as `(r,c)`, with row coordinate `r` and column coordinate `c`. The three line labels are

\[
d=c-r,\qquad s=c+r,\qquad q=c-2r.
\]

Let `B` be any finite set of distinct cells. Let `tau_2(B)` be the minimum number of original cells that must be deleted to leave at most two cells on every `d`- and `s`-line. Let `tau_3(B)` impose all three families. Multiplicities in a line always refer to the original board unless explicitly stated otherwise.

For two families, send a cell to a unit-capacity edge between its `d`- and `s`-vertices, each vertex having retained-degree capacity two. Distinct cells give distinct edges. The map preserves exactly the two-family feasibility predicate. It discards all other collinearities; those require a sidecar or new constraints. The cheapest decisive test for a proposed extension is therefore an exact minimum-deletion calculation, rather than a count of individual bad lines.

For three families, a cell instead becomes a three-partite, three-uniform hyperedge incident with its `d,s,q` vertices. Vertex capacities remain two and cell capacities remain one. This preserves three-family feasibility but does not preserve the integrality of the two-family flow polytope.

## 2. An independent third-direction charge

Define the originally safe subset

\[
U=\{p\in B:\deg_d(p)\le2\ \hbox{and}\ \deg_s(p)\le2\},
\qquad
J_3(B)=\sum_q (|U\cap q|-2)_+.
\]

Then, for every finite board,

\[
\boxed{\tau_3(B)\ge\tau_2(B)+J_3(B).} \tag{1}
\]

Indeed, if a deletion set clears the two old families, removing its intersection with `U` from the deletion set still clears those families: no point of `U` lies on an originally violated old line. Consequently every feasible three-family deletion set deletes at least `tau_2(B)` cells outside `U`. Its deletions inside `U` must leave at most two cells on each `q`-line, so there are at least `J_3(B)` of these. The two charges concern disjoint cells. This proves (1) without an independence assumption and for adaptive repair procedures as well.

Let `I_3(B)` count `q`-lines containing exactly three original cells, each on a singleton `d`-line and singleton `s`-line. These triples are wholly contained in `U`, and distinct `q`-lines are disjoint, so

\[
J_3(B)\ge I_3(B). \tag{2}
\]

The weaker condition defining `U` is intentional: singleton diagonals will pay for a uniform probabilistic event, while degree at most two is the exact deterministic safe condition.

If a final no-three-in-line set `T` is formed by any combination of deletions, insertions, or moves, the retained original cells `B intersection T` must satisfy the three-family constraints. Therefore `|B setminus T| >= tau_3(B)`. Equation (1) charges loss of original cells, not the number of operations of a particular repair implementation.

## 3. Exact two-family dual, strongest Boolean extension, and failure boundary

For a collection `Q` of line vertices, let `N(Q)` be the union of original cells incident with those lines and put

\[
c(Q)=|N(Q)|-2|Q|.
\]

The two-family unit flow gives the exact min-cut formula

\[
\tau_2(B)=\max_{Q\subseteq D\sqcup S}c(Q). \tag{3}
\]

To see the precise correspondence, a cut with source-side old vertices `A subset D` and `H subset S` has capacity

\[
2|D\setminus A|+|E(A,S\setminus H)|+2|H|.
\]

Subtracting this from `|B|` yields `c((D setminus A) union H)`. Every selection of old line vertices occurs in this way. The cell edges have capacity one; omitting that condition would change the problem.

For all three families define the Boolean line-selection bound

\[
\beta_3(B)=\max_{Q\subseteq D\sqcup S\sqcup Q_2}c(Q),
\]

where `Q_2` denotes the set of slope-two line vertices, to distinguish it from a selected collection. Any feasible retained set has at most `2|Q|` cells in `N(Q)`, giving `tau_3 >= beta_3`. The same inequality holds for fractional retained cells, so

\[
\tau_3(B)\ge\tau_{\rm LP}(B)\ge\beta_3(B).
\]

Because the third family partitions the cells, maximizing over its line selection first gives the exact Boolean compiler

\[
\boxed{\beta_3(B)=\max_{Q\subseteq D\sqcup S}
\left(c(Q)+\sum_q(|q\setminus N(Q)|-2)_+\right).} \tag{4}
\]

For a fixed old selection, every newly selected `q`-line contributes its previously uncovered cells minus two; the new lines have no cells in common with one another. This proves (4).

There is an optimum in (3) containing no old line of degree at most two: removing such a vertex loses at most two cells from the union and refunds two. For that optimum, `N(Q)` is disjoint from `U`, and hence (4) also proves

\[
\tau_3\ge\tau_{\rm LP}\ge\beta_3\ge\tau_2+J_3.
\]

The compiler can be stronger than the universal safe-set charge. For example, the two-regular `5 by 5` board with row column-pairs

\[
(01,03,24,34,12)
\]

has `(tau_2,tau_3,J_3,beta_3)=(1,2,0,2)`. Select just the old line `s=6`, with three cells; the three cells of `q=-2` are disjoint from it, so (4) supplies an additional unit. They need not be safe against every originally overfull old line, only disjoint from the particular optimal dual selection.

Neither the Boolean compiler nor the fractional relaxation is generally exact. Consider

\[
(03,12,04,34,12).
\]

Its only overfull lines are the three triples

\[
d=0:\{(0,0),(1,1),(3,3)\},
\quad s=6:\{(2,4),(3,3),(4,2)\},
\quad q=0:\{(0,0),(1,2),(2,4)\}.
\]

Each pair of triples intersects, but the three-way intersection is empty. Thus `tau_2=1` and `tau_3=2`. Selecting one old/new line gives Boolean value one, selecting two gives `5-4=1`, and selecting all gives `6-6=0`; hence `beta_3=1`. Assign fractional deletion `1/2` to each of `(0,0),(2,4),(3,3)` and zero elsewhere. This is feasible with cost `3/2`. Giving each of the three deletion constraints dual weight `1/2` gives the matching lower bound, so

\[
\boxed{\beta_3=1<\tau_{\rm LP}=3/2<\tau_3=2.} \tag{5}
\]

The three constraints restricted to those three shared points have incidence minor

\[
\begin{pmatrix}1&0&1\\0&1&1\\1&1&0\end{pmatrix},
\]

whose determinant is `-2`. This is the first failed implication in extending bipartite flow: hypergraph incidence is not totally unimodular. The strongest survivors are (1), (4), and a genuine fractional optimization if more information is needed. For `n<=4`, every slope-two line has at most two grid cells, so a new-direction gap is impossible. Thus the above `n=5` example is dimension-minimal.

The naive charge `tau_2 + sum_q(deg_q-2)_+` is false even for two-regular boards: `(01,01,23,24,34)` has `tau_2=tau_3=4` and new-direction excess one. The same deletions can clear both obligations. The missing coordinate is which original cells are already charged by the old constraints.

## 4. A uniform positive mean bonus for every fixed row order

Fix any simple bipartite two-regular skeleton `G` on `n+n` vertices and **any** physical row labelling `rho`. Randomize only the column labelling `sigma`, uniformly over its `n!` values. All asymptotics in this section are uniform in `G,rho`. No row order favorable to slope two is assumed.

For an integer `q`, take a triple of distinct rows `r_1,r_2,r_3` for which the targets `c_i=q+2r_i` lie in the grid. Their three `d` labels `q+r_i` and three `s` labels `q+3r_i` are automatically distinct. Each target row has two source-column neighbors, so there are eight choices of three intended source edges.

Call a row triple **eligible** when the following hold. For each selected row `r`, each of its two source-column neighbors has a second row neighbor `v`. No other selected row may be

\[
v,\qquad 2r-v,\qquad (2r+v)/3,
\]

with the last address relevant only when integral. This is at most six excluded other-row addresses for each `r`; repetitions only decrease the count. A slope-two line with `L` possible rows therefore has at most `6L(L-2)` ineligible unordered triples when `L>=3`. Summing over all lines excludes `O(n^3)` triples, with an absolute constant independent of the skeleton and row order.

For an eligible triple, all eight source-edge choices use distinct source columns, because shared source columns would make `v` another selected row. Condition on mapping those columns to `c_i`; the probability of each such assignment is `1/(n)_3`. Its forced companion point `(v,q+2r)` lies on none of the seven target lines:

* it cannot lie on the target `q`-line unless `v=r`;
* it lies on target `d_j` only if `r_j=2r-v`;
* it lies on target `s_j` only if `r_j=(2r+v)/3`.

The case `j=i` in the last two conditions again forces `v=r`, which is impossible in a simple bipartite two-regular skeleton.

The remaining `n-3` source columns are uniformly permuted to the remaining physical column labels. Forbid an assignment if either of its two cells lies on any of the seven target lines. Each forbidden-matrix row has degree at most 14, because a source column has two row neighbors and each physical row meets at most seven lines. Each forbidden-matrix column has degree at most 14, because a physical column meets at most seven lines and each physical row has two source-column neighbors. If the sum of the seven full-grid line lengths is `S`, the total number `m` of forbidden entries is at most `2S`. This bound allows crossings and deliberate overcounting. Since the slope-two line has length at most `ceil(n/2)` and each other line has length at most `n`,

\[
m\le2S\le13n+1.
\]

Apply the supplier's uniform bounded-degree permutation-avoidance lemma with degree bound 14 and `N=n-3`:

\[
\mathbb P(\mathrm{avoid})=e^{-m/N}+o(1)
\ge e^{-2S/n}-o(1). \tag{6}
\]

The last error is still uniform since `S=O(n)`. Avoidance produces exactly three cells on the target `q`-line and singleton opposite diagonals, so it contributes to `I_3`. Different intended assignments, and different triples on the same `q`-line, produce disjoint events. Events on distinct lines need not be independent; expectation is additive.

Ignoring the varying lengths gives the already strict bonus `e^{-13}/4`. The following Jensen computation improves it while retaining uniformity.

## 5. Exact length averaging and the improved coefficient

There are

\[
T_n=\sum_q\binom{L_q}{3}
=\sum_{h=2}^{\lfloor(n-1)/2\rfloor}(n-h)(h-1)(n-2h)
=\frac{n^4}{32}+O(n^3)
\]

potential triples. The middle formula independently counts minimum and maximum row difference `h`, the interior row, and the possible column offset. Thus excluding `O(n^3)` triples leaves the same leading term.

Scale the physical grid to the unit square and use the scaled `q=c-2r` coordinate. Its interval is `[-2,1]`. The possible row interval `I_q` has length

\[
L(q)=\begin{cases}(q+2)/2&-2\le q\le-1,\\1/2&-1\le q\le0,\\(1-q)/2&0\le q\le1.\end{cases}
\]

For a target point with row `r in I_q`, the sum of the lengths of its opposite diagonals is

\[
d(r,q)=2-|r+q|-|3r+q-1|.
\]

Its row integral is

\[
H(q)=\int_{I_q}d(r,q)\,dr
=\begin{cases}(q+2)(q+5)/6&-2\le q\le-1,\\-(q-1)(q+2)/3&-1\le q\le0,\\(q-4)(q-1)/6&0\le q\le1.\end{cases}
\]

The exact integrals, independently checked by polynomial integration and by integrating the absolute-value formula on a rational grid, are

\[
\int L^3=\frac3{16},\qquad
\int L^4=\frac7{80},\qquad
\int L^2H=\frac{187}{720}.
\]

For a potential target triple put `a_T=2S_T/n`, where `S_T` is the sum of its seven line lengths. Each row appears in `binom(L_q-1,2)` triples on its line. Consequently the arithmetic mean of `a_T` over all potential triples converges, by these elementary piecewise-polynomial Riemann sums, to

\[
\overline a
=2\frac{\int L^4+3\int L^2H}{\int L^3}
=\frac{416}{45}. \tag{7}
\]

The costs are bounded by `13+O(1/n)`. Removing `O(n^3)` of `Theta(n^4)` triples therefore changes the average by `o(1)`, uniformly in the fixed skeleton and row order. Jensen's inequality here is applied to the deterministic finite collection of costs, not to purportedly independent events:

\[
\frac1{|\mathcal E|}\sum_{T\in\mathcal E}e^{-a_T}
\ge\exp\!\left(-\frac1{|\mathcal E|}\sum_{T\in\mathcal E}a_T\right)
=e^{-416/45}-o(1).
\]

Summing (6) over the eight choices for every eligible triple yields

\[
\boxed{\liminf_{n\to\infty}\inf_{G,\rho}
\frac{\mathbb E_\sigma I_3(B)}n
\ge\delta_3:=\frac14e^{-416/45}
=0.0000241617729097398\ldots.} \tag{8}
\]

Uniformity of the `o(1)` in (6) matters: after multiplication by `8|E|/(n)_3=Theta(n)`, it contributes `o(n)`, as required. The factor `1/4` is `8` intended edge choices times the triple coefficient `1/32`.

Let the inherited coefficient be

\[
\gamma=1-5e^{-2}
+\frac{e^{-2}}{12}-\frac{e^{-4}}2+\frac{5e^{-6}}4+\frac{e^{-8}}6.
\]

Combining the supplier's uniform mean for `tau_2` with (1), (2), and (8) gives the strictly stronger statement

\[
\boxed{\liminf_{n\to\infty}\inf_{G,\rho}
\frac{\mathbb E_\sigma\tau_3(B)}n\ge\gamma+\delta_3.} \tag{9}
\]

The third direction changes no replacement-Lipschitz constant: if two equal-size boards differ by `s` original points, deleting those `s` points in addition to an optimum repair on the other board gives `|tau_3(B)-tau_3(B')|<=s`. A physical row or column transposition changes at most four points of a two-regular board. Applying the same transposition-exposure concentration argument as the supplier therefore gives, conditionally on `G,rho`,

\[
\mathbb P_\sigma(\tau_3\le k)
\le\exp\!\left[-\frac{(\mathbb E_\sigma\tau_3-k)_+^2}{8(n-1)}\right].
\]

For every fixed `kappa<gamma+delta_3`, the resulting uniform upper exponential rate for retaining all but `kappa n` original cells in any no-three-in-line output is at most

\[
-\frac{(\gamma+\delta_3-\kappa)^2}{8}.
\]

In particular the supplier's uncolored saturated-board mixture inherits the stricter no-three-in-line success rate `-(gamma+delta_3)^2/8`. This remains an exponential density bound within a superexponential-size class, not an extremal nonexistence theorem. No explicit finite threshold for the new asymptotic constant is claimed.

## 6. Finite universe, positive and hostile controls

The standalone producer uses only the Python standard library and explicit exceptions, so `python -O` retains every gate. It enumerates exactly the 2,137 simple two-regular boards on `n by n` grids for `n=2,3,4,5`: respectively 1, 6, 90, and 2,040 boards. All row column-pairs are listed lexicographically, with every column used exactly twice; there is no symmetry quotient or inherited acceptance filter. There are exactly 142 strict `tau_3>tau_2` cases at `n=5`, and no isolated singleton-opposite-diagonal event `I_3>0` in this complete universe.

The source independently compares exact triple-hitting deletion search, unit-capacity flow, the two-family Boolean min-cut value, and the three-family Boolean compiler. It checks (1) on every board in that universe. It does **not** enumerate all `n=6` boards. Instead it verifies the explicit two-regular positive control

\[
(01,13,05,23,24,45)
\]

with `(tau_2,tau_3,J_3,I_3)=(3,4,1,1)`. The isolated triple is `(0,1),(1,3),(2,5)` on `q=1`. Together with the complete smaller universe, this proves that `n=6` is the smallest square dimension supporting this stricter isolated event. The weaker safe-set bonus already occurs at `n=5`: `(01,23,04,12,34)` has `(tau_2,tau_3,J_3)=(0,2,2)`.

The `n=6` minimal isolated example is also a useful eligibility hostile: taking its displayed columns as source columns, its chosen target rows share a source column. The realized isolated event exists, but not all eight intended assignments are valid. Thus the uniform proof must use its explicit eligibility filter; it cannot infer eligibility from a realized event. The stronger filter is sufficient, not necessary.

The finite event bank uses cycle and disjoint-square skeletons, each at `n=6,8,10,12`, with natural row order and the order placing even native rows before odd native rows. These are 16 explicitly declared skeleton/row-order pairs. It checks every one of their 572 eligible triples and all 4,576 intended source-edge choices. For the first at most 16 eligible choices per pair, it constructs the actual residual forbidden matrix, checking all degree and length bounds: 144 full matrices. For the selected matrices at `n<=8`, it exhausts the residual column permutations, totaling 1,920 permutations, and verifies equivalence between matrix avoidance and the literal seven-line event. No finite lower probability bound is inferred from these controls.

Further gates verify every row and column transposition of the five named boards, the exact fractional primal/dual certificate in (5), the nonunimodular minor, both formulas for the potential-triple count for `1<=n<=80`, and the piecewise integral formulas at 181 rational `q` values. The producer passes **20,864 gates**. Its finite universe does not substitute for the analytical proof of (8) or (9).

Reproduction from the repository after installation:

```text
python 04-computation/continuing11_20260908_no3_three_directions.py
python -O 04-computation/continuing11_20260908_no3_three_directions.py
```

Outside the repository, invoke the same source at its current absolute path. The certificate is written beside an outside source, or into `05-knowledge/results` when the source is installed in `04-computation`.

## 7. Research consequence and paid stopping boundary

There is now a strictly stronger uniform repair obstruction, with a mechanism genuinely independent of the already charged original cells. The finite dual compiler often finds further deterministic credit when the universal safe-set charge vanishes. This suggests optimizing an old dual witness jointly with the third-family residual excess, with the triangle obstruction warning that this remains a lower bound.

The next decisive test is a different third slope or a fourth family with an explicitly disjoint original-cell charge; summing raw excesses is already refuted. The route should preserve the residual-cell sidecar and account for the number of forbidden row relations before asserting fixed-row uniformity. For a generic rational slope, the same finite set of forbidden companion-row equations is a plausible reusable mechanism, while the line-length integral and any symmetry reductions must be recomputed. No general-slope theorem is claimed here.

## Accepted independent audits

- [no3 three audit](continuing11_20260908_no3_three_audit.md).
- [no3 three directions analytic audit](continuing11_20260908_no3_three_directions_analytic_audit.md).

The parent accepted these reviews without mathematical repair. Source, output
and certificate bytes remain frozen; this filing changes status and routing only.
