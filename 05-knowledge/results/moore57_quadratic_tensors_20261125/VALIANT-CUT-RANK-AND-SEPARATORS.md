# Moore charts as graph state sums: cuts, width, and arithmetic price

2026-11-26. **PROVED conditional on the graph, not on a choice of automorphism.** The degree-57 graph is not constructed or ruled out. Exact ordered and nonnegative bounds below are independent of the preprint's Euclidean pricing theorem; the unrestricted circuit and ABP bounds are implications of its stated framework.

## 1. What is being transferred

The November 24 [Paths, trees and grids preprint](https://github.com/openai/math/blob/92b516ef8e2abf66c3db7aa31438d3010bc7956f/preprints/Paths-trees-and-grids-separate-Valiants-hierarchy-November-24-2026/paper.pdf) studies the polynomial Z(G,q) formed by labeling every edge of G with one of q states. Each vertex has an independent table indexed by its incident labels. Sum the product of one table entry per vertex over all edge labelings. The coefficient tensor is a product of equality wires, one per edge. The paper combines this fact with slot decompositions for circuits and branching programs and with [Euclidean tensor pricing](https://github.com/openai/math/blob/92b516ef8e2abf66c3db7aa31438d3010bc7956f/preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/paper.pdf). That additional tool addresses summands that use *different* cuts. This note uses the stated pricing theorem, not a replay of the Lean development or a new proof of the hierarchy separation.

Our graph inputs are the [Moore spectral packet](FROZEN-SPECTRUM-AND-ANCHORED-ACTION.md), [Wilson route equations](WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md), and the [shifted mixing bound from the path work](../moore57_graph_spine_tiles_20261124/moore57-hamilton-path.md). The new middle step is a cut bound that sees integral per-vertex contact counts and, if known, integral occupancy of each physical branch. The graph state sums are not the Hamiltonian-path counts of the tournaments encoding those graphs. A specialization can lose all the rank: setting all table variables to one yields just q raised to the number of edges.

## 2. Cut ranks and integral Moore cut certificates

Let A be the adjacency of a d-regular graph on n vertices with all centered eigenvalues in [t,r]. For S of size a, put x=1_S−a1/n and let B(S) be the number of crossing edges. Then

    B(S) = (d−r) a(1−a/n) + xᵀ(rI−A)x.
    ‖(A−rI)x‖² ≤ (r−t) xᵀ(rI−A)x.

The first equality is the Laplacian quadratic form, and the inequality applies to each centered eigenvalue separately. For a full Moore graph, r,t are its two centered eigenvalues, so the second line is equality. Entries of (A−rI)x are integers minus (d−r)a/n when r is integral. Set ρ equal to the fractional part of that rational. Every such entry z satisfies z²≥(1−2ρ)z+ρ(1−ρ), since (z+ρ)(z+ρ−1)≥0. Their sum is zero. Therefore

    B(S) ≥ (d−r) a(1−a/n) + nρ(1−ρ)/(r−t).

Round up to the parity d·a, since the sum of internal degrees counts each internal edge twice. In degree 57, ρ is the fractional part of a/65. This improves the balanced lower bound from 36,123 to **36,170** without assuming additional symmetry. For a half cut the bound is **40,625**, with equality possible only if the cut is equitable: every site has 32 same-side and 25 opposite-side neighbors. No assertion is made that an equitable half cut exists.

The flattening of the coefficient tensor of Z(G,q) on any declared split S has rank exactly q^B(S) over every field: each crossing wire is an identity matrix of rank q and each internal wire has rank one as a factor of this flattening. Thus every fixed-order matrix-product compiler has width at least q^40,625 in the middle of the full graph. Every binary tree tensor network with one vertex table per leaf has a tree edge splitting the leaves into parts both of size ≥n/3 (take a centroid in the unrooted degree-at-most-three tree and its largest branch). The corresponding bond dimension is at least q^36,170. The degree-57 graph has size 3250, not a growing family in graph order here: q is the free parameter. These are exact architectural rank obstructions, not general arithmetic lower bounds by themselves.

## 3. Same calculation with branch occupancies and three frozen stars

On the k retained branches of a rooted Moore chart, there are L sheets per branch; distinct branches are joined by one perfect matching, and no sheet has a same-branch edge. For degree 57 use L=56 and k=57, or k=54 after removing the three frozen leaf branches. The induced graph is (k−1)-regular. On the branch-constant space its centered eigenvalue is −1; the sheet-centered space is invariant. Its adjacency spectrum is in [−8,7], by extending each centered vector by zero to the full graph and bounding its Rayleigh quotient. In the unfrozen leaf graph this subspace in fact retains just the two roots, since extending by zero also intertwines adjacency (the heads see zero branch sums).

For each cut S let a_i be the count of selected sheets in branch i, with sum a and n=kL. Decompose 1_S into its branch means and within-branch deviations. The extra quadratic term from the means is (r+1)/L times the squared deviations of the a_i from a/k. In the sheet-centered component the residual (A−rI) has entries equal to integers minus (a−(r+1)a_i)/L in branch i, because its branch sums are zero and the matching block average in branch i is (a−a_i)/L. Applying the previous integer-square inequality separately in each branch gives

    B(S) ≥ (k−1−r) a(1−a/n)
           + (r+1)/L · Σ_i (a_i−a/k)²
           + L/(r−t) · Σ_i ρ_i(1−ρ_i),

where ρ_i is the fractional part of (a−(r+1)a_i)/L. Again round up to parity (k−1)a, and take the maximum with the global bound of Section 2. For a union of j whole branches the branch variance bound is exact: B=L j(k−j). For example the 28/29 split of the unfrozen chart has 45,472 crossing links, not just the general half-cut floor 39,156.

With only a fixed, the minimizing occupancy vector for the branch-aware bound is balanced to within one. Indeed the separable summand, dropping constants, is (r+1)u²/L+(L/(r−t))·ρ(u)(1−ρ(u)). Its second derivative on each smooth part is 2(r+1)/L·[1−(r+1)/(r−t)], which is positive for the Moore intervals used here, and the first derivative jumps upwards at a wrap of ρ. Discrete marginal prices are therefore nondecreasing. Pairwise transfers balance an unconstrained profile; with lower/upper occupancy caps the same fact gives the exact greedy minimum of the bound under those caps by assigning each free sheet to the smallest available marginal increment. This does not minimize actual cuts; it minimizes the certified lower bound on them.

| Induced graph (all conditional) | sites | edges | circuit-balanced cut floor | half-cut floor |
| --- | ---: | ---: | ---: | ---: |
| Full degree-57 Moore graph | 3250 | 92,625 | 36,170 | 40,625 |
| Rooted leaf graph | 3192 | 89,376 | 34,806 | 39,156 |
| Leaf exterior after three frozen branches | 3024 | 80,136 | 30,958 | 34,776 |

The balanced cut range in the table is n/3<a≤2n/3, as in the source circuit decomposition. Here branch information sharpens the bounds for individual cuts and capacity-restricted families, rather than raising the combined worst-case figures in the last two rows. For example, the 54 exterior branches' split into 27 whole branches on either side has exactly 40,824 crossing edges.

**Search use.** In a proposed partial matching chart, take a fixed S and compute an upper bound on eventual B(S) by independently maximizing crossing entries in each allowed bipartite link (a linear assignment with forbidden positions and pinned contacts). Feasible completions are among those choices, even though the independent optima need not be mutually compatible. If the sum falls below the integer branch-aware bound, reject the node. A contact is forced if forbidding it would yield this strict failure. Taking cuts with uneven branch counts is especially useful; the variance charge explicitly credits them. Passing these tests is not a consistency or completion certificate.

## 4. Large separators cannot be hidden by grouping the sheets

The earlier centered mixing argument works with any d,r,t as above. Set μ=(r+t)/2 and σ=(r−t)/2. Disjoint anticomplete sets of sizes a,b, with remainder c=n−a−b, obey

    (d−μ)² ab ≤ σ²(n−a)(n−b),
    ab ≤ Kc,   K = nσ² / ((d−μ)²−σ²).

For the full Moore graph this is the old mixing inequality written as ab≤(225/4)c. The new consequence is a numerical separator and treewidth bound. A balanced separator C means each component after deleting C has at most n/2 vertices. Write m=n−c for the remaining size. If c<n/4, then every component has size at most 2m/3. One can group components into two anticomplete unions with a,b between m/3 and 2m/3: take a component in that range if possible, or accumulate components smaller than m/3 until crossing m/3. Then ab≥2m²/9>n²/8>Kc for K<n/2. For n/4≤c<n/2 the same grouping works with both sizes in [m−n/2,n/2]: use a component of size ≥m/3 if one exists, otherwise accumulate to a size in [m/3,2m/3]⊆[m−n/2,n/2]. Now ab≥(n/2)(n/2−c). Thus in all cases

    c ≥ ceil((n/2)² / (n/2+K)).

Every tree decomposition of width w supplies a balanced separator of size at most w+1: assign each vertex to one of its bags, use a weighted centroid, and delete that bag. A remaining component stays in one branch because its vertices' bag sets are connected and edges share bags. It follows that the full Moore graph has **treewidth ≥1,570**, its leaf graph ≥1,540, and the three-frozen-branch exterior ≥1,453. These are conditional lower bounds, not exact widths. Small labelled encodings of the Moore graph are not contradicted: a decoder can use global algebra rather than separate the actual sheets into small communicating parts.

## 5. Mixed-cut algebraic sums: an independent monotone test and a priced transfer

**Nonnegative case.** Over the reals, suppose the coefficient tensor is a sum of L nonnegative tensors whose ranks are at most R across possibly different cuts, each crossing at least b edges. Each summand must vanish outside the graph support because there is no cancellation. In its chosen cut flattening, sort by the crossing word. The nonzero blocks lie on the diagonal, so rank ≤R permits at most R nonzero crossing-word blocks. In the target, each such word has exactly q^(e−B) full consistent labelings, where e is the edge count and B is that cut's size. Counting the q^e support entries gives **LR≥q^b**, independent of where each summand cuts. For a nonnegative ABP the source's prescribed-slot decomposition at a half cut gives R≤s (number of ABP vertices) and L≤2^n. For q=2 that gives at least 2^37,375 ABP vertices for the full symbolic tensor and 2^31,752 for the frozen exterior's symbolic tensor. This argument does not cover cancellation or border computations.

**Allow cancellations, using the source price theorem.** Put the two endpoints of each edge on antipodal unit sites of its own axis. Use natural-log prices, w=log q. The graph tensor has price 2ew. For a cut C let a_C be the sum of the C-side terminal vectors; then ‖a_C‖=√B. The preprint's move of one supply site gives a loss estimate; moving both sites gives a slightly sharper finite-graph statement. If their relative displacement is λ a_C/‖a_C‖, choose opposite individual shifts in inverse proportion to the two sides' volumes (numbers of incident terminals). The linear gain is λ w√B. If the volumes are v and 2e−v the quadratic cost is w λ² v(2e−v)/(4e), at most ew λ²/4. The rank wire costs λ log R. Minimizing this quadratic estimate in λ shows that a rank-R term on any cut with B≥b costs at most

    2ew − (w/e)·(max(0, √b − log(R)/w))² .

The same terminals serve all summands. Let Δ be their normalized Gaussian mean width as in the preprint (Δ≤2√(π log(2e)), by the Gaussian exponential-moment bound). Adding L summands costs at most Δ log L beyond the most expensive one, and the price sublevel is Zariski closed. Hence for membership in the closure of these sums one needs

    log R ≥ w√b − √(e Δ w log L).

This finite-graph version does not need the paper's growing-grid sparsity condition e≤2n: that condition controlled the summation penalty in its uniform asymptotics. For graphs in the table n,e are fixed and q→∞. Substituting the preprint's Lemma 5.2 (R≤s, prescribed half-slot cut) and Lemma 4.2 (R≤O((s+1)(n+1)²), balanced slot cut), with L≤2^n, gives the following asymptotic degree-component border lower bounds over any field, under its pricing theorem:

| Graph | Circuit size exponent of q (minus o(1)) | ABP size exponent of q (minus o(1)) |
| --- | ---: | ---: |
| Full Moore | √36170 ≈ 190.184 | √40625 ≈ 201.556 |
| Three-frozen exterior | √30958 ≈ 175.949 | √34776 ≈ 186.483 |

No unrestricted q=2 lower bound is claimed from pricing; its summation penalty can dominate there. The independent q=2 results above concern fixed architecture or nonnegative sums. None of these statements restricts an analytic proof deciding whether the degree-57 graph exists.

Run **python3 moore_state_sum_bounds.py** to replay the rational floors and separator constants, exhaustive Petersen cuts, and 990 exact Hoffman–Singleton tests (including full leaf and three-frozen branch controls).
