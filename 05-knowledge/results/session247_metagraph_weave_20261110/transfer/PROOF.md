# Charge-slice gauge cancellation and token-graph selectors

Session247, 2026-11-10. **PROVED conditional criteria; FINITE-EXACT physical example. Global Chan–Pak coverage remains OPEN.** No old theorem is re-audited here.

## Inheritance and the transferred idea

The inherited protected-ear formula is THM-4114, `01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`. For an actual deletion core C, H(C)=h, its actual two-path-cover coordinates a,W, and a cut S,

    F(S)=h+2 sum_(i in S) a_i - sum_({i,j} subset S) W_ij.

This is the Hamilton-path count after appending a sink and then a root beating S and losing to the sink. W is symmetric, with zero diagonal. All F(S) are odd.

The least-used metagraph sidecar is **conserved incidence margins, with actual port identity retained**. `05-knowledge/results/metagraph_reconstruction_20261023.md`, sections 1 and 4, distinguishes balanced transportation switches from scalar bucket margins and reconstructs actual ports with signed curl. `05-knowledge/results/three_sorted_metagraph_continuation_minimization_codex_S11.out` shows coarse observations can miss coherent phase even after continuation refinement. We transfer conservation of charges, not a claimed identification of metagraph vertices with cuts. The source-to-target map is: bucket occupancy -> selected-vertex occupancy; balanced switch -> token move; discarded row/column potentials -> an explicit slice gauge; physical port -> the literal labelled cut S. Lost information is all other charge slices and future continuation responses, which must remain sidecars.

The nearby THM-5122 shared-core cube theorem kills mixed interactions or bounds every Boolean edge. Here **large dense interactions are allowed** because conserved charges make their common part constant or linear. The general connected-unit-graph intermediate-value principle was already present in THM-4118 (`THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md`); the contribution here is a polynomially testable safe-token graph, its order-statistic certificate, and an exact gauge normal form on charge slices. No claim of inventing that general principle is made.

## 1. Exact slice gauge

Fix |S|=k. For arbitrary rational alpha_i, make the simultaneous coordinate change

    W'_ij=W_ij+alpha_i+alpha_j  (i != j),
    a'_i=a_i+(k-1)alpha_i/2.

Then F'(S)=F(S) for every k-subset. Indeed every selected alpha_i occurs k-1 times in the pair sum. This is an algebraic coordinate change of a response, **not** a claim that W',a' come from another tournament. In particular dense positive W can be removed without changing any observable on this slice.

If W_ij=alpha_i+alpha_j, put c_i=2a_i-(k-1)alpha_i. Then

    F(S)=h+sum_(i in S)c_i.

If c_i=c_0+2z_i with integer z_1<=...<=z_n and z_(i+1)-z_i<=1, the slice image is all odds between h+sum of the k smallest c_i and h+sum of the k largest c_i (provided these endpoints are odd, as automatic for physical responses). Move tokens along the ordered line; every move changes the sum by 0 or 2 and the k-token configuration graph is connected. Selecting a target means following any path between extremizing cuts until that target is reached.

The coefficient condition is sufficient, not necessary. The c_i need not be individually integral before an appropriate common offset is extracted.

## 2. Dense block gauge with multiple conserved charges

Partition the vertex set into nonempty blocks B_r. Fix k_r selected vertices in block r, and K=sum k_r. Suppose, for distinct i,j,

    W_ij=alpha_i+alpha_j+beta_(r(i),r(j)),

with beta symmetric, including within-block beta_rr. Then, exactly,

    F(S)=B(k)+sum_(i in S) [2a_i-(K-1)alpha_i],
    B(k)=h-sum_r beta_rr binom(k_r,2)
              -sum_(r<s) beta_rs k_r k_s.

Thus a fully dense interaction matrix can support a product of independent charge selectors. If within each block the resulting coefficients have a common residue modulo 2 and sorted gaps at most 2, each block supplies an integer interval in that residue, and their sum supplies every odd integer between the actual endpoints. Empty/full blocks contribute constants. This is more flexible than requiring disjoint exchanges to be W-orthogonal: the conserved charges absorb the dense background before selection. Approximate decompositions are handled by the safe-swap criterion below, not by silently deleting the residual matrix.

## 3. Exact robust safe-swap criterion

For p != q and 1<=k<=n-1, let d_j=W_qj-W_pj for j outside {p,q}. Sort these n-2 numbers as d_(1)<=...<=d_(n-2). Put m=k-1 and define empty sums as zero. The swap p->q is safe for **every** k-subset containing p but not q iff

    2(a_q-a_p)-sum_(largest m) d_j >= -2,
    2(a_q-a_p)-sum_(smallest m) d_j <= 2.                 (1)

Proof: write S=R union {p}, |R|=m. Direct subtraction gives

    F(R union {q})-F(R union {p})
       =2(a_q-a_p)-sum_(j in R)(W_qj-W_pj).              (2)

The extreme m-subset sums are precisely the two sorted sums. This also proves invariance of the test under the gauge in section 1: each d_j gains alpha_q-alpha_p, exactly cancelled by the linear term. Reverse swaps obey the same bound.

Let G_k have the core vertices as vertices and its safe pairs as edges. If G_k is connected, the entire k-cut image is a solid odd interval between its actual minimum and maximum. The reason is that the k-token graph of a connected graph is connected. A short proof uses a spanning tree: at a leaf, if two configurations disagree, move the nearest required token (or vacancy) along the unique tree path to settle the leaf, then delete that leaf; induction settles all vertices. Each move exchanges an occupied and unoccupied endpoint. Along this path (1) bounds changes by 2; oddness prevents skipping an intervening odd target.

This gives an O(n^3 log n) arithmetic certificate after a,W are known. It does not make finding endpoints polynomial-time for a general residual quadratic. Witness extrema or any two witness cuts suffice to certify the interval between their counts. For a disconnected G_k, retain the token charge in each connected component; the same proof fills the interval on each fixed component-occupancy sector. It is wrong to merge those sectors merely because their count lattices agree.

The criterion is not just a rank or gcd test: it gives actual permitted moves and a constructive target selector. It also does not require every Boolean cube edge to be small.

## 4. Literal positive witness: a five-vertex core

Take the actual tournament row masks

    [24,17,3,6,12].

Its coordinates are h=15, a=(21,21,21,21,21), and

    W = [0 10  8  8 10
        10  0 10  8  8
         8 10  0 10  8
         8  8 10  0 10
        10  8  8 10  0].

On the k=2 slice, F({i,j})=99-W_ij. Every pair swap is safe, G_2=K5, and the physical image is {89,91}. For example {0,1} gives 89 and {0,2} gives 91. The corresponding literal seven-vertex closures are strong tournaments. The dense mixed Boolean differences are -8 or -10; the one-bit criterion from the full cube is emphatically inapplicable. Subtracting the potential alpha_i=4 removes the common W background 8 and leaves a residual 2 on a five-cycle. The transformed a'_i=19 uses the subtraction gauge; it is only a coordinate representation.

This is a deliberately small witness, not a new numerical coverage bound. The saved check recomputes the ten closures directly from permutations and the coordinates from path-cover tables. It also checks a deliberate unsafe transitive-core edge, so (1) is not vacuous.

## 5. An arbitrarily wide algebraic model, with the precise realization gap

For n>=2, 1<=k<n, and an even positive integer L, set

    W_ij=2L for i!=j,
    a_i=L(k-1)/2+i,      i=0,...,n-1,
    h odd.

Then on |S|=k, F(S)=h+2 sum_(i in S)i. The image is exactly

    [h+k(k-1), h+k(2n-k-1)] intersect (2Z+1),

with width 2k(n-k), although the Boolean interactions have arbitrary magnitude 2L. Adjacent labels are safe and form a spanning path. This proves that the dense-interaction ceiling is absent for conserved-charge selectors. **It does not assert that these arrays are tournament path-cover arrays.**

To advance Chan–Pak with this architecture, one needs actual growing cores whose gauge-reduced rows satisfy (1) on a connected graph, with physical endpoints widening and overlapping across cores/charges. The C5 witness proves the mechanism is physically nonempty. It does not establish that growing realization theorem, absolute endpoint control, or the missing target descent.

## 6. Research leverage

A search should prioritize **small variation of row differences after a charge gauge**, not small entries, small rank, scalar count gcd, or an orthogonal exchange basis. Highly symmetric cores can supply a large removable dense background; carefully broken symmetry must supply nonzero interval width. Perfect symmetry alone gives a singleton (for example the doubly-regular seven-vertex response has constant off-diagonal W). This explains a useful engineering tradeoff without declaring every symmetric family sufficient.

Cheap next test: from a proposed core compute a,W, certify a spanning tree in G_k, and exhibit two far-apart cuts. No audit of the historical metagraph inventory is needed to attempt that construction.
