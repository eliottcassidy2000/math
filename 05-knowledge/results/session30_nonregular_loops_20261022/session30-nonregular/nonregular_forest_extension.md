# Session30: leaf-decorated forest expansion beyond regular tournaments

## Status and inheritance

PROVED identities below concern the exact forest polynomial and individual resummed cluster classes. The complete first- and second-weight expressions are formal cluster identities, with finite computer-assisted local-kernel verification. Their interpretation as uniformly accurate asymptotic approximations to Hamiltonian counts is CONDITIONAL: the long-forest/infinite-cluster remainder is not bounded here. In particular HYP-9365 remains OPEN; the false universal regular-tournament limit e is not inherited. Later source-engine claims are not silently imported.

Closest proved mechanism: audited THM-4757's injective linear forests, partition Möbius inversion, and leaf-free one-/two-loop resummation. Corrected near miss: the source's path-tree penalty alone leaves a mean discrepancy around −0.50 at n=10. Missing sidecar: terminal leaves on cyclic cores. Hostile family: transitive tournaments have score deviations of order one, so the almost-regular power counting does not apply. A sharper hostile family *within* almost-regular tournaments is given in §7: deleted nonvertex-transitive regular blow-ups have a nonzero first-order lollipop contribution.

Notation: S=A−Aᵀ, M=S/n, δ=M1, κ=2^(n−1)H/n!, Φ(z)=Σ_L N_L(S)(z/n)^L. N_L sums signed directed linear forests with L arcs. H counts directed Hamiltonian vertex orders. For all tournaments,

κ=Σ_L N_L/(n)_L.

Let Q=zM and define

G=−Q²(I+Q²)^−1, H=Q²(I−Q²)^−1, O=Q(I−Q²)^−1,
d=diag(H−G), g=G1, h=Q(I+Q²)^−1 1.

The letter H used for the propagator is distinguished by context from the integer Hamiltonian count.

## 1. Exact forest polynomial and parity

Expand ∏_(j=1)^(n−1)(1+S_(σ_j,σ_(j+1))) over all permutations. A chosen signed forest with L arcs is contained in exactly (n−L)! permutations. This proves the displayed κ identity without regularity.

Reversing an odd-length path gives a sign-reversing involution, so only collections of even-length runs contribute. Partition Möbius inversion of the vertex injectivity constraints gives

log Φ(z)=Σ_connected (μ(π)/p!) Z_π(Q),

where p is the number of labelled runs, μ(π)=∏_B(−1)^(|B|−1)(|B|−1)!, and Z_π is the free homomorphism contraction on the quotient multigraph. This includes trees, unlike the regular specialization.

## 2. The correct almost-regular grading

For even n, almost-regularity means δ_i=±1/n and Σ_iδ_i=0. If a connected quotient graph has E arcs, V vertices and ℓ terminal vertices of degree one, define

w=E−V+ℓ=e+ℓ.

Sum all terminal labels first. Each contributes ±δ at its neighbour, of absolute value 1/n. Every remaining edge has absolute weight 1/n. Summing the remaining V−ℓ labels gives the rigorous fixed-graph bound

|Z_π(M)|≤n^(−w).

This argument does not require pseudorandomness or a fitted coefficient. It does not by itself bound the number of graphs or justify summing all weights.

Suppress degree-two vertices. Apart from the bare cycle (weight zero) and the bare path, reduced internal vertices have degree ≥3. If b is their number, then b≤2e+ℓ and total reduced vertices b+ℓ≤2w. Thus fixed weight has finitely many reduced cores.

Weight one has exactly five cores:
- three leaf-free two-loop cores: figure eight, theta, dumbbell;
- one one-leaf unicyclic core: a cycle with an attached path, called the lollipop;
- one two-leaf tree: a path.

Weight two has exactly 26 cores:
- 15 leaf-free three-loop cores;
- seven one-leaf bicyclic cores;
- three two-leaf unicyclic cores;
- one three-leaf tree (a Y).

The core counts and automorphism factors are independently generated in decorated_engine.py and agree with the separate producer engine. One cannot obtain the complete second weight by adding only a vacuum three-loop term and a path-tree cumulant.

## 3. Three closed leaf corrections

### 3.1 Path trees

On a path with 2k edges, a cut into p even runs contributes (−1)^(p−1)2^(p−1), with binomial(k−1,p−1) choices. Summing p gives (−1)^(k−1). Hence the entire path-tree contribution is

T(z)=1ᵀQ²(I+Q²)^−1 1
    =−z²δᵀ(I+z²M²)^−1δ.

### 3.2 Lollipop: a missing first-order term

Fix a reference orientation around a cycle of length a and outward along a tail of length b. Including the order-two core automorphism, the kernel is

4(−1)^(b/2+1), if a≡2 (mod 4) and b is even; zero otherwise.

Thus the complete lollipop contribution is

L(z)=−2 diag(H−G)ᵀG1=−2dᵀg.

Indeed Σ_(a≡2 mod4)Q^a=Q²(I−Q⁴)^−1=(H−G)/2, while the tail series is −4G. This term is identically zero for regular M1=0. It is generally of order 1/n, not uniformly of order 1/n², on almost-regular families.

The kernel is obtained by local half-edge matchings at the trivalent vertex and degree-two subdivision vertices. Degree-two path generating functions have only fourth-root-of-unity poles, so their coefficients are four-periodic. Explicit residue representatives are checked by core_kernel.py, without the previous source engine. The bivariate series is 4x²y²/((1−x⁴)(1+y²)).

### 3.3 Three-leaf tree

The complete Y-tree contribution, of weight two, is

Y(z)=−2Σ_v g_v h_v²−(2/3)Σ_v g_v³.

An analytic derivation: the signed chain sum for length a is F_a=2(−1)^(a/2−1) for positive even a, and zero for odd a. At the trivalent vertex the no-pair Möbius factor is 2; each of the three one-pair configurations has factor −1. If two arms a,b are joined, their unbroken-center sum is J_ab=(−1)^a F_(a+b)+F_aF_b. Consequently the kernel is

−F_aF_bF_c−Σ_cyclic (−1)^a F_(a+b)F_c.

Divide by the order-six automorphism group and sum arm lengths. Write r=R_i1=g+ih and its conjugate rbar=g−ih. Then Σ F_aQ^a1=−2g, and the two-arm term is −2r rbar. The result simplifies to the displayed Y formula.

## 4. Exact low-degree score identities

For every tournament, with δ=M1,

[z²]logΦ=(n−1)/n−||δ||²,

[z⁴]logΦ=(n−1)(15−7n)/(2n³)
 +(8n−13)||δ||²/n²+δᵀM²δ
 +2〈δ∘δ,Mδ〉−(1/2)Σ_iδ_i⁴.

For almost-regular tournaments this reduces to

[z⁴]logΦ=(n−1)(15−7n)/(2n³)+δᵀM²δ+(16n−31)/(2n³).

Thus the path-only correction misses an additional +8/n² at leading order in this coefficient. This is not a fitted correction. It is an exact coefficient identity.

cluster_shapes.py independently enumerates connected run partitions at L=2,4 and aggregates signed quotient graphs. The separate audit reproduced this identity on every labelled tournament through n=6 and 128 sampled order-seven tournaments, using direct permutation products rather than producer formulas.

## 5. Finite parameter-free engine

For each reduced multigraph, choose a matching of half-edges at each vertex. A vertex of degree d with p through-pairs contributes μ=(−1)^(d−p−1)(d−p−1)!. For every reduced edge, choose an uncut path or a path containing an internal junction. The resulting components are macro-runs. Reject closed runs and orient each open macro-run in both ways.

For endpoint signs s,t∈{+1,−1}, each reduced edge receives one of the exact transfer matrices

R_s=sQ(I−sQ)^−1,
C_(s,t)=−st Q²(I−Q²)(I−sQ)^−1(I−tQ)^−1(I+Q²)^−1.

The first sums uncut positive-length paths. The second sums paths with at least one internal junction. Contract the endpoint indices of the reduced graph, multiply its local Möbius factors, and divide by its multigraph automorphism order, including permutations of parallel edges and reversal of loops. This gives the exact resummed contribution of that core. All matrices are real for s,t=±1.

Independent implementation: decorated_engine.py. It was written from the above local rule, without importing the earlier or separately adapted engine. It reproduces the closed T,L,Y and audited vacuum X1 on random nonregular orders 4,5,8 to 2.3e−16. All 124 per-core comparisons (31 cores times four fixtures: transitive4, transitive6, R7, deleted R7) agree with the separate implementation to at most 5.33e−15. See engine_comparison.json. These numerical comparisons are verification, not a substitute for the local decomposition proof.

## 6. First and second formal predictions, and the open remainder

Let F(z)=exp X0(z), with

X0(z)=(1/2)log det(I−z²M²)−(1/2)log det(I+z²M²), ρ=F(1).

Let A1(z) be the sum of the five weight-one cores:

A1=X1_vacuum+T+L.

Let A2(z) be the sum of all 26 weight-two cores. Define D=z∂_z, P(D)=(D²−D)/2, Q(D)=(3D⁴−2D³−3D²+2D)/24. These coefficients come from

n^L/(n)_L=1+L(L−1)/(2n)+L(L−1)(3L²+L−2)/(24n²)+….

At z=1, the full first-weight formal predictor is

κ1/ρ=1+A1+(F^−1P(D)F)/n,
y1=−nA1−F^−1P(D)F.

The full second-weight formal predictor is

κ2/ρ=1+A1+(F^−1P(D)F)/n
 +A2+A1²/2+(F^−1P(D)(FA1))/n+(F^−1Q(D)F)/n²,
y2=n(1−κ2/ρ).

Every coefficient is prescribed. No regression is involved. Derivatives of T and L must be included in the mixed P(D)(FA1) term.

For implementation, if f=X0 and f_j=D^jf, then p0=(f_2+f_1²−f_1)/2 and

F^−1P(D)(FA1)=(D²A1)/2+(f_1−1/2)DA1+p0A1.

The Q term follows from the standard Bell polynomials in f_1,...,f_4.

Formal identities are not an unconditional asymptotic theorem. A proof of κ−κ2=O(n^−3), uniformly over almost-regular tournaments, still requires control of all omitted cluster weights and the falling-factorial-weighted long-forest tail. Bounds on a fixed core or analyticity near z=1 alone do not give that result: n^L/(n)_L grows substantially for L comparable to n. The source's uniform regular remainder was already open, and this work does not close it by relabelling a finite expansion.

## 7. A genuine leading-order obstruction to path-only predictions

Let B be the explicit regular order-seven skew matrix in lollipop_obstruction.json. Six vertices have one local spectral diagonal and one vertex has a distinct one. For odd t, replace every base vertex by the same rotational regular tournament on t vertices, then delete one vertex from base block a. The resulting order n=7t−1 tournament is almost-regular.

Writing G_B,H_B for the propagators of B/7 and d_B=diag(H_B−G_B), the lollipop contribution obeys

lim_(t→∞) nL(T_t)=14(G_Bd_B)_a.

Exact rational evaluation gives

12607/496782 for the exceptional block,
−12607/2980692 for each of the other six blocks.

The limits are nonzero and sum to zero. They certify that almost-regular score balance does not remove the need for vertex-resolved spectral data. A scalar irregularity norm cannot see which base block was deleted.

Proof details and the deletion error estimate are supplied separately in lollipop_limit_proof.md. The independent testing lane recovered the same rational vector on its differently labelled nonvertex-transitive base, and its nL scaling data approaches these limits. This obstruction concerns the missing formal first-weight term; it does not assume or prove the global asymptotic remainder for Hamiltonian counts.

## 8. Unconditional coefficientwise remainder theorem

For every fixed L, uniformly over all almost-regular tournaments of even order n,

[z^L](log Φ−X0−A1−A2)=O_L(n^−3).

For odd L both sides vanish. For L=2k>0, an explicit (deliberately loose) admissible constant is

C_L=Σ_(p=1)^k binomial(k−1,p−1) Bell(L+p)(L+p−1)!/p!.

Proof: there are binomial(k−1,p−1) ordered positive even run compositions, with q=L+p run positions. Partition Möbius coefficients have absolute value at most (q−1)!, and there are Bell(q) partitions. All connected graphs of weights zero, one and two are precisely accounted for by X0,A1,A2. Every remaining nonzero quotient has weight at least three, so its leaf-summed homomorphism contraction is bounded by n^−3. Sum the finite set of partitions. A connected graph consisting of one edge would need a separate two-leaf convention, but it never occurs here because odd runs cancel; self-loops carrying a single M diagonal entry vanish.

This is a genuine uniform theorem at each fixed coefficient, not merely tested evidence. It does not permit exchanging the limit with evaluation at z=1 or the falling-factorial summation. That distinction is the remaining analytic obligation.

## 9. Exact symbolic audit upgrade

Beyond the floating-point comparisons in §5, compare_symbolic.py independently rebuilt all 31 local kernel tensors and compared every rational coefficient after the exact basis reduction

C_(−,−)=−R_+−R_−−C_(+,+).

All 31 tensors are exactly equal. symbolic_comparison.json and symbolic_comparison.out record the per-core sizes and verdicts. This extends the verification to the formal tensor mechanism of the 15 vacuum three-loop cores as well as the newly decorated cores. It does not assert any uniform asymptotic error estimate.

prove_closed_kernels.py additionally establishes exact scalar rational identities for the entire path-length generating functions, using symbolic algebra rather than a finite residue sample:

path: x²/(1+x²),
lollipop: 4x²y²/((1−x⁴)(1+y²)),
Y: 2xyz(xyz+x+y+z)/(3(1+x²)(1+y²)(1+z²)).

Every symbolic difference from the closed derivation is exactly zero. These functions are generating kernels, not scalar eigenvalue models for an arbitrary branched graph: each edge variable is replaced by its corresponding matrix propagator and endpoint indices are contracted as specified by the graph.
