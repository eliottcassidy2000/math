# Three repo-derived bridges: nonlinear gluing, critical characters, and finite moments

Date: 2026-10-19. Status: PROVED in the stated abstract scopes, with FINITE-EXACT companion controls. This is a research artifact, not a claim of world-first priority. The general ingredients have classical antecedents; the explicit transfers, adversaries, and repair calculations below are the session's repository contribution.

Repository baseline inspected: f96ed29dfbc8. Inherited computations were not assumed to have been replayed; the literal theorem data used below were independently reconstructed.

## Results first

1. **A nonlinear conductor obstruction invisible to the tangent matroid.** For every r≥6 there are two rational arrangements of r lines with distinct tangent directions through the origin in affine three-space, with the same linear tangent matroid U(3,r), whose seminormal defect lengths differ by floor((r−4)/2). Their complete Hilbert-rank and branch-conductor formulas are proved below. At r=6 there is a complete classification: the defects are 4 and 3, and conductor exponents are respectively (3,3,3,3,3,3) and (2,2,2,2,2,2), according as the six projective directions lie on a conic or not. Six is the minimal uniform-rank-three boundary. The missing information is quadratic and higher evaluation on the projective directions.
2. **An exact price for turning the repo's origin bit into a critical-group character.** On the 22-edge core of THM-3260/3278, four edge duplications are necessary and sufficient, and there are exactly two optimum choices. Requiring invariance under the intrinsic swap (17 21) raises the minimum to six. There are exactly 27 invariant six-edge choices. This changes the integral graph, rather than claiming that the original forbidden character already exists.
3. **A bounded-run replacement for the one-ramp decoder.** A finite subset of the integers with at most k nonempty consecutive runs is uniquely recoverable from its first 2k power moments, by a triangular finite-difference transform followed by an exact rational decoder. Zero endpoints are retained by homogeneous endpoint polynomials. The moment budget is sharp for the corresponding continuous interval problem, and explicitly sharp for discrete k=1 and k=2. No all-k discrete lower bound is asserted.

The themes were recovered independently from the repository's critical-group, singularity, and graceful-gap files. None is a rephrasing of the current Collatz/Gilbreath/tournament/period targets.

## 1. Inheritance and live concept board

### A. Conductor and nonlinear tangent data

- Closest mechanism: THM-4067, `01-canon/theorems/THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction.md`, identifies the gap between graph-continuous functions and ambient functions as the seminormal defect. THM-4381, `01-canon/theorems/THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification.md`, isolates a missing first-jet direction at one ordinary plane triple point.
- Canonical hostile/corrected near miss: THM-4067's figure-eight has zero graph periods without ambient descent; ordinary plane multiple points are not automatically seminormal. Graph incidence does not see contact and jet congruences.
- Least-used useful sidecar: the actual degree-by-degree branch-evaluation spaces, rather than just the number of branches, first-jet rank, or tangent dependence matroid.
- Cheap decisive test: six tangent directions on a conic versus five on that conic and one displaced; compute the degree-two evaluation rank.

### B. Critical group and a selector cut

- Closest mechanism: THM-3273, `01-canon/theorems/THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary.md`, computes the cyclic critical group of the core. THM-3278, `01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md`, identifies the selector origin bit with a bipartition cut and proves it is not an original critical-group character.
- Hostile: row 11 has origin bit 1 and even critical coordinate. The first failed implication would be “intrinsic cut implies harmonic character.”
- Useful sidecar: divergence of the cut modulo two, equivalently the degree-parity syndrome; keep the intrinsic involution as a separate constraint.
- Cheap decisive test: compute the eight odd-degree vertices and look for a matching on them.

### C. Moments and ordered masks

- Closest mechanism: THM-2793, `01-canon/theorems/THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction.md`, reconstructs a single interval from length and its first moment. THM-2534, `01-canon/theorems/THM-2534-positive-kakeya-boundary-transform-and-crofton-reconstruction.md`, distinguishes terminal-support recovery from the lost lengths of runs.
- Hostile: a run-terminal set alone forgets run lengths. Three moments already collide for two discrete runs, as shown below.
- Useful sidecar: the signed entrance/exit boundary, together with its full rational endpoint ratio.
- Cheap decisive test: construct two two-run sets having equal cardinality, sum, and squared sum.

The live concepts are: graph equalizer; graded branch evaluations; critical-group character; degree-parity repair; bounded-run moment transform. Their common methodological lesson is precise, but does not identify their mathematical types: one must retain the exact obstruction that the chosen quotient discards. This is not a common numerical “entropy invariant.”

## 2. All-line-arrangement formula

Work over an algebraically closed field k of characteristic zero. Let v_1,…,v_r be nonzero, pairwise nonproportional vectors in k^d. Let X be the reduced union of their lines through zero. Its coordinate ring A embeds in

B = direct sum over i of k[t_i],

by x_j ↦ (v_{1j}t_1,…,v_{rj}t_r). Put

E_q = image of the homogeneous degree-q evaluation map F ↦ (F(v_1),…,F(v_r)),
h_q = dim E_q.

Then the degree-q part of A is exactly E_q t^q. The graph value-equalizer is

S = { (f_i) in B : f_1(0)=…=f_r(0) }.

Consequently

length_k(S/A) = sum over q≥1 of (r−h_q),
length_k(B/A) = (r−1) + sum over q≥1 of (r−h_q).                 (2.1)

These sums are finite. Indeed, for each i and each j≠i choose a linear form vanishing at v_j but not at v_i. Their product is a degree r−1 separator. Thus E_(r−1)=k^r. A linear form nonzero at every v_i exists because k is infinite; multiplication by its invertible diagonal value vector shows E_q=k^r for every q≥r−1.

For clarity, S is the seminormalization in this situation. It is finite over A by the preceding conductor bound. It changes neither residue fields nor the point set: away from the origin the lines are disjoint and already normal, and the r origins remain a single point. It is seminormal: if an element of B has square and cube in S, its branch constants have identical squares and cubes, hence are identical. This is the standard seminormal description by equality of values on normalization fibres. Formula (2.1), however, also follows directly as a statement about the explicitly defined value-equalizer, without needing its terminology.

### Branchwise conductor, not just its total length

Rescale each v_i so that a fixed linear form takes value 1 on every vector. Let

d_i = least q for which the coordinate vector e_i belongs to E_q.

Then

(A:B) = direct sum over i of t_i^(d_i) k[t_i].                 (2.2)

To prove this, a function supported on branch i and homogeneous of degree q lies in A exactly when e_i∈E_q. Once a separator occurs, multiplication by the chosen linear form gives separators at every later degree. An ideal of the product ring B decomposes branchwise, giving exactly (2.2). The same formulas hold after completion at the common origin, replacing polynomials by formal series; the quotients have finite length.

These statements are specific to unions of straight lines. For arbitrary curved branches, higher jets can carry further information not determined by tangent directions.

## 3. Same uniform matroid, unboundedly different conductor defects

Fix r≥6. Consider the two configurations in k^3:

C_r: v_t=(1,t,t²), for t=0,…,r−1;
D_r: v_t=(1,t,t²), for t=0,…,r−2, and w=(1,r−1,(r−1)²+1).

### 3.1 The linear matroids agree

Every triple from C_r has nonzero Vandermonde determinant. For D_r, triples of old points have the same property. A projective line through two old points with parameters a<b has affine equation z=(a+b)y−ab. At y=R=r−1, its value is

R²−(R−a)(R−b) < R² < R²+1.

It therefore misses w. Thus every triple in D_r is independent too. Both configurations realize the same uniform matroid U(3,r), with the same r smooth branches, common origin, and pairwise scheme-theoretic intersection equal to the reduced origin. Their actual numeric tangent configurations are different; the claim concerns the loss under passage to the linear dependence matroid, not equality of the full tangent-coordinate matrices.

### 3.2 Exact Hilbert ranks

On C_r, degree-q homogeneous monomials restrict to polynomials in t of every degree 0,…,2q. Therefore

h_q(C_r)=min(r,2q+1).                                         (3.1)

On D_r, h_0=1 and h_1=3. For q≥2, restriction to the first r−1 points has rank min(r−1,2q+1). The degree-q form

x^(q−2)(xz−y²)

vanishes on every old point and has value 1 at w. It supplies the independent last-coordinate direction. Hence

h_q(D_r)=min(r,2q+2), for q≥2.                               (3.2)

No generic-position assertion is needed.

### 3.3 Defect formulas

Substituting (3.1) and (3.2) in (2.1) gives

length(S_C/A_C) = floor((r−2)²/4),
length(S_D/A_D) = floor((r−2)²/4) − floor((r−4)/2).             (3.3)

The difference is unbounded as r grows. Thus neither the branch-incidence model nor the complete linear tangent matroid determines the seminormal defect, even within this explicit rational family of smooth branches whose pairwise scheme-theoretic intersections are the reduced origin.

The corresponding normalization defects are obtained by adding r−1. These two lengths must not be conflated.

### 3.4 Exact conductor exponents

For C_r, a degree-q separator restricts to a polynomial of degree at most 2q that vanishes at r−1 distinct parameters. Such a separator exists exactly when 2q≥r−1. Therefore every branch has

d_i(C_r)=ceil((r−1)/2).                                      (3.4)

For D_r, the displaced branch has separator xz−y² and cannot have a linear separator because the old directions span k^3. Thus d_w=2.

For an old branch, restriction to the old conic must vanish at r−2 other old parameters, so q≥ceil((r−2)/2). Conversely choose the required polynomial on the old conic, lift it to a homogeneous form, then subtract an appropriate multiple of x^(q−2)(xz−y²) to make its value at w zero. Since r≥6, q≥2 and this correction is available. Hence

d_i(D_r)=ceil((r−2)/2) on old branches, and d_w(D_r)=2.        (3.5)

The full coordinate-space sidecar E_1 does determine all E_q by coordinatewise products. The failure is precisely that the abstract matroid retains zero/nonzero minors but not these numerical products.

### 3.5 The six-line witness in one calculation

For r=6:

- C_6 has Hilbert ranks 1,3,5,6,…, seminormal defect 4, normalization defect 9, and conductor exponent 3 on all six branches.
- D_6 has ranks 1,3,6,6,…, seminormal defect 3, normalization defect 8, and conductor exponent 2 on all six branches.

On C_6, the vector (−1,5,−10,10,−5,1) annihilates all quadratic evaluations, because it is a fifth finite difference on polynomials of degree at most four. A tuple supported on the last branch with coefficient t_5² violates that relation and represents a missing quadratic class. On D_6 the same tuple is realized by xz−y², so that class disappears.

### 3.6 Six is the first possible uniform-rank-three boundary

For 3≤r≤5, any realization of U(3,r) has h_1=3 and h_2=r. To prove the latter, fix one of the projective points. Partition the other at most four points into two groups of at most two points. For a two-point group, take its joining line; for a singleton, choose a line through it avoiding the fixed point. No three points are collinear, so each chosen line avoids the fixed point. Their product is a quadratic separator. Hence E_2=k^r.

It follows that the seminormal defect is forced to be r−3 for r=3,4,5. At r=3 every conductor exponent is one. At r=4,5 every exponent is two: a linear separator would vanish at three independent remaining directions, which is impossible. Thus six branches are the first possible value at which the uniform rank-three tangent matroid fails to determine either the defect or the branchwise conductor exponents, and C_6,D_6 attain that boundary.

### 3.7 Complete six-direction classification

Let six projective directions in P²(k) realize U(3,6), over the same algebraically closed characteristic-zero field. Every five have independent quadratic evaluations by the separator proof in Section 3.6. The six-by-six quadratic evaluation matrix therefore has rank either five or six. Its rank is five exactly when a nonzero homogeneous quadratic vanishes at all six directions, that is, when they lie on a conic.

Such a conic is irreducible and smooth. A reducible conic is two projective lines, with repetition allowed, but no line contains three of the six points; two lines can cover at most four of them. Over an algebraically closed field a singular plane conic is reducible, so singularity is excluded too.

If the quadratic rank is five, its one-dimensional left annihilator has every coefficient nonzero: a zero coefficient would give a linear dependence among the other five evaluation rows. Hence no coordinate vector e_i belongs to E_2. Cubic separators do exist: cover the other five directions by two joining lines and one additional line, all avoiding the selected direction, and multiply their equations. Thus E_3=k^6 and every conductor exponent is exactly three. If the quadratic rank is six, every exponent is exactly two, since the other five directions span k³ and prevent a linear separator.

Therefore the full classification is:

- six directions on a conic: ranks 1,3,5,6,…; seminormal defect 4; normalization defect 9; conductor exponents all 3;
- six directions not on a conic: ranks 1,3,6,6,…; seminormal defect 3; normalization defect 8; conductor exponents all 2.

This exhausts straight-line cones with tangent matroid U(3,6). It does not classify general six-branch curve singularities with the same tangent matroid.

**Connection contract.** Source: the cone of projective branch directions. Target: graph-equalizer quotient and conductor ideal. Map: homogeneous evaluation followed by forgetting to branch incidence or the linear tangent matroid. Preserved: branch count, smoothness, reduced pairwise intersection, and linear dependence data. Lost: higher evaluation ranks and branch separator degrees. Sidecar: E_q, or equivalently its Veronese evaluation matrices with coordinates. Decisive test: the r=6 quadratic rank 5 versus 6. No global manifold filling, S6 complex structure, Keller map, or Jacobian conjecture consequence follows.

## 4. When a cut is a critical-group character

Let G be a connected finite graph with positive integer edge weights w, interpreted as parallel-edge multiplicities. Let L_w be its weighted Laplacian, and let c:V→Z/m have c(root)=0. Pair a degree-zero divisor D with c by

χ_c(D)=sum_v c(v)D(v) mod m.

This descends to a character Jac(G_w)→Z/m if and only if

L_w c = 0 mod m.                                             (4.1)

Indeed, principal divisors are L_w z and the pairing with every such divisor vanishes exactly when c^T L_w=0; symmetry of the Laplacian gives (4.1). This is the usual harmonic-character description.

If G is bipartite and c is its 0/1 bipartition, then for m=2,

(L_w c)(v)=weighted_degree(v) mod 2.                         (4.2)

Thus the origin cut becomes a genuine order-two character exactly when all weighted degrees are even. Being a coboundary on edges is not sufficient for being harmonic on vertices.

### 4.1 Minimum duplication is the parity-join problem

Start with unit weights, let T be the odd-degree vertices, and allow additional copies only of existing edges, at unit cost per copy. In a minimum solution every edge receives either zero or one additional copy: removing two copies preserves all degree parities and reduces cost. The selected extra-edge set F must satisfy

odd_boundary(F)=T.                                           (4.3)

Thus the optimization is precisely a minimum T-join in the existing graph. Requiring invariance under a graph automorphism adds that requirement to F. This is a structural connection to the classical parity-join/Eulerization problem, not an identification of the original and repaired critical groups.

### 4.2 Exact core and two optimal repairs

Reconstruct the core as the disjoint union of the two trees printed in THM-3269, `01-canon/theorems/THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization.md`:

T = {(2,13),(2,17),(2,19),(2,21),(3,11),(3,16),(7,22),(10,11),(16,17),(18,19),(21,22)},
T'={(2,10),(3,22),(7,18),(10,16),(10,22),(11,13),(13,18),(13,22),(16,21),(17,22),(19,22)}.

The bipartition carrying the origin bit is

{2,11,16,18,22} | {3,7,10,13,17,19,21}.

Its odd-degree vertices are

T_odd={2,3,11,17,18,19,21,22}.                                (4.4)

Every extra edge can repair at most two of these vertices, so the cost is at least four. At cost four, every selected edge must lie in a perfect matching on T_odd. Vertex 11 forces (3,11); vertex 18 then forces (18,19). On the remaining vertices {2,17,21,22} the induced graph is K_(2,2), giving precisely two choices:

F_1={(2,17),(3,11),(18,19),(21,22)},
F_2={(2,21),(3,11),(18,19),(17,22)}.                         (4.5)

Thus four is the exact minimum and these are all minimizers.

### 4.3 An exact symmetry premium

The intrinsic involution σ=(17 21) exchanges F_1 and F_2, so no cost-four repair is invariant. Every admissible F has even cardinality: sum its degree parities over the first bipartition class; that class contains four vertices of T_odd, so |F|≡4≡0 mod 2. Hence an invariant repair costs at least six.

The invariant six-edge set

F_sym={(2,10),(2,17),(2,21),(3,11),(10,22),(18,19)}            (4.6)

has boundary T_odd. Therefore six is the exact invariant minimum, a symmetry premium of two edges. Exhaustive enumeration of all subsets of the 22 core edges through cardinality six finds exactly 27 invariant minimizers. The bound and the witness are proved structurally; the number 27 is a FINITE-EXACT census.

### 4.4 Integral check and boundary

The exact reduced-Laplacian Smith groups are:

- original: C_74748, matching the inherited theorem;
- repair F_1: C_2 ⊕ C_4 ⊕ C_49156, order 393248;
- repair F_sym: C_4 ⊕ C_4 ⊕ C_46620, order 745920.

The character c is harmonic in both repaired graphs and not in the original graph. These changed groups make the boundary explicit: duplication changes the integral lattice. The analytic clutch strengths of THM-3269 are not the integer weights used here. Nothing proves that these extra copies are lawful physical response transitions or that the repair solves the selector problem dynamically.

**Connection contract.** Source: the core bipartition cut. Target: an order-two character of a modified graph Jacobian. Map: pair the same vertex 0/1 function with divisors after integer edge duplication. Preserved: vertices, underlying simple graph, cut, and optionally involution. Lost: original Laplacian, spanning-tree counts, critical group, and any unproved physical interpretation. Sidecar: explicit edge multiplicities and parity syndrome. Decisive test: (4.4) and the two perfect matchings.

## 5. Bounded-run reconstruction from exactly specified moments

Let S be a finite subset of Z with s≤k nonempty maximal runs

S = disjoint union over i=1,…,s of {a_i,a_i+1,…,b_i−1},

a_i<b_i<a_(i+1).

The strict last inequality says the runs are separated by at least one missing integer; touching runs have already been merged. All 2s endpoints are distinct. Define

M_j = sum over x in S of x^j, for j=0,…,2k−1.

By telescoping,

D_n := sum_i (b_i^n−a_i^n)
     = sum_(j=0)^(n−1) binom(n,j) M_j, for n=1,…,2k.         (5.1)

Define the monic endpoint polynomials

A(T)=product_i(T−a_i),  B(T)=product_i(T−b_i).

They have the same degree s and no common factor. Their normalized ratio at infinity satisfies

A(T)/B(T) = exp(sum_(n≥1) D_n/(n T^n)).                       (5.2)

Thus the given moments determine the expansion through T^(−2k).

### 5.1 Uniqueness

Suppose another set with s'≤k runs has the same moments. Its monic endpoint polynomials A',B' have the same ratio expansion through order 2k. Therefore

AB'−A'B = O(T^(s+s'−2k−1)) as T→∞.

The exponent is at most −1. A polynomial that is O(T^(−1)) is zero, so AB'=A'B. Since both ratios are reduced and their numerator and denominator polynomials are monic, A=A' and B=B'. Sorting the endpoints reconstructs S. The empty set is included separately by M_0=0.

### 5.2 Constructive decoder and zero endpoints

Set z=1/T and compute coefficients c_n of

R(z)=exp(sum_(n=1)^(2k) D_n z^n/n) modulo z^(2k+1).

They satisfy c_0=1 and

n c_n = sum_(j=1)^n D_j c_(n−j).                            (5.3)

For d=0,…,k, solve for q_1,…,q_d in

c_n + sum_(j=1)^d q_j c_(n−j)=0, for n=d+1,…,2k.             (5.4)

The first reduced solution gives Q(z)=1+sum q_j z^j and P(z)=(Q R) through degree d. Their common homogeneous degree is

d=max(deg P,deg Q).

Return the monic polynomials

A(T)=T^d P(1/T),  B(T)=T^d Q(1/T).                           (5.5)

If an endpoint is zero, one of P or Q has degree d−1. Formula (5.5) restores its missing factor T. Because the endpoints are distinct, zero cannot be a root of both A and B; therefore max(deg P,deg Q)=s and no run is lost. It is incorrect to read both endpoint lists solely from the finite reciprocal roots while dropping this common-degree record.

If using this as a recognition procedure on arbitrary supplied data, also require the recovered endpoints to be distinct integers, ordered a_1<b_1<a_2<…<b_s, s≤k, and recheck the original moments. An arbitrary rational Padé pair need not describe a Boolean set.

### 5.3 What is sharp, and in which category

For one run, cardinality alone plainly misses position, so two moments are necessary.

For two runs, take

S={1,5,6},  S'={2,3,7}.

Their moment lists through degree three are

(3,12,62,342),  (3,12,62,378).

Thus three moments do not suffice and four do. This is an exact integer hostile, not a parameter-count heuristic.

For the continuous analogue E=union_i[a_i,b_i] with 2k strictly ordered real endpoints, define C_j=integral_E x^j dx. The Jacobian of the map

(a_1,b_1,…,a_k,b_k) ↦ (C_0,…,C_(2k−1))

has columns −(1,a_i,…,a_i^(2k−1)) and +(1,b_i,…,b_i^(2k−1)). Its determinant is a nonzero signed Vandermonde product. The first 2k−1 rows have rank 2k−1, so by the implicit-function theorem their level sets locally contain a one-dimensional family inside the strict endpoint chamber. Hence 2k moments are necessary locally in the real endpoint problem and sufficient by the same ratio argument, using D_n=n C_(n−1).

This continuous proof does not establish distinct integer solutions for every k. The general all-k discrete sharpness question is left OPEN here.

### 5.4 Exact transfer to the repo's interval Gram decoder

For rows C_i=1_(S_i), with each S_i having at most k runs on an ordered integer line, G_ii gives M_0 and cross-products against the 2k−1 reference vectors

(x), (x²), …, (x^(2k−1))

give the remaining moments. Therefore the full ordered row masks are recovered. Off-diagonal Gram entries remain compatibility tests. At k=1 this is exactly the one-ramp mechanism of THM-2793.

The extension reconstructs masks. It does not create a graceful labeling on a new tree, prove that arbitrary cut lifts form a nested ranked flag, or recover an origin or orientation that was never supplied.

**Connection contract.** Source: bounded-run Boolean masks with ordered integer coordinates. Target: monic entrance/exit polynomials. Map: moments → telescoping boundary powers → rational endpoint ratio. Preserved: exact occupied set and ordered run endpoints, assuming the run budget. Lost without the sidecar: coordinates, Boolean weights, and a valid bound on the number of runs. Sidecar: k, orientation/coordinate frame, homogeneous degree, and endpoint admissibility. Decisive hostile: the two three-element sets above.

## 6. Reproduction and finite universes

Run in this artifact directory:

python3 -B verify_unprompted_bridges.py
python3 -O -B verify_unprompted_bridges.py

The companion uses exact rational arithmetic and explicit exceptions, with no `assert` gates or floating-point arithmetic. It checks:

- all subsets of {1,…,12} with at most k runs for k=1,2,3, including full rational decoding and signature uniqueness;
- all subsets of {−5,…,4} with at most three runs, covering negative and zero endpoints and reciprocal degree drop;
- both cone families for every r=6,…,18, all triple determinants, all required homogeneous ranks, branch separator degrees, and defect formulas; ranks are independently recomputed by Fraction elimination and SymPy;
- the explicit six-branch quadratic finite-difference annihilator;
- all subsets of the 22-edge core of cardinality at most six, with and without involution invariance;
- exact Laplacian harmonicity and Smith/determinant controls for the original and two repaired graphs.

The infinite claims are proved in the text. Finite windows corroborate rather than replace those proofs. Output and SHA-256 manifest accompany this note.

## 7. Classical antecedents and attribution boundary

- The rational transform behind interval reconstruction belongs to the classical Markov moment tradition. Gabriela and Mihai Putinar, *Reconstruction of algebraic sets from dynamic moments* (2007), introduction, explicitly discuss the exponential transform and Padé recovery of unions of intervals: https://web.math.ucsb.edu/~mputinar/mom1.pdf and DOI https://doi.org/10.5802/afst.1163. Our discrete finite-difference translation and repo cross-Gram application are given in full rather than presented as a new general moment theory.
- Harmonic functions and graph critical groups are standard. David Jekel and Avi Levy, *Torsion of the graph Laplacian: sandpiles, electrical networks, and homological algebra* (the linked April 2016 v1), https://sites.math.washington.edu/~reu/papers/2017/avi/torsion.pdf, explains the Q/Z-valued harmonic description. The elementary proof (4.1) suffices here. Minimum parity-join repair is likewise classical; the new scoped calculation is its exact application to this repo core, including the two-edge symmetry premium.
- Graded evaluation on concurrent line arrangements is classical. Leslie G. Roberts, *Ordinary singularities with decreasing Hilbert function*, Canadian Journal of Mathematics 34 (1982), DOI https://doi.org/10.4153/CJM-1982-010-1, Section 2, describes straight-line directions and coordinatewise monomial evaluations: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E9D1AB1F50426FB4891D7CD245750555/S0008414X00026596a.pdf/ordinary-singularities-with-decreasing-hilbert-function.pdf. The note's explicit same-U(3,r) hostile, conductor formulas, and transfer to the repo's lost-gluing coordinate are independently proved here. No literature-wide novelty claim has been established.

## 8. Stopping conditions and useful next gates

- Conductor lane reached an all-r theorem, exact hostile, and branchwise repair coordinate. A next question worth pursuing is the sharp range of seminormal defects among all rational realizations of U(3,r); the present construction only supplies one unbounded difference family. No claim of extremality is made.

  **Current status, 2026-10-19:** The [sharp minimum and maximum, their rational realizations, and the maximum-equality classification](../../conductor_extrema_20261019/conductor_extremal_classification_session7.md) are now proved. The intermediate attainable defect spectrum remains open. The earlier stopping point above is retained as provenance.
- Critical lane reached a structural minimum and exact symmetry premium. It stops before asking whether edge multiplicities correspond to lawful response dynamics, which requires information the abstract core does not carry.
- Moment lane reached a constructive all-k sufficiency theorem, sharp continuous boundary, and exact discrete hostiles. It stops before an all-k discrete sharpness claim or a noise-stability claim. No uniform conditioning is implied when endpoints approach collision.

These boundaries keep the three results usable without silently transporting them to the session's main conjectures.
