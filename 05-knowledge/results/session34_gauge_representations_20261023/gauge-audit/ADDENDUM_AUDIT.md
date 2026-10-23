# Independent audit of final character, logical-distance, and continuum-boundary additions

Verdict: PASS. No substantive mathematical correction required in the files reviewed below. This report is separate from the main audit and does not modify its files.

Reviewed:
- session34-gauge/TETRAHEDRON_CHARACTER_ADDENDUM.md
- session34-gauge/LOGICAL_EDGE_DISTANCE_ADDENDUM.md
- session34-root/FOURTH_ORDER_AND_CONTINUUM_BOUNDARY.md

## 1. Finite-group tetrahedral formula

The normalized character convention w=Σ d_r a_r χ_r is correct. For a real, inversion-symmetric class function, the coefficients a_r are real and agree for dual representations. Normalized convolution gives χ_r*χ_s=δ_rs χ_r/d_r. Consequently convolution gluing supplies Σ d_r² a_r⁴ at the identity.

Root gauge fixing has no missing factor with normalized Haar/counting measure. In the unreduced count, each rooted configuration has exactly |G|³ representatives. Centrality permits cyclic reordering inside a class function; inversion symmetry handles the face orientations. The stated rooted four-face integral is therefore valid for nonabelian groups as well.

Independent checker method: enumerated ALL46,656 six-edge S3 assignments, formed their four actual oriented triangular holonomies, and accumulated the symbolic polynomial for three arbitrary class weights x,y,z. This bypasses the producer's rooted three-link formula. The resulting normalized polynomial agrees identically with

    ((x+3y+2z)/6)^4 + ((x-3y+2z)/6)^4 + 4((x-z)/6)^4.

Thus agreement holds for arbitrary real class weights, not merely the producer's three numerical tests. It finds216 full flat configurations. Setting y=z=1 gives the rooted polynomial x⁴+30x²+80x+105 after multiplication by216. A second nonabelian check with the order-eight dihedral group independently confirms the universal identity-weight polynomial. Symbolic expansion also verifies the formula for general group order q.

The term “partition function” here denotes the stated algebraic integral. For signed w, it need not define a positive probability measure. The note's signed tests are valid tests of the identity, not claims of probabilistic positivity.

This is a two-dimensional sphere-gluing formula. It does not establish a four-dimensional bulk theory. Extension to infinite compact-group character expansions properly requires convergence justification, which the note explicitly does not assume.

## 2. Logical edge distances

The independent checker constructs the rainbow evaluation kernel directly by GF(2) nullspace elimination on the geometric21 cycle vectors. It does NOT assume the kernel is cuts plus a Fano code. It then partitions all1024 kernel words using canonical root gauge fixing.

The resulting16 cosets each have64 elements. Minimum distances and numbers of minimizers agree exactly:

- One coset: minimum0, one minimizer.
- Seven cosets: minimum3, one minimizer each.
- Seven cosets: minimum6, three minimizers each.
- One coset: minimum9, thirty-five minimizers.

The complete kernel weight enumerator agrees coefficient by coefficient. The matching/cut proof in the addendum is valid. In particular a pencil's nine-edge representative differs from a three-edge near-perfect matching by a star cut. It is important that this representative is a matching, not a triangle. The seven matchings partition all21 edges.

The crown graph description follows from the index-two even/odd decomposition of C and exclusion of one odd generator. Its graph-distance layers1,7,7,1 agree with the exact edge-distance layers0,3,6,9, proving the stated factor-three metric relation. Distance alone is not a tunneling coefficient; the addendum correctly reserves that claim for a specified kinetic operator and intermediate energies.

## 3. Classical scaling and oscillator/flat-direction examples

The scaling A_lambda(x)=lambda A(lambda x) gives F_lambda=lambda²F(lambda x) at fixed coupling, and the squared-curvature integral scales as lambda^(4-d). The d=4 classical scale invariance conclusion is exact and does not imply quantum scale invariance or a mass scale.

For the quartic oscillator the unitary dilation can be made explicit:

    (U_g f)(y)=g^(-1/12) f(g^(-1/6)y).

It maps L²(dx) to L²(dy), and U_g H_g U_g^-1=g^(1/3)H_1. Compact resolvent follows from local compactness and the confining tail bound ∫_(|x|>R)|f|² <= (gR⁴)^(-1)∫g x⁴|f|² for form-bounded families. One-dimensional ground-state simplicity then gives a strictly positive first gap. This supplies the stated conclusion without inferring a gap from a Taylor Hessian.

For K_g, orthogonal u,v coordinates preserve the Laplacian and give the coefficient4g correctly. Expansion in the u-oscillator eigenbasis and Fourier transformation in v yield a direct sum of multiplication operators E_n(4g)+k². Hence the spectrum is exactly [E_0(4g),infinity). The threshold has no L² eigenfunction because k=0 has zero Lebesgue measure. There is neither a normalizable full ground state nor a spectral interval above the threshold free of spectrum. These statements concern the standard self-adjoint/Friedrichs realizations on the full Euclidean spaces, not finite boxes or altered boundary conditions.

The examples therefore support exactly the claimed distinction: quartic positivity and a vanishing Hessian alone determine neither confinement nor a quantum gap.

The topology example is also correct: full clique completion of a three-vertex directed cycle is a filled triangle, whereas its transitive-face complex is its unfilled boundary. The construction must be specified before assigning topology. No vertex algebra follows merely from the oriented one-skeleton.

## 4. Primary-source boundary verification

Checked the following primary pages on2026-10-23:

- Christiansen and Halvorsen, [A simplicial gauge theory](https://arxiv.org/abs/1006.2059), arXiv:1006.2059v4, published J. Math. Phys.53,033501 (2012), DOI10.1063/1.3692167. The authors' abstract states discrete gauge invariance and consistency of their specified simplicial action. It does not claim a continuum mass-gap theorem. Only the public abstract was used here.
- [Clay Mathematics Institute: Yang–Mills & the Mass Gap](https://www.claymath.org/millennium/yang-mills-the-maths-gap/). The official page labels the problem unsolved and distinguishes mathematical proof from evidence by experiments and simulations.

These primary sources support the boundary claims. The elementary finite-group, scaling, and oscillator calculations above are independently derived; neither source is cited as proving those unrelated calculations. No attached manuscript or private paper content was used.

## Reproduction

Run addendum_checker.py normally and with python3 -O. The resulting addendum_certificate.json and addendum_certificate.optimized.json agree. The checker verifies symbolic S3 class weights by full-edge enumeration, an independent D8 control, direct rainbow nullspace distances and multiplicities, and algebraic scaling identities. Analytical spectral compactness and continuum boundary statements are supported by the reasoning above, not presented as consequences of the finite checker.

## 5. Final actual-edge kinetic theorem: fourth-order audit

Also reviewed session34-gauge/EDGE_KINETIC_FOURTH_ORDER_ADDENDUM.md. Verdict: PASS, including the improved remainder O(h⁵/J⁴) for the gap.

Independently checked every ordered two-, three-, and FOUR-step edge sequence, using the auditor's geometrically constructed kernel and exact rational potential values. The four-step enumeration includes all21⁴=194481 words and enforces excited intermediate states rather than assuming the producer's path classification. The results are:

    B0 = B_lambda = 21I,
    B1 = 6 A_crown,
    B2 = 518I.

The420 ordered distinct edge pairs have energies J,3J/2,2J with multiplicities42,168,210 respectively. Their exact rational denominators are essential. The independent16×16 logical adjacency has eigenvalues7,1,−1,−7 with multiplicities1,7,7,1.

After shifting21h and dividing by J, the Schur equation uses lambda=O(t²). Expanding its inverse introduces the term −lambda t² B_lambda. Substitution of lambda=−21t²+O(t³) contributes +441t⁴, giving the scalar quartic correction −518+441=−77. All fourth-order off-diagonal entries vanish, as confirmed by the complete four-step enumeration. Therefore the fourth-order correction cancels in differences of low eigenvalues, leaving

    Delta_edge = 36h³/J² + O(h⁵/J⁴).

The remainder argument is valid for this fixed finite matrix: D_Q>=I, ||A||<=21, and small-t spectral separation allow a locally uniform Neumann expansion. Counting eigenvalues by the isolated low projection or scalar code-character sectors justifies the stated multiplicities through the displayed order. The seven first-excited branches need not be declared exactly degenerate to all orders merely on this calculation; the theorem only needs their shared expansion.

The producer explicitly states the necessary scope restriction: this is the32768-dimensional B-invariant Hilbert space, not the unrestricted2²¹-dimensional space. The gauge transformations B are a declared state-space constraint of this finite model. Code translations remain symmetries and are not projected away. The result is neither the root's128-state line-flip Hamiltonian nor a continuum physical mass-gap theorem.

The independent addendum certificates now contain20 explicit aggregate checks, encompassing the exhaustive symbolic group, nullspace, and virtual-path enumerations. Normal and optimized runs remain identical. No previous main audit artifact was modified.
