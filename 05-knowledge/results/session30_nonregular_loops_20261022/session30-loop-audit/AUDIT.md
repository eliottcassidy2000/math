# Independent audit of the nonregular forest expansion

Session30, 2026-10-22. This packet changes no repository source.

## Conclusions and scope

1. **PROVED:** the proposed general-tournament coefficients of z² and z⁴ in log Φ. An independent elementary proof is below, using simple paths and intersecting triples rather than the producer's partition enumeration.
2. **FINITE-EXACT:** independent permutation-product verification on all 33,867 labelled tournaments of orders 1 through 6, and 128 seeded order-7 examples. The almost-regular simplification passes all 2, 24, 2640 labelled examples at orders 2, 4, 6.
3. **VERIFIED:** the candidate decorated-core engine agrees with an independent exact connected-partition reference through z⁶, separately for every excess/leaf type at weights 1 and 2, on a weighted skew three-state matrix and a hostile transitive four-tournament. The engine coefficients were extracted numerically by a Cauchy sum. This is a numerical comparison with exact references, not an exact-arithmetic audit of the engine.
4. **FINITE-EXACT:** the lowest lollipop coefficient is −4, independently obtained by multivariate clone-multiplicity extraction.
5. **PROOF AUDIT PASS:** root's score-flow centering, uniform resolvent bounds, Schur determinant ratio, unique rooted regular completion, and deleted-doubly-regular identities are correct under their stated hypotheses.
6. **FINITE-EXACT:** the order-6 isospectral H=41/45 witness and order-8 isospectral H=653/651 witness were independently checked using every vertex ordering and exact characteristic polynomials. The order-8 pair also has the same path-tree scalar, 11169728/10230017. Spectrum plus that scalar therefore does not determine H.

No uniform Hamiltonian asymptotic remainder is proved here. The higher-order core resummation and numerical predictor need the scopes stated in their own reports; agreement at finitely many coefficients cannot establish them globally.

## Exact conventions

S is a skew tournament sign matrix, d=S1, M=S/n, δ=M1. Let N_L be the sum of signed directed linear forests with L arcs, and Φ(z)=Σ N_L(z/n)^L. Reversing an odd-length path cancels all forests containing one. Hence N_odd=0, N₂ is the simple directed two-path sum P₂, and N₄ is P₄ plus the unordered pair of disjoint directed two-paths.

The identities are

    N₂ = n(n−1) − ||d||²,

    2N₄−N₂² = n(n−1)(15−7n) + 2(8n−13)||d||²
               −2||Sd||² +4Σ_i d_i²(Sd)_i −Σ_i d_i⁴.

Dividing by n² and 2n⁴ respectively gives

    [z²]logΦ = (n−1)/n − ||δ||²,

    [z⁴]logΦ = (n−1)(15−7n)/(2n³)
               +(8n−13)||δ||²/n² +δᵀM²δ
               +2〈δ∘δ,Mδ〉 −(1/2)Σ_iδ_i⁴.

For an almost-regular even tournament,

    [z⁴]logΦ = (n−1)(15−7n)/(2n³)
               +δᵀM²δ +(16n−31)/(2n³).

These tournament-specialized constants must not be substituted unchanged for a general weighted skew matrix.

## Independent elementary proof

Put a=n−1 and e=Sd. Unrestricted two-walks have sum 1ᵀS²1=−||d||². The backtracks have sum −n(n−1), proving P₂=n(n−1)−||d||².

For each three-element vertex set U, let Q_U be the sum of its six directed two-path signs. It is 6 for a cyclic triple and −2 for a transitive triple. Define Q_v=Σ_{U∋v}Q_U and Q_uv=Σ_{U⊃{u,v}}Q_U. Summing paths according to whether v is an end or the middle gives

    Q_v=2e_v+3a−d_v²,
    Q_uv=2[S_uv(d_v−d_u)+2+(S²)_uv].

The cyclic-triple formula c₃=n(n²−1)/24−||d||²/8 gives

    Σ_U Q_U² = 2n²a−4||d||².

Inclusion-exclusion on the intersection of two triples gives exactly

    2N₄−N₂² = 2P₄−Σ_vQ_v²+Σ_{u<v}Q_uv²−Σ_UQ_U².

The needed four-path identity is

    P₄=||Sd||²−tr(S⁴)+3(n−2)P₂.

Here is a direct derivation. Sum a simple three-path (a,b,c,d), weighted by the row sum d_d at its terminal vertex; call this A. The unrestricted sum is ||Sd||². Removing repeated positions a=c and b=d subtracts (n−1)||d||² each, their intersection adds ||d||², and a=d contributes zero because (S³)_aa=0. Thus A=||Sd||²−(2n−3)||d||². On extending by a final vertex, the forbidden return to a contributes the simple closed four-walk sum tr(S⁴)−n(n−1)(2n−3). Return to b cancels by reversing its triangle; return to c contributes −(n−3)P₂. Substitution proves the displayed P₄ identity.

Finally,

    Σ_vQ_v² = 4||e||²+9na²+Σd_i⁴−18a||d||²−4Σd_i²e_i,
    Σ_{u<v}Q_uv² = 4(n−6)||d||²+16na+2tr(S⁴)−2na².

The potentially mixed term Σ_{u<v}S_uv(d_v−d_u)(S²)_uv vanishes: it equals a sum of d_v times a diagonal entry of S³. Combining all four displayed identities proves the quartic coefficient. For almost-regular d_i²=1, use Σ_i(Sd)_i=−||d||²=−n to get the simplification.

## Independent computation methods

check_low_degree.cpp computes the coefficient of t^L in Σ_permutations ∏(1+tS_next), divides by (n−L)!, and checks the integer identities above. It does not call any producer code.

check_lollipop_low.py extracts the coefficient of q₀q₁q₂q₃ x₀₁²x₀₂x₂₃ using independently generated base walks and falling clone multiplicities. The single-four-run contribution is +2; the two-two-run contribution is −6. Their sum −4 agrees with the proposed −2 diag(H−G)ᵀG1 formula at its first nonzero coefficient. This is not an all-length proof.

check_partition_coefficients.py independently generates set partitions as subsets containing the least unassigned position. It applies the partition Möbius weights, tests quotient connectivity, counts leaves, and evaluates free homomorphism sums exactly. It imports the producer only to obtain the separately compared numerical engine values. All exact references and numerical errors are stored in partition-engine-results.json. The test uses 32 Cauchy points at radius 0.35 and tolerance 2e−8. Reference coefficients are exact rational numbers; the Cauchy error is not enclosed by certified interval arithmetic.

check_spectral_witness.py independently verifies the first two producer witnesses, using all permutations and SymPy exact determinants/characteristic polynomials. It also independently computes the tree scalar.

## Reproduction

    g++ -O3 -std=c++17 check_low_degree.cpp -o /tmp/session30_check_low_degree
    /tmp/session30_check_low_degree
    python3 check_lollipop_low.py
    OPENBLAS_NUM_THREADS=1 python3 check_partition_coefficients.py
    python3 check_spectral_witness.py

The two cross-packet scripts expect sibling directories session30-decorated-engine and session30-almost-regular in this staging layout. Publication should preserve or explicitly adapt these relative paths and rerun the tests.

## Additional completed audits

### Full weight-one/two core inventory

check_core_inventory.py generates nondecreasing degree sequences first, then fills symmetric multigraph multiplicities under the remaining-degree constraints. This is separate from the producer's edge-multiset-first generator. The handshake inequality V≤2w makes the search complete. All 5 weight-one and 26 weight-two cores, and every automorphism factor including repeated edges and loop flips, agree exactly. See core-inventory-results.json.

### Independent connected-partition extension through length eight

partition_reference.cpp is a second reference implementation. It enumerates restricted-growth partitions of all run positions, rejects zero skew loops, tests connectivity, and evaluates homomorphism sums by exact integer color assignments. At L=8 there are 1,585,529 connected partitions, with 205,929 in the requested weight-one/two groups, for each control. The weighted skew3 and transitive4 references are in partition-reference-L8.jsonl. check_engine_L8.py compares 48-point Cauchy extractions at radius0.45, with maximum errors respectively 6.21e−15 and 4.26e−14. This is still a numerical engine comparison to exact references, not interval certification. Every (excess,leaf) group is checked separately, including the seven-bicyclic-one-leaf group.

### Inverse-factorial and predictor algebra

For D=z∂z, the first two exact coefficientwise operators in the expansion of n^L/(n)_L are

    A=(D²−D)/2,
    B=(3D⁴−2D³−3D²+2D)/24.

Writing F₁ for the weight-one core sum, F₂ for weight two, and A₀=e^(−X₀)Ae^X₀, B₀=e^(−X₀)Be^X₀, the normalized second-order formal expression is

    1 + F₁ + A₀/n
      + F₂ + F₁²/2
      + [F₁A₀ +(D²F₁ +(2DX₀−1)DF₁)/2]/n + B₀/n².

This independently derived expression agrees with predict.py, including its mixed Euler-derivative term and spectral Bell polynomials. The finite spectral sum to odd powers through149 has an elementary exponentially small tail from ||M||²<1/2; the Cauchy derivative implementation is floating-point, not certified interval arithmetic. No uniform expansion remainder follows from the algebra.

A wording repair was requested for the reduced-core power count: handle the isolated two-terminal path core separately by its exact δ-resolvent quadratic form. All other leaf-bearing reduced cores attach each leaf to an internal vertex, so their terminal contractions may be bounded independently. The unreduced single-edge graph is identically zero by skew symmetry.

### Complete order-eight almost-regular universe

check_census8.py independently examines all 4!4! score-preserving vertex orders of each of the 85 supplied tournaments. Their canonical integers are distinct, every stabilizer count agrees, and Σ8!/|Aut|=3,230,080. The producer's degree-sequence recurrence is logically sound: it chooses the last vertex's labelled outneighbors, subtracts exactly those remaining scores corresponding to an edge into that last vertex, and exploits permutation-invariance only after this multiplicity has been retained. The sum equals the labelled regular-nine count by unique completion at a fixed ninth label. Thus random discovery of representatives does not weaken this exact mass certificate; no uniform sampling claim is required. See census8-audit.json.

### Root's random-orientation obstruction and the lollipop limit

The score-scale bound D^ell n^(−e−ell), the exact identities Eκ=1, E||d||²=n(n−1), Var||d||²=2n(n−1)(n−2), and the bounds exp((n−1)/n)≤ρ≤(3n−1)/(n+1) are correct. E tr(M⁴)=(n−1)(2n−3)/n³ and the stated atanh bound give Eρ→e. Hence liminf E|κ−ρ|≥e−1 is a valid mean-absolute-error obstruction. It asserts neither κ concentration nor independence.

The theory lane's lollipop_limit_proof.md passes independent algebraic review. In the regular completion normalized by the original n, the Schur resolvent identity has the stated minus sign; its old-old correction is entrywise O(n⁻²). The lexicographic-block functional calculus then gives lim nL=2m[G(B/m)d(B/m)]_b. Its nonzero rational seven-vertex example proves that this formal lollipop term is genuinely first-order. Without a uniform bound on all omitted clusters, it does not by itself prove an asymptotic error theorem for the full Hamiltonian-count predictor.

Additional reproduction commands:

    python3 check_core_inventory.py
    g++ -O3 -std=c++17 partition_reference.cpp -o /tmp/session30_partition_reference
    /tmp/session30_partition_reference
    OPENBLAS_NUM_THREADS=1 python3 check_engine_L8.py
    python3 check_census8.py
