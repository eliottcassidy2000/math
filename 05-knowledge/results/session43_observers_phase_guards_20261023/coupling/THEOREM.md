# When the 8/9 coupling gap survives restrictions

Session43, 23 October 2026. Proven statements below concern the graph/tournament quotient coupling, not the native Hamiltonian-path metagraph walk. No literature-priority claim is made. The starting theorem is Session42's signed-Walsh diagonalization: on the uniform edge cube, K=P_G P_T P_G has nonconstant eigenvalues mu(F)^2 <= 1/9, with equality precisely for P3 and C4 shapes after deleting isolates. The signed basis itself is inherited from repository THM-4721 / THM-4531.

## 1. Three operations that must not be conflated

(A) Restrict the original linear operator to a reducing subspace. Every surviving nonconstant eigenvalue stays at most 1/9. The gap is exactly 8/9 if and only if the subspace contains a nonzero vector in the P3/C4 eigenspace; otherwise the gap is strictly larger (when the mean-zero space is nontrivial).

(B) Coarsen either observed partition while keeping the full uniform edge cube. This also preserves the lower bound 8/9, without requiring the new graph-statistic subspace to reduce K. Indeed let U' be a subspace of graph-class functions and V' a subspace of tournament-class functions, both containing constants, with orthogonal projections Q_U,Q_V. For mean-zero f in U',

    <f,Q_U Q_V Q_U f> = ||Q_V f||^2
                             <= ||P_T f||^2 <= (1/9)||f||^2.

Equality is possible exactly when f lies in the original top eigenspace and P_T f lies in V'. Thus retaining an arbitrary coarse statistic cannot slow this coupling. Complement/converse pairing is a special case retaining the even-Walsh-degree nonzero modes. Summarizing a tournament by its Hamiltonian count modulo 8 is another legal coarsening and inherits the bound. CONDITIONING on a particular residue is not this operation.

(C) Condition the underlying words to a subset A and rebuild conditional expectations. The original actions need not preserve A. Write C_ij=|A intersect graph class i intersect tournament class j|, r_i=sum_j C_ij, c_j=sum_i C_ij. Then the honest restricted chain is

    K_A = D_r^{-1} C D_c^{-1} C^T.

It is still positive semidefinite, reversible and stochastic. But no 8/9 lower bound follows. Its eigenvalue-one multiplicity is exactly the number of connected components of the bipartite incidence support. Proof: a function survives both orthogonal projections unchanged exactly when it is constant on each support component. In particular common information, not a failed numerical approximation, causes disconnectedness.

There is no nonempty proper subset invariant under BOTH full original actions. Their generated group contains every single-edge translation (compose the two versions of an adjacent transposition, then conjugate by graph permutations), so is transitive on the cube. A purported nontrivial common-invariant conditioning must therefore change the actions or the model.

## 2. Hamiltonian observer conditioning destroys the global gap at every size

Use x_ij=1 for backward tournament arcs and graph edges. Fix the Hamiltonian observer path 1 -> 2 -> ... -> n. Its tile cube is

    A_H = {x : x_(i,i+1)=0 for every i},
    |A_H| = 2^((n-1)(n-2)/2).

**Theorem.** For every n>=3 the conditioned coupling K_AH has spectral gap zero, both unpaired and complement/converse-paired.

**Proof.** A transitive tournament has a unique Hamiltonian vertex order. If its observer path is 1,...,n, every arc points forward, so its backward-edge graph is empty. Conversely the empty graph is exactly this word. Hence the empty-graph row and transitive-tournament column form an isolated 1-by-1 incidence block. Other words exist for n>=3, for example reversing only edge {1,3} gives a directed triangle while preserving the observer path. Pairing does not remove the block: the only graph complement mate of the empty graph is the complete graph, which is absent from A_H, and the transitive tournament is self-converse. Therefore there is a nonconstant common invariant, namely the indicator of this block. QED.

This is the smallest nontrivial counterexample already at n=3, where C is the 2-by-2 identity and K is the identity. n=2 has only one tile and no nontrivial mean-zero space.

Removing the transitive tile is not a universal repair. At n=4, using lexicographic edge bits 12,13,14,23,24,34, the five graph representatives are masks 0,1,3,12,13 (empty, single edge, P3, two disjoint edges, P4). With tournament representatives 0,2,4,5 the incidence matrix is

    [1 0 0 0]
    [0 1 1 1]
    [0 0 2 0]
    [0 0 1 0]
    [0 0 1 0].

The seven-word nontransitive component has transition matrix

    [11/15  2/15  1/15  1/15]
    [ 1/5   2/5   1/5   1/5 ]
    [ 1/5   2/5   1/5   1/5 ]
    [ 1/5   2/5   1/5   1/5 ].

Its eigenvalues are 1,8/15,0,0, hence gap 7/15. Converse pairing merges the two columns with singleton mass, but produces exactly the same K here. The stationary weights of this component are 3,2,1,1 divided by 7, not uniform graph-class weights.

### The observer changes the measure, not just labels

For an unlabelled tournament T, the number of A_H words of type T is h(T)/|Aut(T)|. Indeed labelings whose label order is a Hamiltonian path correspond to directed Hamiltonian vertex sequences, and automorphisms act freely on such sequences. Summing these ratios gives |A_H|. This differs from the full-cube mass n!/|Aut(T)|. An observer selection is a Hamiltonian-count bias, with correlations between graph and tournament observations; the Walsh orthogonality argument cannot be imported unchanged.

Converse has an honest lift preserving this observer cube: reverse the vertex order AND reverse every tournament arc. On backward-edge bits this is plain reversal of coordinate indices, so it preserves the graph class and takes T to its converse. Merely complementing every bit does not preserve A_H. This distinction is essential for native ports.

## 3. Actual modulo-8 conditioning can both slow and accelerate

Condition the FULL edge cube on h(T) congruent to a chosen residue modulo 8, retaining the original unpaired class partitions. Exhaustive exact checks through order 5 give:

- At n=3 and n=4 every nonempty residue-conditioned chain has only eigenvalues 1 and 0, and gap 1.
- At n=5, the h(T)=15 (residue 7) subset has 64 words and 22 graph states. Its spectrum is 1,89/315, and twenty zero eigenvalues. Thus its exact gap is 226/315, strictly below 8/9. This is the smallest order where any nonempty modulo-8 Hamiltonian-count filter violates the 8/9 bound, among n>=2; n=2 is trivial.
- At n=5, residue 5 has 360 words, 28 graph states, and spectrum 1,143903/1801800, and twenty-six zeros. Its exact gap is 1657897/1801800, strictly ABOVE 8/9.

Thus both ends of the discrepancy are meaningful: a residue as an OBSERVED SUMMARY preserves the bound by coarsening; a residue as a CONDITION changes weights and can speed or slow the chain. Exact characteristic polynomials for every nonempty residue through n=5 are recorded by the checker. None of this licenses discarding mod-8 guards in a native compiler.

## 4. Endpoint marks and shared ports

Marks typically refine, rather than coarsen, the isomorphism partitions. The uniform-cube comparison proof does not extend to arbitrary refinements. For example replace S_n on both sides by its subgroup fixing two named endpoints a,b individually. The edge coordinate x_ab is then fixed under BOTH actions: the affine inversion correction at that pair is zero. Its nonconstant indicator gives eigenvalue 1, so the unpaired marked coupling has gap zero. This is an actual legal subgroup action, not a guessed group restriction.

Converse averaging may remove this particular odd-degree witness, so this endpoint example alone does not establish a paired obstruction. With three individually fixed vertices, a product of two distinct edge characters among them is even under global complement and fixed by both actions; it gives a paired obstruction. Fixing the port values and working inside one common-information component removes that particular witness but requires a fresh incidence/spectral calculation. A native compiler must retain the coordinates rather than infer compatibility from a gap on an unmarked quotient.

## 5. Sufficient criteria and exact certificates

For arbitrary A, the gap is at least gamma if and only if the following rational symmetric matrix is positive semidefinite:

    D_r - C D_c^{-1} C^T
        - gamma (D_r - r r^T/|A|).

This follows directly from the reversible Rayleigh quotient and is a complete finite certificate, retaining actual incidence and weights. It does not reduce to vertex/edge counts.

A strong sufficient combinatorial criterion for preserving the ENTIRE original chain is balanced thinning: C_A = alpha C_full for a common positive alpha. In that case both marginal diagonals scale by alpha and K_A=K_full. Separate row/column balance alone is not asserted sufficient; the joint incidence is what matters.

A robust nonzero-weight variant is also available. Let w be a strictly positive normalized density on the full cube, with 0<a<=w<=b, and keep the original partitions. Then gap(K_w)>= (a/b)*(8/9). For graph-class f,

    E_w Var_w(f|T) = inf_(h tournament-invariant) E_w (f-h)^2
       >= a inf_h E (f-h)^2 >= a*(8/9) Var(f),
    Var_w(f) <= b Var(f).

Their ratio is the gap Rayleigh quotient. This proves that bounded reweighting cannot disconnect the coupling, while hard conditioning has a=0 and evades the comparison exactly as the observer example does. In particular weights 1 on accepted words and epsilon on rejected words guarantee gap >= 8 epsilon/9. Exact finite controls verify epsilon=1/2 for the order-four observer guard and order-five residue-seven guard. This is a different weighted model, not a license to relax a required hard compatibility condition.

## 6. What is proved and what remains

Proved: coarsening stability, reducing-subspace criterion, exact incidence certificate, all-size observer isolation, legal endpoint-mark obstructions, bounded-density comparison, and the stated exact small restricted spectra. Inherited: full-cube signed-Walsh spectrum and its equality classification. Not proved: a native port-preserving metagraph gap, an all-size nontransitive observer-core spectrum, or a characterization of all hard filters retaining 8/9. A useful next problem is to characterize common information and connected components BEFORE estimating mixing inside each component.
