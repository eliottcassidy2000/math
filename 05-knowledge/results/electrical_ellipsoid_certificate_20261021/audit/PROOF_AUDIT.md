# Independent audit: electrical-network ellipsoid certificate

Verdict: PASS for the pinned producer packet. No blocking mathematical error found.

Audit date: 2026-10-21 UTC. The audit owns only this directory, staged outside the repository. No producer file or Git state was changed. The four producer files were pinned before review. Root subsequently made only the three readability clarifications recommended below; the exact diff was reviewed. The final four producer files are copied into `frozen-inputs/`, with final hashes in `SOURCE_PINS.json`; the initial README and pins remain in `initial-snapshot/`. Any further producer edit is not automatically covered by this audit.

## 1. Scope and independence

The reviewed result is a standard maximum-determinant/contact-duality specialization to the projected cube, expressed through a grounded graph Laplacian and effective resistances. It is not a new general duality theorem. Neither this result nor this audit provides an arithmetic orbit hit, a lonely time, a lattice-incidence theorem, a uniform LRC improvement, or a literature-priority claim.

The base packet `session19-ellipsoid-optimum/README.md` and the triangle packet `session19-triangle-ellipsoid-repair/README.md` were read and pinned as context. The carrier, gauge, certificate, existence, matrix-tree identity, tree obstruction, general one-pair construction, and two literal examples are rederived below rather than accepted as unproved inputs. The triangle packet's additional diagonal-lift, comparison, and complete-triple results are outside this audit.

The independent executable imports only the Python standard library and neither reads nor imports producer files. Its determinants use the Leibniz formula and inverses use adjugate cofactors, unlike the producer's Gaussian elimination. Weighted tree counts are independently enumerated by graph connectivity, not by incidence minors. Every executable check raises an explicit exception; none depends on `assert` or `__debug__`.

## 2. Carrier and grounded-coordinate audit

Fix the normalization D = P[-1,1]^k in H = v-perp, with all v_i positive and r_i = 1/v_i. The unambiguous grounded map is

    Phi(x)_i = x_i/v_i - x_k/v_k,       i < k.

Writing y_k = 0 is a choice of representative in the quotient by common constants. It does not mean that the original representative x in H has x_k = 0. This distinction is respected by the producer's quotient language.

For any grounded y, put ytilde = (y_1,...,y_(k-1),0), s = sum_i v_i^2, and

    alpha = sum_(i<k) v_i^2 y_i / s,
    x_i = v_i (ytilde_i - alpha).

Then x lies in H and Phi(x) = y. This proves that Phi is a linear isomorphism.

Membership in the projected cube is equivalent to existence of a common shift t with

    |v_i (ytilde_i + t)| <= 1 for every i.

Equivalently, the intervals [-r_i-ytilde_i, r_i-ytilde_i] have a common point. A finite collection of real intervals has a common point exactly when every pair intersects, giving precisely

    |ytilde_i-ytilde_j| <= r_i+r_j.

This directly proves the carrier without requiring polarity. The producer's cap normalization is correct for [-1,1]^k; a unit-side cube would require the corresponding factor of one half.

For additional normalization control, the inverse map G has Gram matrix

    G^T G = diag(v_1^2,...,v_(k-1)^2) - aa^T/s,
    a = (v_1^2,...,v_(k-1)^2).

The determinant lemma gives det(G^T G) = product_i v_i^2 / s. Thus the squared Euclidean-volume determinant of an ellipsoid is its grounded covariance determinant times this factor. This confirms that the coordinate determinants in the producer receipts are not being mistaken for Euclidean determinants.

## 3. Feasibility, electrical interpretation, and sufficiency

For covariance Sigma > 0, the support of Sigma^(1/2)B in the b_e direction is sqrt(b_e^T Sigma b_e). This proves the pair feasibility conditions. An orientation change b_e -> -b_e has no effect.

With lambda_e >= 0, the quadratic form of the grounded Laplacian is

    z^T L z = sum_e lambda_e (z_i-z_j)^2, with z_k = 0.

Consequently L is positive definite exactly when the positive-weight graph is connected. If disconnected, the indicator potential of a component not containing the ground is a nonzero kernel vector. This is also why connectivity must concern positive weights rather than the full list of possible pairs.

For unit injection b_e, the grounded potential is p = L^(-1)b_e. Its terminal voltage b_e^T p and energy p^T Lp both equal b_e^T L^(-1)b_e. This verifies the effective-resistance interpretation, including nonedges.

If Sigma = L^(-1) is feasible and every positive-weight edge is tight, then

    sum_e lambda_e c_e^2 = tr(L Sigma) = d.

Every feasible T satisfies tr(LT) <= d. The positive eigenvalues of L^(1/2)TL^(1/2) satisfy AM-GM, giving det(T)det(L) <= 1. Equality forces that matrix to be the identity, hence T = Sigma. This proves the unique covariance claim.

For a translated ellipsoid, the slab condition is

    |b_e^T h| + sqrt(b_e^T T b_e) <= c_e.

Its centered covariance is feasible. Equality in volume first forces T = Sigma; each positive contact then forces b_e^T h = 0. Connectivity and grounding imply h = 0. Singular covariances have zero volume and cannot maximize because strictly feasible positive-definite shapes exist. Therefore uniqueness includes the center, not merely the covariance.

For arbitrary connected nonnegative weights, exactly the same argument with C = sum_e lambda_e c_e^2 gives

    det(T) <= (C/d)^d / det(L).

No feasibility of L^(-1) is needed for this upper bound. A failed resistance check is not evidence that the original convex body is empty.

## 4. Necessity and duality audit

For full rigor, take the closed feasible covariance set in the positive-semidefinite cone. The grounded-pair constraints bound each diagonal entry by (r_i+r_k)^2. Positive semidefiniteness gives |Sigma_ij|^2 <= Sigma_ii Sigma_jj, so all entries are bounded. The set is compact.

For sufficiently small epsilon > 0, epsilon I satisfies every constraint strictly. Hence the continuous determinant attains a positive maximum on the compact set, and any maximizer is positive definite. Maximizing its logarithm is equivalent. Strict concavity of log det gives uniqueness on this convex feasible set. This is the compact-closure detail implicit in the producer's existence paragraph; there is no boundary-attainment gap.

Let A_e = b_e b_e^T and G = Sigma^(-1). The cone generated by active A_e is closed. Besides the standard finite-generated-cone fact, this follows here directly because every generator has positive trace, so coefficients of a convergent nonnegative combination are bounded. If G is not in this cone, separation gives a symmetric U with

    tr(GU) > 0 and tr(A_e U) <= 0 for every active e.

A sufficiently small positive step Sigma+tU remains positive definite. Active inequalities cannot increase, and finitely many slack inequalities remain slack for a sufficiently small step. The derivative of log det at t=0 is tr(GU)>0, contradicting optimality. Thus G is a nonnegative active combination. Positive definiteness forces its support graph connected. The effective-resistance and complementary-contact conditions follow immediately.

The dual objective is C - log det L - d on the domain lambda >= 0, L > 0. The certificate attains equality with the primal log determinant, so both the dual formula and equality claim are correct. The final producer revision explicitly states this positive-definite domain.

## 5. Matrix-tree identity

With reduced incidence columns B and W = diag(lambda), L = BWB^T. Cauchy-Binet gives a sum of squared square incidence minors times edge-weight products. A k-1-edge subset has nonzero reduced incidence determinant exactly when connected, hence exactly when a tree. Leaf removal proves its determinant has absolute value one. Disconnected subsets have deficient rank. This proves the stated weighted spanning-tree expansion with no orientation issue.

An unweighted tree count omits essential weight information and cannot replace either the determinant or resistance inequalities. The producer explicitly avoids that invalid inference.

## 6. Every tree certificate is a dominant star

In a tree, unit injection at the endpoints of an edge sends unit current through that edge and no current through a side branch. Hence its resistance is 1/lambda_e. For general terminals, conservation forces unit current along their unique path, so the resistance is the sum of reciprocal conductances on that path.

Contact equality on every positive tree edge forces lambda_e = 1/c_e^2. Feasibility is therefore exactly the stated squared path criterion.

For a length-two path a-b-c, expansion of that criterion gives

    r_b^2 + r_b(r_a+r_c) <= r_a r_c.

If r_b >= r_a, the left side already strictly exceeds r_a r_c; similarly r_b cannot be at least r_c. Thus r_b < min(r_a,r_c). A length-three path a-b-c-d would imply both r_b < r_c and r_c < r_b. Every nonstar tree has a length-three path. Thus a certifying tree must be a star.

If its center has speed M and two leaves have speeds a,b, the two-edge inequality is, after multiplying by the positive denominator a b M^2,

    M^2 >= M(a+b)+ab.

This implies M>a+b, so the center is the unique fastest vertex. The right side is strictly increasing in each leaf speed, making the two largest remaining speeds decisive. Conversely those inequalities make the star with lambda_i = 1/(r_i+r_center)^2 a certificate.

Outside this cone, connectivity requires at least k-1 positive edges and the k-1 case has been excluded. Every certificate therefore has at least k positive edges and contains an undirected cycle. This is not a temporal, directed, or tournament-cycle statement. Boundary equality is allowed inside the cone and does not force a cycle: an additional tight facet may have zero multiplier.

## 7. Triangle construction and literal controls

The general one-pair construction can be verified directly in grounded coordinates. Put s_i = r_i+r_k, where v_k is largest. If only the pair p,q violates the diagonal-star support test, set

    rho = (s_p^2+s_q^2-c_pq^2)/(2 s_p s_q),
    Sigma_ii = s_i^2,
    Sigma_pq = Sigma_qp = rho s_p s_q,

and all other off-diagonal entries zero. Since c_pq > |s_p-s_q| and the numerator is positive, 0 < rho < 1. The exceptional pair becomes tight, all other supports are unchanged, and Sigma is positive definite.

The inverse has one negative off-diagonal pair, giving the positive conductance rho/[s_p s_q(1-rho^2)]. The two affected grounded-edge weights are positive: for example the relevant numerator is proportional to

    s_q^2 - rho s_p s_q = (s_q^2-s_p^2+c_pq^2)/2 > 0,

because v_k >= v_q implies s_p <= c_pq. The other endpoint is symmetric. Unaffected grounded edges have weight 1/s_i^2. Every positive edge is tight, so the certificate theorem applies. This establishes the triangle cross-reference without depending on the companion's broader results.

The independently reconstructed literal examples are:

- Speeds (1,2,4): grounded determinant 225/256; weighted tree sum 256/225; Euclidean determinant 75/28.
- Speeds (1,2,3,5): grounded determinant 124/625; weighted tree sum 625/124; Euclidean determinant 1488/325. The only nonzero leaf covariance entry is 1/25, giving normalized correlation 3/28.

All conductances and all pair resistances match the pinned receipts exactly. Deleting the triangle's extra edge leaves a connected tree but fails a resistance constraint, confirming that connectivity alone is insufficient.

## 8. Reproduction, coverage, and nonblocking clarity notes

Run from any working directory, using absolute paths as needed:

    python3 -B independent_check.py
    python3 -B -O independent_check.py

The two outputs must match `independent-normal.json` and `independent-optimized.json`, which are byte-identical. The producer was separately replayed in normal and optimized modes from its frozen input copy; both outputs match the corresponding pinned receipt byte-for-byte.

The independent run contains 4,965 explicit checks, including:

- 348 quotient, projected-vertex, and interval-lift cases.
- 24 determinant-bound samples with arbitrary connected conductances and nonhomothetic feasible covariances, not just scaled copies of a certified shape.
- All labelled spanning trees in declared nondecreasing positive-integer speed boxes: k=3 with speeds 1..6; k=4 with speeds 1..5; k=5 with speeds 1..3. This is 3,913 tree/speed cases.
- Every k-1-column incidence minor for the complete graphs through k=5, unequal-weight matrix-tree controls, a disconnected-support singularity hostile, and translation/contact controls for both literal certificates.

Finite controls do not replace any universal proof above. The producer's own feasible covariance samples are homothetic, so its AM-GM sampling is narrow; the independent checker fills that testing gap without changing the proof verdict.

Resolved nonblocking clarity suggestions: root explicitly named D=P[-1,1]^k, added Phi(x)_i=x_i/v_i-x_k/v_k and its inverse in the opening paragraph, and stated lambda>=0, L>0 beside the dual objective. The reviewed final README hash is 493e454602518eb85ab458f6d866f73bbb34a84e066c95cc3abc4a1a37cc513f. The checker and both producer receipts are unchanged. These clarifications do not alter the theorem or proof verdict.
