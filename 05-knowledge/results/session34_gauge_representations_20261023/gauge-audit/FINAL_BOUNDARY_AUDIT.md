# Final independent audit: closed amplitudes, boundary positivity, and quartic comparison

Verdict: PASS. No substantive correction required.

Reviewed session34-root/CLOSED_AMPLITUDES_BOUNDARY_POSITIVITY.md and its checker/certificate, plus check_line_model_fourth.py and line_model_fourth.json. This separate report leaves earlier audit reports/checkers unchanged. No attached manuscript was used, and no named-paper critique is made.

## 1. General genus and normalization, independently derived

For normalized Haar counting, convolution of central characters is χ_r*χ_s=δ_rs χ_r/d_r. Face gluing along a dual spanning tree therefore replaces the F face densities by w^{*F}, whose expansion is Σ_r d_r a_r^F χ_r. Tree gauge fixing loses no factor under normalized measure. The remaining polygon boundary is the product of g commutators.

There is a second direct derivation of the genus exponent. For an irreducible unitary representation rho of dimension d, Schur averaging and matrix-element orthogonality give

    ∫_(a,b) rho(aba^-1b^-1) da db = I/d².

Indeed averaging a first gives χ(b)I/d before multiplication by rho(b)^-1, and integrating χ(b)rho(b)^-1 gives I/d. For g independent handles, the average matrix is I/d^(2g), hence the averaged character is d^(1−2g). Multiplication by the face coefficient d a^F gives d^(2−2g)a^F. This verifies the normalization and Euler exponent without relying on the tetrahedron or merely inspecting the proposed formula.

Thus for every connected closed orientable triangulated surface,

    Z = Σ_r d_r^(2−2g) a_r^F.

The incidence identity3F=2E forces F even. Reversing the sign coefficient therefore leaves every such closed amplitude unchanged. The claim correctly does not extend to arbitrary cell decompositions with odd face count, or assert triangulation/subdivision independence. For a fixed mesh it is exact equality of the two models.

The independent checker computes the full S3 commutator distribution directly from36 pairs, then convolves handle distributions and face weights. It verifies the formula for genus0,1,2,3 and face powers2,4,6 for BOTH signs. The lower face counts here test the polygon/handle algebra; the general triangulated-surface conclusion is proved by gluing and parity, not by claiming that each tested pair is an abstract simplicial triangulation.

## 2. Boundary kernel and real-time refinement

The class values, strict positivity, and Haar normalization of both w_plus and w_minus are correct. The standard character coefficient1/4 corresponds to normalized representation coefficient1/8, since d_standard=2. No dimension factor is missing.

For the given six-state convolution kernel, the eigenvalues are1,±1/4,1/8, the last with multiplicity4. Its self-adjointness follows from inversion symmetry. The normalized sign quadratic form is exactly±1/4. Therefore the negative model fails precisely the stated one-step boundary positivity test despite its positive entrywise density. Entrywise positivity and positive-semidefinite operator positivity are different statements.

The independent checker constructs the constant and sign orthogonal projectors P0 and Ps, and Pt=I−P0−Ps. The generator divided by log2 is L=−2Ps−3Pt. Every transposition off-diagonal rate is1/3, every nonidentity three-cycle rate is2/3, and row sums are zero. The time-one matrix is

    P0 + Ps/4 + Pt/8 = K_plus.

Thus the claimed central symmetric continuous-time Markov realization is exact.

For K_minus,

    det K_minus = −1/16384.

Every finite real matrix exponential has determinant exp(tr Q)>0. The claimed obstruction to a real time-homogeneous continuous semigroup on this same six-state space is therefore valid. It does not exclude enlarged hidden-state models, changed observations, different step conventions, or positivity after squaring. The source note states these limits appropriately.

The boundary operator is explicitly specified. The finite counterexample proves that these closed scalar amplitudes do not determine that operator's positivity or same-space continuous-time refinement. It does not establish or refute a continuum reconstruction theorem.

## 3. Independent fourth-order line-model calculation

The independent checker does NOT reuse the producer's weight-enumerator polynomial expansion. It reconstructs all16 real-space signed syndrome adjacencies from their code characters, as in the earlier independent spectral audit. Write each adjacency as (0,v^T;v,Q). The exact Schur equation for the dimensionless low eigenvalue is

    z = 7t − t² v^T [I+t(7I−Q)−zI]^-1 v.

Using z=7t−7t²+O(t³), its inverse through t² is

    I+tQ+t²(Q²−7I)+O(t³).

Since v^T v=7, this gives independently

    c3=−v^T Qv,    c4=49−v^T Q²v.

The exact pairs (c3,c4) and sector multiplicities are

    (−42,−203):1; (−6,−11):7; (+6,−11):7; (+42,−203):1.

Consequently the line-model gap is

    36h³/J² +192h⁴/J³ +O(h⁵/J⁴).

The separately audited actual-edge model instead has common quartic scalar−77 in its low band, so its gap has no fourth-order term. Agreement of cubic crown adjacency does not make the kinetic models equivalent. Their first-order scalar shifts, state spaces, elementary moves, and fourth-order band corrections differ.

## 4. Sources and scope

The relevant historical primary paper is Edward Witten, [On quantum gauge theories in two dimensions](https://doi.org/10.1007/BF02100009), Communications in Mathematical Physics141(1),153–209 (1991). The paper's institutional [Princeton publication record](https://collaborate.princeton.edu/en/publications/on-quantum-gauge-theories-in-two-dimensions/) confirms the bibliographic details and two-dimensional scope. Checked on2026-10-23. The publisher PDF link surfaced only an iframe to this audit's browser, so no uninspected paper passage is attributed as proof. The finite-group formula and normalization above have their own derivation and exact controls.

All conclusions concern the stated finite groups, finite matrices, explicit boundary kernel, and fixed-model perturbative limit. No continuum Yang–Mills existence, physical mass gap, volume-uniform estimate, or unnamed physical interpretation follows.

## Reproduction

New files only:
- final_boundary_checker.py
- final_boundary_certificate.json
- final_boundary_certificate.optimized.json

Normal and python3 -O runs agree. Checks use explicit failures, not removable assertions. The certificate records all tested genus/face amplitudes, the negative determinant, and the independently derived cubic/quartic band coefficients.
