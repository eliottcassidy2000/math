# Independent review: critical injectivity and incidence trees

Date: 2026-10-19. Verdict: **accepted; no mathematical gap found; no correction requested**.

Reviewed artifact: `packet138-stage/05-knowledge/results/polynomial_critical_injectivity_incidence_tree_20261019.md`.

Reviewed final SHA-256: `6b5f4921ff90c0af539af978085f0ba66a5739b6bf628a914da7281d51520ffe`.

The mathematical text has SHA-256 `ecf104b00410e320468137b996f4475bd0a5f0d4698ef230cc703968de492487`; the final artifact only appends the producer’s verification-scope section. This was confirmed by hashing the prefix before that appended section. The reported 13-check receipt is not an independent replay by this reviewer.

## Focused assessment

1. **Simultaneous cycles and the tree.** Standard based loops around all finite branch values give a simultaneous generating tuple in one labeling; the infinity permutation is already a product of these generators. The punctured polynomial cover is connected. Consequently the support-incidence graph is connected, including every sheet vertex. Each generator is one cycle because critical injectivity leaves only one ramified point over that finite branch value. The edge count is sum(d_i+1)=n−1+q, exactly one less than n+q vertices. It is therefore a tree, and with at least two cycle vertices some pair of cycle supports intersects in exactly one sheet. This avoids the invalid alternative of choosing independent inertia conjugates and presuming they generate transitively.

2. **Commutator and elementary alternating-group argument.** For cycles meeting only at c, with rightmost action first, the commutator sends c to sigma(c), sigma(c) to tau(c), and tau(c) to c, fixing the other sheets. Thus the displayed 3-cycle and its orientation are correct. The connected-support lemma is also correct: A_U acts transitively on unordered pairs for |U|>=3, including |U|=3; taking inverses supplies either orientation. A support meeting U in one point can be joined using (a x y)(b x y)^(-1)=(a x b), after which two successive two-point overlaps enlarge the alternating group. The normal closure of the initial 3-cycle is nontrivial and normal in the primitive group, hence transitive. Its connected support hypergraph then yields A_n without a classification theorem.

3. **Primitivity and parity.** Critical injectivity rules out a nontrivial polynomial composition by forcing the inner factor to be a translated pure power over any critical point of the outer factor. The outer factor then has only one critical point, which forces the whole composition to have only one critical point, contrary to hypothesis. The standard intermediate-field/Luroth argument gives primitivity. Since the simultaneous finite cycles generate the group, their parity proves the alternating upper bound when all d_i are even, and gives the symmetric group when some d_i is odd.

4. **Extension of constants.** The reduction through a finitely generated characteristic-zero coefficient field and its algebraic closure does not require the original algebraically closed field to embed in C. A field L over the algebraically closed field k is geometrically integral as a k-algebra. Thus L tensor_k Omega is a domain. Localizing at the nonzero polynomials in Omega[Y] gives precisely L tensor_(k(Y)) Omega(Y); those polynomials remain nonzero because Y remains transcendental over Omega. This localization is a domain, and its finite dimension over Omega(Y) makes it a field of the unchanged extension degree. Automorphisms and root generators therefore preserve the full permutation representation of the splitting-field group. The stated base-change proof is sufficient even when the ambient algebraically closed field has large transcendence degree.

5. **Arithmetic and scope.** If all derivative multiplicities are even, the unique monic square root of f'/(n*leading_coefficient(f)) descends to K in characteristic zero. The monic discriminant formula has square class (-1)^((n−1)/2)n, exactly as stated. The arithmetic group lies between the geometric alternating group and S_n, and a constant becomes square in K(Y) only if it was square in K. The theorem remains generic-target; the fixed-time transfer is correctly kept as a separate input from packets 135 and 137.

## Primary-source checks and limits

Independently inspected the simultaneous branch-cycle description in Peter Muller, *Primitive Monodromy Groups of Polynomials*, section 2.1: https://ypfmde.github.io/Papers/mon.pdf . Only its covering setup is used; its classification results are unnecessary here.

Independently inspected Stacks Project Lemma 10.49.3, which gives the stated geometric-integrality criterion for arbitrary k-algebras: https://stacks.math.columbia.edu/tag/0FWF .

No broad computation or Git action was performed. The argument supplies a valid CFSG-free replacement for the relevant group-identification steps in packets 135 and 137, without altering their separately established source-transfer, coordinate-field, or specialization scope. This separate review may be copied into the producer packet.
