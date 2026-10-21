# Independent review: Fano supports and integral flow series

Reviewed frozen producer packet `../session23b-flows/`, 2026-10-21.

## Verdict

PASS in the stated scope, conditional on the explicitly cited inherited all-order finite-strong-alphabet and center laws. No mathematical formula correction was needed. One harmless degenerate-case proof wording was identified and has now been corrected in the final producer source: in the anti-coset proof, use the **signed/tagged extension** of automorphisms and anti-isomorphisms, or treat the zero-edge singleton separately. For K1 the ordinary and anti vertex-permutation sets are literally the same singleton, so an untagged union does not have twice the order. Formula (6) is nevertheless correct, with both series equal to 1.

The producer's normal and optimized controls were rerun independently and match its frozen receipts byte for byte. Independent controls are in `independent_check.py`; their normal and optimized receipts are `normal.json` and `optimized.json`, also byte-identical (42,029 always-active checks). Their embedded proof hash identifies the initially frozen text. `proof_precision_edit.json` separately pins the final proof and verifies that reversing only the tagged-group prose correction exactly recovers that initial hash. No numerical code changed.

## Scope actually checked

1. All 64 labelled K4 orientations, their cycle supports, and their opposite-orientation fibers.
2. All 7 Fano points, 7 lines and 168 incidence-preserving permutations, with the actual faithful tetrahedral S4 action and its 24-element partition stabilizer.
3. T4's integral cone by directly solving incidence equations on a spanning tree, through total arc use 24, independently of its advertised cycle parameters. This gives 119 flows, exactly the free cone on its two triangles and quadrilateral.
4. R5 directly from incidence equations through total arc use 24: 7,315 flows, 1,475 automorphism orbits, and 309 self-converse flow-orbits. Both ordinary and anti-coset Burnside agree with explicit orbit representatives.
5. All strong ordinary tournament types through order 5, reconstructed from labelled tournaments, plus a real non-self-converse order-6 pair with h=19. The order-6 types are specified by triangular masks 20 and 26 in the script's documented encoding. This is not a complete order-6 census.
6. All 838,695 decorated SCC words on that finite alphabet with at most 12 vertices and total arc use at most 12. Literal converse-orbit counts agree with formulas (9)–(15), across 28 underlying-h slices and 33 underlying-degree slices. In particular the degree-zero test includes the nontrivial R5 center, and the test alphabet contains genuine non-self-converse strong atoms, unlike the producer's smaller word alphabet.
7. The explicit R5 anti-fixed parameter sum (18), not only abstract anti-coset Burnside, checked against direct incidence-kernel enumeration through total use 24.

These are controls of the proofs, not claims that the finite samples prove the all-order conclusions or exhaust strong types beyond the stated alphabet.

## Main proof audit

### Canonical K4/Fano structure

The even-degree edge kernel has binary dimension three. Its nonzero elements are exactly four faces and three quadrilaterals. The four face vectors have only their full-sum dependence, so arbitrary face permutations extend uniquely to the 24 tetrahedral linear symmetries. The larger Fano automorphism group has order 168 and does not preserve the extra tetrahedral partition in general.

The tournament census is correct: 24 transitive, 16 with one directed triangle, and 24 strong with two triangles and one quadrilateral. The strong three-point sets are noncollinear. Each of 12 such support triples has exactly the two opposite tournament orientations. Thus support-stabilizer order two must not be substituted for the trivial tournament automorphism group.

### Integral T4 cone

The displayed balance equations give the claimed vector (a+b+c,a+c,a,b,b+c,c), with a,b,c literal arc coordinates. The corresponding identity minor proves lattice unimodularity, not merely real simpliciality. The anti-isomorphism swaps the two triangle coordinates and fixes the quadrilateral coordinate. Therefore ordinary and self-converse series have the claimed respective denominators (1-q^3)^2(1-q^4) and (1-q^6)(1-q^4).

### Fixed rational cones and Burnside

The fixed cones are rational and pointed. A rational simplicial subdivision plus fundamental-parallelepiped numerators proves rationality; the strictly positive total-coordinate grading makes the one-variable specialization well-defined. No unimodularity assumption enters the general argument. Fixed-weight slices are finite, so ordinary and anti-coset Burnside are legitimate coefficientwise. The converse-fixed count is an average over the whole anti-coset, not generally the fixed-point count for one arbitrarily chosen anti-isomorphism. The signed-extension wording above handles K1 without changing any formula.

### SCC factorization and three genuinely different carriers

Summing conservation over forward SCC cuts proves all cross-SCC flow coordinates vanish. Ordered SCCs have a direct-product automorphism group with no permutation factor even for repeated equal atoms. This proves the ordinary word carrier.

Three constructions must remain separate:

- Ordinary half U: its usual flow-orbit series f_U(q).
- Mirrored decorated halves: f_U(q^2), since the opposite decoration is forced.
- Independently decorated halves with the same underlying U: f_U(q)^2, since only the underlying opposite type is forced.

Consequently the independent neutral carrier has denominator 1-z^2-z^6/(1-q^3)^2, while the mirrored carrier has denominator 1-z^2-z^6/(1-q^6). They are not interchangeable. The direct neutral hostile C3 followed by C3 gives coefficient 2 versus 0 at z^6 q^3. The producer's larger K1,C3,T4 alphabet gives 4 versus 0 there, correctly explained by the extra K1,T4,K1 word.

For fixed underlying h, the factorization sum has finitely many nonrigid factors because their odd h-values exceed one. The P_j may have unknown coefficients but are finite rational sums by the inherited order bound. This is an effective rationality argument, not a computed full atom classification.

For fixed underlying degree, L counts all ordinary decorations over self-converse underlying types; M counts the actually self-converse decorations; F-U removes all underlying self-converse types from the non-self-converse branch. Formula (14) is exactly the disjoint two-case Burnside quotient. In particular substituting M for L or dividing L alone by two is wrong. Formula (15) uses the complete inherited degree-zero center alphabet, includes R5, and retains both carriers. The all-order formulas concern underlying h and degree, not invariants recomputed with reduced decorated stabilizers.

### R5 cone and anti-fixed sum

Writing t=x_i+y_i+y_(i-1), the balance equations say all five such values agree. This gives the advertised integral inverse and weight 5t-sum(y_i). At t=1 the inequalities are those of the fractional stable-set polytope of the pentagon. Its vertices are the zero vector, five singleton stable sets, five two-point stable sets, and the all-half vector. For the last assertion, after any zero coordinate the tight positive edge equations lie on paths and force integral endpoint values; without a zero, a nonintegral vertex requires all five edge inequalities tight, giving the all-half vector. Thus the twelve stated primitive rays are complete, including the second Hamiltonian ray from homogenizing the all-half vector.

Negation maps y_i to y_(-i-2), giving (a,b,b,a,c). Its five inequalities reduce exactly to t>=max(a+b,2b,a+c). The weight is 5t-2a-2b-c. Summing the five inequalities gives 2 sum(y_i)<=5t, so any flow of weight at most W has t<=floor(2W/5). This makes the independent finite evaluation of (18) exhaustive, not a parameter cutoff guess.

Nontrivial rotations fix constant x and y, which proves their series 1/(1-q^5)^2. The five reflections are conjugate under rotations, so their fixed-flow series equals the anti-coset average in this particular tournament. This does not justify using one anti-isomorphism for an arbitrary tournament. The checks also exhibit different sums of two primitive cycles producing the same R5 flow, rejecting a free twelve-generator interpretation.

No reduced rational numerator is needed or claimed by the frozen packet; any subsequent closed-form addendum requires its own audit.

## Inherited-input inspection

I inspected the stated source statements, conventions and supersession status rather than purporting to reprove the whole repository:

- `fixed_degree_metagraph_rationality_20261021/session22-fixed-degree/PROOF.md` states the fundamental paired-black degree convention, the finite multiplicative factorization, the exact center-deficit law, and the full degree-zero formula.
- `effective_black_defect_order_bound_20261020/PROOF.txt` states the uniform strong-order bound for bounded Delta=h-b. Because 0<=b<=h, fixed h gives a bounded deficit and hence a finite strong alphabet. K1 is handled separately from the bound's n>=3 statement.
- `metagraph_classification_20261021a/README.md` supplies the unique strong h=3 atom R5 and identifies the convention that its H=15, |Aut|=5 and b=3.
- The current supersession at `THM-2454-pure-blue-alphabet-and-center-law.md` explicitly makes the strong pure-blue alphabet K1,C3,R5 an all-order theorem and points to the separate rigid-strong and prefix-rerouting proofs. Historical text below that supersession is not the current theorem status.

The new proof's use of these inputs is faithful. The q=0 specialization retains precisely the zero decoration on each underlying type and therefore recovers the inherited formulas. q=1 is not a generally valid forgetful specialization.

## Nonclaims retained

A circulation is an arc-multiplicity vector, not an Euler traversal or a Hamilton word. Cross-SCC flow can have disconnected support and no single Euler tour. Neither Fano incidence nor rational flow counts establish metagraph adjacency, connectivity, bridge structure, or any unrelated open conjecture. The packet is careful about these distinctions.

## Reproduction and frozen-source identity

Run from any directory:

    python -B session23b-flows-review/independent_check.py
    python -B -O session23b-flows-review/independent_check.py

No assertions disabled by Python optimization are used. The source hashes are embedded in both receipts. Frozen producer SHA256 values:

- Initial PROOF.md: f65e28015417ad36c638416a08606bcb3e75a30ed99b4aca49fecaec611ecc14
- Final PROOF.md after reviewed precision edit: 253c95e582ab94ce4490e31ab115363ce41b74e96050d6fee1ba1241410b8ad0
- check.py: 7116893ed7e5fd94e9a4b39403ca33ce5a9d2a8f166dbc89168e17806e41ebaa
- normal.json and optimized.json: b032bc4e790516d37b072d24d2bf833c805d19ba87bc860a899b6a326d7ed25a

Only this review directory was written; no Git actions were performed.
