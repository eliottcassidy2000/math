# Independent audit: finite inverse-graph selector and Euler gate

20 October 2026. **PASS, under the manuscript's stated global hypotheses.**

Reviewed manuscript: `session11-jacobian/FINITE_SELECTOR_EULER_GATE.md`.

Reviewed proof SHA-256: `29a2b95b756e1e7ed9895641344993ea5b7bbc445b8bdb84368afff528c6986a`.

Reviewed corrected producer verifier SHA-256: `4d2511b21016e71dbc3e328a30aec52ac17fb9a1bf3a1380e2c9add9cb5121b3`.

## Scope

The conclusions are conditional on an actual complex polynomial Keller map of generic degree d≥3, an affine source coordinate x primitive in the source function field, P monic in the companion source coordinate y, and polynomial traces Tr(x^k) for 1≤k≤d−2. The theorem does not establish those trace hypotheses for all Keller maps. No arbitrary local root packet or rational symplectic map is substituted for an actual affine-plane source.

The audit covers the finite graph, common denominator, omitted-point description, boundary fiber argument, Euler identity, Bézout bound, and selector poles. The independent local trace-lattice review is separate. Any subsequent source-line or boundary-classification note is outside this frozen audit.

## Proof obligations checked

1. **Inherited finite-sheet bridge.** The elementary argument in the prior `TRACE_POLE_BOUNDARY.md` was inspected: a target prime divisor cannot be wholly omitted, because its equation would pull back to a nonzero constant on affine space, contradicting dominance. Étale openness then gives an actual finite sheet generically. Monicity in y means an integral x value forces integral y. The low-moment maximum-pole argument therefore leaves precisely d−1 escaping sheets and one finite sheet on each nonproper component. This review does not depend on the prior packet's obsolete planar-degree-four frontier wording or its broader monodromy claims.

2. **Coefficient denominators.** Newton identities give polynomial c₁,…,c_(d−2). At a boundary divisor, c_(d−1) has the first polar trace's order, while substitution of the finite root bounds the pole of c_d by that same order. Thus c_(d−1)=A/H in reduced form has denominator support exactly the entire nonproper curve, and Hc_d=B is polynomial by normality/UFD intersection over affine prime divisors. Equality of the actual reduced denominators of the last two coefficients is neither needed nor asserted.

3. **Finite graph, without missing saturation.** With S=C[P,Q,x], the equation P(x,y)−u is monic in y, so C[x,y]=S[y] is finite integral over S. The primitive cleared equation R=Hχ is irreducible by Gauss's lemma; gcd(H,A)=1 suffices for primitivity. Its principal prime is the whole height-one kernel. Lying over makes the actual map onto V(R) finite and surjective.

4. **Exact omitted locus.** Specializing R at a target point gives a degree-d polynomial, a linear polynomial, a nonzero constant, or the zero polynomial, exactly according to H,A,B. The zero-polynomial case would make the target fiber contain a finite surjective preimage of an affine line, contradicting Keller quasi-finiteness. Hence O=V(H,A), the triple locus is empty, and the Nullstellensatz makes B a unit even in the possibly nonreduced finite intersection algebra. No generic-only conclusion is being promoted to an exact fiber statement here.

5. **Singleton fibers everywhere on the boundary.** Two distinct preimages of any boundary point produce disjoint local inverse neighborhoods. Their target neighborhoods overlap a generic smooth point of a boundary component, contradicting the generic singleton property. This uses global étaleness essentially. Restricting an étale local homeomorphism to the boundary gives a bijective open map onto C minus O, hence a homeomorphism. Nodal normalization alone is not a counterexample because it is not an étale map between ambient smooth surfaces.

6. **Euler and intersection budget.** Off C, F is a finite étale degree-d covering; on C minus O it is a homeomorphism. Compactly supported Euler additivity yields 1=d(1−χ_c(C))+χ_c(C)−#O, so #O=(d−1)(1−χ_c(C)). Since h=rad(H) and A have no common curve component, reduced affine cardinality is at most deg(h)deg(A). The normalization formula counts branch identifications b_p−1, rather than delta invariants, and is asserted only when the normalization is a disjoint union of affine lines.

7. **Selector poles.** At an omitted point B is nonzero and A vanishes along every normalization branch of C; coprimality prevents A from vanishing identically on a component. Thus −B/A has a genuine pole there, of the branch order of A, and no other finite poles. This does not produce a global polynomial inverse or close the higher-degree sector.

## Correction made during review

The original verifier's nonmonicity hostile used the companion equation uZ−v while naming the chosen finite primitive observable q for F(q,p)=(q,qp). The correct minimal equation of that observable is Z−u. The producer corrected the control: at (u,v)=(0,1), Z−u has the bounded root zero, while the companion relation up−v has no affine solution. The manuscript proof was unchanged. The control now demonstrates precisely why monicity/finite graph surjectivity cannot be omitted.

## Computational controls

The corrected producer verifier was copied into this audit directory and run under both normal Python and Python -O. Both runs passed 1,143 checks: 750 specialized coefficient fibers and 36 excluded coefficient atlases, plus other controls. Normal and optimized stdout and JSON are byte-identical, and the produced JSON matches the producer's final receipt.

The independent script `independent_controls.py` passes 64 exact controls in normal and optimized Python, with identical receipts. These include a nonreduced intersection-algebra unit inverse, the difference between reduced and multiplicity-weighted intersection budgets, actual non-Keller polynomial graphs, the corrected nonmonicity hostile, nodal normalization, and compactly supported Euler accounting. No large historical census was rerun. These controls test representative algebra and failure modes; they do not replace the universal proofs or construct any noninvertible Keller map.

## Literature and limits

The primary publisher page for Domrina's 2000 paper was independently checked during this audit: it states the nonexistence of four-sheeted polynomial maps C²→C² with nonzero constant Jacobian. The manuscript's degree-four correction is therefore appropriate. No claim that degree six is the first open generic degree is made.

Source: [A. V. Domrina, On four-sheeted polynomial mappings of C². II. The general case](https://www.mathnet.ru/eng/im273), Izvestiya: Mathematics 64 (2000), 1–33; DOI 10.1070/im2000v064n01ABEH000273.

Only the review staging directory was written. No repository, checkout, Git index, branch, or remote was changed.
