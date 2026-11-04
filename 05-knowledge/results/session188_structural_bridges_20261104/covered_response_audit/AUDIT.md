# Independent audit of the equal-profile cyclic-quotient separator

**PASS.** Reviewed the final unconditional proof in `../covered_response/PROOF.md`, its Python verifier, C++ full-response engine, and retained data. The exact uniform query formula and its strict sign are valid. The common-profile hypothesis remains necessary to this derivation; arbitrary unequal profiles and arbitrary strong outer module contexts are not claimed.

## Formal run calculation

The old-run matrix N/(1-XYZ) has the correct start/end indexing. Every entry enumerates precisely the cyclic color words with the specified endpoints. Transposing N gives the converse quotient with the same physical child labels. Prefix and suffix factors allow an empty segment, while the middle segment is optional only when x precedes y. When y precedes x it is nonempty, because the fixed port edge is x->y. Thus the numerator accounts for all placements of both ports, including separated placements and all higher-degree attachment effects.

The six-permutation identity has the correct normalization: symmetrizing the numerator gives 2(1-XYZ)K, and the common-profile functional turns division by six into the factor one third in K/(1-XYZ)^2. The four monomial classes in K produce, respectively, q_k q_(k+2)^2, twice q_k q_(k+1)q_(k+2), minus q_(k+1)^2 q_(k+2), and minus twice q_(k+1)^3. Their factorization is exactly the stated negative cubic functional, with the series coefficient k+1.

The alternate aggregate cover skew has the right row/column orientation. Its factors are the column-1 times row-2 sum minus column-2 times row-1 sum of the old-run matrix. Its equality to the negative query gap requires the same common-profile symmetry; the proof does not extend it to unequal profiles.

## Positivity and physical scope

The standalone positive-square identity removes the log-concavity hypothesis. Its required convention is q_0=0, since the child blocks are nonempty. The upper boundary j=m+1 is retained. The independent audit of that algebra is in `../covered_positivity_review`; the present audit does not substitute finite sign tests for that proof.

The query is genuinely uniform on each child, so conversing any child leaves its integer value unchanged. It distinguishes the two quotient orientations even when each child has its own converse bit. Once the quotient orientation is known, the incoming disjoint-isolation theorem supplies one child decoder per block. The resulting four adaptive or seven fixed queries distinguish all sixteen states when the common block order is at least two. Equality of complete child profiles suffices; child isomorphism is not required.

The restriction on arbitrary ambient transfer is correct. A gap divisible by four cannot be lifted merely by the modulo-four substitution theorem. The result is for the cyclic substitution itself, with the stated ordered SCC padding. No full arbitrary covered-laminar theorem is claimed.

## Independent finite verification

`verify.py` imports only the independent primitives from `../covered_quotient/verify.py`, not code from the producer being audited.

- It directly enumerates bounded old-run words with the two ports each used once, for run cutoffs one through seven. After collecting equal-profile monomials, it verifies the entire cubic polynomial formula, not selected numerical values.
- Literal full attachments verify transitive profiles of orders one, two and three, plus a directed-triangle profile. Children with conversed interiors are included, so the implementation does not silently assume identical labelled orientations. The gaps are -2, -36, -2384 and -2052.
- It independently reinterpolates all 6,348 retained coefficients of the three nine-name response tables. Whole-converse quartics cancel; the quotient-converse difference has exactly 54 quadratic, 432 cubic and 540 quartic nonzero coordinates.
- Twelve dense physical queries are recomputed by an independently written recursive subset DP and agree with interpolation. The retained tables include named zero coefficients.
- Normal and Python -O outputs are byte-identical. Every check raises an explicit exception on failure.

The C++ forward recurrence and Python predecessor recurrence were also inspected. At eleven augmented vertices the exact counts and signed interpolation intermediates fit their integer types. These finite controls validate the implementations and indexing; the all-order proof is the run identity and positive-square algebra.

The receipt pins the final producer proof, verifier, C++ source, full retained response table, and independent helper. No correction to the final strengthened source was needed.
