# Adversarial review of the session 9 plane-coloring extensions

Reviewed 2026-10-20 UTC. Scope: the staged binary-transfer note, finite-norm note (including the subsequently added eight-point forcing lemma), circle-antiperiodicity note, and their executable certificates. No checkout or Git writes were made. The supplied blueprint's major theorems and claimed Lean verification remain source claims; this review does not validate them.

## Verdict

The principal standalone arguments pass review. The q17 refutation, q23 rational enclosure, eight-point forcing, and finite circle checks all replay successfully. Two open-question descriptions needed correction, and one elementary counting consequence gives a useful strengthening and an explicit global separating inequality. These findings were sent immediately to the parent.

## Corrections and strengthening

1. The parent corrected the q23 placement total before completion of this review: 529 × 24 = 12,696. The executable calculation already used the correct value. The stale script comment was also flagged.
2. The initial H-F3 sentence claiming the single-circle test is inconclusive only for some q17 shapes was incorrect: F3 is shape-independent. The subsequent eight-point argument resolves forcing for every distinct q17 template, so the original open question is superseded.
3. In H-B2, the deterministic attainable defect vectors need not form a convex set under arbitrary placement measures. Atomic measures already give finite, generally nonconvex attainable sets. Discuss the attainable set and its convex hull separately. For exact terminology, its *closed* convex hull is the closure of the randomized defect region; closure should not silently assert exact attainment.
4. Let M(s) be the number of monochromatic unordered norm-one triangles at q17 and C(s) the number of cut edges in the norm-one graph. Every graph edge belongs to exactly two triangles, giving the exact identity

       M(s) = 1734 − C(s).

   Every degree is 18, so every cut has even cardinality, and M(s) is even. The checked refutation thus implies M(s) ≥ 2 and C(s) ≤ 1732. Consequently the explicit cut/correlation-polytope inequality

       sum over norm-one edges {i,j} of (1 − rho(i,j))/2 ≤ 1732

   is valid for every genuine sign correlation. The radial PSD hostile has rho(i,j) = −1/3 on all 2601 such edges, so its left side is 1734: a violation by 2. An explicit global separating inequality is therefore already available. The remaining interesting H-F2 problem is smaller support, sharper separation, or a simpler proof of the inequality.
5. With the preceding parity strengthening, the eight-point comparison yields

       D_T ≥ D_eq(a)/6 ≥ (2/1734)/6 = 1/5202.

   This is one ordered-placement quantum for the 5202 finite placements.

## Binary transfer and degree boundary

- B1 is valid for the stated jointly measurable bounded field and specified placement measure. A common positive separation permits a countable square partition whose cells separate every triple. Conditional independence gives the exact expected monochromatic defect Q/4; Tonelli and realization selection give the claimed loss-one bound. No stationarity is needed for this analytic rounding statement.
- The growing complete three-uniform hypergraph on 2m fixed points proves asymptotic sharpness of coefficient one in the arbitrary-placement formulation. It does not establish sharpness for Haar-uniform plane placements, and the note correctly separates those claims.
- The source-to-plane application correctly remains conditional on the projected field and its nonzero pair-correlation identity. Preserving a PSD kernel alone is insufficient. The higher-transcendence carrier caveat is necessary and is present.
- The countable zero-defect threshold argument, sharp factor-four threshold bound, exact Boolean Fourier degree, and four-bit parity obstruction are correct.
- The separation-free epsilon-loss claim is valid: the small-distance collision probability tends to zero for a probability measure supported on distinct triples. Exact attainment is not asserted.
- A concrete simultaneous-coordinate obstruction is available: use five fixed points, ten atomic triple objectives, and F identically zero. Each projected defect is 1/4, but every two-coloring makes at least one triple monochromatic, so that coordinate is 1. Hence a single deterministic coloring need not improve every coordinate.

## Finite norm geometry and certificates

- The trace-pairing Fourier normalization, the minus-Kloosterman multiplier, the common dual-norm parameter, F2, and the mean-refined lower bound are correct. The root-count derivation works over every odd finite field, with a fixed nontrivial additive character.
- The q23 enclosure is rigorous rational arithmetic: alternating arctangent bounds for Machin's identity, interval evaluation of the cosine polynomial, and a valid order-50 remainder (the order-49 coefficient vanishes). It certifies the unique maximizer c=11, the displayed rational upper envelope, and a lower bound exceeding 15/12696. Integrality gives at least 16 ordered monochromatic placements; equilateral multiplicity three gives at least six unordered triples.
- Independent adjacency reconstruction verifies the q5 and q11 witnesses and the q17 graph/hypergraph counts. The q17 verifier checks a complete binary refutation tree, every propagation reason, every conflict, and both branch assignments, with global sign complementation handling the root fix.
- The radial PSD hostile is a genuine real positive-definite rotation-invariant kernel normalized at zero. It is not globally Boolean-realizable. The finite refutation establishes this without any implication for a real-plane lower bound.
- The eight-point forcing lemma is valid over the claimed formal ring, including auxiliary coincidences. The six multipliers have norm one. Averaging is justified because composition by each fixed motion permutes the finite translation-rotation group. Scaling the coloring by a nonzero side vector transfers the equilateral obstruction to that side, hence to every distinct template.

## Circle and measurable-plane argument

- The common-frequency criterion is correct, including an arbitrary angle-index set and its per-angle exceptional null sets. One nonzero Fourier coefficient serves every angle; no uncountable null-set intersection is used.
- The rational 2-adic classification, odd/odd classification with pi/3, explicit odd closed walks, weighted defect inequalities, exact finite maximum-cut reduction, pinned-angle optimum, and irrational zero-infimum construction are correct.
- The measurable-plane proof is sound. The nested balanced-square argument supplies a nonempty null density boundary. The apex lemma gives a genuinely positive-measure family of monochromatic placements from two density-one base vertices. The two-circle map has nonzero Jacobian off a null critical set, so the self-supporting-null-set contradiction is valid.
- These statements give qualitative measurable-plane exclusion only. Compatibility does not construct a plane coloring, the circle defect infimum does not give a plane quantitative bound, and no nonmeasurable-plane conclusion is asserted.

## Replayed outputs

- `verify_q17.json`: PASS; 289 vertices, 2601 graph edges, 1734 triangles, 3451 nodes, 1726 conflict leaves, 125099 forced assignments; SHA-256 `3980c13e87296452da52c9d0f8a6da1a5bfbbf1130e9d8a19f1b00ec0052d9e6`.
- `kloosterman_bounds.json`: exact q17 and q23 enclosures; q23 maximizer 11 and placement count 12696.
- `circle_finite_cases.txt`: PASS for 92 reduced rational angle pairs with denominators 1 through 12, plus four full finite cut examples.
- `binary_controls.txt`: PASS for 729 rational Q checks, Boolean degrees through arity 8, four-bit parity controls, and q5/q11 positive controls.
- `small_coloring_replay.txt`: separate graph-based replay of both smaller-field avoiding colorings.
- `eight_point_forcing.json`: PASS for six formal motion identities, six norm-one multipliers, and all 256 eight-label assignments.
- `triangle_cut_identity.txt`: exact incidence and degree checks supporting the global cut inequality, plus 100 seeded finite sanity checks of the counting identity. The general counting/parity argument is algebraic, not inferred from the samples.
