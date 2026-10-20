# Independent review: rooted graph costs and finite probability patching

2026-10-20. **PASS within the stated finite scopes.** No blocking mathematical defect found. This review does not certify the supplied infinitary blueprint, a Lean build, geometric realizability of arbitrary star banks, or a new chromatic bound.

Reviewed payload: `session9-linear-graphs/README.md` and `check_linear_graph_tools.py`. The relevant inherited scopes in THM-3396 (four-bit Fourier cone), THM-2221 (pinned context cut response), and `boundary_palette_compression_20261011.md` were checked against their actual statements. The current note uses them consistently.

## Proof checks

1. **Star coefficient and triangular tomography: PASS.** For an even star with 2m leaves, its top Walsh coefficient is (-1)^m binom(2m−2,m−1)/2^(2m−1). The sign and factor of two agree with independent finite differences through 20 leaves. A star on T depends only on T; its Fourier support consists of even subsets of T, and its top coefficient is nonzero. Size-ordered elimination therefore gives a basis of all globally sign-invariant functions. The numerical measurements recover the symmetrized law, exactly losing odd moments. The count 2^(r−1)−1 applies for r≥1; the empty boundary is correctly handled separately.

2. **Three-root boundary: PASS.** Global sign reversal removes all odd Fourier degrees. On up to three roots only constant and pair terms remain. Separately minimizing each interior component proves the sufficient three-port-component condition. The four-star has a nonzero quartic coefficient, and the even/odd parity laws have identical proper marginals but different expected costs 3/2 and 1. This is a genuine minimal-root obstruction. The precise pairwise-sufficiency criterion is for all input laws, as stated.

3. **Min-plus gluing: PASS with the stated ownership guards.** Disjoint interior variable sets make minimization additive. A shared boundary edge must be counted once; a shared interior sign must be promoted to a named port. The text includes both hostiles. For simple graphs, the repeated-copy sentence must continue to be read under this preceding disjoint-edge condition: a root-root edge cannot become m distinct edges merely by giving copies different owners. In the blueprint's convention every edge meets the private interior, so this caveat disappears. No correction is required if the existing ownership/disjointness language is preserved.

4. **Finite replication/minimax theorem: PASS.** For every simplex weight rho, m_i=ceil(k rho_i) satisfies m_i≤k and sum m_i≤k+r−1. Nonnegative costs give the pointwise weighted comparison before minimization. The rounding residues lie in [0,1), have integer sum ≤r−1, and consequently their size-weighted sum is at most the sum of all sizes except the smallest. Finite minimax applies to deterministic kernels, because root words can choose their signs independently before averaging. Caratheodory gives at most r+1 global kernels in r-dimensional cost space. The proof keeps the shared-vertex charge C|W| and supplies an explicit slack, rather than claiming exact successor closure.

5. **Five-path hostile: PASS.** Opposite endpoint signs obstruct a proper coloring only when the complete four-edge path is present; deleting a root or edge removes that obstruction. Thus the old-root C=1/3 bound holds for every rooted edge subpattern. After promoting the middle vertex to a root, the two length-two costs are complementary, so both cannot have expectation at most 1/3. All one-copy inequalities pass and the (2,2) replication fails. The latter is a hypothetical extra-copy test, not an assertion that the original five-vertex ambient graph already contains these copies.

6. **One-sided bounds versus complete measurement: PASS.** The two-root even path and odd path have complementary costs. Thus star tomography as an invertible numerical instrument does not turn one-sided star inequalities into arbitrary cost inequalities. This is an important retained scope boundary.

## Computation independently replayed

- The author's checker was rerun normally: **672,505 checks**, output byte-identical to the supplied report. Its explicit `check` function remains active under `-O`; no verification relies on an optimized-away assertion.
- `independent_review.py` tests the integrated finite-patching theorem by an independent exact rational linear-program vertex solver. It uses two root words with masses 1/3 and 2/3, one shared sign, three constraints, and 600 reproducible sampled nonnegative-cost problems (seed 20261020). Of these, 119 meet every required bounded-copy aggregate hypothesis. All 119 satisfy the claimed tariff; 44 have strictly positive optimal excess, so the test includes genuine slack-dependent cases. The other 481 fail the hypotheses and are not counted as tests of the theorem's conclusion.
- The independent checker computes top star coefficients by successive finite differences for even sizes 2 through 20, separately from both of the author's coefficient methods.
- It checks all 89 combinations of endpoint-root subsets, interior subsets, and allowed edge subpatterns of the five-path. This independently covers the omitted-root statement, beyond the author's 34 full-endpoint-root tests.
- The independent checker is run both normally and with `-O`; its JSON reports agree byte-for-byte.

These finite checks support, and do not replace, the proofs. The sample of 600 general cost problems is explicitly a reproducible finite sample, not an exhaustive universe of all three-constraint games.

## Receipt files

- `author_checker_rerun.out`
- `independent_review.py`
- `independent_review.json`
- `independent_review.optimized.json`

No author payload was edited by this review.
