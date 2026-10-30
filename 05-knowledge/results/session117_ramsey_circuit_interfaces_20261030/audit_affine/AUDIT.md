# Independent affine-Ramsey geometry audit — Session 117

**Verdict: PASS for the stated n ≥ 3 geometry claims.** This is an independent proof and literal-enumeration audit, not a literature novelty assessment. The exact B₂ gate theorem is imported and must pass its separate circuit audit.

## Proof review

On a flat avoiding the defect, the function is affine and is constant exactly when its direction annihilates r. On a flat through the defect, flipping one point of a constant affine function destroys constancy for every positive dimension. Flipping one point of a nonconstant affine function removes one of its two colors exactly when each color previously had one point, i.e. dimension one. Such a through-defect line has r·u=1 and resulting color c+1. Singletons are a separate, always-monochromatic case. These cases are exhaustive.

For each d-dimensional direction in ker(r), exactly one of its 2^(n−d) cosets contains the defect. When r is nonzero, the cosets split equally by affine color. Deleting the defect coset subtracts one from color c, not c+1. Adding the through-defect lines gives exactly the total and color-resolved formulas in the producer's note. For r=0 all positive-dimensional monochromatic flats have color epsilon. The singleton color formulas also hold, including n=0 for r=0.

For n≥3, a through-defect hyperplane cannot be monochromatic. A nonzero functional has only one kernel hyperplane direction, with exactly one coset avoiding the defect. A zero functional admits all 2^n−1 hyperplane directions, each with one avoiding coset. Thus the maximum monochromatic dimension is n−1 and the hyperplane counts are correct.

Translate the defect to zero: the function becomes c+r·x+delta_0(x). The invertible linear group is transitive on nonzero functionals, proving sufficiency of c as the orbit invariant. Its two values give distinct truth-table weights 2^(n−1)+1 and 2^(n−1)−1, proving necessity. No output complementation is required inside either orbit. For n≥3 the affine-plus-single-defect representation is unique because distinct affine functions disagree on at least four points, whereas two defect indicators disagree on at most two.

An invertible input map acts bijectively on points and affine flats, preserves inclusion and dimension, and transports the function values. Therefore the **abstract colored incidence structures** are isomorphic, not necessarily identical as coordinate-labelled structures. In epsilon=v=0, Walsh summation gives one coefficient 2^n−2 at r and −2 everywhere else; the claimed multiset equality follows.

## Low-dimensional exceptions

- n=0: only one point and one singleton; f=epsilon+1. There is no nonzero-r family and no hyperplane. The positive-dimensional formulas have empty scope.
- n=1, r=0: one monochromatic singleton of each color; the full line is nonconstant.
- n=1, r≠0: the defect makes the whole space constant, of color c+1. Its maximum monochromatic dimension is 1, not n−1. The total-flat formula including its line correction still gives the correct answer; the claim that full-dimensional color counts vanish must retain n≥2 (the producer assumes n≥3 throughout).
- n=2, r≠0: there are three monochromatic lines, all color c+1: the avoiding kernel coset and two exceptional through-defect lines. Hence uniqueness of the monochromatic hyperplane needs n≥3 as explicitly stated. Unique affine-plus-defect representation also need not hold in this dimension.

## Circuit bridge scope

Conditional on the independently audited companion theorem C(f)=n−1 for weight 0 or 1 and C(f)=n+wt(r)−2 otherwise, the nonzero-r normalized GL orbit has costs attaining every integer from n−1 to 2n−2. Consequently any orbit-only quantity that lower-bounds every member's ordinary B₂ complexity is bounded above by the minimum n−1. The maximum shortfall is n−1 gates and the maximum ratio is 2. This is an obstruction for this particular orbit, not a universal impossibility theorem about affine invariants or arbitrary Boolean families. The producer's note correctly limits the claim. No new explicit NP circuit lower-bound record follows.

## Independent executable evidence

`check_independent.py` imports no producer code. It constructs linear spaces by closing sets under adjoining vectors, constructs affine flats as literal cosets, and exhaustively tests every r, v, epsilon for n=0,1,2,3,4. It checks:

- 164,090 literal flat classifications;
- total and color-resolved counts, including singleton and low-dimensional full-space exceptions;
- all Walsh coefficients for every parameter triple;
- the normalized-family action of every invertible linear map in dimensions 0 through 3 (including all 168 maps in dimension 3).

Every check uses explicit exceptions, not assert. Normal and optimized Python both passed; `normal.json` and `optimized.json` are byte-identical. The proof, not finite testing, establishes the unrestricted-n statements.
