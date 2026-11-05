# Independent review: root-preserving context composition

**PASS.** Frozen source in `../session197-odd-cycle-path-polynomials/`:

- `ROOTED_COMPOSITION.md`: `a8069ad433107862d73bb042757a50f8d2a170fb70c1db3b13b63151e46b06ed`
- `verify_rooted_composition.py`: `c768ccdcf6a66c88f8034e09897368ccbc5a4c0c35d75944a7de5bee736697d8`
- Copied normal/optimized stdout, both identical: `fd834fc61b31371fa555aafb3aa5d06d4feaf8c50c3ebff78cfb77ff43bc0c6e`.

No source changes were made. Dependencies are the previously reviewed universal context, matching-amplitude, and marked-context graph helpers.

## Proof and scope

Rooted module substitution is physically associative with the stated operator order: the inner Q operator acts before the outer P operator. At the retained original clone site, the two off-diagonal root-cross blocks vanish for opposite reasons. An outer in-neighbor dominates every inner vertex, while every inner vertex dominates an outer out-neighbor. The remaining root-cross graph is exactly the disjoint union of the two old graphs, so both matching-rank addition and maximum-matching multiplication hold, including rank-zero components.

For left-normal shift coefficients, S^i q_j(X)=q_j(X+i)S^i. The unique pair of shifts attaining the maximum sum is i=nu_P,j=nu_Q. Hence the top coefficient is p_top(x) q_top(x+nu_P), with the positive shift on the inner coefficient. Its degree and leading coefficient agree with the independently proved matching theorem. Positivity prevents cancellation of this term.

At scalar readout x=1, the highest exponential base is nu+1 with strictly positive coefficient; negative bases have smaller absolute value. An active outer template increases nu and therefore creates a genuinely larger asymptotic base. A passive template has its root as a singleton SCC and multiplies the full profile by B_(P-p)(x), hence scales every scalar clone coefficient by the same positive H(P-p). The negative-offset obstruction follows for a fixed finite word of outer whole-module templates at the same original clone site.

The proof does not rule out rerooting, splitting the original input's exterior neighborhood types, changing another vertex, or a target-dependent family of templates. The example's positive offset after E is correctly distinguished from a pure base-three return: the positive 4^m mode is still present. The formal m=0 zero-base term is explicitly retained where needed.

## Independent controls

`rooted-composition-audit/independent.py` checks the exact disjoint cross-edge sets and literal graph associativity on all 841 rooted-template pairs through order three. Matching counts are computed by edge-subset enumeration.

For 25 additional template pairs, including a rank-three composite, it obtains actual clone bin polynomials by literal directed odd-cycle enumeration and disjoint cycle packings. It recovers shift coefficients at each scalar x by exact rational Vandermonde inversion of those actual monomial images. This uses no source skeleton-operator compiler or noncommutative multiplication routine. All 175 top-coefficient evaluation comparisons satisfy the shifted product law.

For the eight-vertex clock, actual clone profiles at m=0,...,4 independently recover:

- At x=1 and bases -1,0,1,2,3: coefficients 0,0,-1,0,60.
- At x=2 and bases 0,1,2,3,4: coefficients 6,0,-176,480,630.
- Eight independent evaluations of the top coefficient agree with 3x+14x^2+23x^3+16x^4+4x^5.

A separate unrestricted-walk inclusion-exclusion counter verifies all twelve actual clock/D/E graphs for m=1,...,4. It shares neither source endpoint DP nor the literal cycle-packing algorithm. All twelve Hamilton counts match the displayed formulas.

No mathematical, orientation, boundary, or scope correction is required.
