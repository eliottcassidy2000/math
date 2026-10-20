# Independent review: generic Apollonius squareclass code

Status: **PASS**, with a minor graph-domain clarification.

Reviewed `APOLLONIUS_SQUARECLASS_CODE.md` from the 2026-10-19 Carmichael session. The argument is valid over the explicitly stated rational function field for labelled circles whose radii are part of the base field.

1. Independent symbolic calculation uses centres `(0,0),(a,b),(c,d)` directly, with no rotation or square-root base extension. If `Delta=ad-bc`, the discriminant is exactly the product of the three defects divided by `Delta^2`. This confirms both the squareclass and the absence of a hidden nonsquare constant.
2. Each defect is a rank-three homogeneous quadratic after an invertible linear coordinate change, hence absolutely irreducible. The six factors are pairwise nonassociate. Their distinct prime-divisor valuations prove independent squareclasses in the rational function field.
3. Independently constructing sign-class incidence gives rows `21,41,38,26`, with exactly the zero relation and the all-four relation. Thus rank is three. The induced action image is exactly the even-weight four-bit code, of order eight.
4. The four two-point orbits have fixed-point histogram `{0:1,4:6,8:1}` under the uniform sign group. This finite statement does not establish an error term or a uniform finite-field specialization distribution.
5. Generic signed radius is recoverable from the affine centre coordinates. Different sign classes cannot describe the same generic unoriented circle: after fixing the first input sign, changing the solution-radius sign forces the solution radius or first input radius to vanish; a change of any other input sign similarly forces a zero radius. Thus no generic pair merger invalidates the eight-point action.
6. Clarify that the graph generalization is for finite **loopless** graphs (simple graphs suffice). A loop-only vertex would invalidate the statement if a loop were counted once in its incident product. The supplied verifier already uses simple graphs.

Reproduce with `python3 check_generic_discriminant.py`; exact output is `apollonius_review.json`. All gates use explicit exceptions; normal and `python3 -O` outputs are byte-identical (`apollonius_review_optimized.json`). No external analytic input is used.
