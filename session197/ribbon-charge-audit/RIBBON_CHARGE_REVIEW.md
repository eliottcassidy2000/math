# Independent review: the width-three charge-minus-two obstruction

**PASS.** This report is bound to the following frozen source snapshot in `../session197/ribbon-charge/`:

- `PROOF.md`: `9cf2326af8d4881ef777d5eb5e84f020fa565800cdf05ee34112d3adcd5ac88d`
- `verify.py`: `bfcebc3255fe65043e4ef36e10c25490535f3e788d4416df1df0c03c369248e2`
- `closure.json`: `376bcc84f28b996c69aac9fac86bedc3086c784f80a7d4934ef514965f09d622`
- `VALIDATION.json`: `2dbd10813130f668ca81d431a1395c60b9cc66471d8ecd3eb2b515d26c3b11eb`

No source files were edited. Copies of the source checker and its pinned dependencies were run under normal Python and `-O`; both stdout files exactly reproduce the frozen validation receipt.

## Proof review

The three readouts are correct. For R, the component beginning at the last root is uniquely the second path Q. The other path P is allowed to be empty. Hence there is no two-component ordering multiplier. No desired cover is lost by the transfer's restriction to at most one forgotten start, even at intermediate prefixes: Q's eventual start is the final newly added root; the only possible forgotten start belongs to P.

Coordinatewise saturation is exact as an abstraction because transitions have nonnegative integer coefficients. Retention by H<=13 guarantees that every H coordinate is below the cap and hence exact; e is an H sub-sum and is also exact. R is not presumed bounded in the retained state space. The tested identity requires R<=24, while one saturated R summand alone contributes at least 27. Therefore an overflow cannot produce a false charge equality. Conversely a true equality has all its R summands below the cap, so it is retained exactly. The only negative contribution, -2H, is already exact. This is the necessary headroom argument; subtracting two unbounded saturated observables would not have sufficed.

The finite certificate is a closed reachable state space, not a graph-order census. Initial-state inclusion, closure under every retained transition, and monotonicity of H under vertex insertion prove all-order coverage. Empty initial H=0 is properly excluded from the claimed nonempty spectrum.

The exact scope is a rooted tournament admitting a span-three ordering **with that root last**. The proof does not claim an exclusion for arbitrary marked roots or for scalar H=13. The displayed span-two H=13 graph correctly has charge -6.

## Independent finite controls

`charge-audit/independent.py` uses the previously reviewed transfer-state/transition generator, but none of the source checker's advance, closure, graph-profile, or readout functions. It converts the transfer to dense integer matrices, independently derives the three observable columns, and regenerates capped reachability by dense multiplication.

- Cap 27: exactly 2,382 vectors and 5,738 retained edges; the full vector set agrees with the supplied certificate.
- Cap 25: 2,382 vectors and 5,738 retained edges; the same compatible count set.
- Cap 31: 2,392 vectors and 5,768 retained edges; the same compatible count set.
- In all three cases, the charge-compatible positive counts are exactly 3,5,9,11.
- Independently generated all 4,683 ordered span-three graphs of orders 1 through 6. Literal vertex permutations and a cut immediately before the designated root independently count H, e, and R. Every result agrees with the frontier readout. This method shares neither source subset dynamic program.

The alternative cap checks are robustness controls, not the reason for all-order validity. The proof and the closed cap-27 certificate supply that validity. No mathematical or implementation correction is required.
