# Primitive cancellation, exact targets and leaf rigidity

This three-packet release is independently reviewed. Run `python3 audit_primitive.py` to replay every verifier under normal and optimized Python and compare the frozen receipts byte-for-byte.

- `PRIMITIVE_CANCELLATION.md`: two explicit physical modules give the ray 1998447616*3^m-5304875. It contains a certified prime count 5990037973 and a certified 7p count 161868952021=7*23124136003, with graph orders 26 and 29. The entire fixed-bank cancellation language is classified by three nonnegative copy parameters, and every ray has exact content 1. The large graph counts are derived from the audited arbitrary-module identity, not claimed as large subset-DP computations.
- `PRIMITIVE_TARGET_INVERSE.md`: a finite exact inverse decides any target count inside this complete language and returns every parameter recipe. The nearby prime 5990037979 is rigorously absent from the language.
- `PRIMITIVE_LEAF_RIGIDITY.md`: even when the transitive leaf length is unrestricted, this bank permits cancellation only at length 7. An elementary order 36 obstruction removes the other valuation mode at every leaf length.

The independent physical and arithmetic reports are included verbatim under `primitive-reviews/`. Their SHA256 hashes are respectively `dcd885e4a84fa845fef7790516a7dd462f96421473fec4259a0278375010fbb0` and `59dca4541772d50c6fe2a289d3d488ccb238501955261130e4b003f76984cba2`.

`PRIMITIVE_DEPENDENCIES.json` records the exact pre-existing proof and executable dependencies. This release does not change those files, any earlier frozen packet or the navigation guide.

The family is a constructive escape from a fixed-divisor obstruction. It is not an assertion of infinitely many primes, arbitrary prime targeting or full {7,21} spectrum coverage.
