# Independent audit of the same-parity Stirling-cone obstruction

**PASS.** Both displayed cubic coefficient vectors and both positive evaluations are exact. The examples invalidate the proposed higher-cover domination argument based only on positive same-parity transitive-profile cone membership. They do not invalidate domination for actual tournament profiles, establish an ambient cancellation, or give full named-response collisions.

## Proof and source review

The primary Grinberg–Stanley detailed PDF, dated June 13, 2026, was checked directly at https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf. Its Theorem 1.39 is on printed page 14 and supplies the stated positive odd-cycle expansion. Principal specialization gives powers indexed by the number of cycles; expansion in binomial polynomials gives the stated Stirling-profile combination. Converse invariance of ordered path-cover counts handles the directed-complement convention of the symmetric function. Same parity, unit highest coefficient, and even lower coefficients are valid necessary conditions, but not sufficient for realizability.

The ordered-cover functional is used within its proved equal-profile scope. Its polarization is algebraic: the three possible placements of one large transitive child yield the mixed coefficient after summing, not three individually equal literal contrasts. The proof makes that distinction correctly. Both artificial profiles violate the universal bound H<=n!, so neither can be mistaken for a tournament profile.

For the odd family, the cubic coefficient vector is (-2910041940,1746120,60,0). The value at 1582 is 2483340, while the value at 1580 is negative. Positive linear and quadratic coefficients prove strict increase for nonnegative weights, hence the stated first even crossing. For the even family, the vector is (-223898369213472,41383812864,834144,-192) and its value at 8192 is 65543165754336. No first-crossing claim is made for that family.

## Independent computation

The self-contained audit imports no producer code. It builds the symmetric numerator and old-slot polynomial independently, obtains transitive cover profiles by inclusion-exclusion for surjections, and reconstructs each entire cubic by summing its eight trilinear generator choices with exact Fraction arithmetic.

It independently builds all six physical cyclic substitutions, including the fixed port names and masks. A backward physical Hamilton-path recurrence computes every induced-subset path count. A canonical unordered-partition recurrence, whose first block contains the least remaining vertex, then computes complete unordered path-cover profiles; multiplication by k! gives ordered profiles. All twelve complete profiles agree with the source records.

The three odd polarized differences are 5940, -6120, 240, summing to 60. The three even differences are 69991968, -70991232, 1833408, summing to 834144. Thus the mixed coefficients are checked by actual physical examples, not inferred solely from the artificial profiles.

Independent normal and optimized outputs agree byte for byte. Fresh producer normal and optimized outputs also match their frozen receipts byte for byte. Every audit condition uses explicit errors and survives optimization.

## Reproduction and allowlist

Run `python3 verify.py > normal.json`, `python3 -O verify.py > optimized.json`, then compare the outputs. Only the standard library is required. Frozen publication allowlist: `REVIEW.md`, `verify.py`, `normal.json`, `optimized.json`, `SHA256SUMS`. Redundant producer-replay copies need not be published.

The full producer pins, all audited cubic coefficients and all twelve cover profiles are retained in the audit receipts. Reviewed proof hash: ee2c2be62f355ba6dae151a30b8943c4314936d2fac787e37a0c2ede0f012705.
