# Independent audit of the cubic positive-square identity

**PASS.** The exact identity, zero-extension convention, strict positivity, quantitative boundary lower bound, and verifier all check. The review independently expands the squares and reindexes every term; it does not infer an unbounded result from the finite runs.

## Mathematical review

Write a=q_(j-1), b=q_j, c=q_(j+1). The first square expands to j(a^3+3a^2b-6abc-3ac^2+3bc^2+2c^3). The second expands to (j+2)b^3-3b^2c-3jbc^2+(2j+1)c^3. The bc^2 terms cancel. On summing, the cube coefficient is 6j for q_j^3, with the sole exception of q_1^3, whose coefficient five is repaired by the separately displayed q_1^3. The remaining adjacent terms satisfy sum(3j q_(j-1)^2 q_j-3q_j^2 q_(j+1))=sum 3j q_j^2 q_(j+1). The remaining mixed terms are -6jabc-3jac^2. This is precisely three times the stated cubic functional.

The first square sum MUST include j=m+1. Its contribution is (m+1)q_m^3. The second sum at j=m contributes (m+2)q_m^3. Together with the separate q_1^3, these prove the claimed bound 3P >= q_1^3+(2m+3)q_m^3, where m is the last positive index. This also establishes strict positivity when q_1=0 but a later coordinate is positive. For m=1, equality holds and P=2q_1^3. Interior or trailing zeros cause no difficulty. The zero sequence has P=0. Signed sequences retain the polynomial identity but not the nonnegative-term argument.

The exact response application has the correct sign: F_T-F_U=-P for the physical query whose x mask is the first child and y mask is the second child. Equality of complete child profiles, not child isomorphism, is the required symmetry. Arbitrary unequal profiles and an arbitrary strong outer ambient module remain outside this positivity proof.

## Independent verifier review

The source checker has no disabled assertions: every condition uses an explicit exception. Its formal coefficient dictionaries include every boundary monomial. Its exhaustive sequence loop includes interior zeros and its signed controls test the identity independently of the positivity condition. Its normal and optimized outputs reproduce the stored receipts byte for byte.

The independent audit script imports no producer code. It uses exponent-vector polynomial multiplication rather than the source's sorted-variable-index coefficient insertion, verifies every cubic coefficient at cutoffs 1 through 30, and checks all sequences of lengths one through six over {0,1,2,3}, plus signed and boundary controls. It separately reruns the producer normally and with -O. The independent receipts agree byte for byte. Source hashes are pinned in those receipts.

## Reproduction

Run `python audit.py > receipt.json` and `python -O audit.py > receipt.optimized.json` from this directory. The adjacent `covered_positivity` packet must be present. This audit is distinct from the separate audit of the cyclic physical-query generating-function identity.
