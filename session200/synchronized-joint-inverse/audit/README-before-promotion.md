# Synchronized joint source inverse

Proof candidate and executable verifier for an exact square-root-height
inverse of `(G(s),2U-d)` on physical sources, with a physical Mersenne
boundary criterion for prescribed C1 targets. Universal coverage remains open.

Run serially:

    python3 verify.py > VALIDATION-normal.json
    python3 -O verify.py > VALIDATION-optimized.json
    cmp VALIDATION-normal.json VALIDATION-optimized.json

The default target cap is 10,001. The full forward bank through U=10,000
is complete for every C1 target at most that cap. All checks use explicit
exceptions. See PROOF.md for scope and constructive recipes.
