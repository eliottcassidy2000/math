# Construct eleven adjacent counts

Explicit strong contexts produce all eleven odd integers from 36*2^k−7 through 36*2^k+13, for every k>=1, on at most k+7 vertices. More generally their counts are 36Q_U(2)+dH(U), for d=-7,-5,...,13 and any nonempty tournament U.

This is a proved fixed-core certificate plus an all-order bijection, independently audited. It does not settle the Chan–Pak conjecture or extend an unrestricted coverage record. It supplies a context-selection family, not an edge-flip staircase.

Read PROOF.md. Reproduce:

    python3 verify.py > replay.json
    cmp replay.json VERIFY.json
    python3 -O verify.py > replay_opt.json
    cmp replay_opt.json VERIFY_OPT.json

The verifier checks complete core-gap coefficients, displayed Hamilton cycles, and 891 literal composites with an independent subset-DP count. No inherited files are needed.
