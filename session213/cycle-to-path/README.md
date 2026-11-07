# Exact path response of the cycle comb

The Session212 universal cycle-count family has an explicit quadratic Hamiltonian-path readout. Its negative pair coefficients count ordered two-path covers and depend only on bit-distance. This provides a fast exact inverse for the family, not a solution of Chan–Pak.

- `PROOF.md`: attachment identity, closed comb formula, complete length bound, and the induced-core-retaining compiler obstruction.
- `formula.py`: self-contained constructor, closed formula, and Gray-code support/inverse routine.
- `verify.py`, `RESULT.json`: independent subset-DP/permutation/cover controls and complete all-length family support through65535.
- `SOURCES.json`: pinned inherited source paths.

The direct comb family misses H=11 and35 at every order. Every odd h>=23 remains the ambient target; only this family has been ruled out as a complete direct readout.

Reproduce with Python3, no dependencies:

    python3 verify.py > RESULT.json
    python3 -O verify.py > RESULT.optimized.json
    cmp RESULT.json RESULT.optimized.json
