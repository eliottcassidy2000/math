# Restricted theta optimum addendum

[PROOF.md](PROOF.md) certifies the known n=12 value for theta(5,5,2) embeddings with the two-edge path common to the two fields and the straight terminal chord inside the outer polygon. Convex outer boundaries satisfy the chord hypothesis. Other embeddings and nonconvex chord placements remain outside scope.

Reproduce the independent algebraic derivative and area calculation:

    PYTHONDONTWRITEBYTECODE=1 python3 verify.py verification.json
    PYTHONDONTWRITEBYTECODE=1 python3 -O verify.py verification.optimized.json
    cmp verification.json verification.optimized.json

This addendum does not modify or replace the base packet's exact geometry checker. The fixed-side cyclic-polygon maximum is classical cited input; the all-configuration conclusion is restricted by the explicit embedding hypotheses.
