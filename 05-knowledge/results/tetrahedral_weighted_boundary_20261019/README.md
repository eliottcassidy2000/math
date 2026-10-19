# A4 weighted-boundary packet

The main artifact is [A4_WEIGHTED_BOUNDARY.md](A4_WEIGHTED_BOUNDARY.md), independently accepted with its stated scope.

It constructs the exact polynomial map

    (h,e) -> (9 h^4 e^2 - 8 h e, -27 h^6 e^3 + 36 h^3 e^2 - 8 e)

with Jacobian 64e and generic A4 monodromy. On A1 x Gm it is etale; its
finite normalization is b^2=ac, and its cubic ramification curve is rational,
unibranch, and forest-compatible. Two inverse trace moments vanish exactly.
The cone's Z/2 class group nevertheless excludes every same-field affine-plane
Keller chart, and hence every common polynomial source ring stable under the
two target-dual derivations.

The packet also proves the classical/scoped finite-factorial-domination gate
and gives a global-standard-cusp normalization classification. It explicitly
credits Chau's 1999 stronger exclusion when the entire nonproper-value set is
homeomorphic to C. No general A4 exclusion or Jacobian-conjecture solution is
claimed. Reflection and determinant-one binary-tetrahedral representations
are kept distinct.

Run:

    python verify_a4_weighted_boundary.py
    python -O verify_a4_weighted_boundary.py

Both runs pass 130 exact checks and have byte-identical stdout. The verifier
writes `verification-results.json`; it checks finite algebra and groups,
not a substitute for the separately written geometric proofs.

Read the [independent audit](review/REVIEW.md) and the separate
[inverse-trace predecessor](../inverse_trace_pole_boundary_20261019/README.md).

## Independent reproduction

Use Python 3 with SymPy (verified here with SymPy 1.14.0). From this package directory in a full repository checkout, run:

    python3 -B review/independent_checks.py
    python3 -B -O review/independent_checks.py

Both modes perform 32 explicit independent checks, including the seven repository dependency hashes in `SOURCE-DEPENDENCIES.json`. Compare stdout with `review/independent.out` and `review/independent.optimized.out`, and the generated JSON with `review/independent-results.json`. For the 130-check author verifier, both commands above in the Run section remain effective under optimization.

The main paper is byte-identical to the final accepted version. `SOURCE_PROVENANCE.json` records original and published hashes. Presentation changes remove internal staging prose and add live links and reproduction instructions. The independent checker's provenance paths are package-relative and its source-file list is explicit; its mathematical checks are unchanged.
