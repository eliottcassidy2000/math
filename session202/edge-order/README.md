# Two-port Hamilton-path order sectors

An independently audited actual C1 construction has arbitrarily small forward
boundary-order fraction, even though the boundary arc points forward. Its exact
sectors are f=2^p+2^q+1 and r=2^(p+q), for p,q>=1. The more general private-module
formula retains a product for the total count, rather than additive closure.

The count f concerns x before y, not consecutive use of x -> y. No eventual
odd-spectrum threshold, universal gluing inequality, or sharp global bound on
r in terms of f is claimed. The arbitrary-module formula does not assert C1;
unique Hamilton cycles are proved for the transitive-module family only.

Read `PROOF.md` and `audit/AUDIT.md`. The original reviewed inputs are preserved
byte for byte under `originals/`. Only status prose, the explicit sector warning,
and file-location plumbing changed during release assembly; see `PROMOTION.json`.
No outside theorem dependency is required for the case decomposition.

From this directory, run:

    python3 verify.py > VALIDATION_NORMAL.json
    python3 -O verify.py > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json
    python3 audit/check.py > audit/RELEASE_NORMAL.json
    python3 -O audit/check.py > audit/RELEASE_OPTIMIZED.json
    cmp audit/RELEASE_NORMAL.json audit/RELEASE_OPTIMIZED.json
    cmp audit/RESULT.json audit/RELEASE_NORMAL.json
    sha256sum -c RELEASE_FILES.sha256

The author DP covers 16 transitive pairs through order 11 and all 121 pairs of
labeled modules of orders 1 through 3. The independent literal-permutation check
covers those same 121 module pairs and 9 transitive pairs through order 9.
