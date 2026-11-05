# Session 200: canonical physical deficit cores

Proved in the stated scope, with an independent audit and reproducible
implementation. Universal coverage remains open.

- `PROOF.md`: sharp endpoint/height dichotomy; canonical sink-core projection;
  exact physical dyadic-ray joint fibers; constructive C1 target compiler;
  quantitative obstruction to uniformly small deficit.
- `verify.py`: complete finite inventories, minimum J costs, formula and
  hostile controls, literal graph tests.
- Frozen dependencies: `CARMICHAEL-DEPENDENCY.md`, `JOINT-DEPENDENCY.md`.

Reproduce:

    python3 verify.py --cap 5001 > VALIDATION_NORMAL.json
    python3 -O verify.py --cap 5001 > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The full {7,21} coverage theorem remains open. No external conjecture or
Diophantine-definability theorem is assumed.

## Audit and status promotion

The independent audit is `audit/REPORT.md`; its original inputs are recorded
in `audit/INPUTS.json`, and the pre-promotion proof and README are archived
beside it. Only status and audit-navigation prose was changed for this
release. The mathematics, verifier, and original receipts are unchanged.
The author receipts retain the scope string recorded before the audit;
`audit/REPORT.md` and `audit/PROMOTION.md` give the completed review status.

From the audit directory, reproduce the independent check:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

From this packet's root, `sha256sum -c MANIFEST.sha256` checks every release
input. All verifiers use the Python standard library and need no repository
checkout. Supplemental frozen physical-source dependencies are included in
`dependencies/` with their original bytes.

The sibling `../gap-one-core-continuation/` is explicitly unreviewed research.
It is not part of this audit or a universal coverage claim.
