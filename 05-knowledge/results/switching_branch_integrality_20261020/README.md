# Switching branch packet

Staged outside the repository for centralized publication and independent review. No repository file or Git state was changed by this worker.

- `QUOTIENT_WALSH_AUDIT.md`: standalone switching-quotient determinant identity needed by the uniform 3/5 equivariance obstruction.
- `BRANCH_INTEGRALITY.md`: all-order branch integrality, sharp n = 2 exception, cited Frobenius dependency, and general p-regular restricted-Burnside criterion.
- `PAIRING_BOUNDARY.md`: exact n = 4 failure of matching flip, n = 6 positive control, and inherited n = 8 impossibility of an involution on all 27 unadorned even-Aut types.
- `check_branches.py`: 158,744 odd-partition terms through n = 64, independently cross-checked by binary linear algebra on all 138 cycle types through n = 10.
- `check_matching_pairing.py`: direct labelled quotient/permutation enumeration at n = 4,6; separately verifies arithmetic of the inherited n = 8 census.

Reproduction from this directory:

    python3 -B check_branches.py --max-n 64 --linear-n 10 --out branch_checks.json > branch_checks.out
    python3 -B -O check_branches.py --max-n 64 --linear-n 10 --out branch_checks.optimized.json > branch_checks.optimized.out
    cmp branch_checks.json branch_checks.optimized.json
    cmp branch_checks.out branch_checks.optimized.out
    python3 -B check_matching_pairing.py > matching_pairing_probe.out
    cp matching_pairing_probe.json matching_pairing_probe.normal.json
    python3 -B -O check_matching_pairing.py > matching_pairing_probe.optimized.out
    cmp matching_pairing_probe.normal.json matching_pairing_probe.json
    cmp matching_pairing_probe.out matching_pairing_probe.optimized.out

All commands passed. Mathematical checks use explicit conditional exceptions, not assertions that disappear under -O. Finite controls are not the proof of the all-order theorem. No all-order structural pairing is claimed.
