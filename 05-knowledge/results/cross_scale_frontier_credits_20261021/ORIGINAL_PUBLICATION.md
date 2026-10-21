# Current audited cross-scale credit packet

2026-10-21. See [the mathematical proof](README.md) and [independent final audit](review/REVIEW.md). The final producer and review files are retained unchanged. Superseded drafts, exploratory probes, and obsolete output snapshots are excluded.

The all-scale inverse compiler, stable three-phase scalar law, separate mod105 period-three/six strata, and permanent-capacity criterion have proved scoped statements. The minimum capacities 273 and 474, with 79 and 154 distinct interaction owners, are finite exact results for the two frozen sibling banks. Uniform native 153N coverage and an unbounded-scale remaining-remainder bound remain open.

## Portable replay

From this directory, replacing /path/to/math with the repository root:

    python3 -B check.py --repo-root /path/to/math --verify
    python3 -B -O check.py --repo-root /path/to/math --verify
    python3 -B review/check_independent.py --repo-root /path/to/math
    python3 -B -O review/check_independent.py --repo-root /path/to/math
    python3 -B review/check_handoffs.py --repo-root /path/to/math --producer-root .
    python3 -B -O review/check_handoffs.py --repo-root /path/to/math --producer-root .

Generation output normal.json records 4,143,337 producer checks. Verification includes one additional equality check against the saved receipt, so verification output matches optimized.json with 4,143,338 checks in either Python mode. No mathematical receipt changes between the modes. Audit generation and verification similarly distinguish the final saved-receipt equality check.

All original input dependencies are repository files authenticated by the checkers or their current receipts. The independent handoff checker imports only its companion independent core, not the producer implementation.
