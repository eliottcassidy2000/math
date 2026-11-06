# Frozen audit and portable replay

`AUDIT.md`, `RESULT.json` and the two AUTHOR_RERUN receipts are copied without
changes from the independent auditor's packet. `original-check.py` is its exact
original script. The portable `check.py` changes only BASE to point to the
release's frozen `originals/` directory. It imports no author functions and uses
literal permutations, separately encoded edge predicates, and direct path cuts.
The release replay must match the original RESULT.json byte for byte.

The audit's historical workspace-root command is preserved in its original
report; the portable command for this package is `python3 audit/check.py`.
