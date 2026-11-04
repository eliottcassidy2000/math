# Independent nested-reversal audit

[AUDIT.md](AUDIT.md) verifies the exact weighted run-reversal bijection,
both parities of the anchor battery, and the nested/exposed-laminar scope.

    python3 verify_independent.py > receipt.json
    python3 -O verify_independent.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker imports no producer code and retains literal attachment counts,
so cubic and quartic effects cannot be accidentally discarded.

[UNIFORM_AUDIT.md](UNIFORM_AUDIT.md) independently verifies the stronger
uniform-module query compression, including a literal support-level involution,
compressed order two, and exact restored-module query counts.

    python3 verify_uniform.py > uniform_receipt.json
    python3 -O verify_uniform.py > uniform_receipt_optimized.json
    cmp uniform_receipt.json uniform_receipt_optimized.json
