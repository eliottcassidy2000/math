# Equal-profile cyclic-quotient audit

See AUDIT.md. Reproduce with the adjacent covered_response packet and the
independent covered_quotient/verify.py helper present:

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker imports no covered_response producer code. It checks independent
run-word coefficients, literal attachments, complete retained interpolation,
and dense physical queries. Its receipt pins every relevant source.

Publication allowlist: AUDIT.md, README.md, verify.py, receipt.json,
receipt_optimized.json, SHA256SUMS. Preserve the adjacent helper dependency.
