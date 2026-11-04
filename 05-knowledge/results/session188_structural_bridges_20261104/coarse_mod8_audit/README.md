# Coarse modulo-eight audit

See AUDIT.md. Run:

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

Dependencies: adjacent coarse_mod8_boundary source and receipt, mod8_transfer
proof, covered_quotient/verify.py independent primitives, and g++ for the local
literal.cpp physical-count control. No installation is performed.

Publication files: AUDIT.md, README.md, verify.py, literal.cpp, receipt.json,
receipt_optimized.json, SHA256SUMS.
