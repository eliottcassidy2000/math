# Rational cyclotomic CM observer: portable text-only supplement

THEOREM.md proves a coefficient-independent necessary finite-conductor bound for nonconstant rational observers of real cyclotomic generators whose output is a singular modulus. This supplement does not change the tangent-Legendre classification packets.

Reproduce:

    python check.py > results.normal.json
    python -O check.py > results.optimized.json
    cmp results.normal.json results.optimized.json
    python test_rejections.py > rejections.normal.json
    python -O test_rejections.py > rejections.optimized.json
    cmp rejections.normal.json rejections.optimized.json
    python verify_manifest.py

Only Python's standard library is required. All arithmetic is exact. The output includes N_d for d1..8, local unit-group formula checks, a d1 positive CM control, and four hostile boundary examples. The n3..500 group scan is explicitly only a sanity test; it does not justify the all-n theorem.

SOURCES.json is bibliographic provenance, with hash receipts for privately inspected articles. External PDFs are not included or required to execute the certificate. The manifest enforces the exact recursive file set and rejects missing, altered, or unlisted files. Independent mathematical and relocation audits are kept separately.
