# Exact inverse tournament sizes through count 661

Read PROOF.md for the complete finite classification, construction bounds, and open strong-first-occurrence question. minimum_orders.json contains all329 attainable counts in[1,661] and their exact first orders; each then occurs at every larger order. REVIEW.md and PACKAGING_REVIEW.md record independent review.

This is a targeted nine-witness extension of the exhaustive order-eight census, not an exhaustive order-nine census. census.json is copied unchanged from the reviewed size-spectrum release at commit d3b91e52c6b756d32ab90fe7f1a406f23633f60d. No new global completeness claim is made.

Reproduce from this directory with assertions enabled:

    g++ -O3 -std=c++17 fill.cpp -o /tmp/session207-prefix-fill
    /tmp/session207-prefix-fill < inputs.txt > /tmp/session207-prefix-fill.out
    cmp fill.out /tmp/session207-prefix-fill.out
    python3 check.py > /tmp/session207-prefix-check.out
    cmp check.out /tmp/session207-prefix-check.out
    python3 independent_audit.py > /tmp/session207-prefix-independent.out
    cmp independent_audit.out /tmp/session207-prefix-independent.out

check.py regenerates minimum_orders.json deterministically. It independently counts each new witness by all9! permutations. independent_audit.py uses a separately written memoized endpoint recurrence. manifest.json records SHA-256 digests for every other release file.
