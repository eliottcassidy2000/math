# Equal-ten independent audit

Frozen publication allowlist: `REVIEW.md`, `README.md`, `audit.py`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Run `python audit.py > normal.json`, `python -O audit.py > optimized.json`, then compare. The independent mathematics uses only the Python standard library. Its producer replay invokes the existing C++17 compiler through the producer driver. No producer implementation is imported.

PASS: independent permutation-completion counts of all 64 order-four and 32,768 order-six tournaments agree with the entire source histograms, and every one of 1,085 order-ten simplex coefficients matches an independent physical-word/simplex-vertex reconstruction. The analytic review checks the 6+4 bound, profile domain, strict scalar term and actual physical arbitrary-ambient lift. Earlier orders are not recounted here.
