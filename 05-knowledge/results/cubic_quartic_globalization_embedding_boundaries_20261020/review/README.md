# Independent review

**Verdict: PASS within the stated scope.**

- `AUDIT.md`: proof, quantifier, embedding, and primary-source audit.
- `REVIEW_SUMMARY.json`: reviewed producer hashes, inherited-source hash
  checks, exact-control counts, and output comparisons.
- `independent_checks.py`: independent finite-exact controls.
- `independent-normal.json`, `independent-optimized.json`: identical
  independent receipts, 56 checks.
- `producer-normal.json`, `producer-optimized.json`: identical producer
  replays, 74 checks; they equal the producer's frozen result.
- `producer-hash-verification.txt`: final producer-manifest verification.

`startup.txt` and `messages.txt` are local read-only startup receipts, not
part of the proposed research publication payload.

Reproduce from the packet directory:

    python3 review/independent_checks.py > /tmp/cq-review.json
    python3 -O review/independent_checks.py > /tmp/cq-review-O.json
    cmp /tmp/cq-review.json /tmp/cq-review-O.json
    python3 verify_cubic_quartic_boundaries.py > /tmp/cq-producer.json
    cmp /tmp/cq-producer.json results.json
    sha256sum -c SHA256SUMS.txt

All review writes were confined to this directory. No producer file or
shared repository was edited or published by the reviewer.
