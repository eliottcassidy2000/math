# Final publication review and portable verification

**Independently audited and approved for publication.** All six producer
files are unchanged. Their README's request for review records the earlier
submission state; the completed [independent audit](independent_audit/AUDIT.md)
and [homotopy supplement](independent_audit/HOMOTOPY_SUPPLEMENT.md) provide the
final reviewed status and additional proved scope. No source proof is rewritten.

The graph carrier is classical. The surface and normalization statements
retain their exact odd-characteristic hypotheses and the proved no-torus
corollary. No global tournament or Collatz conclusion follows. This is a separate continuation of the
[marked F5 model](../session17_prime_shapes_20261020/F5_PATH_PAIR_ICOSAHEDRON.md).

Historical frozen-directory commands in the audit describe its original
review environment. From this published directory, run:

    python3 -B verify_surfaces.py > /tmp/surfaces-producer.json
    cmp /tmp/surfaces-producer.json receipt.json
    python3 -B -O verify_surfaces.py > /tmp/surfaces-producer-optimized.json
    cmp /tmp/surfaces-producer-optimized.json receipt.json
    python3 -B independent_audit/independent_surfaces.py > /tmp/surfaces-independent.json
    cmp /tmp/surfaces-independent.json independent_audit/independent-normal.json
    python3 -B -O independent_audit/independent_surfaces.py > /tmp/surfaces-independent-optimized.json
    cmp /tmp/surfaces-independent-optimized.json independent_audit/independent-normal.json
    sha256sum -c SHA256SUMS

All four copied-layout runs match the frozen receipts byte-for-byte. Neither
checker imports producer or third-party code; finite checks corroborate the
printed structural proofs rather than extrapolating their all-order claims.
