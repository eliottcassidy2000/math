# Fixed-presentation tournament renewal limit laws

Status: PROVED from the explicitly inherited finite-alphabet results; independent review PASS. FINITE-EXACT controls support the proof but do not establish asymptotic limits by extrapolation.

The frozen [proof](session23-root/RENEWAL_LIMIT_LAWS.md) and [independent review](session23-renewal-review/REVIEW.md) establish renewal and sampling results for fixed-presentation strata, not uniformly random unrestricted tournaments. Inherited SCC and finite-atom laws are dependencies, not new classifications. Synthetic controls are explicitly distinguished from actual labelled tournament enumeration.

Original source and review files are preserved byte-for-byte in sibling directories. Source-pin paths describe their original repository/staging provenance; no duplicated inherited source bank is required. There are no navigation or canonical changes in this packet.

## Replay

From this packet directory, using Python 3 with SymPy and mpmath available:

```sh
python3 session23-root/check_renewal.py > /tmp/renewal-producer.json
python3 -O session23-root/check_renewal.py > /tmp/renewal-producer-O.json
python3 session23-renewal-review/check.py > /tmp/renewal-review.json
python3 -O session23-renewal-review/check.py > /tmp/renewal-review-O.json
cmp /tmp/renewal-producer.json session23-root/normal.json
cmp /tmp/renewal-producer-O.json session23-root/optimized.json
cmp /tmp/renewal-review.json session23-renewal-review/normal.json
cmp /tmp/renewal-review-O.json session23-renewal-review/optimized.json
```

The producer records 2,164 controls. The independent checker records 8,928 controls, including all 33,867 labelled tournaments through order 6, 11,747 synthetic words, and 4,786 sampler words. The packet manifest covers all published payload files with exact byte lengths and SHA256 hashes.
