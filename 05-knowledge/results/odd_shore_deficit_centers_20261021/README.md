# Odd-shore bounds for deficit-six centers

Current analytic result: degree R6 <=47. This improves inherited bounds but does not compute the full center polynomial or metagraph wiring. The proof retains exact marked-cycle rotation index, shore scores and source-bound prefix receipts. Earlier 56 and schema-defect stages are retained only as explicitly superseded audit lineage.

- [Current proof](session22-centers/PROOF.md) and [producer guide](session22-centers/README.md)
- [Second review: repaired checker and schema](session22-centers-review/SECOND_REVIEW.md)
- [Third review: current degree47 bound](session22-centers-review/THIRD_REVIEW.md)
- [Non-self-converse quadratic corollary](session22-root/NONSELFCONVERSE_QUADRATIC.md) and [separate review](session22-centers-review/NONSELFCONVERSE_REVIEW.md)
- [Independent trivalent recurrence audit](session22-fixed-degree-review/TRIVALENT_RECURRENCE_REVIEW.md)

The trivalent rational series has minimal eventual recurrence order18 and safe onset n>=60 under degree R6<=47. No earliest-onset or complete coefficient census claim follows. Small C3/R5/Paley7 controls and source-bound witness tests are not all-class enumeration.

## Reproduce

From this packet directory:

```sh
(cd session22-centers && python3 -B check.py)
(cd session22-centers && python3 -B -O check.py)
python3 session22-centers-review/repaired_audit.py
python3 session22-centers-review/repaired_schema.py
python3 -O session22-centers-review/repaired_schema.py
python3 session22-centers-review/score_run_control.py
python3 session22-fixed-degree-review/trivalent_recurrence_check.py
python3 -O session22-fixed-degree-review/trivalent_recurrence_check.py
```

Producer runs overwrite receipt.json with fresh actual clock fields; run in a temporary copy to preserve the published receipt. Compare scientific fields after excluding start_utc/end_utc/wall_seconds. The independent repaired_audit.py uses assertions and must run without -O. The strict-schema and symbolic recurrence controls use always-active checks.

## Packaging provenance

Original final source and review pins remain unchanged. Only two review drivers have publication adaptations: receipt paths now read the final producer receipt_normal.json/receipt.json rather than duplicated review snapshots. portable-paths.patch records the exact path-only diff, and portable-source-pins.json records original and published SHA-256 hashes. These adapted drivers are not claimed byte-identical to their staged originals. No copied checker snapshots, old normal/optimized directories, startup logs or scratch files are included.

Inherited proofs are linked from the producer proof and source_pins.json. This packet's source baseline is commit777c7240f75d3c574c6d0b91e8f6e65a70b2fde6. The top manifest pins every included file except itself. New proofs and controls do not claim recovery of lost Session21 audits.
