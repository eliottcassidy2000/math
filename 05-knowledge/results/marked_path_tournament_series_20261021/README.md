# Marked-path tournament series and rooted Burnside repair

Status: PROVED in the stated marked-path domains; independent review PASS. FINITE-EXACT checks are controls, not substitutes for the all-order proofs.

The frozen [research proof](session23-symmetry/README.md) and [independent review](session23-symmetry-review/REVIEW.md) concern ordered directed paths on k distinct marked vertices, with n >= k >= 1. Converse reverses both arcs and marked order. The packet distinguishes this action from pointwise fixed roots and distinguishes normalized rational series from raw non-D-finite series.

All producer and review files are unchanged. The original correction preparation is retained as historical provenance: its PREPARED / NOT APPLIED status describes the preparation stage only. In this checkpoint the exact patch is applied to THM-260 and the mistakes ledger after checking its preimages. See [patch](session23-symmetry-correction/THM260-correction.patch) and [original/proposed hashes](session23-symmetry-correction/correction-receipt.json). The corrected canonical theorem is authoritative; historical source pins refer to its pre-correction version. Correct tables, numerical prefixes, banners, and unrelated interpretations are preserved. This is not a re-audit of those interpretations.

## Replay

From this packet directory with Python 3:

```sh
python3 session23-symmetry/check.py --receipt /tmp/marked-producer.json
python3 -O session23-symmetry/check.py --receipt /tmp/marked-producer-O.json
python3 session23-symmetry-review/audit.py --receipt /tmp/marked-review.json
python3 -O session23-symmetry-review/audit.py --receipt /tmp/marked-review-O.json
cmp /tmp/marked-producer.json session23-symmetry/receipt.json
cmp /tmp/marked-producer-O.json session23-symmetry/optimized-receipt.json
cmp /tmp/marked-review.json session23-symmetry-review/audit-normal.json
cmp /tmp/marked-review-O.json session23-symmetry-review/audit-optimized.json
```

The producer records 8,954 checks and the independent reviewer 499,356. The top-level manifest covers all exact payload bytes. No navigation changes accompany this packet.
