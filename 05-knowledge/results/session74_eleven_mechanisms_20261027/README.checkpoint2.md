# Session 74 — final checkpoint: useful elevens and their missing coordinates

**Read the [final findings](session74-root/FINAL_FINDINGS.md) first.**
Research source: `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`.
This is the current guide; the original [checkpoint-one guide](README.md) and
all C1 native files are preserved byte-for-byte as a historical checkpoint.
Its then-open transitive distance is now exactly fifteen.

## Results and evidence

| Result | Proof / certificate | Independent review |
|---|---|---|
| Exact Paley11 distance 6 to forest frames, 15 to the transitive class | [Score proof](session74-paley/REPORT.md), [six moves](session74-paley/certificate.json), [fifteen moves](session74-paley/supplement/certificate15.json), [extension](session74-paley/supplement/REPORT.md) | [Forest audit](session74-score-audit/REPORT.md), [fifteen-move audit](session74-score-audit-supplement/REPORT.md) |
| All-period golden fibres and source-preserving ideal lifts | [Golden proof](session74-clocks/RESULT.md) | [Independent carry proof](session74-golden-audit/AUDIT.md), [prior-art boundary](session74-clocks/PRIMARY_LITERATURE.md), [final-byte review](session74-golden-final-hash-review/FINAL_HASH_REVIEW.json) |
| Weighted bracket counts and prime-power lifts | [Arithmetic proof](session74-arithmetic/README.md), [marked Paley interface](session74-root/WEIGHTED_PALEY_INTERFACE.md) | [Arithmetic audit](session74-root/ARITHMETIC_AUDIT.md) |
| Return multipliers 2, 3 and 11 swap, rotate or annihilate projected labels | [Three operators](session74-root/THREE_RETURN_OPERATORS.md) | [Golden audit, operator review](session74-golden-audit/AUDIT.md) |
| `R_15 = F_4 x Z/341`, conductor-11 axis `124a`, four-source kernel | [CRT proof](session74-clocks-supplement/RESULT.md) | [Independent CRT audit](session74-crt-audit/AUDIT.md), [final-byte review](session74-crt-final-review/FINAL-REVIEW.json) |

The bracket's original doubling prime set `{2,3,11}` is inherited, not newly
claimed. Weighted arithmetic is not prime counting. Golden word length is a
presentation period, not always a source's least period. No positive-integer
Collatz convergence, LRC(14) closure, all-order tournament accessibility or
historical-priority claim is asserted. The exact tournament endpoints are
classes, not every prescribed labelled target. The source11 axis is not the
75-element unit subgroup; dropping both nonunit axes loses correction36,
with12 supplied by the11-axis.

[The inheritance map](session74-inheritance/INHERITANCE.md) separately recovers
square-sum Q25's marked ports, the finite valuation-potential obstruction,
Terras period11 versus standard period18, Mills' pair and the spectrum count.
The [session plan](session74-root/SESSION_PLAN.md) is historical planning, not
a claim that every candidate was completed.

## Reproduce and inspect

From this directory, with Python3 and g++ (C++17):

```bash
python3 -B verify_packet.py --replay
python3 -O -B verify_packet.py --replay
python3 -B verify_provenance.py
python3 -O -B verify_provenance.py
python3 -B test_packet_corruption.py .
python3 -O -B test_packet_corruption.py .
python3 -B test_provenance_corruption.py .
python3 -O -B test_provenance_corruption.py .
```

The verifier chooses the latest additive manifest. All15 replay jobs run both
normally and with `-O` in a temporary portable packet; compiled binaries remain
temporary. Search programs are retained as provenance, not substituted for the
independent literal certificate checks. See [C2 verification scope](VERIFICATION.checkpoint2.md),
[dependency origins](DEPENDENCY_ORIGINS.checkpoint2.json), and
[audit bindings](AUDIT_BINDINGS.checkpoint2.json).

The [current canonical correction route](../../../00-navigation/S74.md) links
these proofs to THM-4592, THM-4810 and HYP-9420 without rewriting their historical
proofs, finite censuses, or original close-out letters.
