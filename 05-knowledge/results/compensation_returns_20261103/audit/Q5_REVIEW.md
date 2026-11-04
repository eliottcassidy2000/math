# Independent review of the single q5 carrier

2026-11-03. **VERIFIED fixed-carrier formula and literal witness audit.**
The claims in [`star/q5/PROOF.md`](../star/q5/PROOF.md) pass with their
declared scope: `TT8 -> C4`, the specified C4 rows `[8,1,3,6]`, the same13
physical names, and all4096 cuts around mark0. This is not an arbitrary-q5
classification and supplies no native program.

The independently compiled [`exact_counter.cpp`](exact_counter.cpp) gives
five full literal profiles: the deleted carrier, the four-vertex cap, the
unique strong target-count witness, its single-edge root repair, and its
unchanged actual source child. The C4 count is5 and P2 is32. The complete
fresh deleted E/S/J profile agrees with the producer's coefficient data.

Using those independently counted coefficients, the exact apex-last,
apex-first, and interior-join identity replays every one of the4096 saved
cuts. Independent reachability also checks strongness on every cut.
This is a formula replay, not a claim that4096 full profiles were computed
independently. The only strong target hit is mask3021, H1409; H1427 has no
strong hit. There are1920 strong cuts and960 with root99 end defect0.

For the actual target witness, every name, row, incoming/outgoing cut,
source-relative changed pair, E1/S1/E2/S2/J coordinate, deletion count and
fresh phase observer was independently checked. It has q0=5 and root99
end defect5. No vertex has end defect0 or start defect1/2. Thus selecting
another port cannot rescue this witness. Reversing {0,10} gives a strong
graph with H975, q0=5, and root99 end defect0: the guard repair is real,
and its count change is real. The original H1409 child remains the positive
control with q0=3 and root99 end defect0.

The final README's q1/q3 impossibility scopes and CAP_WIDTH's exact cap
feedback number, wrapper bounds, and growing-width consequence agree with
the audited proofs. One prose issue was reported for correction: complete
stored J matrices belong to the source/fixed-carrier target profiles;
the relocation receipt stores the two current ordered join values per
flip and its summary guards, rather than a complete J matrix at every
relocation stage. That storage distinction does not change any count or
proof.

Reproduce from the repository root:

```
python 05-knowledge/results/compensation_returns_20261103/audit/check_q5_review.py
python -O 05-knowledge/results/compensation_returns_20261103/audit/check_q5_review.py
```

The receipts `q5_review.json` and `q5_review_optimized.json` match
byte-for-byte and pin the reviewed producer proof/data and independent
checker sources. No producer file was changed, no new carrier was searched,
and no broad q5 obstruction follows from this review.
