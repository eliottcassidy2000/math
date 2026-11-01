# Verification receipt

- **PASS:** `PYTHONDONTWRITEBYTECODE=1 python3 verify.py > normal.json`.
- **PASS:** `PYTHONDONTWRITEBYTECODE=1 python3 -O verify.py > optimized.json`.
- **PASS:** normal and optimized receipts are byte-identical.
- **PASS:** all 22 generated indices/shards/tables/IR/reverse lists/metadata/fraction files rebuild byte-identically (`REBUILD.json`).
- **PASS:** both large tables are stored as authenticated plaintext indices and seven total row shards; every publication file is <=90,000 bytes. Shards concatenate ordinary transitions without expansion.
- Literal tables were exercised, not interpreted through high-level matrix instructions.
- 729 pure-blue windows in `{0,1,2}^6`, scratch 0, clock 3; continuous physical tape through checkpoint 24. Final output 60, 10,686 original fraction steps, **1,149,296 actual one-tape transitions**.
- Six actual compact-tail rounds (orders 6→12), one full-migration round, and two invariant-face rounds. The latter outputs 53 then 4213; round times 2,176 and **11,449,520 actual five-tape transitions**.
- Eight additional positive-matrix cases (zero, scalar, permutation, nonsymmetric, multiplicity and nilpotent controls).
- 1,764 reverse-compiler instruction cases, eight self-loop cases, halting control, five malformed optimized-compiler inputs, and nine general-list cases covering first-eligible priority, rollback, cancellation, inert factors, and halting.

The independent auditor additionally reported its own serialized-table implementation passes 1,764 instruction cases, 1,764 reverse-fraction cases, 210 original-fraction boundary checks and 1,000 arbitrary-list cases. Its separate source and receipts are maintained in the session independent-audit lane; this note does not replace them.

The tests are FINITE-EXACT. The unrestricted nonnegative-input and exact-cost theorems are proved in `PROOF.md`. No Lean or other proof-assistant certification is claimed. Hashes certify delivered bytes, not proofs.

Research-method cards used: “Verify the consequence, not the model’s own assumptions” (literal tables and original fraction integers) and “Use redundant paths as detectors” (separate interpreter, arithmetic formulas, normal/optimized replay and external independent auditor).
