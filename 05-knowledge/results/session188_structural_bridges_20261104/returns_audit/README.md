# Independent return-family audit

Both the transitive-family theorem and its arbitrary-tournament-module
extension pass after correcting the F/G prose definitions. See
[AUDIT.md](AUDIT.md) for the model, mask, context, and clock boundaries.
This is separate from the producer's verification packet.

Run from this directory:

    python3 verify_independent.py > receipt.json
    python3 -O verify_independent.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json
    python3 verify_modules_independent.py > modules_receipt.json
    python3 -O verify_modules_independent.py > modules_receipt_optimized.json
    cmp modules_receipt.json modules_receipt_optimized.json

The checker imports no producer module and enumerates literal directed paths.
It reads the sibling `returns` packet only to compare its published small
cases, source hashes, and stored normal/optimized outputs. The all-order
proof audit is mathematical; the finite experiment is not its substitute.
