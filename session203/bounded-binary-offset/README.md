# Bounded binary offsets

**PROVED / independently audited in the stated scope.** For every fixed C,
only O_C(log^2 N) odd parent counts through N admit a physical append move
2n+r with |r|<=C, even after arbitrary target-preserving source rebasing.
An exact Fibonacci core bound gives a necessary polynomial-sized offset
for a scheme covering every parent in an interval.

This is an obstruction to an append-only binary construction. It does not
prove an eventual Hamilton-path spectrum threshold, establish compact-form
occupancy, or prove any target absent from unrestricted tournaments.

- `PROOF.md`: all-order proof and exact scope.
- `verify.py`: standalone author controls.
- `audit/AUDIT.md`: independent mathematical and computational review.
- `audit/check.py`: independent verifier, importing no author code.
- `dependencies/PHYSICAL_COMPILER.md`: reviewed physical construction.
- `dependencies/EXACT_BINARY_CHILDREN.md`: earlier +/-1 classification.
- `originals/`: frozen unpromoted author packet.
- `PROMOTION.diff`: status-only promotion after independent audit.

Reproduce from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/check.py
    python3 -O audit/check.py

The independent verifier uses `audit/PROOF-frozen.md` as its frozen proof
reference. Promotion changes only the status line. Author output hashes
therefore change; mathematical output does not.
