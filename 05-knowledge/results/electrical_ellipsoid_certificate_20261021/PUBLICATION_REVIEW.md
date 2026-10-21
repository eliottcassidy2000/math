# Publication review: electrical ellipsoid certificate

**Independent review: PASS, 2026-10-21.** The [proof](README.md) is preserved
byte-for-byte. Its opening “awaiting independent review” describes the
producer-stage snapshot, not the completed review recorded here. See the
[full independent audit](audit/PROOF_AUDIT.md) and its [frozen manifest](audit/MANIFEST.json).

The contact-network proof is a standard maximum-determinant duality
specialized to this projected-cube carrier, not a new general duality theorem.
It establishes exact optimality certificates, a matrix-tree determinant
bound, and the necessity of a cycle outside the dominant-speed cone. It does
not establish an arithmetic lattice hit or a new Lonely Runner case.

Companion proofs: [projected-cube/star optimum](../projected_cube_ellipsoid_20261021/README.md)
and [triangle contact repair](../triangle_ellipsoid_repair_20261021/README.md).
The grounding gauge, contact slacks, and arithmetic lattice are distinct data.
The graph orientation used in incidence matrices is a gauge, not a tournament.

Reproduce `python3 -B check.py` and `python3 -B -O check.py`, comparing both
with `receipt.json`. Independently run `python3 -B audit/independent_check.py`
and its `-O` counterpart, comparing with `audit/independent-normal.json`.
All four runs use only the Python standard library. Both retained audit
receipts contain 4,965 exact checks; the universal proof is audited separately.
