# Publication review: one-pair triangle ellipsoid repair

**Independent proof and frozen-artifact review: PASS, 2026-10-21.**
The [triangle proof](README.md), [producer manifest](MANIFEST.json), and
its [hash sidecar](MANIFEST.sha256) are preserved byte-for-byte. Their
producer-stage status is historical. The current review is the
[combined independent audit](../projected_cube_ellipsoid_20261021/audit/README.md),
whose [manifest](../projected_cube_ellipsoid_20261021/audit/MANIFEST.json)
and [source pins](../projected_cube_ellipsoid_20261021/audit/SOURCE-PINS.md)
cover both this exact packet and the [base proof](../projected_cube_ellipsoid_20261021/README.md).
The base packet's [review wrapper](../projected_cube_ellipsoid_20261021/PUBLICATION_REVIEW.md)
explains its frozen historical status and the newly settled one-pair frontier.

Exactly one violating leaf pair admits the unique unrestricted maximum-volume
ellipsoid, including translated competitors. Together with the dominant
cone, this gives the complete positive three-speed classification. In
dimension at least four in the one-pair regime, no real diagonal lift
represents the optimum, and the full feasible nonnegative diagonal class
has a strict volume gap. Multiple violating pairs, a quantitative gap,
uniform LRC improvement, lattice incidence, and literature priority are
not claimed. Volume gain is not set containment or pointwise norm dominance.

## Reproduction

From this directory run `python3 -B check_exact.py` and
`python3 -B -O check_exact.py`; both outputs must equal [check_exact.out](check_exact.out),
with 31,344 exact gates. For the independent full replay run
`python3 -B ../projected_cube_ellipsoid_20261021/audit/replay.py` and its
`-O` counterpart. That audit contains both complete frozen source snapshots
and compares all normal/optimized runs with its retained receipts.
Its independent checker has 7,914 gates and imports no producer code.
All checks are standard-library-only; do not use `--record` for verification.
The written proof audit, not a finite enumeration or proof assistant,
establishes the universal claims for positive real speeds.
