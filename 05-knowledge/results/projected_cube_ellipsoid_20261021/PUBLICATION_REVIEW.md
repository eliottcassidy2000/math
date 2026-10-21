# Publication review: projected-cube ellipsoid optimum

**Independent proof and frozen-artifact review: PASS, 2026-10-21.**
The [base proof](README.md), [source pins](SOURCE-PINS.md), and
[producer manifest](MANIFEST.json) are preserved byte-for-byte. References
in those frozen records to pending review describe the producer-stage
snapshot. The completed review is the [combined independent audit](audit/README.md),
pinned by its [manifest](audit/MANIFEST.json) and [hash sidecar](audit/MANIFEST.sha256).
The audit's [independent source pins](audit/SOURCE-PINS.md) record the bounded
external attributions and historical repository observations.

The base theorem gives the unique maximum-volume ellipsoid, including all
translated competitors, throughout the exact dominant-speed cone. It also
classifies nonnegative diagonal boundary optimizers, records the two-speed
exception, and excludes every nonstar pair-facet spanning-tree certificate.
The source comparison concerns fixed-shape volume, not containment, uniform
norm dominance, a global Lonely Runner bound, or a literature-priority claim.

The retained [triangle addendum snapshot](audit/frozen_sources/session19-triangle-ellipsoid-repair/README.md)
settles the exactly-one-violating-pair regime and hence all positive triples;
it does not settle multiple violating pairs. That is the precise extension
of the older frontier statement in the frozen base proof. This audit stores
both complete source snapshots, so verification does not depend on the
separate triangle publication having happened yet.

## Reproduction

From this directory run `python3 -B check_exact.py` and
`python3 -B -O check_exact.py`; both outputs must equal [check_exact.out](check_exact.out).
The base companion has 71,504 exact gates. Run `python3 -B audit/replay.py`
and `python3 -B -O audit/replay.py` for the complete independent replay.
Each command verifies every frozen audit/source byte and runs the base,
triangle, and independent checkers in both normal and optimized modes.
The triangle companion has 31,344 gates and the independent checker 7,914.
Retained outputs are summarized in the [replay receipt](audit/receipts/replay.json).
All scripts use only the Python standard library and require no network,
repository, or original scratch directory. Do not use `--record` to verify.

Finite checks corroborate, rather than formalize, the universal real-speed
proofs. The audit explicitly retains the nonblocking fact that the base
producer's `grid_unique_diagonal_lift` branch has zero executions; analytic
uniqueness and independent equality controls provide the relevant evidence.
No lattice-point incidence, lonely-time witness, external computational-proof
certification, or improved global LRC constant is asserted.
