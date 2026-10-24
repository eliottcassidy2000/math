# Quantitative blocks, with their resource boundary

**PROVED, scoped:** [Balanced charge blocks and the cost of retained precedence](BOUNDED_CIRCUITS_AND_PRECEDENCE.md).

Every prescribed balanced endpoint difference in a jointly compatible binary
U/V supply decomposes monotonically into blocks of at most 2D(G) changed slots.
Olson's classical p-group formula gives 16 for the original paired mod-8 code
and 264 for the full dyadic code. This is a bound on authenticated compatible
supplies, not a claim that the actual atlas realizes the abstract worst case.

A unit-capacity vacancy path shows why endpoint feasibility alone is too weak:
singleton implication components and the tiny group Z/2 can still force an
arbitrarily large simultaneous guarded move. The root
[endpoint-filter theorem](../interface/ENDPOINT_FILTERS.md) retains precisely
that missing precedence. The [native supply](../native/README.md) studies a
different escape: helper sources which start and finish unchanged can leave the
original endpoint interval and lower its bottleneck.

`check.py` verifies 210,034 conditions before frozen-receipt replay, using
explicit exceptions so the same checks run under `python -O`. It includes
exhaustive small group/supply universes, sharp abstract controls, literal
resource occupancy and hashes of the inherited native charge dependencies.

Reproduce from the repository root:

    python 05-knowledge/results/native_packet_decomposition_20261023/algebra/check.py
    python -O 05-knowledge/results/native_packet_decomposition_20261023/algebra/check.py

The independent root-interface audit is [AUDIT_ENDPOINT_FILTERS.md](AUDIT_ENDPOINT_FILTERS.md).
