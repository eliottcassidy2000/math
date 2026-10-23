# Shared endpoint pins: at most three further adaptive owner queries

**FINITE-EXACT supplement.** Retain the same four original source matrices, full actual path inventories, original electrical character, and original target atlas as ADAPTIVE_OBSERVERS.md. Suppose the two real endpoint owners `P[0],P[7]` are already known and shared with other consumers. They are fixed observations, not discarded labels and not newly charged queries.

For every nonempty endpoint fibre of all four sources, at most three further adaptive queries among clocks1 through6 identify the frozen phase, the actual merged target, or the whole literal Hamiltonian path. The exact distributions of minimum extra phase-query depth are:

| Source | Depth0 | Depth1 | Depth2 | Depth3 | Endpoint fibres |
|---:|---:|---:|---:|---:|---:|
|1238|1|13|29|6|49|
|2242|1|15|27|6|49|
|2246|1|11|40|4|56|
|2315|3|10|42|1|56|

The maximum three is attained in every source, so the bound is sharp for these individual complete source inventories. It is not an all-source or all-order theorem.

Knowing the endpoints does not make adaptivity redundant. At source1238, endpoint fibres `(7,3),(7,4),(7,6)` have23,22,25 paths respectively; their phase needs exactly four additional nonadaptive clocks, but only three adaptive clocks. There are further three-to-two savings, including11 endpoint fibres of source2315. At source2242, the phase's adaptive and fixed extra-query costs coincide on every endpoint fibre even though unrestricted full-source adaptivity saves two queries.

The endpoint pins remain literal ordered owners throughout. Any shared incoming constraint referring to those same owners is preserved by conditioning on its actual values. This is a read-only observer statement, not a modification of the path or a coupled-move theorem. An electrical phase guard may select a smaller fibre afterward; these query costs do not imply that its accepted move graph is connected.

endpoint_observers.py exhausts all210 nonempty endpoint fibres and all three observables. It uses independent tuple-state and bitset recurrences, and exhausts all64 fixed internal-clock banks for comparison. The638 checks survive optimized Python; the two receipt outputs agree byte-for-byte. The full per-endpoint answers are retained in endpoint.normal.json.

    python3 -B endpoint_observers.py > endpoint.normal.json
    python3 -O -B endpoint_observers.py > endpoint.optimized.json
    cmp endpoint.normal.json endpoint.optimized.json

This supplement does not alter the frozen main packet or its MANIFEST.json. Its own hashes are in ENDPOINT_MANIFEST.json.
