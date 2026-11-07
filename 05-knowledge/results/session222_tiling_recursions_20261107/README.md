# Session222: tiling growth, doubling and the odd-chain plane

**IN PROGRESS: checkpoint 1 establishes the shell comparison and its exact Hamiltonian-count boundary. Chan–Pak remains OPEN.**

The user's starting picture joins an L along the two legs to two new hypotenuse diagonals, then asks how tournament +2 growth meets doubling and the plane of odd chains. The session starts from main `aa977c095d7651e1f27c3b626fa61bed75cf4ae9` and recovers the original coordinate definitions before comparing the operations.

## Same enlarged staircase, different graph operation

For the n-vertex fixed-path tournament carrier, both constructions add **2n−1 cells** and reach the (n+2)-vertex carrier. An explicit cell permutation maps the leg shell to the two diagonal layers. It preserves the Boolean cube, not tournament incidence or Hamiltonian-path count H.

The smallest decisive example is the directed triangle, H=3. Zero leg growth keeps H=3 by adjoining a source and a sink. Zero diagonal growth places its single reversed arc across a longer interval and gives H=9. Thus vertex order +2, tile count +2n−1, and a change in H are different quantities.

[Full proof, recovered sources and loss ledger](shells/REPORT.md) · [Coordinate diagram](shells/shells.svg).

On the one-reversed-arc family, adding one diagonal doubles H−1 and adding two quadruples it. This is a genuine doubling law with a stated carrier, not H→2H for arbitrary tournaments.

## Exact experimental question

Zero two-diagonal anchoring strictly increases H for all 32,767 nonzero framed order-seven tilings and for all smaller orders checked. Three thousand seeded nonzero order-eight samples also increase it. These are **FINITE-EXACT** bounded checks (the order-eight sample is not exhaustive). An all-order strict-monotonicity theorem or counterexample remains **OPEN**; no path injection is supplied.

Run `python3 shells/verify.py` and `python3 shells/probe_monotonicity.py`. Explicit checks remain active under `python3 -O`. The independent audit uses a separate subset-DP implementation. Later checkpoints will integrate the odd-chain plane, the physically nested doubling operation and corrected doubling-gauge statements.

## Research posture

Anchor: interpret the user's tiling recursions without conflating observables. Niche: recover the native odd-core/dyadic grid and its compatible tournament family. Wildcard: connect discrete cell traversal to a separately defined continuous plane-filling curve.

The live comparisons are cell geometry versus vertex incidence, count versus count difference, odd-core growth versus literal integer addition, path order versus integer labels, and finite evidence versus all-order construction. The inherited near miss is the corrected Mode-B boundary count in THM-291; the closest mechanism is source/sink preservation. Every result retains its missing coordinate and cheapest hostile example. No literature-priority claim is made.
