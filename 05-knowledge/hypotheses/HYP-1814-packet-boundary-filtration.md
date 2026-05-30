# HYP-1814: Packet-Boundary Filtration Of Hamiltonian Paths

**Status:** EXPLORATORY formalization target.  
**Source:** codex-2026-05-30-S6.  
**Related:** `TournamentH7.RootSigns`, `TournamentH7.RootPackets`, OCF,
HYP-1806, HYP-1808.

## Statement

Hamiltonian path counting should lift from the scalar OCF identity to a
filtration by endpoint-root fibers and compatible closed root packets:

```text
Hamiltonian paths
  -> endpoint-root fibers
  -> closed odd root packets inside each fiber
  -> compatible packet collections
  -> scalar evaluation at fugacity 2.
```

The scalar identity

```text
H(T) = I(Omega(T),2)
```

should be the decategorified evaluation of this packet-boundary filtration,
not the primitive object.

## Evidence

- `RootSigns.lean` proves finite root walks telescope to endpoint boundary and
  closed walks have zero total root.
- `RootPackets.lean` packages open walks and closed packets, and converts every
  existing `DirectedCycle T k` into a zero-root packet.
- OCF already says `H` is an independence-polynomial evaluation over
  vertex-disjoint odd-cycle packets.
- Endpoint-fiber theorems in Lean show Hamiltonian paths already decompose by
  start and end vertices.
- Paley/Interval data shows packet abundance and packet compatibility differ:
  more odd cycles can lose to more packable cycles.

## Predictions

1. Endpoint-root Hamiltonian path distributions should separate examples that
   have similar score or total `H` but different packet geometry.
2. A concrete odd-root-packet module should make parts of `alphaCount`
   de-opaquable before OCF itself is formalized.
3. Character-resolved packet ledgers for circulants should explain the
   Paley/Interval crossover better than root-sign spectra alone.
4. Real-root and fugacity phenomena of `I(Omega,x)` should become stability
   statements for the packet algebra.

## Next Tests

- Add `RootSupport.lean`: support of roots, walks, packets, and directed-cycle
  packets.
- Add an endpoint-root Hamiltonian path count extractor and compare transitive,
  Paley, interval, H=63, H=37, and THM-025 examples.
- Build a small packet algebra data structure in Python before attempting a
  Lean formalization of support disjointness.
