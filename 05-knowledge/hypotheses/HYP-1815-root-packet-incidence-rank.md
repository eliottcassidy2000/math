# HYP-1815: Root-Packet Incidence Rank Controls Transfer

**Status:** PARTLY CONFIRMED proof technology; two exact cross-domain controls
now realize opposite incidence-rank extremes.
**Source:** codex-2026-05-30-S6.  
**Related:** HYP-1796, HYP-1806, HYP-2430, HYP-2448, endpoint transfer, path
homology, `TournamentH7.RootPackets`.

## Statement

When support-level packet compatibility has the expected shape but a theorem
or reconstruction still fails, the missing invariant is rank or torsion in a
packet incidence matrix.

In short:

```text
support explains disjointness;
incidence rank explains effective transfer.
```

The relevant matrices may have rows indexed by vertices, arcs, endpoint
fibers, deletion fibers, or quotient buckets, and columns indexed by closed
root packets or problem-specific analogues.

## Evidence

- Endpoint-transfer work showed that private child witnesses imply rank, but
  merged quotients can keep rank after private witnesses disappear.
- Even-graph endpoint transfer can have plausible support matching while
  failing over `F_2`, revealing cancellation invisible to adjacency alone.
- Smith normal form sessions found torsion structure in transfer matrices.
- `RootPackets.lean` gives a formal closed-packet object whose total type-A
  boundary is zero; support and incidence can now be separated cleanly.
- Path-homology computations already live in boundary matrices and cokernels.
- Lonely Runner endpoint protection and Caccetta-Haggkvist return residue both
  turn into boundary-erasure problems when translated away from tournaments.

## Predictions

1. Packet incidence rank will separate THM-025-style near-kill examples from
   exact-kill examples better than support loss alone.
2. Real-root failures or small Newton margins should be enriched in packet
   incidence rank/torsion anomalies.
3. Endpoint-transfer private-pivot phenomena should have an Omega-packet
   analogue: some packet columns are private in a quotient even when support
   shadows collide.
4. Any artificial all-protected Lonely Runner endpoint core should show a
   nontrivial incidence-rank obstruction or quotient torsion profile.
5. Active ranking acquisition should be improved by expected packet-incidence
   rank collapse, not only expected `H` drop.

## 2026-07-22 exact realization

The GMC factor-recombination transfer and the `[72,36,16]` design gate now
realize the slogan with opposite exact matrices.

For separable `f` of degree `d`, HYP-2448's root-subset packet

```text
S |-> (sum_(alpha in S) alpha^k/f'(alpha))_(0<=k<d)
```

is a Vandermonde matrix times an invertible diagonal matrix, hence injective.
In characteristic zero its first `d-1` rows have only the full-root indicator
as kernel, so they still detect every proper Galois-stable factor packet.

By contrast, HYP-2430's five-incidence map from formal `16`-blocks on `72`
points to `5`-subsets has enormous kernel. Even after imposing the order-five
action of HYP-2441 and freezing every fixed block, its moving invariant kernel
has dimension at least `823261001634556`.

Thus support agreement transfers factor recombination because the weighted
root incidence has full rank, while the same style of marginal data cannot
transfer code-support realizability because its incidence map is massively
noninjective. The missing code sidecar is binary linear/matroid realization,
not another scalar or low-marginal identity.

## Next Tests

- Build vertex-by-packet and arc-by-packet incidence matrices for all small
  tournaments where odd cycles are explicitly enumerable.
- Compare incidence ranks with Omega support graph features, endpoint-transfer
  ranks, path-homology phases, and real-root margins.
- Compute Smith normal forms for packet incidence matrices in the standard
  contrast set: transitive, Paley, interval, H=63 single-core classes, H=37,
  and THM-025.
