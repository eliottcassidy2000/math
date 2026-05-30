# Root Packets Feedback Loop

Session: codex-2026-05-30-S6  
Status: formalization plus exploratory synthesis  
Anchor: `TournamentH7.RootPackets`

## Executive Thesis

The root-sign formalization became more useful once it stopped being only a
metaphor for triangles.  The new formal object is smaller and more structural:

```text
open walk    -> endpoint boundary
closed walk  -> zero-root packet
directed cycle -> closed root packet
```

The session's useful abstraction is therefore:

```text
Do not start with the scalar invariant.
First expose the boundary, then the closed packets, then compatibility,
then incidence rank, and only at the end evaluate a scalar.
```

For tournaments this ladder looks like:

```text
root signs
  -> endpoint-root fibers of Hamiltonian paths
  -> closed odd-cycle packets
  -> packet conflict graph Omega(T)
  -> independence/packet algebra
  -> H(T)=I(Omega(T),2)
```

The new Lean module does not prove OCF.  It makes the first non-slogan bridge:
the existing `DirectedCycle` record now has a canonical zero-root packet
shadow.  That is the right kind of small formal step because later support,
disjointness, and packet algebra can attach to it without redoing the type-A
telescoping arithmetic.

## What Was Formalized

Added `TournamentH7.RootPackets`.

In namespace `TypeA`:

```lean
structure RootWalk (n : Nat) where
  source : Fin n
  target : Fin n
  middle : List (Fin n)
```

with vertex list `source :: middle ++ [target]`.  Its theorem is:

```lean
RootWalk.rootTotal_eq_boundary :
  W.rootTotal = root W.source W.target
```

So an open walk carries exactly its endpoint charge.

Also:

```lean
structure RootPacket (n : Nat) where
  base : Fin n
  middle : List (Fin n)
```

with vertex list `base :: middle ++ [base]`.  Its theorem is:

```lean
RootPacket.rootTotal_eq_zero :
  P.rootTotal = 0
```

So a closed packet carries no endpoint charge.

Finally, in namespace `Tournament.DirectedCycle`, the module defines:

```lean
rootPacketBase
rootPacketMiddle
toRootPacket
toRootPacket_rootTotal
```

The last theorem says every existing `DirectedCycle T k` induces a closed
type-A root packet with zero total root.  This is not yet a support theorem,
and not yet an odd-cycle theorem.  It is the exact bridge needed before those
can be meaningful.

## Why This Is The Right Boundary

The previous root-sign theorem proved:

```text
(e_a-e_b) + (e_b-e_c) + ... + (e_y-e_z) = e_a-e_z.
```

That is the incidence boundary map of the complete directed graph.  The new
module packages it in the language the rest of the repo already wants:

- open objects have endpoint residue;
- closed objects are kernel packets;
- cycles are not opaque objects, but certified kernel packets;
- support/disjointness is now a later layer, not mixed into the algebra.

This distinction matters.  A directed 3-cycle is not fundamental because it is
the smallest odd cycle.  It is fundamental because it is the smallest closed
root packet that can be odd and tournament-compatible.

## Repo Threads Re-read Through The Packet Boundary

### OCF And Real-Rootedness

OCF says:

```text
H(T) = I(Omega(T),2).
```

The old reading was "Hamiltonian paths equal a weighted count of compatible
odd cycles."  The packet reading is sharper:

```text
Hamiltonian path fibers are evaluated by a hard-core gas on closed
root-kernel packets.
```

This makes real-rootedness of `I(Omega,x)` feel less like a graph coincidence
and more like a stability question for a packet algebra.  The root-spectrum
of `I(Omega,x)` may be recording modes of packet compatibility, not merely
roots of an auxiliary polynomial.

### Endpoint Transfer And Good Cuts

Endpoint-transfer sessions repeatedly found that support matching is weaker
than incidence rank.  Private child witnesses imply rank; merged quotients can
lose private witnesses but keep rank; even graphs can have matching-looking
support while failing over `F_2`.

The root-packet reading says why: support is the shadow of a boundary map, but
rank lives in the map itself.  A support-only Omega cannot see all transfer
phenomena.  The next packet object should therefore have both:

```text
support(P)      vertices touched by a packet
incidence(P, v) signed or parity incidence at vertices/fibers
```

Disjointness is a support fact.  Cancellation, torsion, and private pivots are
incidence facts.

### Path Homology

The path-homology scripts already talk about boundary matrices, cokernels,
and surviving cycles.  Root packets are the tournament-side low-dimensional
boundary primitive: closed directed walks land in the kernel of the type-A
endpoint map.

This suggests a disciplined path-homology bridge:

```text
root packet kernel
  -> allowed-path boundary matrices
  -> Omega packet compatibility
  -> homology phase / torsion profile
```

The point is not that root packets solve homology.  They provide the first
shared chain-language object, so path homology and OCF no longer have to be
compared only through scalar features like `H`, `alpha`, or Betti numbers.

### Fejer, Paley, And Circulant Phase

In prime circulants, root signs collapse to signs on cyclic pairs `{d,-d}`.
Fourier characters read those signs by sine projection.  The all-one chamber
is the Fejer kernel; Paley is a multiplicative-character chamber.

The packet upgrade says that a good `phase_profile` should eventually be
packet-resolved:

```text
character channel of root signs
character channel of closed packets
character channel of packet compatibility
```

This may explain why Fejer concentration alone is not monotone for `H`.
Spectral concentration sees the sign field.  `H` sees compatible closed
packets.  The missing bridge is the exclusion diagram from spectral walk
packets to simple odd-cycle packets.

### Lonely Runner And Caccetta-Haggkvist

Lonely Runner endpoint protection and Caccetta-Haggkvist return residue are not
tournaments, but they share the same boundary grammar.

For Lonely Runner:

```text
open interval cover leaves endpoint boundary
all-protected endpoint graph erases every boundary witness
```

For Caccetta-Haggkvist:

```text
rooted expansion avoids return residue
first short cycle is the forced closed packet
```

The transfer is not literal type-A roots.  The transferable idea is:

```text
find the endpoint boundary map first;
then classify closed packets and their compatibility.
```

That is why endpoint-protection graphs and return-residue layer graphs feel
similar to Omega even though their vertices are not tournament odd cycles.

### Active Ranking And Reconstruction

A partial ranking is a partially filled chamber in the type-A root arrangement.
Querying a comparison chooses one root sign.  The operational objective should
not be "reduce uncertainty" in the abstract.  It should be:

```text
reduce endpoint-fiber ambiguity and collapse incompatible closed packets.
```

This sharpens the active-ranking hypothesis.  A useful product output would
not only be a ranking; it would report:

- score/divergence pressure;
- endpoint-root fiber uncertainty;
- closed packet witnesses of ambiguity;
- expected packet collapse after candidate queries;
- completeness defect when pair data is missing, tied, or weighted.

## A Boundary-Fiber Ladder

The feedback loop produced a ladder that should guide both formalization and
computation:

```text
0. root sign atom
1. open-walk endpoint boundary
2. closed root packet
3. packet support
4. packet disjointness/conflict
5. packet incidence rank/torsion
6. packet character/phase decomposition
7. scalar evaluation
```

Current Lean now covers levels 0--2.  Existing OCF axiomatics live mostly at
levels 4 and 7.  Endpoint-transfer work lives at level 5.  Fejer/circulant
work lives at level 6.  A lot of confusion in the repo comes from jumping
between non-adjacent levels.

The practical rule:

```text
When a scalar theorem is stuck, move one level down.
When a support theorem is false, move one level up to incidence rank.
When a residue theorem is too local, move sideways to character phase.
```

## Concrete Next Problems

### 1. `RootSupport.lean`

Define support for roots, walks, and packets:

```text
support(root i j)
support(RootPacket)
```

Prove basic facts:

- degenerate roots have empty or singleton-degenerate support;
- nondegenerate root support is exactly `{i,j}`;
- `RootPacket` support contains its base and middle vertices;
- `DirectedCycle.toRootPacket` support agrees with the cycle vertex image.

This is the missing bridge to Omega disjointness.

### 2. `OddRootPacket.lean`

Wrap a directed cycle with:

```text
length
oddness
root packet
support
zero-root certificate
```

Keep this separate from `alphaCount`, which is still opaque.  The goal is to
make "odd directed cycle" into a concrete packet object even if OCF remains an
external theorem.

### 3. Endpoint-Root Hamiltonian Fibers

Add a Python extractor:

```text
endpoint_root_counts(T)[a,b] = number of Hamiltonian paths from a to b.
```

Compare transitive, Paley, interval, H=63 single-core, THM-025, and H=37.
Predictions:

- transitive is concentrated on one endpoint root;
- regular/Paley is endpoint-balanced but packet-rich;
- H=63 single-core examples have concentrated packet residue even if endpoint
  fibers look less extreme;
- THM-025 should show a tiny endpoint residue that still has nontrivial
  packet compatibility.

### 4. Packet Incidence Matrices

Build matrices:

```text
vertices x packets
arcs x packets
endpoint fibers x packets
```

Compute ranks and Smith normal forms for small tournaments.  Compare:

- Omega support graph;
- packet incidence rank;
- real-root failure / Newton margin;
- endpoint-transfer private pivots;
- path-homology phase.

The working suspicion is that support explains disjointness, but incidence
rank explains when disjointness is algebraically effective.

### 5. Character-Resolved Packet Ledgers

For circulants, decompose packet counts by Fourier character orbit.  The
Paley/Interval crossover should then be phrased as:

```text
multiplicative flat sign phase
versus additive Fejer sign phase
after packet-exclusion corrections.
```

This is more precise than comparing spectra alone.

## New Hypotheses

### HYP-1814: Packet-Boundary Filtration Of Hamiltonian Paths

Hamiltonian path counting should lift from the scalar OCF identity to a
filtration by endpoint-root fibers and compatible closed root packets.  The
associated graded object should specialize to `I(Omega(T),2)`.

### HYP-1815: Root-Packet Incidence Rank Controls Transfer

When support-level packet compatibility has the expected shape but a theorem
still fails, the missing invariant should be rank or torsion in a packet
incidence matrix.  Endpoint transfer, path homology, real-root failures, LRC
endpoint protection, and CH return layers are all versions of this warning.

## Summary

The formal gain is small but real: directed cycles now have a certified
zero-root packet shadow in Lean.  The creative gain is larger: many unrelated
threads become instances of one boundary process.

```text
open object -> endpoint residue
closed object -> packet
many packets -> compatibility
compatibility with signs -> incidence/phase
scalar theorem -> final evaluation
```

The next formal target should stay modest: support and disjointness for
`DirectedCycle.toRootPacket`.  That would be enough to begin replacing opaque
cycle language with concrete packet language without pretending OCF itself has
been formalized.
