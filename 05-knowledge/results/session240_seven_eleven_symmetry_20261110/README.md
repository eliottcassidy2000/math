# Session240: seven and eleven, with the symmetry coordinate restored

2026-11-10. **COMPLETE. Independently audited; normal and optimized final replays PASS.**
Chan–Pak completeness remains **OPEN**. H counts ordered directed Hamiltonian
paths, not isomorphism classes or paths modulo automorphisms.

## Three precise advances

1. **Sharp strong elevenfold gap.** A strong tournament with an order-11
   automorphism either has exactly 11 vertices and H in
   {92411,93027,93467,95095}, or has H at least **1,848,759**.
   The latter bound is attained on 13 vertices. This extends Session239's
   sharp sevenfold threshold2163 using its inherited prime-action reduction,
   plus a new complete two-path-cover census. [Proof](gaps/PROOF.md) and
   [independent audit](audit_gap/AUDIT.md). The core minimum92411 is inherited.
2. **Exactly what the exceptional design symmetry forgets.** The Paley
   neighborhood designs admit exactly **8** compatible tournament matchings
   at seven and **12** at eleven. Each fibre is a single faithful
   2-transitive design-group orbit; stabilizers recover the actual tournament
   groups21 and55 from168 and660. The correct induced action on all marked
   Hamiltonian-path sets has the same quotient9 or1729 as the original
   tournament action. Larger symmetry buys no stronger raw-H divisibility.
   [Proof and complete order11 controls](paley/REPORT.md),
   [independent audit](audit_design/AUDIT.md).
3. **First joint carrier of seven and eleven.** An actual order-77
   automorphism in a strong tournament requires at least **19 vertices**.
   At that minimum order the sharp path floor is **3,049,061,631**,
   attained by C3[Q7,R11,1]. For a specified (7)(11)(1) action there are
   exactly120 minimizing orientations. [Proof](mixed/PROOF.md),
   [independent direct 19-vertex DP](audit_gap/AUDIT.md).
   This is a minimum-order-stratum result, not an all-order77 path floor.

## Recovered mathematics, not a numerical analogy

[Provenance and explicit transfer contracts](provenance/README.md) merge
Fano/biplane designs, projective versus affine actions, golden period-five
Collatz labels, binary codes, octonionic pairs, and marked norm tori.
Every connection names its map, preserved predicate, lost information,
sidecar, and hostile test. Sources are pinned by exact Git blobs.

A particularly useful inherited contrast is independently reproduced:
over F2, P7 adjacency is a rank4 projection, while P11 adjacency is
invertible. Their parity-extended image codes are [8,4,4] and [12,11,2].
Thus shared affine-group structure does not give the same coding mechanism.
The latter code is not binary Golay. Physical eleven-square marks and
Collatz arithmetic predicates also cannot be discarded.

## What remains open

- Arbitrary odd target realization in the Chan–Pak problem.
- A complete ambient equality classification for the elevenfold gap.
- The minimum H over **all** orders carrying an order77 automorphism.
- A general classification of compatible matching fibres beyond the tested
  Paley3/7/11/19 and all32 order11 circulants.
- Global Paley Hamiltonian extremality: no claim is promoted from HYP-9357.

All finite tables distinguish labelled orientations, matching fibres,
actions, and raw Hamiltonian-path counts. No literature-priority claim.
[Session plan](PLAN.md) records inherited mechanisms and the initial
interruption. [Reproduction](REPRODUCE.md) gives commands and verification scope.
