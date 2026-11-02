# Session137: overnight proof frontiers

## Checkpoint1: the Hamiltonian determinant bound is proved

**PROVED, independently audited.** For every nonempty tournament with adjacency matrix A,

    H(T) >= det(I+A)+|det A|+C_nonHamilton(T)
         >= max(1,s(T)) disc(T) >= disc(T).

Here C_nonHamilton counts directed cycles that do not span the tournament. Equality in the strongest bound holds exactly when no two directed triangles are vertex-disjoint. In particular H=disc exactly for transitive tournaments. This closes the formerly open target in THM1950 and THM3729.

[Complete proof](discriminant/GLOBAL_PROOF.md) · [Independent mathematical audit](independent_global/GLOBAL_DETERMINANT_AUDIT.md) · [Reproduction and checkpoint scope](CHECKPOINT1.md).

The proof combines the inherited odd-cycle collection formula with determinant cycle expansion and a weaker induction whose empty base is valid. It does not infer an all-order statement from a finite census. The empty tournament is deliberately excluded from the strong inequality.

The weaker bound H>=det(I+A)+|det A| also holds for real weighted tournament matrices with Aij in[0,1] and Aij+Aji=1. Its equality classification includes every weighted order-two matrix; for larger orders equality requires a deterministic transitive tournament or C3. The proper-cycle refinement uses integrality and is not exported to fractional matrices.

## Counting certificates

The positive recurrences also place four specific deficit functions, including H-det(I+A)-|det A| and its proper-cycle refinement, in #P through canonical polynomial-time certificate verifiers. [Counting construction and scope](global_referee/COUNTING_GAPS.md). This does not claim that every nonnegative GapP function lies in #P.

## Corrections and independent controls

The exact SCC discriminant product law remains true, but its old unconditional supermultiplicativity interpretation was false: a strong factor with disc6 and s2/3 has ordinal-square disc26, below36. Reducibility also does not force s below1. The canonical addenda correct those descriptions while preserving their provenance; the new global proof does not depend on either error.

Root checks cover33,867 ordinary and4,325 weighted matrices, with explicit wrong-sign, missing-factor, empty-domain and weighted-two-cycle controls. Separate implementations audit34,380 ordinary tournaments,42,532 root insertions and4,748,109 paths; rational weighted audits and a proper-cycle audit are also included. These supplement independent proof review, not replace it.

## Other overnight lanes

Sharp contextual modular signatures, source-cone/constructor rigidity, natural cycle-carrier classification, cross-parent editing obstructions, and harmonic owner-window extensions are being integrated in later checkpoints. They are not assumptions of this proof.

**The full odd-value spectrum conjecture and LRC(14) remain OPEN.** This checkpoint asserts no efficient exact H algorithm, formal proof-assistant certification, or external literature-priority claim.

[Concept board](BOARD.md).
