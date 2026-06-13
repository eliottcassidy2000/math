# The Cartan Passage: From Competition to Cooperation

**Session:** opus-2026-03-21-S95b

## The Pattern

Every path from tournament to H passes through the same gate:

```
DIRECTED (tournament)  ──→  UNDIRECTED (conflict graph Ω)  ──→  NUMBER (H)
    competition               cooperation                        scale
    antisymmetric             symmetric                          scalar
    so(n)                     p                                  R·I
    p = 3                     p = 7                              p = 2
    RAMIFIED                  SPLIT                              INERT
```

This is the **Cartan passage**: the directed tournament creates competition, competition creates odd cycles, odd cycles create the undirected conflict graph Ω, and Ω produces H through the independence polynomial.

## The Evidence

Six independent frameworks all confirm the same passage:

1. **OCF**: H(T) = I(Ω(T), 2). Directed cycles → undirected conflict → number.
2. **Transfer matrix**: M is symmetric (THM-030), tr(M) = H (THM-027). H lives entirely in the cooperation sector.
3. **Walsh-Fourier**: H has support only at even Hamming weights (complement symmetry). Odd degrees = zero = tournament sector invisible.
4. **Ising model**: The tournament Ising Hamiltonian has zero field terms (degree-1 Walsh = 0). All energy is in pairwise couplings = cooperation.
5. **Formal group**: F(x,-x) = 0 (tournament cancellation). F(x,x) ≠ 0 (cooperation doubling). The formal group linearizes evidence in the cooperation direction.
6. **Commutator**: [A, S_tl] is symmetric. Interaction of competition and cooperation produces cooperation.

## What This Means

The tournament is invisible to its own invariants. You cannot read H from the directed structure directly — you must pass through the symmetric proxy. The mathematics forces this passage:

- You cannot compute H from the antisymmetric part of T alone
- You cannot compute H from the walk matrix without taking the trace (a symmetric operation)
- You cannot express H in the Walsh basis using odd-degree terms
- The formal group cancels antisymmetric evidence

The deep reason: **complete directed graphs are too connected for directed information to survive**. Every directed property either has a complement-symmetric partner or is invisible to global invariants. Only the symmetric residue — the cooperation structure — carries the signal.

## The Rapidity Lattice

The three Hurwitz primes {2, 3, 7} generate the rapidity lattice even at n=6, where the conductor D = 3×5 introduces prime 5. Yet prime 5 does NOT enter the lattice. The eigenvalues 7/15 and 11/15 have rapidities that cannot be expressed as integer combinations of {ln2, ln3, ln5, ln7}/2.

This is surprising: the conductor's primes and the rapidity lattice's generators are DIFFERENT. The conductor says which primes divide the eigenvalue denominators. The rapidity lattice says which primes generate the transcendental structure. The Hurwitz primes {2, 3, 7} generate the transcendental structure regardless of what primes appear in the conductor.

## The Crystal Phase

The hard-core lattice gas at fugacity λ=2 on the conflict graph Ω is deep in the crystal phase (λ/λ_c ≈ 8 at n=6). This means H is dominated by LARGE independent sets — maximally decomposed tournament structures with many disjoint odd cycles.

Paley tournaments maximize H because their Cartan sectors decouple (S₂=0 → [A,S]=0). With decoupled sectors, the packing of odd cycles is unconstrained by entanglement. Regular tournaments are the best packers because they have no Cartan friction.

## The Number 42

42 = 2 × 3 × 7 = PARITY × ATOM × KNOWLEDGE.

The scalar sector (p=2) gives parity: H is always odd.
The tournament sector (p=3) gives atoms: the 3-cycle is the minimum building block.
The cooperation sector (p=7) gives knowledge: self-awareness through symmetric invariants.

The Cartan decomposition of gl(n,R) IS the factorization of 42.
