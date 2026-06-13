# Waggly as a Schurian Coherent Configuration

**Session:** kind-pasteur-2026-03-24-S20fu
**Date:** 2026-03-24

## The Framework

The waggly structure is precisely captured by a chain of three nested mathematical objects:

### Layer 1: The Hamming Scheme H(m, 2)
The ambient space Q_m = {0,1}^m with Hamming distance relations R_0, R_1, ..., R_m. This is a classical association scheme with Krawtchouk polynomial eigenvalues. Its automorphism group is S_2 wr S_m (hyperoctahedral group).

### Layer 2: The Schurian Coherent Configuration CC(S_n, Q_m)
S_n embeds into Aut(H(m,2)) via its action on 2-element subsets of {1,...,n}. The orbits of S_n on Q_m x Q_m form a coherent configuration — generally non-commutative (not an association scheme). Each Hamming class R_d splits into multiple S_n-orbitals.

### Layer 3: The Waggly Layers (Fusion)
Fusing all orbitals within each Hamming class R_d back together gives the waggly layers. Each waggly layer d is the d-th Hamming class projected to the quotient Q_m/S_n. This fusion is NOT an association scheme (verified computationally: residual 11-35% at n=4..7).

## The Central Paradox

Two opposing trends govern the broken Hamming scheme:

| n | m | V | Scheme breaking | Markov accuracy |
|---|---|---|---|---|
| 4 | 3 | 4 | 10.6% residual | P2 coeff = 0.63 |
| 5 | 6 | 12 | 22.9% | 0.83 |
| 6 | 10 | 56 | 31.8% | 0.90 |
| 7 | 15 | 456 | 35.5% | 0.94 |

**Structure becomes LESS regular** (residual grows) but **dynamics become MORE predictable** (P2 coefficient approaches 1).

This is the law of large numbers applied to the metagraph: more classes = more averaging = better Markov approximation, even as the individual class structure becomes more complex.

## The n=8 Catastrophe

Five independent structures simultaneously break at n=8:
1. SC backbone (spine) **fragments** (was connected through n=7)
2. Conflict graph **perfectness fails** (53.8% failure)
3. **beta_3 ≤ 1 breaks** (beta_3=2 found, beta_4 appears)
4. **Seesaw breaks** (beta_1 and beta_3 coexist)
5. Complement-flip overlap **returns** (was absent at n=7)

These are all manifestations of the SAME phenomenon: at n=8, the S_n quotient crosses a complexity threshold where the approximate regularities that held for n ≤ 7 no longer suffice.

In the Hamming scheme language: n=8 is where the Schurian coherent configuration becomes sufficiently non-commutative that the fusion to waggly layers loses qualitative features, not just quantitative accuracy.

## Connection to the Cayley-Dickson Tower

The n=8 catastrophe aligns with the Cayley-Dickson tower: n=9 = 2^3+1 is the octonion threshold where associativity is lost. The n=8 transition is the "pre-octonion" phase where the quaternionic structure (n=5) begins to fail.

The tower: R(n=2) → C(n=3) → H(n=5) → O(n=9) → S(n=17)

Each transition loses a property. The waggly scheme breaking parallels this: the Hamming scheme (perfect commutativity) progressively loses its algebraic structure as n increases through these thresholds.

## Key Theorems from the Framework

1. **Krawtchouk eigenvalues bound waggly layer spectra**: The quotient eigenvalues are averages of Krawtchouk polynomials over S_n-orbitals.

2. **Completeness order = diameter** (proved): k* = minimum d for full coverage equals the metagraph diameter. This follows from the triangle inequality on Hamming distance between orbits.

3. **MacWilliams-type identity**: The complement involution (d=m waggly) acts as a MacWilliams transform on the weight distribution of orbits.

4. **Three-phase filling**: The Schur-Weyl decomposition W = V_(n) + V_(n-1,1) + V_(n-2,2) predicts that the filling function F(k) has three phases: the trivial rep (captured at d=1), the standard rep (d~2), and the two-row rep (d~diameter).

## The Broken Scheme as Effective Theory

The waggly structure is best understood not as an exact mathematical object (it's not an association scheme) but as an **effective theory**: the Hamming scheme provides the approximate framework, and the S_n quotient introduces corrections that grow with n but become statistically irrelevant for the Markov dynamics.

This parallels how effective field theories work in physics: the "exact" theory (the Hamming scheme) is modified by a symmetry-breaking group action (S_n), producing a "broken" theory that nonetheless makes accurate predictions for macroscopic quantities (mixing times, spectral gaps, coverage rates).
