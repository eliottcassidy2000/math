# Tournaments at n = ∞

**opus-2026-03-25-S339**

## The Three Faces of Infinity

Lachlan (1984) proved: there are exactly **three** countable homogeneous tournaments.

| | (Q, <) | S(2) | T_ω |
|---|--------|------|-----|
| Name | Rational order | Dense local order | Random tournament |
| Structure | Transitive (acyclic) | Circular | Universal |
| Hamiltonian paths | 1 (the order itself) | None (all cycles) | Uncountably many |
| 3-cycles | 0 | Dense | Dense |
| Automorphism group | Aut(Q,<) | Related to PSL₂(R) | Simple, oligomorphic |
| Orbit count | 1 per n | Depends on n | A000568(n) |
| Limit of | Transitive T_n | Paley T_p (circular) | Random T_n |
| In metagraph | Source (H=1) | Regular (H max) | Generic point |

Our project lives in the space connecting these three.

## A000568 = Orbit Count of Aut(T_ω)

The deepest interpretation of our central sequence:

**A000568(n) = the number of orbits of Aut(T_ω) on n-element subsets**

Since T_ω contains ALL finite tournaments as induced subtournaments, the number of "types" of n-element subtournaments equals the number of isomorphism classes. This is Cameron's oligomorphic framework: Aut(T_ω) is an oligomorphic permutation group whose orbit-counting sequence is A000568.

Our entire project studies the structure of this automorphism group.

## Tournament Space → Cantor Set

| n | Tiling space | Dimension | Size |
|---|-------------|-----------|------|
| 3 | Q₁ = {0,1}¹ | 1 | 2 |
| 5 | Q₆ = {0,1}⁶ | 6 | 64 |
| 7 | Q₁₅ = {0,1}¹⁵ | 15 | 32,768 |
| ∞ | 2^ω = {0,1}^ℕ | ∞ | Cantor set |

The projective limit of Q_m as m → ∞ is the **Cantor space** 2^ω. Under the product topology:
- T_ω is a **comeager** (residual, topologically generic) point
- A random point is T_ω with probability 1 (the 0-1 law)
- The metagraph G_n becomes the quotient 2^ω / Aut(T_ω)

## What Survives at n = ∞

- **Rédei**: every countable tournament has a spanning ω-path
- **Erdős-Moser**: every infinite tournament contains an infinite transitive subtournament
- **Score sequence**: all vertices of T_ω have out-degree ω
- **Universality**: T_ω contains every finite tournament

## What Dies at n = ∞

- **Parity**: the odd-parity theorem is inherently finite
- **Burnside counting**: A000568(n) diverges; classes form a continuous space
- **Chromatic number**: χ = n-1 → ∞
- **Spectral gap**: 2/n → 0; the spectrum becomes continuous
- **Band-limitedness**: Walsh degree → ∞; all frequencies used
- **Fiber fraction**: f(n) → 0; no class has positive measure

## Tournamentons

The graph limit of a convergent sequence of tournaments is a **tournamenton**: a measurable function W: [0,1]² → [0,1] with W(x,y) + W(y,x) = 1 a.e.

The three homogeneous tournaments correspond to:
- W(x,y) = 1_{x<y} → (Q, <) (transitive)
- W(x,y) = 1_{0<y-x<1/2} → S(2) (circular)
- W(x,y) = 1/2 → T_ω (random)

## The Central Insight

Our project studies **finite-dimensional shadows of T_ω**. Every theorem we prove (band-limitedness, χ = n-1, fiber fraction, Krawtchouk structure) is a statement about the finite-dimensional approximations to Aut(T_ω). The limits as n → ∞ characterize the infinite group itself.

The metagraph G_n is a finite quotient of Q_m by S_n. At n = ∞, this becomes the quotient of the Cantor set by the oligomorphic group Aut(T_ω). The geometry of this quotient — its topology, its measure theory, its spectral theory — is what our sequences encode.
