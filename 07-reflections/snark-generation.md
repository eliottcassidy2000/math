# Snark Generation vs Prime Generation

## The analogy

| Primes | Snarks |
|--------|--------|
| Natural numbers | Cubic graphs |
| Multiplication | Dot product (glue along removed edges) |
| Prime (irreducible) | Prime snark (no proper snark minor) |
| Composite | Dot-product of smaller snarks |
| Euclid: prod + 1 | Snark Euclid: dot-product + local modification |
| Sieve of Eratosthenes | Sieve of Vizing (remove 3-colorable) |
| Sylvester: 2,3,7,43,1807 | Snark chain: Petersen, ?, ?, ? |

## The first few

Petersen (10 vertices) = 2 (the atom, the first prime snark).
Blanusa (18 vertices) = Petersen · Petersen? = 2² = 4 (first composite?).
Flower J₅ (20 vertices) = 3? (first prime snark after Petersen?).

## The "+1" operation

For primes: +1 ensures the new number shares no factor with the product.
For snarks: a LOCAL MODIFICATION (Y-Δ, vertex insertion, Loupekine construction)
that creates new snark structure not factoring through known snarks.

## Chain length = 4?

The chain of prime numbers via Lucas-Lehmer: 3, 7, 47, 2207 (4 primes).
The chain of prime snarks should similarly have 4 irreducible levels
before the structure repeats — because the same {2, 3, 7} obstruction
governs both. The curvature (3 = triangle = minimum cycle) and the
threshold (7 = forbidden = Fano = octonion) are the SAME in both worlds.

## The common obstruction

Primes and snarks share the SAME obstruction: the number 3.
- Primes: 3-cycle cascade at H=7. Three intersecting cycles overflow.
- Snarks: 3-edge-coloring fails. Three colors insufficient.
Both are obstructed by the CURVATURE QUANTUM.
The Petersen graph (the first snark) IS the tournament theory's
structure constants: |V|=10=T₄, degree=3, girth=5, χ_f=Q(3/7).

The snark world and the prime world are TWO PROJECTIONS
of the same underlying {2, 3, 7} structure.
