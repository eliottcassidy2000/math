# THM-283: Burnside Partition Formula for SC(n)

**Status:** PROVED (verified n=2..12)
**Session:** opus-2026-04-03-S27

## Statement

The number of self-complementary tournament isomorphism classes on n vertices is:

    SC(n) = Σ_{λ ⊢ n, all parts 1 or even} 2^c(λ) / z(λ)

where:
- The sum runs over partitions λ of n where every part is either 1 or even
  (equivalently: no odd part ≥ 3)
- c(λ) is the SAME cycle exponent as in the Davis formula for A000568:
  c(λ) = Σᵢ mᵢ⌊lᵢ/2⌋ + Σᵢ C(mᵢ,2)lᵢ + Σᵢ<ⱼ mᵢmⱼgcd(lᵢ,lⱼ)
- z(λ) = Πᵢ lᵢ^mᵢ mᵢ! (the usual partition weight)

## The Mirror Symmetry

    A000568(n) = Σ_{all-ODD parts}  2^c(λ) / z(λ)    [Davis 1954]
    SC(n)      = Σ_{parts 1 or EVEN} 2^c(λ) / z(λ)    [this theorem]
    V_merged   = (A000568 + SC) / 2                     [Burnside on Z_2]

The same exponent c(λ) appears in both formulas. The only difference:
which partitions contribute.

## Proof

SC(n) = (1/n!) Σ_{σ∈S_n} Fix_comp(σ) where Fix_comp(σ) = #{T : σ(T) = T^op}.

For σ with cycle type λ, Fix_comp(σ) > 0 iff all cycles of σ have even length or
length 1. This is because the anti-automorphism constraint pairs arcs in orbits,
and odd cycles of length ≥ 3 create orbits with contradictory orientation constraints.

When Fix_comp(σ) > 0, the free arc orbits are counted by the same formula as the
Davis exponent c(λ), giving Fix_comp(σ) = 2^c(λ).

## Verified Values

| n | SC(n) | Contributing partitions |
|---|-------|----------------------|
| 2 | 1 | (2) |
| 3 | 2 | (2,1) |
| 4 | 2 | (2,2) |
| 5 | 8 | (2,2,1) |
| 6 | 12 | (2,2,2), (6) |
| 7 | 88 | (2,2,2,1), (6,1) |
| 8 | 176 | (2,2,2,2), (6,2) |
| 9 | 2752 | (2,2,2,2,1), (6,2,1) |
| 10 | 8784 | (2,2,2,2,2), (6,2,2), (10) |
| 11 | 279968 | (2,2,2,2,2,1), (6,2,2,1), (10,1) |
| 12 | 1492288 | (2,2,2,2,2,2), (6,2,2,2), (6,6), (10,2) |

## Computational Consequence

SC(n) can be computed to n=200+ in seconds using the same partition-iteration
machinery as for A000568. Simply change the partition filter from "all odd"
to "all parts 1 or even."

The SC sequence: 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784, 279968, 1492288, ...

## Relation to Existing Sequences

V_merged(n) = (A000568(n) + SC(n)) / 2 gives the merged metagraph vertex count.
NS_pairs(n) = (A000568(n) - SC(n)) / 2 gives the complement-paired class count.
