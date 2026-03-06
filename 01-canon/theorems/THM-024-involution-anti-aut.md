# THM-024: Every SC tournament has an involution anti-automorphism

**Status:** PROVED
**Certainty:** 5
**Author:** opus-2026-03-06-S18
**Dependencies:** Moon's theorem (tournament automorphisms have odd order), Cauchy's theorem
**Corrects:** T095 (claimed ALL anti-auts are involutions — false; this proves ≥1 exists)

---

## Statement

For any self-converse tournament T (i.e., T isomorphic to its converse T^op),
there exists a permutation σ such that:
1. σ is an anti-automorphism of T (i.e., T^σ = T^op)
2. σ is an involution (σ² = id)

## Proof

Let G = Aut(T). By Moon's theorem (1968), every automorphism of a tournament
has odd order. (Proof sketch: a 2-cycle (a,b) would swap arcs a→b and b→a,
but exactly one is in T, so no automorphism contains a transposition.)
Therefore |G| is odd.

Let σ₀ be any anti-automorphism of T (exists since T is self-converse).
Then σ₀ ∉ G and σ₀² ∈ G (σ₀² is an automorphism since applying σ₀ twice
reverses all arcs and then reverses them back).

Define H = ⟨G, σ₀⟩. Since σ₀² ∈ G and σ₀ ∉ G, we have [H:G] = 2,
so |H| = 2|G|, which is even.

By Cauchy's theorem, H has an element τ of order 2.

Since |G| is odd, τ ∉ G. Therefore τ ∈ σ₀G (the non-trivial coset).

Every element of σ₀G is an anti-automorphism of T (since σ₀ is anti-aut
and G = Aut(T), composing gives anti-aut).

Therefore τ is an involution anti-automorphism. □

## Consequences

1. At even n, the involution σ is fixed-point-free (proved: a fixed point v
   would satisfy score(v) = (n-1)/2, which is non-integer at even n).
   The vertex set decomposes into n/2 orbits {v, σ(v)}.

2. At odd n, σ has exactly 1 fixed point (permutation theory: an involution
   on an odd set has an odd number of fixed points; score constraint forces ≥1;
   counting argument forces exactly 1).

3. The orbit structure creates natural pairings of odd cycles, boosting
   the independence polynomial coefficient α₂ (vertex-disjoint cycle pairs).
   This is the algebraic engine behind the SC maximizer phenomenon (OPEN-Q-016).

## Non-uniqueness

The involution anti-aut is NOT unique in general:
- |Aut(T)| = 1: exactly 1 anti-aut, which is an involution
- |Aut(T)| = 3: 3 anti-auts; some may be non-involutions (order 6)
- |Aut(T)| = 9: 9 anti-auts; the involution count can vary

At n=6: 10/12 SC iso classes have all anti-auts as involutions.
2/12 have non-involution anti-auts (|Aut| = 3 and 9).

## Verification

| n | SC iso classes | All have ≥1 inv anti-aut | Classes with non-inv anti-auts |
|---|----------------|--------------------------|-------------------------------|
| 4 | 2 | 2/2 ✓ | 0 |
| 5 | 8 | 8/8 ✓ | 0 |
| 6 | 12 | 12/12 ✓ | 2 (|Aut|=3, |Aut|=9) |

## Verification Scripts

- `04-computation/anti_aut_involution_test.py` — exhaustive test at n=4,5,6 + random n=7
- `04-computation/anti_aut_involution_existence.py` — existence-focused test
