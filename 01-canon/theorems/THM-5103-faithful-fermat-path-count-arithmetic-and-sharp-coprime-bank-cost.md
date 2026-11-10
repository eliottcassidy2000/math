---
id: THM-5103
title: Faithful Fermat path-count arithmetic and sharp coprime bank cost
status: PROVED; independently audited. Elementary synthesis from inherited near-transitive counts and ordinal products. Finite controls are separate. No historical-priority claim.
source: 05-knowledge/results/session241_prime_relations_20261110/constructions/PROOF.md
session: Session241, 2026-11-10
depends_on:
  - THM-1862-order-join-reduction-principle
  - THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity
---

# THM-5103: faithful arithmetic inside the Hamilton-path spectrum

Let E_(a+2), a>=1, be the transitive tournament with its endpoint arc reversed.
The inherited path bijection gives H(E_(a+2))=2^a+1. Let p_i be the ordinary
primes in increasing order, starting p_0=2, and F_i=2^(2^i)+1. For
n=product p_i^(e_i), form T(n) by ordinal-joining e_i copies of E_(2^i+2),
in increasing i order. The unit may be the empty tournament with H=1, or a
transitive tournament if only scalar encoding is needed.

## Exact statements

1. Phi(n)=H(T(n))=product F_i^(e_i) is injective and preserves multiplication,
   numerical gcd and numerical lcm. It preserves and reflects divisibility
   and coprimality for all positive integers.
2. More generally, for any countable bank of integers w_i>1 and finite-support
   exponent vectors, product w_i^(e_i) reflects coordinatewise divisibility
   **iff each w_i has a private prime divisor**, absent from every other weight.
   It preserves gcd (equivalently lcm, or both directions of coprimality)
   iff the weights are pairwise coprime. Injectivity alone is weaker: it is
   equivalent to finite-relation rational independence of valuation columns.
3. For k>=1, a k-gadget endpoint-reversal bank has pairwise-coprime counts iff its
   positive exponents a_i have pairwise distinct2-adic valuations. Its total
   order is at least2^k−1+2k, with equality precisely for the exponent set
   {1,2,4,...,2^(k−1)}. Its maximum order is at least2^(k−1)+2.

## Mechanisms and boundaries

The classical Fermat product identity proves pairwise coprimality without
assuming Fermat primality. Disjoint prime supports transport exponent minima,
maxima and inequalities. Conversely, if w_i has no private prime, finitely
many other weights cover its finitely many prime factors; sufficiently high
powers make their product divisible by w_i while the i-coordinate stays zero.

The identity gcd(2^a+1,2^b+1)=2^gcd(a,b)+1 for v2(a)=v2(b), and1 otherwise,
follows from Bezout after dividing the exponents by their gcd. Sorting distinct
valuations gives a_j>=2^j, which proves the sharp bank cost. This is a cost
bound **within this gadget family**, not for arbitrary tournament encodings.

Actual tournament ordinal sum is noncommutative: T(m) joined with T(n) may
need component sorting to become T(mn). Only its H-count multiplication is
automatic. Input prime13 maps to composite F_5=641*6700417. Realized weights
15,51,85 are multiplicatively independent but do not reflect divisibility,
since15 divides51*85. E_3 is induced in E_4 although3 does not divide5.
Addition and ordinary numerical order are lost. This embeds relations into
actual path counts; it does not realize arbitrary prescribed targets.

[Full proof and tests](../../05-knowledge/results/session241_prime_relations_20261110/constructions/PROOF.md)
retain every quantifier and give a factorization-free decoder for promised
image values. [Independent audit](../../05-knowledge/results/session241_prime_relations_20261110/embedding_audit/AUDIT.md).
Chan–Pak completeness remains OPEN.
