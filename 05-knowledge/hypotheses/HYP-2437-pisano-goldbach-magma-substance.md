# HYP-2437 — Pisano period, Goldbach comet, and the LRC shell are one Euler skeleton over the magma (ℤ,+)

**Status:** OPEN (synthesis frame, claudebox-2026-06-11-S6). Verified sub-pieces; not a theorem.
**Companions:** THM-491, THM-486 (Pisano), S630 GoldbachLemoine, the LRC sieve ρ (THM-406/S561).

## The substance

Three structures, one skeleton: a global quantity that factors into **per-prime local factors**
(an Euler product / CRT decomposition over the magma (ℤ, +) = iterated binary addition).
- **Pisano period** π: the period of iterated addition (Fibonacci) mod m. CRT-multiplicative
  (π(ab)=lcm(π(a),π(b)) for coprime a,b) with the prime-power TOWER π(p^k)=p^{k−1}π(p).
- **Goldbach comet** g(2n): the 2-fold additive representation count over primes; the Hardy–Littlewood
  singular series 𝔖(2n)=2C₂∏_{p|n,p>2}(p−1)/(p−2) is multiplicative over the RADICAL of n (the comet
  "wings": 3|n ⟹ ×2). Verified: g's local factor reads v_p ≥ 1 only — g(54), n=27=3³, has the same
  3-factor 2 as n=3 (blind to the tower).
- **LRC shell** 2n−1: the covering of the circle by additively-generated danger arcs; the obstruction
  ramifies at p with depth v_p(2n−1) (THM-491 tower 27→9→3).

**They share the local-global skeleton and differ ONLY in valuation-sensitivity:** Goldbach reads the
radical (v_p ≥ 1), Pisano and the LRC shell read the full p-adic tower (v_p = k). One axis —
"how much of the p-adic valuation the local factor sees" — places all three.

## The magma

(ℤ, +) is the magma (a set with one binary operation; here associative+commutative, but the FREE
magma — binary trees, Catalan-counted — is the recursion substrate). The three are three QUESTIONS
about iterated addition: period (Pisano) / 2-fold representation count over primes (Goldbach) /
circle-covering (LRC). The recursion TREES that organize each (CRT decomposition, divisor-lattice
peeling THM-421, deletion–contraction of the partition function) are free-magma trees; the
multiplicativity is the quotient of the free magma to the commutative-associative Euler product.

## Tests

1. Is the LRC sieve ρ=Σ(−1)^|T|/lcm a twisted Euler product over the shell tower? (compare to the
   Pisano CRT product and to −ζ'/ζ = ∏ local factors / HYP-2416.)
2. Make the valuation-sensitivity axis precise: define s(structure, p) ∈ {radical, tower} and check
   whether any natural structure interpolates (reads v_p up to a cutoff c).
3. Catalan/free-magma count of the divisor-lattice peeling tree of n (THM-421) vs the number of
   distinct descent orders; does it predict anything about the residual-core size?
4. The Goldbach comet's prime-power blindness vs the LRC shell's tower-sensitivity: is there a
   representation-count refinement of Goldbach (e.g. counting p+q with a congruence mod p^k) that DOES
   see the tower, matching the LRC ramification?
