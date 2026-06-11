# The ramification tower, the freed clock, and the magma substance (claudebox-2026-06-11-S6)

User dispatch: work LRC back and forth between n=14 and n=19; see everything recursively; separate
off small aspects predictable by formula; find how the Pisano period, the Goldbach comet, and magmas
get at one underlying substance.

## n=14 vs n=19, made recursive

The two frontier cases differ in exactly one coordinate — the **prime-power structure of the shell
2n−1** — and that coordinate is recursive:
- **n=14, shell 27 = 3³ (ramified, depth 3).** The runners split into units of ℤ/27 (coprime to 3,
  dodged by the doubling orbit, ord₂₇(2)=18=φ(27)) and non-units (the multiples of 3). The non-unit
  core, divided by 3, IS the shell-9 problem; its core ÷3 is shell-3. So the difficulty descends a
  3-adic tower 27 → 9 → 3 — and this is the SAME recursion as the Pisano period: π(3ᵏ)=8·3ᵏ⁻¹
  (8,24,72). n=14 is the smallest n with depth ≥ 3; it is also composite (2·7), so it is two-headed
  (the divisor-clock peeling THM-421 on one side, the 3-adic shell tower on the other).
- **n=19, shell 37 (prime, depth 1).** No tower. Both n and 2n−1 are prime, so there is neither a CRT
  fiber nor a ramified core. The recursion is cyclotomic: (ℤ/37)* is a single 36-cycle under doubling,
  and the only hard configs are ±-transversals (THM-420) — of which a random sample has zero. Clean.

The lesson: **the LRC frontier is indexed by ramification depth v_p(2n−1)**, and n=14's notoriety is
not bad luck — it is the first shell deep enough (3³) to need more than one descent.

## A small aspect, separated off and predicted exactly (the freed clock)

Perturb the tight AP {1,…,n−1} by dropping runner j (and adding n). If no other speed is a multiple
of j, the j-clock t=1/j is a witness and M jumps to exactly 1/j. But if j | n, the runner n stays a
multiple of j and BLOCKS its own clock. So the composite/prime nature of n is visible in the very
smallest perturbations: n=14 has self-blocked clocks at its divisors {2,7} (drop 7 → 1/11, not 1/7);
n=19, prime, frees every large clock cleanly (drop j → 1/j for all j=10..18). An exact rational,
predicted from j and the divisibility of n — the kind of formula-predictable sub-structure the brief
asked for, and it re-expresses the same composite-vs-prime split as the tower story.

## The common substance: one Euler skeleton, three valuation-sensitivities, over the magma

Pisano period, Goldbach comet, and the LRC shell are three readings of ONE local-global skeleton: a
global quantity factoring into per-prime local factors. They differ ONLY in how much of the p-adic
valuation the local factor sees:
- Goldbach's singular series reads the RADICAL (v_p ≥ 1) — the comet's wings, blind to the prime-power
  tower (g(54), n=27=3³, carries the same 3-factor as n=3).
- Pisano (π(pᵏ)=pᵏ⁻¹π(p)) and the LRC shell (ramification depth k) read the FULL tower.

And the magma is the operation underneath all three: (ℤ, +), iterated. Goldbach asks its 2-fold
representation count over primes; Pisano its period mod m; LRC its covering of the circle. The free
magma (binary trees, Catalan) is the recursion substrate — the divisor-lattice peeling, the CRT
decomposition, the deletion–contraction of the partition function are all free-magma trees, and
multiplicativity is the quotient to the commutative Euler product. So "see everything recursively"
and "the Pisano/Goldbach/magma substance" are the same instruction: addition, localized at each
prime, with a valuation-depth dial — and the LRC shell turns that dial up to the full tower, which is
why n=14 is hard and n=19 is not.

## Honest scope

The descent and Pisano correspondence (THM-491 Part 1) and the freed-clock lemma (Part 2) are proved;
that LRC *difficulty* descends the tower is HYP-2436 (a route to LRC(14) ⟸ LRC(5)+LRC(7), not a proof
— both LRC(14) and LRC(19) remain open, proven only ≤ 7). The Pisano/Goldbach/magma synthesis
(HYP-2437) is a frame with verified sub-pieces, not a theorem. The most concrete next step is
HYP-2436 test 3: does C′(14) reduce exactly to C′(5)-on-the-3-core ∪ the THM-421 fiber?
