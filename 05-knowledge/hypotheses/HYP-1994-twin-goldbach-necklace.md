---
id: HYP-1994
status: OPEN
source: oracle-2026-06-01-S521
related:
  - HYP-1965
  - HYP-1966
  - HYP-1987
---

# HYP-1994: twin-prime Goldbach reduces to binary Goldbach on the twin-center necklace; last exception is 4208

**Setup.** Twin prime = prime in a twin pair; T={3,5,7,11,13,17,19,29,31,...}.
Even numbers NOT of the form t1+t2 (t1,t2 in T) = OEIS A007534.

**Verified (exhaustive to 2,000,000):** exactly **35** such even numbers, largest
**4208**. Apart from {2,4} they are **11 triples** of consecutive evens
{6m-2,6m,6m+2}.

**Structural claim (VERIFIED exactly on even n in (8,6000], 0 mismatches).**
Twin pairs are centered at multiples of 6 (pair {6k-1,6k+1} has center 6k; only
{3,5}, center 4, is anomalous). Let C = twin-pair centers, K = C/6 =
{k : 6k±1 both prime} = {1,2,3,5,7,10,12,17,18,...} (the "twin-center necklace").
A sum of two twin primes lands in {c1+c2-2, c1+c2, c1+c2+2}; since C⊂6Z the
covering sum is a multiple of 6, so each mod-6 triple {6m-2,6m,6m+2} is
representable as a UNIT iff 6m∈C+C iff m∈K+K. Hence:

  **even n>8 is a sum of two twin primes  ⇔  round(n/6) ∈ K+K.**

The 11 exception-triples ⇔ the 11 integers not in K+K:
m = 16, 67, 86, 131, 151, 186, 191, 211, 226, 541, **701** (×6 → 96,402,516,
786,906,1116,1146,1266,1356,3246,4206). The triples come in threes purely as the
±2 thickening of one center-sum 6m; the genuine content is the 11-hole comb of
K+K.

**Conjecture (OPEN):** m=701 is the LAST hole of K+K — i.e. every m>701 is a sum
of two twin-center necklace elements, so twin-Goldbach has no exception above
4208. K has density ~ 2·C2·n/log^2 n (sparser than primes); the holes are a
finite startup deficit before the representation count r_K(m) stays positive.

**To do:** (a) verify r_K(m)>0 for m well past 701 (push the bound to 10^7+);
(b) lower-bound r_K(m) via Hardy–Littlewood twin-prime density to make the
"no holes past 701" claim conditional-provable; (c) explain the specific 11
m-values (do they sit in local thinnings of K?).

**Evidence:** `04-computation/twin_goldbach_exceptions_s521.py`,
`04-computation/twin_goldbach_necklace_reduction_s521.py` (+ `.out` in
`05-knowledge/results/`). **Reflection:**
`07-reflections/twin-goldbach-necklace-triples-s521.md`.
**Lens link:** the center coordinate 6k is exactly the pair-lens coordinate of
HYP-1965/1966 (`pair-first-twin-prime-lens-s502`); twin-Goldbach is Goldbach
"on pairs."
