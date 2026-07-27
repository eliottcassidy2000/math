---
id: HYP-9028
title: "The circulant Hamiltonian excess tends to e"
status: >
  OPEN (exact data n <= 17 rotational, p <= 19 Paley; both excess
  sequences increase toward ~2.52 with O(1/n)-shaped corrections
  consistent with limit e). Includes the proved free-rotation lemma
  n | H and the honest kill of the 1729/Carmichael numerology.
source: kind-pasteur-2026-07-26-S134
related:
  - THM-2454-pure-blue-alphabet-and-center-law
  - THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum
script: 04-computation/ham_excess_spectra_kps_S134.py
output: 05-knowledge/results/ham_excess_spectra_kps_S134.out
---

# HYP-9028 -- regular circulant tournaments carry e times the mean Hamiltonian mass

Define the **Hamiltonian excess** of an n-tournament as
`exc(T) = H(T) / (n!/2^{n-1})` (denominator = the mean over uniform
random tournaments). Exact data:

```text
rotational R_n = circulant{1..(n-1)/2}:
n    3     5     7      9      11      13       15         17
H    3     15    175    3267   93027   3711175  198464295  13689269499
exc  2.000 2.000 2.222  2.305  2.386   2.441    2.487      2.522

Paley QR_p:
p    7      11      19
H    189    95095   1172695746915
exc  2.400  2.4395  2.5272
```

**Conjecture.** For both families (and plausibly for every
quasirandom regular tournament family), `exc -> e`, with
`exc(n) = e (1 - alpha/n + o(1/n))`, fitted `alpha ~ 1.2 - 1.3`.

Mechanism sketches:

1. (Paley, harmonics form) `H = 2^{-(p-1)} sum_sigma prod_i
   (1 + chi(sigma_{i+1} - sigma_i))`; the empty-set term of the
   product expansion is exactly the random mean, and every
   correction is a Weil-type character sum over difference patterns
   -- the excess is a convergent series of prime harmonics, the
   third appearance of the THM-2447 grammar this week.
2. (regularity form) `e` is the standard second-order factor for
   counting long paths under exact degree conditioning (permanent-
   style corrections); a proof should go through Cuckler-style
   counting for regular tournaments.

Proved along the way: the rotation group `Z/n` acts freely on
directed Hamiltonian paths, so `n | H(R_n)` and `p | H(QR_p)`
(verified on all rows).

Numerology control (recorded to stop a false lead): `tc(QR_11) =
H/|Aut| = 95095/55 = 1729 = 7 * 13 * 19` is Chernick's first
universal-form Carmichael number and contains the LRC stalk modulus
`91 = 7 * 13`; but `tc(QR_19) = 6,857,869,865` FAILS Korselt, so
the Carmichael property is a p = 11 coincidence, not a law.

Cheapest decisive tests: (i) compute `exc` for R_19, R_21 (bit-DP,
C++) and check the `e(1 - alpha/n)` fit tightens; (ii) compute the
first character-harmonic correction analytically at p -> infinity
and compare with `e - exc(p)`; (iii) OEIS: the H(R_n) sequence's
novelty is UNVERIFIED (lookups rate-limited); check before any
submission claim.
