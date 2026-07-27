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

Same-session extension (04-computation/ham_excess_ext_kps_S134.cpp,
output ham_excess_ext_kps_S134.out):

```text
R_19  H=1184212824763       exc 2.551967
R_21  H=125547534942879     exc 2.576702
R_23  H=16011537490557279   exc 2.597757
QR_23 H=15760206976379349   exc 2.556980   (tc = 62,293,308,207,033)
```

Consecutive rotational differences `0.0297, 0.0247, 0.0211` shrink
geometrically (~0.85 ratio), extrapolating the limit to `~2.72`.
The `e(1 - alpha/n)` fit gives `alpha = 1.226, 1.163, 1.094, 1.020`
at `n = 17, 19, 21, 23` -- drifting toward 1.

**Conjecture.** For both families (and plausibly for every
quasirandom regular tournament family), `exc -> e`; sharpened by
the extension data to `exc(n) = e (1 - (1 + o(1))/n)` -- the fitted
`alpha` passes through `~1.02` at `n = 23` and is still falling, so
`alpha = 1` exactly is the natural target (a second-order term of
the permanent-style regularity correction).

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

**Universality upgrade (same session, hostile control PASSED
toward the conjecture).** The maximally structured regular family
`C3[R_k, R_k, R_k]` (tower-composed, the opposite of quasirandom)
also climbs to `e`:

```text
n    9        15        21
exc  2.228571 2.420817  2.516651     (rotational at same n:
                                      2.3048, 2.4866, 2.5767)
```

with the family gap shrinking `0.076 -> 0.066 -> 0.060`. So `e` is
NOT a quasirandomness constant; the sharpened conjecture is
**universality over regular tournaments**: every sequence of regular
tournaments has `exc -> e`, with family-dependent `O(1/n)`
coefficient. (The `n = 9` composed value `H = 3159` independently
matches the S132 hand computation via the THM-2453 run-transfer
expansion -- which is also the proof route: the run-transfer formula
should PROPAGATE the `e`-limit through `C3`-composition, making the
tower grammar an induction scheme for the conjecture.)

Cheapest decisive tests: (i) DONE (R_19..23, alpha -> 1);
(ii) compute the first character-harmonic correction analytically
at p -> infinity and compare with `e - exc(p)`; (iii) OEIS: the
H(R_n) sequence's novelty is UNVERIFIED (lookups rate-limited);
(iv) prove the C3-composition transfer: if exc(A_k) -> e then
exc(C3[A_k,A_k,A_k]) -> e via the run-transfer expansion --
the first candidate PROOF route into the conjecture;
(v) sample non-circulant regular tournaments at n = 13..17 for the
spread's shrink rate.
