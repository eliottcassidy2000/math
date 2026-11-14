---
id: HYP-9682
title: "Invertible root-change conjecture: (A) a one-apex source followed by invertible blocks only (J and J O v Z, any length; no singular block, so no reset) attains every odd N > 6223; (B) without a source, the orbit of the root state (1,1,1,1) under the invertible blocks attains every odd N > 179,367; (C) its multiplicities obey an L^2 local limit law with the Haar singular series of THM-5162, so density one holds; (D) a finite block set (|w| <= 5) already suffices beyond 12,811"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-14-density, wave 24 of chanpak_waves_20261113.
  FINITE-EXACT: (A) exactly 43 odd numbers below 2^26 are missed, the largest 6223 (scripts/invorb.c, mode 1);
  (B) exactly 288 below 2^26, the largest 179,367 (mode 2);
  (D) |w| <= 6 gives the same 43 misses below 2^20; |w| <= 5 gives exactly 48 below 2^26 (the largest 12,811);
  |w| <= 4 gives 71 below 2^20 (the largest 968,947).
  NUMERICAL, pre-registered (P51-P58, P62-P65 of wave 24, all CONFIRMED): the local law of THM-5162 (D)-(E) in the census;
  multiplicity exponent 0.2777 for all N and 0.2771 for primes; the L^2 discrepancy V(X) of (C) falls like
  8/mean(X).
  PROVED inputs: THM-5162 (strong approximation at every odd modulus, S = {2}; Haar law of long random words;
  singular series (7/8) zeta(3) prod_{p | N} (1 - p^-2)). NOT PROVED: (C) in any form, uniformity of mixing in the modulus, the archimedean local limit, minor arcs.
  Relation to earlier conjectures: (A) is stronger than HYP-9494 (D) (the full O/Z/J grammar misses exactly 29 odd
  numbers) because the language is a sub-language; it is independent of HYP-9677 (bounded runs, singular blocks
  allowed). (A) together with the census below 6223 implies Chan-Pak.
session: mac-mini-2026-11-14-density
date: 2026-11-14
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-24-the-invertible-root-change-orbit--reset-free-haar-distributed-and-complete-beyond-6223
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/invorb.c (census; modes 1, 2; block cap via CAP)
  - 05-knowledge/results/chanpak_waves_20261113/scripts/srcdiv.c (source p-divisible share, P58)
  - 05-knowledge/results/chanpak_waves_20261113/scripts/l2test.py (L^2 test, P55)
related:
  - THM-5162 (local half, proved), THM-5120 (resets need singular blocks), HYP-9494, HYP-9677, HYP-9660 (B)
---

# HYP-9682 - the invertible root-change conjecture

**Language.** Letters O, Z, J act on (H, E, L, h) as in THM-5162, starting from (1, 1, 1, 1). The *invertible
language* consists of any O/Z word (the one-apex source), followed by any sequence of invertible blocks: J alone, or
J O v Z with v in {O, Z}* of any length.

In letters, once the first J has occurred:
- no Z directly after a J;
- no J directly after an O.

No singular block ever occurs. So by THM-5120 no rank-one reset can occur, and by THM-5162 (E) no state can become
p-divisible after the source.

**(A) Free source.** Every odd N > 6223 is a count of the invertible language.
- FINITE-EXACT to 2^26: the misses are exactly 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251,
  277, 281, 299, 313, 335, 421, 449, 503, 533, 551, 569, 583, 637, 667, 821, 1067, 1339, 1387, 1661, 1811, 1861,
  1877, 1957, 2647, 3523, 3551, 5231, 6223.
- The full grammar's 29 misses are among them.
- Wave 20 reported that invertible blocks alone miss 34% of the odd numbers below 2^21. That census capped the
  source at four letters; with a free source the picture changes completely.

**(B) No source.** Every odd N > 179,367 is a count of a word that starts like a block (J, or O v Z then J, ...).
This is a pure semigroup orbit of one vector in GL_3(Z[1/2]). FINITE-EXACT: exactly 288 misses below 2^26, the
largest 179,367.

**(C) L^2 local limit law.** Let:
- r(N) be the number of words with count N;
- mean(X) the average of r over odd N in [X, 2X);
- S(N) the class factor of THM-5162 (E), the Haar law nu_p corrected by the source share pi_p.

Conjecture: (1/#) sum (r(N)/(mean(X) S(N)) - 1)^2 -> 0. By Chebyshev's inequality this gives density one.

Data on [2^23, 2^24) (and [2^25, 2^26) where stated):
- *Free source:* mean multiplicity 369, which grows like X^0.2777, and primes at X^0.2771. V = 0.023 = 8.4/mean (0.0166 = 8.96/mean on [2^25, 2^26), minimum multiplicity 196);
  the share of N with r < model/2 is 1.1·10^-4.
- *No source:* the class factors at p = 3, 5, 7, 11 are 1.040, 1.008, 1.003, 1.001, against the Haar values
  p^3/(p^3 - 1) = 1.0385, 1.0081, 1.0029, 1.0008. The prime mean is 1.0526 times the mean, against
  (7/8) zeta(3) = 1.0518.
- *Prime powers (P65):* the density depends only on which primes divide N. N ≡ 0 vs 3, 6 (mod 9) gives 0.998;
  N ≡ 0 vs the other multiples of 5 (mod 25) gives 1.000; the ratio (3 | N)/(3 ∤ N) is 0.8845, against 8/9.

**(D) Finitely many blocks.** Blocks with |w| <= 5 (J and 15 others) attain every odd N > 12,811. FINITE-EXACT to
2^20. With |w| <= 6 the census below 2^20 already coincides with the unbounded one.

**Why it matters.** This is the cleanest form of Chan-Pak found so far. The orbit lives in a semigroup of
GL_3(Z[1/2]), Zariski dense (wave 23), and THM-5162 settles its local behaviour completely: strong approximation at
every odd prime and the Haar law. Primes are not penalized; their singular series exceeds 1. So a counterexample
above 6223 would have to be a purely archimedean, sporadic failure of a supercritical orbit with critical exponent
about 1.28. That is the regime of Apollonian packings (critical exponent 1.3056, density one by Bourgain–Kontorovich).
