# Small-number coincidences: Chan-Pak, the Fermat primes, the square-sum problem, and the root-change exceptions

kind-pasteur-2026-11-07-S15. This note answers the owner's request to compare small-number coincidences whose
structure we may be able to explain.

Facts are cited to their sources in the repository; the survey was assembled by a read-only sweep. The repository's
standing rule (MISTAKE-559) applies: a numerical coincidence is numerology unless a map carries the mechanism. So every
row below states its mechanism.

## 1. The exceptional sets side by side

| Phenomenon | Exceptions (largest) | "Only these" status | Mechanism of the exceptions | How witnesses grow |
| --- | --- | --- | --- | --- |
| Chan-Pak: odd Hamiltonian-path counts | {7, 21} (21) | Exclusion PROVED for every order (THM-1370). Completeness verified through 2^32 - 1 (S78, THM-4875 (D); recertified by Session203) and 2^35 (S82/S83, one run); OPEN beyond | **Structural.** H multiplies over strong components; Moon's minimum over strong tournaments is non-decreasing (3, 5, 9, 15, 25, 45 for n = 3..8). 7 falls in the gap (5, 9). 21 needs a strong component with H = 21, which would need n <= 6 because f(7) = 25; THM-1370's exhaustive spectra (file THM-1370-h-spectrum-omits-7-21-all-n.md) show none exists, and 21 = 3 * 7 would need a strong 7 | m(h) >= 6 for 63 < h < 625 and >= 73 on [375, 623] (THM-4884); the counts at fixed H grow polynomially in n (Session219) |
| Strong C = 1 (uniquely Hamiltonian) counts | 29 values (2647) | FINITE-EXACT below 10^9 (THM-4904 / HYP-9494) | Statistical within an infinite grammar: the root-change depth is at most 7 | mean M_2 about ln^2 N |
| Two root changes, F_<=2 | 140 values (3,780,967) | FINITE-EXACT below 10^9; R1 extends to 2^32 (section 4 of PROOF) | Statistical, with a divisor-type local factor (PROOF section 2): large misses are 1 (mod 3) (14 of 15 above 10^5), never 0 (mod 3) | mean M_2 about ln^2 N; heavy lower tail (min 8 at 2.8 * 10^8) |
| Fermat primes 2^(2^n) + 1 | {3, 5, 17, 257, 65537} (65537) | Composite for 5 <= n <= 32; OPEN, with barriers (THM-4847: no covering congruence, Pi^0_1) | **Statistical in reverse.** P(F_n prime) decays like 2^-n. Boklan-Conway (arXiv:1605.01371) expect under 10^-9 new primes; the repo's sieve model (HYP-9410) gives 6.6 * 10^-9 | none: the "witness count" is the primality probability, which shrinks |
| Square-sum paths on {1..n} | n in {2..14, 18..22, 24} (24) | PROVED (Gerbicz 2018: finite base 41..2032 plus a 49-fold blow-up) | Structural for small n: leaves, degree law, cut certificates (THM-4505). Pell forces Q_15's unique path | paths: 10, 12, 35, 52, 19, 20, 349, 392, ... for n = 25.. (A071983) |
| Square-sum cycles | n <= 31 (31) | PROVED (Gerbicz) | Q_31 has 349 paths and no cycle; min degree >= 2 only from 31 | cycles: 1, 1, 11, 57, 31, 20, ... for n = 32.. |
| Terminal-form square-sum cycles | holes up to 33, and 39 | PROVED (THM-4851 (c)) | finite certificate + threshold theorem | 50 cycles at 39, none terminal |
| Nine cubes | {23, 239} | PROVED (Dickson 1939; the attribution is the survey's, not the repository's) | the number of representations by 8 cubes grows, so only small numbers fail | DP to 400,000 (HYP-2490) |
| Untouchable primes | {2, 5} | conditional on distinct-prime Goldbach (THM-4805) | s(pq) = p + q + 1 | Goldbach representations grow |
| Three-core grammar | {47, 63} (63) | holes PROVED permanent in the grammar; coverage of N >= 65 OPEN (Session216) | positive-cover padding bound | covers [23, 10^6] otherwise |
| Primes in square-sum order (HYP-9451) | paths first exist at k = 50,396 | the minimum is PROVED | arrangements need Chebyshev's mod-4 bias in {0, 1, 2} and six mod-36 conditions, which below 5 * 10^6 hold on only two windows; whether infinitely many bounds admit one is OPEN | a cautionary case |

## 2. A common framework, and where it stops

Each statistical row has a witness count r(N) whose mean mu(N) grows (or, for Fermat, a success probability that
shrinks).

**Heuristic.** The expected number of exceptions beyond X is sum over N > X of P(r(N) = 0). With a Poisson lower tail
this is about sum exp(-L(N) mu(N)), where L(N) is the local factor.

The last exception then sits where L(N) mu(N) is about ln N. In practice the lower tail is heavier than Poisson, so the
last exception comes later than this predicts.

The two-root-change data show both effects:
- The local factor is divisor-type. In the letter-count limit it is exactly L_p(0) = 1 + (3p-2)(p-1)/(2p^3) at each
  prime p (THM-4916), about 1 + 3/(2p). At census scale there is a strong finite-size boost at the Fermat-type
  primes 3, 5, 17 (and at 11, which divides 33).
  14 of the 15 misses above 10^5 are 1 (mod 3), none is 0, and 6 have a prime factor below 20.
- Even so, mu(N) is about 227 at the last miss, 3,780,967 (exact counts on [3.5e6, 4e6); the 180.8 in `rc_1e9.txt` is an
  8-bit saturated lower bound), while Poisson with L = 0.8 would make a miss there impossible (e^-182). The
  lower tail is structural: the near-misses come from correlated families, not independent trials.

Hence the census cannot be extrapolated with Fermat-style confidence. The repository has no calibrated expected count
of misses beyond 2^32 for F_<=2. The analogous computation for Fermat primes is credible only because primality
events are nearly independent once the sieve is accounted for.

**Proved rows.** They close the statistical gap in one of three ways:
1. **Structural obstruction plus finite search:** Chan-Pak's 7 and 21 (Moon's minimum gap), and the small square-sum
   cases.
2. **Constructive self-similar induction:** Gerbicz's square-sum blow-up n -> 49n + res, from a finite base.
3. **Analysis with explicit constants plus finite computation:** nine cubes, and in the same spirit Helfgott's ternary
   Goldbach.

The two-root-change cutoff currently has neither 2 nor 3 (PROOF section 3). Its family is bilinear in two states
(THM-4916 (1d)). That points to route 3, a second-moment bound for K(v,w) over independent states, or to a state-cone
version of route 2.

## 3. Fermat primes inside the root-change grammar

Fermat numbers enter the grammar in three proved ways:
- H(t . Z^r) = 2^(r+1) + 1, so every Fermat number is a root-tail value (Session200);
- the multiplier identity (THM-4916 (1c)): value(p J Z O^(j-1) Z) = (2^(j-1) + 1) G(p);
- the local-factor boost at the multiplier primes (PROOF section 2).

The Fermat primes 3, 5, 17, 257, 65537 are exactly the primes that are themselves multipliers 2^i + 1. Every prime
with even ord_p(2) divides some multiplier (11 | 33, 13 | 65, 19 | 513).

Fermat-type numbers enter the census through two exact families:
- the (1c) multiplier words;
- the words that begin with the root rectangle Z^a O^b J, all of whose values are divisible by 2^a + 1, because
  t.Z^a O^b J = (h, g, h, g) with g = 2^a + 1 dividing h.

At [2^23, 2^24) these carry about 23% and 21% of the zero-class excess at p = 5, and 22% and 19% at p = 17. Seed
confinement at the last root change carries the rest, about half. The R3 ordering persists with both families
removed: rho = 1.363, 1.339, 1.324, 1.247 for p = 5, 17, 11, 7. So R3 confirms the ordering but not this mechanism
(audit M3). In the letter-count limit the whole excess is seed confinement (THM-4916, closed form).

Session219 records further Fermat structure in the grammar:
- **five-run escape:** every O/Z recipe with at most four runs obeys 6700417 | H => 641 | H, from F_5 = 641 * 6700417;
- **Fermat-hole totient law:** #{x : phi(x) = 2^a} = a + 2 for 0 <= a <= 31, and 32 for 32 <= a <= 63 (Session219;
  beyond 63 it depends on further Fermat numbers).

## 4. Comparison with the square-sum problem

- The square-sum graph on 2^a + 1 has only the edge 3-33, and prime factors carry no internal edge (Session223,
  first checkpoint).
- Session211 proves that no map from square-sum graphs to tournaments can have H = 2p(Q_n) + 1: p(Q_23) = 3 and
  p(Q_25) = 10 would force H = 7 and 21. The two exceptional sets do not share a law.

The structural lesson is still useful. The square-sum problem was settled by a self-similar construction from a finite
base. That is exactly the missing piece for the two-root-change cutoff (PROOF section 3.3-3.4).
