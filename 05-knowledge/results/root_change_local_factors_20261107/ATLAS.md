# Small-number coincidences: Chan-Pak, the Fermat primes, the square-sum problem, and the root-change exceptions

kind-pasteur-2026-11-07-S15. This note answers the owner's request to compare small-number coincidences whose
structure we may be able to explain.

Facts are cited to their sources in the repository; the survey was assembled by a read-only sweep. The repository's
standing rule (MISTAKE-559) applies: a numerical coincidence is numerology unless a map carries the mechanism. So every
row below states its mechanism.

## 1. The exceptional sets side by side

| Phenomenon | Exceptions (largest) | "Only these" status | Mechanism of the exceptions | How witnesses grow |
| --- | --- | --- | --- | --- |
| Chan-Pak: odd Hamiltonian-path counts | {7, 21} (21) | Exclusion PROVED for every order (THM-1370). Completeness verified through 2^32 - 1 (Session203) and 2^35 (S82/S83, one run); OPEN beyond | **Structural.** H multiplies over strong components; Moon's minimum over strong tournaments is non-decreasing (3, 5, 9, 15, 25, 45 for n = 3..8). 7 falls in the gap (5, 9); 21 = 3 * 7 or falls in (15, 25) | m(h) >= 6 for 63 < h < 625 and >= 73 on [375, 623] (THM-4884); the counts at fixed H grow polynomially in n (Session219) |
| Strong C = 1 (uniquely Hamiltonian) counts | 29 values (2647) | FINITE-EXACT below 10^9 (THM-4904 / HYP-9494) | Statistical within an infinite grammar: the root-change depth is at most 7 | M_2 mean about 0.5 ln^2 N |
| Two root changes, F_<=2 | 140 values (3,780,967) | FINITE-EXACT below 10^9; R1 extends to 2^32 (section 4 of PROOF) | Statistical, with a divisor-type local factor (PROOF section 2): misses are rough, 1 (mod 3), never 0 (mod 3) | M_2 mean about 0.5 ln^2 N; heavy lower tail (min 8 at 2.8 * 10^8) |
| Fermat primes 2^(2^n) + 1 | {3, 5, 17, 257, 65537} (65537) | Composite for 5 <= n <= 32; OPEN, with barriers (THM-4847: no covering congruence, Pi^0_1) | **Statistical in reverse.** P(F_n prime) decays like 2^-n. Boklan-Conway (arXiv:1605.01371) expect under 10^-9 new primes; the repo's sieve model (HYP-9410) gives 6.6 * 10^-9 | none: the "witness count" is the primality probability, which shrinks |
| Square-sum paths on {1..n} | n in {2..14, 18..22, 24} (24) | PROVED (Gerbicz 2018: finite base 41..2032 plus a 49-fold blow-up) | Structural for small n: leaves, degree law, cut certificates (THM-4505). Pell forces Q_15's unique path | paths: 10, 12, 35, 52, 19, 20, 349, 392, ... for n = 25.. (A071983) |
| Square-sum cycles | n <= 31 (31) | PROVED (Gerbicz) | Q_31 has 349 paths and no cycle; min degree >= 2 only from 31 | cycles: 1, 1, 11, 57, 31, 20, ... for n = 32.. |
| Terminal-form square-sum cycles | holes up to 33, and 39 | PROVED (THM-4851 (c)) | finite certificate + threshold theorem | 50 cycles at 39, none terminal |
| Nine cubes | {23, 239} | PROVED (Dickson 1939; the attribution is the survey's, not the repository's) | the number of representations by 8 cubes grows, so only small numbers fail | DP to 400,000 (HYP-2490) |
| Untouchable primes | {2, 5} | conditional on distinct-prime Goldbach (THM-4805) | s(pq) = p + q + 1 | Goldbach representations grow |
| Three-core grammar | {47, 63} (63) | holes PROVED permanent in the grammar; coverage of N >= 65 OPEN (Session216) | positive-cover padding bound | covers [23, 10^6] otherwise |
| Primes in square-sum order (HYP-9451) | paths first exist at k = 50,396 | the minimum is PROVED | Chebyshev's mod-4 bias must lie in {0, 1, 2} | a cautionary case: "eventually always" fails for a long time |

## 2. A common framework, and where it stops

Each statistical row has a witness count r(N) whose mean mu(N) grows (or, for Fermat, a success probability that
shrinks).

**Heuristic.** The expected number of exceptions beyond X is sum over N > X of P(r(N) = 0). With a Poisson lower tail
this is about sum exp(-L(N) mu(N)), where L(N) is the local factor.

The last exception then sits where L(N) mu(N) is about ln N. In practice the lower tail is heavier than Poisson, so the
last exception comes later than this predicts.

The two-root-change data show both effects:
- The local factor L(N) is divisor-type (PROOF section 2), about 1 + 1.5/p for each prime p dividing N
  asymptotically, with a strong finite-size boost at the Fermat-type primes 3, 5, 17 (and at 11, which divides 33).
  Every large miss is rough and 1 (mod 3).
- Even so, mu(N) is about 180 at the last miss, 3,780,967, while Poisson with L = 0.8 would make a miss there
  impossible (e^-144). The lower tail is structural: the near-misses come from correlated families, not independent
  trials.

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

So the Fermat primes 3, 5, 17, 257 are exactly the primes p for which a single cheap tail makes every state's balanced
readout a multiple of p. This is why 3, 5 and 17 dominate the local factors at census scale (R3 CONFIRMED). This is a
mechanism, not a coincidence: it is an exact identity plus its measured consequence.

Session219 records further Fermat structure in the grammar:
- **five-run escape:** every O/Z recipe with at most four runs obeys 6700417 | H => 641 | H, from F_5 = 641 * 6700417;
- **Fermat-hole totient law:** #{x : phi(x) = 2^a} = a + 2 for a <= 31, then exactly 32.

## 4. Comparison with the square-sum problem

- The square-sum graph on 2^a + 1 has only the edge 3-33, and prime factors carry no internal edge (Session223,
  first checkpoint).
- Session211's port bridge sends square-sum path counts to tournaments by H = 2p + 1. That would force H = 7 and 21,
  so no such map exists, and Session211 says the two exceptional sets do not share a law.

The structural lesson is still useful. The square-sum problem was settled by a self-similar construction from a finite
base. That is exactly the missing piece for the two-root-change cutoff (PROOF section 3.3-3.4).
