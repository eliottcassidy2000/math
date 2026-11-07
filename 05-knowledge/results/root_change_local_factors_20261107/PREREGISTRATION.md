# Pre-registration: local factors of the two-root-change multiplicity, and the census to 2^32

kind-pasteur-2026-11-07-S15. Reserves **THM-4916** and **HYP-9506**.

Owner prompt (2026-11-07, second of the day): compare small-number coincidences whose structure we may be able to
explain (Chan-Pak, the Fermat primes, the square-sum conjecture); derive the 3-adic local factor of the
two-root-change multiplicity; creatively construct a proven cutoff beyond which two root changes always suffice.

**Concurrent same-prompt work.** Session223 (first checkpoint 8865bb666b, 13:20) has:
- the exact two-run product rectangle and its collision defect;
- the square-sum graph on 2^a + 1, which has only the edge 3-33;
- full-state cancellation and square-reflection closure.

It lists the two-root-change cutoff and a "height-weighted 3-adic factor" as under investigation. Credited; the method
below is different.

Method credited as KNOWN:
- opus-S82 predicted a mod-3 local factor for its order-<=7 bank by seed confinement
  (`04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py`);
- this uses Session193 modular_sync Theorems 3-4 and Session206's weighted capacities.

## Theory computed before any registered test (exploratory)

**Confinement lemma (exact, trivial).** A seed s = (h,R,e,d) gives H_s(m) = h - e (m mod 2) (mod gcd(R,d)). For a
root-change seed, gcd(R,d) = gcd(D,U) of the pre-J state.

A small-scale exact decomposition (`decomp.c`, X = 2^18, window [2^17, 2^18)) confirms the confined odd-m mass lands
only on N = 0 (mod 3^j). It also shows that unconfined seeds are not class-uniform at that scale (max/min 1.87 over
the classes mod 729). Hence the method below.

**Residue automaton (`resdp.c`).** The moves are integer-linear, so (D,E,U,d) mod q is a finite automaton. Count every
prefix with exactly n letters O/Z and at most 2 letters J, weight 2 if it has at most 1 J and 1 otherwise, by the
class of F = D - E + U mod q. The local factor is L_n(r) = q * share of class r. No census data are used.

The asymptotic values at n = 1000, for N = 0 (mod p), with nonzero classes equal to within 0.001:

| p | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L(0) | 1.263 | 1.2125 | 1.1687 | 1.1219 | 1.1023 | 1.0877 | 1.0754 | 1.0644 | 1.0511 | 1.0479 |

For p = 3 the value at n = 2000 is 1.2612, and the nonzero classes are 0.8692 and 0.8696. Conjecture: the limits are
34/27 and 47/54. Asymptotically L(0) - 1 is about 1.5/p. At n near 20 the Fermat-type primes 3, 5 and 17 are
strongly boosted, a finite-size effect that decays like 1/n.

The bands over n = 18..23 (`dp_bands.txt`) are:
- the ratio rho_p = mean(N = 0 mod p) / mean(N != 0 mod p);
- for q = 27, the factors by v3(N).

## Registered tests (window: odd N in [2^23, 2^24), census with X = 2^24, exact 16-bit counters)

- **R1 (census to 2^32).** F_<=2 has no odd miss in (10^9, 2^32]: the misses below 2^32 are exactly the 140 below
  3,780,967. It will be run with the auditor's independent segmented enumerator, `audit/census_audit.c` of the S14
  packet.
- **R2 (fresh primes).** The census ratio rho_p lies in the band:

  | p | 5 | 7 | 11 | 13 | 17 | 23 | 31 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | band | [1.4375, 1.5674] | [1.1093, 1.2409] | [1.2604, 1.3569] | [1.1001, 1.1707] | [1.3379, 1.4857] | [1.0650, 1.1486] | [1.0514, 1.1227] |

  Each band is the DP band for n = 18..23, widened by 3% on each side. All seven must hold.
- **R3 (Fermat boost ordering).** rho_5 > rho_17 > rho_11 > max(rho_7, rho_13, rho_23, rho_31).
- **R4 (nonzero classes).** For each p in R2, the means over the p - 1 nonzero classes have max/min <= 1.10.
- **R5 (classes mod 27).** For at least 25 of the 27 classes r mod 27, the census mean in class r divided by the
  overall mean lies in the DP per-class band (n = 18..23) widened by 5%. The v3 aggregates of the census at this
  window are already known from S14 (205.8, 332.0, 457.3, ...). So R5 is blind only at the level of individual
  classes.

Any failure is reported as a failure.

## Not registered

- The proven-cutoff attack: bilinear product families from Session199's K-form, and barrier lemmas.
- The small-number coincidence atlas.
- The exact limits (34/27, 47/54).
