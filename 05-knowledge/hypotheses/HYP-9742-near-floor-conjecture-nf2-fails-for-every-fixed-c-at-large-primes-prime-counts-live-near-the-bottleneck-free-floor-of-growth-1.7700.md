---
id: HYP-9742
title: "The near-floor conjecture NF2(c) fails for every fixed c at all sufficiently large primes: a prime count of a strong n-tournament is at least about rho*^n with rho* = ((9 + sqrt 113)/2)^(1/4) = 1.7700005 (the bottleneck-free growth of THM-5202 (D)), exponentially above the Moon-Busch floor f(n) ~ 5^(n/3); extrapolated failure near 1.5e8 for c = 3 and 1.5e10 for c = 4. Corrected form: primes live near the bottleneck-free floor"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-bottlenecks, wave 28 of chanpak_waves_20261113.
  PROVED ingredient: THM-5202 (C), bulk floor runs force exact factorizations H = L·R in bounded-span ribbons.
  FINITE-EXACT ingredient: THM-5202 (D), every periodic block cheaper than 1.7700 is a bottleneck (periods ≤ 10 in
  span 2, ≤ 6 in span 3, plus 1,926 floor-plus-spoiler words).
  NUMERICAL: the smallest prime count over f(k) rises 1.336 → 1.713 from k = 10 to k = 18, i.e. 3.2% per vertex,
  against the predicted 3.51%.
  Concurrent and independent: the 7-21 session (notes/06 E2–E3) conjectures the same failure from p-adic statistics.
session: mac-mini-2026-11-15-bottlenecks
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 28)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w28_nearfloor_ribbon_dfs.c, w28_bottleneck_free_premium.py,
    w28_bottleneck_free_charpoly.py
related:
  - THM-5202; 7-21 repository notes/05 D9 (NF2(c)) and notes/06 (two regimes); THM-4744 (Busch floor)
---

# HYP-9742 - why near-floor witnesses run out at primes

**Naming.** "Near-floor conjecture" has two meanings in the current packets.
- Here it is the *coverage* statement NF2 of the 7-21 repository.
- The *lower bound* H(T) ≥ f(n) + 3(hc(T) − 1) is a different statement. It is proved in
  chanpak_near_floor_crossing_20261115 (NEAR_FLOOR.md, GLOBAL.md), which also computes the exact span-2 ribbon
  minimum m(16 + 2k) = 3375·3^k. That agrees with this session's span-2 census (minimum 3375 at order 16, 10125 at
  order 18).

**NF2(c)** (7-21 notes/05 D9) says that every odd N ≥ 25 is the count of a strong n-tournament with H < c·f(n).
It holds for c = 4 below 12,779 (7-21 notes/06 E2).

**Conjecture.**
1. Let π(n) be the smallest prime H(T) over strong n-tournaments. Then (1/n) ln π(n) → ln ρ*, with
   ρ* = ((9 + √113)/2)^(1/4).
2. Hence π(n)/f(n) grows like 1.0351^n, and for every c there is p_c beyond which no prime satisfies NF2(c).

**Mechanism (THM-5202).**
- A word within a bounded ratio of the floor consists mostly of floor units.
- A run of 3 (5-units) or 4 (3-units) consecutive units is a rank-one bottleneck, so the count factors exactly as
  L(prefix)·R(suffix).
- A prime count must therefore break every run before it reaches bottleneck length.
- The cheapest way to do that has growth ρ*, the dominant root of x² - 9x - 8 per four vertices.
- Composites survive near the floor as products. Their misses are the numbers whose factorizations do not split into
  available end values: every window miss at orders 12–14 is prime or has a prime factor ≥ 13.

**Extrapolated thresholds (HEURISTIC).** Take π(k)/f(k) ≈ 1.713 · 1.0351^(k-18).
- The ratio passes 3 at k ≈ 34 and 4 at k ≈ 42.
- A prime p > c·f(n_c) must use an order above n_c, where no prime count has ratio below c.
- So NF2(3) should fail beyond about 3·f(34) ≈ 1.5·10^8, and NF2(4) beyond about 4·f(42) ≈ 1.5·10^10.

These numbers depend on extrapolating eight orders of data and should be read as orders of magnitude.

**What survives (corrected statement, CONJECTURE).** Every prime p ≥ 11 is the count of a strong n-tournament with
H < C·ρ*^n, for an absolute C. This is the near-floor idea with the floor replaced by the bottleneck-free floor.
Composite N then follow by multiplicativity (ordinal sums), except for the multiples of 7 that need direct witnesses.

**Caveats.**
- The lower-growth claim is proved only for periodic words of the stated lengths, not for all aperiodic words.
- Non-ribbon tournaments are not covered. Near-floor classes have span ≤ 4 (7-21 census), and through order 10 the
  smallest prime counts of all strong tournaments equal the span-2 ones.
- This does not touch Chan–Pak itself. Primes are realized above the floor by additive families (recursive trees,
  carrier forms).

**Predictions.**
- *P100.* Using `w28_nearfloor_ribbon_dfs.c` with span 2, the smallest prime count at order 19 exceeds 1.70·f(19), and
  at order 20 it exceeds 1.75·f(20).
- *P101.* No bulk-strong span-3 periodic block of period 7 or 8 is bottleneck-free with growth below 1.7700. The test
  needs a C port of `w28_bottleneck_free_premium.py` to cover all 8^7 + 8^8 blocks.
