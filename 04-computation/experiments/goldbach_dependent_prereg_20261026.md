# Pre-registered predictions for the Goldbach-dependent census at N = 10^8 (opus-2026-10-26-S61)

Written and committed before `./gdep 100000000` was run. Baseline (`goldbach_dependent_odd_values_20261026.cpp` at
N = 10^7, matching S57's counts plus the untouchable 5): dependent share 0.09255 overall; per decade 0.1525 (10^4),
0.1067 (10^5), 0.09033 (10^6); the fit share(decade) = c / ln(5 x 10^d) gives c = 1.400 and 1.393 for the last two
decades. Within v - 1 = 2 x prime the dependent share is 0.14657; the pattern 2^1 x p q holds 0.2153 of the dependent
set; the enrichment of the pattern 2^1 x p (dependent share / overall share) is 0.1104 / 0.0697 = 1.58.

Predictions for N = 10^8:
1. Dependent share in the decade [10^7, 10^8): in [0.074, 0.083] (c in [1.31, 1.47] at ln(5 x 10^7) = 17.73).
2. Cumulative dependent count at 10^8 equals S57's 4,117,839 plus one (the untouchable 5): 4,117,840 (consistency, not a
   prediction).
3. Within v - 1 = 2 x prime, the cumulative dependent share at 10^8 lies in [0.125, 0.145].
4. The pattern 2^1 x p q holds a share of the dependent set in [0.200, 0.240].
5. The enrichment of 2^1 x p is in [1.45, 1.75].
