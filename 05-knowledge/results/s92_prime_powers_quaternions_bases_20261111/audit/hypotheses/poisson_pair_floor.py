#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (B). Benchmark for the L1 remainder mean|r - 2 kappa S ln N| / mean r measured in
# series_counts*.out.  r counts ORDERED triples, so off the diagonal r = 2u with u the number of unordered pairs;
# if u were Poisson with mean lambda/2 the relative mean absolute deviation of r would be
#   E|2u - lambda| / lambda = 2 E|u - lambda/2| / lambda       ("Poisson-pair floor").
# Listed against the measured values (sheet 0 / sheet 1) on the windows [2^(E-1), 2^E).
import math


def mad_poisson(mu):
    k = int(math.floor(mu))
    return 2 * math.exp(-mu + (k + 1) * math.log(mu) - math.lgamma(k + 1))


rows = [
    (17, 4.1278, 5.1261, 0.6051, 0.5700, 0.6359, 0.6303),
    (19, 4.5752, 5.7139, 0.5705, 0.5389, 0.6134, 0.5820),
    (21, 5.0222, 6.3034, 0.5396, 0.5079, 0.5894, 0.5467),
    (23, 5.4700, 6.8925, 0.5152, 0.4803, 0.5625, 0.5271),
    (25, 5.9180, 7.4817, 0.4967, 0.4556, 0.5317, 0.5046),
]
print("# poisson_pair_floor.py (auditor hypotheses). Inputs copied from series_counts_trend.out and series_counts_25.out")
print("window            mean r (s0/s1)    measured L1/mean with S (s0/s1)   with a constant predictor (s0/s1)   Poisson-pair floor (s0/s1)")
for E, m0, m1, a0, a1, c0, c1 in rows:
    f0 = 2 * mad_poisson(m0 / 2) / m0
    f1 = 2 * mad_poisson(m1 / 2) / m1
    print("[2^%d, 2^%d)   %.3f / %.3f      %.4f / %.4f                   %.4f / %.4f                     %.4f / %.4f" % (
        E - 1, E, m0, m1, a0, a1, c0, c1, f0, f1))
