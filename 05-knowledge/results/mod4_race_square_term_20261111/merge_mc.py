#!/usr/bin/env python3
"""Merge the one-line outputs of biasmc.c (independent seeds) and print the estimates with standard errors:
D_(m+1) = E[eps_R b_R], rho_m = 2^m E[b_R^2], E|b_R|, the rank law of beta_R and the non-zero-bias frequency by rank.
Standard errors: for D from the pooled second moment (Var(eps b) = E[b^2] - D^2); for fractions, binomial.
Usage: python3 merge_mc.py file [file ...]"""
import sys, re, math
N = 0; S1 = S2 = SA = SU = 0.0; fails = [0, 0, 0, 0]; r = {}; m = None; per = []
for path in sys.argv[1:]:
    for line in open(path):
        g = re.match(r'm (\d+) N (\d+) failures\(mu odd, a odd, C3, edges\) (\d+) (\d+) (\d+) (\d+) S1 (\S+) S2 (\S+) SA (\S+) SU (\S+) npos (-?\d+)(.*)', line)
        if not g: continue
        m = int(g.group(1)); n_ = int(g.group(2)); N += n_
        for k in range(4): fails[k] += int(g.group(3 + k))
        s1 = float(g.group(7)); S1 += s1; S2 += float(g.group(8)); SA += float(g.group(9)); SU += float(g.group(10))
        per.append(s1 / n_)
        for h in re.finditer(r' r(\d+) (\d+) (\d+)', g.group(12)):
            k = int(h.group(1)); a, b = int(h.group(2)), int(h.group(3))
            if a or b:
                r.setdefault(k, [0, 0]); r[k][0] += a; r[k][1] += b
D = S1 / N; E2 = S2 / N; se = math.sqrt(max(E2 - D * D, 0) / N)
n = m + 1
print("m = %d, N = %d samples in %d runs; control failures (mu odd, a odd, Theorem C (3), even number of edges) = %s" % (m, N, len(per), fails))
print("D_%d = %.6e +- %.2e (standard error);  D_%d * 4^%d = %.1f +- %.1f;  z = %.1f" % (n, D, se, n, n, D * 4 ** n, se * 4 ** n, D / se))
print("per-run estimates of D_%d * 4^%d: %s" % (n, n, ", ".join("%.0f" % (x * 4 ** n) for x in per)))
print("rho_%d = 2^m E[b^2] = %.5f ;  E|b_R| = %.6f ;  E[b_R] = %.6e (forced part 2^(1-m) = %.6e)" % (m, E2 * 2 ** m, SA / N, SU / N, 2.0 ** (1 - m)))
for k in sorted(r):
    a, b = r[k]; fr = (a + b) / N
    print("  rank %2d (corank %d): fraction %.6f +- %.6f ; non-zero-bias frequency %.5f +- %.5f (uniform-linear-part model %.5f)" % (
        k, m - k, fr, math.sqrt(fr * (1 - fr) / N), b / (a + b), math.sqrt((b / (a + b)) * (1 - b / (a + b)) / (a + b)), 2.0 ** (k - m)))
