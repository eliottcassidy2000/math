#!/usr/bin/env python3
"""Run C: estimates of E(-1)^c3, E(-1)^c5, E(-1)^(c3+c5) from the outputs of c35mc.c (one line per order:
n N s3 s5 s35), with standard errors, and the verdicts of the pre-registered R7, R8 and the control.
Usage: python3 c35_tables.py out/c35_*.out"""
import sys, math
rows = {}
for fn in sys.argv[1:]:
    for line in open(fn):
        n, N, s3, s5, s35 = map(int, line.split()); rows[n] = (N, s3 / N, s5 / N, s35 / N)
print(" n      N        E[{3}] (exact)            E[{5}]            E[{3,5}]         E[{3,5}] 2^floor((n-1)/2)   s.e.")
for n in sorted(rows):
    N, e3, e5, e35 = rows[n]; se = 1 / math.sqrt(N); ex = 2.0 ** (-((n - 1) // 2))
    print("%2d  %9d  %+.6f (%.6f, z=%+.1f)  %+.6f  %+.6f      %.4f                    %.6f" % (n, N, e3, ex, (e3 - ex) / se, e5, e35, e35 / ex, se))
print()
def v(ok): return "CONFIRMED" if ok else "FAILED"
for n in (14, 18):
    if n in rows:
        N, e3, e5, e35 = rows[n]; se = 1 / math.sqrt(N); pred = 2.0 ** (-(n - 2) // 2)
        print("R7 at n = %d: estimate %.6f, prediction 2^(-(n-2)/2) = %.6f, deviation %+.2f s.e.: %s" % (n, e35, pred, (e35 - pred) / se, v(abs(e35 - pred) <= 3 * se)))
for a, b in ((11, 12), (15, 16)):
    if a in rows and b in rows:
        d = rows[b][3] - rows[a][3]; se = math.sqrt(1 / rows[a][0] + 1 / rows[b][0])
        print("R8 pair (%d, %d): E_%d - E_%d = %+.6f = %+.2f s.e.: %s" % (a, b, b, a, d, d / se, v(abs(d) <= 3 * se)))
bad = [n for n in rows if abs(rows[n][1] - 2.0 ** (-((n - 1) // 2))) > 3 / math.sqrt(rows[n][0])]
print("Control (Theorem T within 3 s.e. at every order):", "OK" if not bad else "FAILS at %s" % bad)
print()
print("Comparison with the exact values of Theorem P (pentagon_closed_forms.py; derived after this run):")
from pentagon_closed_forms import E5, E35
print(" n   E[{5}] exact   estimate   z       E[{3,5}] exact   estimate   z")
zs = []
for n in sorted(rows):
    N, e3, e5, e35 = rows[n]; se = 1 / math.sqrt(N); x5 = float(E5(n)); x35 = float(E35(n))
    zs += [(e5 - x5) / se, (e35 - x35) / se]
    print("%2d   %.7f      %+.7f  %+.2f    %.7f        %+.7f  %+.2f" % (n, x5, e5, (e5 - x5) / se, x35, e35, (e35 - x35) / se))
print("largest |z| among the %d comparisons: %.2f" % (len(zs), max(abs(z) for z in zs)))
