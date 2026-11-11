#!/usr/bin/env python3
"""Verifier of THM-5134: the comparison of P3 at the multiples of 7, read from the session outputs.

For every group g0 > 1 the outputs print: measured factor | (a-priori prediction |) first-order value Phi | limit L.
P3-type comparison: is |measured - Phi| < |measured - L| ?   "Remainder against correction": is
|measured - Phi| > |Phi - L| ?
Files: cum_22.out, cum_24.out, cum_26.out, cum_28.out (cumulative), window_insample_*.out (four in-sample windows),
window_prereg_28_27.out (registered window), all in 04-computation/experiments/s92_prime_powers_20261111/.
"""
import os
import re
import sys

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "..",
                    "04-computation", "experiments", "s92_prime_powers_20261111")
FILES = ["cum_22.out", "cum_24.out", "cum_26.out", "cum_28.out", "window_insample_24_23.out",
         "window_insample_25_24.out", "window_insample_26_25.out", "window_insample_27_26.out",
         "window_prereg_28_27.out"]
ok_all = True


def groups(path, m):
    out = []
    for line in open(path, encoding="utf-8"):
        mm = re.match(r"\s*m =\s*(\d+) sheet (\d): (.*)", line)
        if not mm or int(mm.group(1)) != m:
            continue
        sheet = int(mm.group(2))
        for part in mm.group(3).split(";"):
            g = re.match(r"\s*g0=(\d+): ([\d.]+) \[[^\]]*\] \| (.*)", part)
            vals = [float(v) for v in g.group(3).split("|")]
            # cumulative: Phi | limit ; windows: a-priori | Phi (measured types) | limit
            phi, lim = vals[-2], vals[-1]
            out.append((sheet, int(g.group(1)), float(g.group(2)), phi, lim))
    return out


print("m = 49, groups g0 > 1 (sheet, g0): measured, Phi, limit -> closer to Phi than to the limit? remainder > correction?")
for f in FILES:
    gs = [g for g in groups(os.path.join(BASE, f), 49) if g[1] > 1]
    closer = sum(1 for (_, _, me, phi, lim) in gs if abs(me - phi) < abs(me - lim))
    larger = sum(1 for (_, _, me, phi, lim) in gs if abs(me - phi) > abs(phi - lim))
    print("  %-28s closer to Phi: %d of %d;  |measured - Phi| > |Phi - limit|: %d of %d;  %s"
          % (f, closer, len(gs), larger, len(gs),
             "  ".join("(%d, %d): %.4f %.4f %.4f" % g for g in gs)))
    if closer != 0 or larger != len(gs) or len(gs) != 4:
        ok_all = False
print("  -> at 49 the comparison fails in every group (0 of 4) in the nine outputs, and the remainder exceeds the "
      "first-order correction in every group: %s" % ok_all)

print("X = 2^28 (cum_28.out), m = 63 and m = 21")
want = {63: (6, 10), 21: (4, 6)}
for m in (63, 21):
    gs = [g for g in groups(os.path.join(BASE, "cum_28.out"), m) if g[1] > 1]
    closer = [g for g in gs if abs(g[2] - g[3]) < abs(g[2] - g[4])]
    failing = [g for g in gs if not abs(g[2] - g[3]) < abs(g[2] - g[4])]
    print("  m = %d: closer to Phi in %d of %d groups; failing groups (sheet, g0): %s; in the failing groups "
          "|measured - Phi| > |Phi - limit|: %s"
          % (m, len(closer), len(gs), [(g[0], g[1]) for g in failing],
             all(abs(g[2] - g[3]) > abs(g[3] - g[4]) for g in failing)))
    if (len(closer), len(gs)) != want[m] or not all(abs(g[2] - g[3]) > abs(g[3] - g[4]) for g in failing):
        ok_all = False
    if sorted(set(g[1] for g in failing)) != ([7, 21] if m == 63 else [7]):
        ok_all = False
print("  -> 6 of 10 at 63 and 4 of 6 at 21; the failing groups are g0 = 7 (63 and 21) and g0 = 21 (63): %s" % ok_all)

print("X = 2^28: remainders at 9, 27, 81 on the two sheets, 2^24 against 2^28 (cum_24.out, cum_28.out)")


def rems(path):
    out = {}
    for line in open(path, encoding="utf-8"):
        mm = re.match(r"\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", line)
        if mm:
            out[int(mm.group(1))] = (float(mm.group(2)), float(mm.group(3)))
    return out


r24, r28 = rems(os.path.join(BASE, "cum_24.out")), rems(os.path.join(BASE, "cum_28.out"))
for m in (9, 27, 81):
    print("  m = %2d: sheet 0 %.5f -> %.5f (%s), sheet 1 %.5f -> %.5f (%s)"
          % (m, r24[m][0], r28[m][0], "smaller" if r28[m][0] < r24[m][0] else "NOT smaller",
             r24[m][1], r28[m][1], "smaller" if r28[m][1] < r24[m][1] else "NOT smaller"))
    if not (r28[m][0] > r24[m][0] and r28[m][1] < r24[m][1]):
        ok_all = False
print("  -> from 2^24 to 2^28 the remainder at 9, 27, 81 is not smaller on sheet 0 and is smaller on sheet 1: %s" % ok_all)
print()
print("ALL CHECKS PASSED" if ok_all else "FAILURE")
sys.exit(0 if ok_all else 1)
