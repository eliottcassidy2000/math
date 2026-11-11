#!/usr/bin/env python3
"""Verifier checks for HYP-9661 (A), (B), (C), (D): numbers recomputed from the committed receipts.

Reads only committed text outputs (slope_lists_le73.txt, source_supply.out, line_density_classes.out,
line_density_order10.out, head_values_check.out). No session code is imported.
"""
import math
import os
import re
import sys
from fractions import Fraction

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))
EXP = os.path.join(ROOT, "04-computation", "experiments", "s91_f1_local_densities_20261110", "rooted_spectrum")
FAILS = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" : " + detail if detail else ""))
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------- slope lists
lists = {}
with open(os.path.join(EXP, "slope_lists_le73.txt")) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        head, rest = line.split(":", 1)
        D = int(head)
        lists[D] = [(int(a), int(b)) for a, b in (tok.split(":") for tok in rest.split())]
check("35 attained offsets <= 73", len(lists) == 35, str(len(lists)))
check("1066 odd parts in all", sum(len(v) for v in lists.values()) == 1066)


def is_slope(D, U):
    if D not in lists or U < 1:
        return False
    j = 0
    d = U
    while d % 2 == 0:
        d //= 2
        j += 1
    for dd, jj in lists[D]:
        if dd == d:
            return j >= jj
    return False


def count(D, lo, hi):
    """number of slopes U of offset D with lo <= U <= hi (integers)"""
    return sum(1 for U in range(lo, hi + 1) if is_slope(D, U))


# windows [h, 3h/2]
win = {}
for h in sorted(lists):
    n_int = (3 * h) // 2 - h + 1
    win[h] = (count(h, h, (3 * h) // 2), n_int)
full = [h for h in win if win[h][0] == win[h][1]]
check("full windows only at 1 and 5", full == [1, 5], str(full))
excl = [h for h in sorted(lists) if len(lists[h]) < win[h][1]]
check("(A3) bound excludes a full window exactly at 3, 13, 35, 39, 43, 47, 63, 67",
      excl == [3, 13, 35, 39, 43, 47, 63, 67], str(excl))
check("odd parts at 33, 65 are 34, 66 against 17, 33 integers",
      (len(lists[33]), len(lists[65]), win[33][1], win[65][1]) == (34, 66, 17, 33))
sh = [win[h][0] / win[h][1] for h in sorted(lists) if 37 <= h <= 73]
check("mean covered share on 37..73 is 0.500 (min 0.156, max 0.743)",
      abs(sum(sh) / len(sh) - 0.500) < 5e-4 and abs(min(sh) - 0.156) < 5e-4 and abs(max(sh) - 0.743) < 5e-4,
      "%.4f %.4f %.4f" % (sum(sh) / len(sh), min(sh), max(sh)))

# shares of [h, 9h/8], [h, 2h), [2h, 4h) at h = 17, 33, 65 (HYP-9661 (A), "At a single offset")
quoted = {17: (0.67, 0.41, 0.29), 33: (0.80, 0.67, 0.44), 65: (0.78, 0.57, 0.38)}
for h, q in quoted.items():
    a = count(h, h, (9 * h) // 8) / ((9 * h) // 8 - h + 1)
    b = count(h, h, 2 * h - 1) / h
    c = count(h, 2 * h, 4 * h - 1) / (2 * h)
    check("single-offset shares at h = %d" % h,
          all(abs(x - y) < 0.005 + 1e-9 for x, y in zip((a, b, c), q)), "%.3f %.3f %.3f" % (a, b, c))

# r(N) on [3, 225): D odd, D <= N/3, (N - D)/2 a slope of D
zero = []
rvals = {}
for N in range(3, 225, 2):
    r = sum(1 for D in range(1, N // 3 + 1, 2) if is_slope(D, (N - D) // 2))
    rvals[N] = r
    if r == 0:
        zero.append(N)
check("r(N) = 0 on [3, 225) exactly at 7, 11, 13, 21, 47, 67", zero == [7, 11, 13, 21, 47, 67], str(zero))
check("r(N) >= 5 on [129, 225)", min(rvals[N] for N in range(129, 225, 2)) == 5)
check("r(N) >= N/64 on [129, 225)", all(64 * rvals[N] >= N for N in range(129, 225, 2)))

# heads of (C): the four numbers are D + 2 U_even with (D, U) sources
for N, D, Ue, U in ((11, 5, 3, 7), (13, 5, 4, 9), (47, 17, 15, 31), (67, 33, 17, 41)):
    check("head %d = %d + 2 * %d, (D, U) = (%d, %d) a source, D > N/3" % (N, D, Ue, D, U),
          D + 2 * Ue == N and is_slope(D, U) and 3 * D > N)

# (D) Consequence: offset 3 has no dyadic slope; cone 5 + 2^m supplied from m = 6; head 13, 21, 37
check("offset 3: only odd part 3", lists[3] == [(3, 0)])
check("offset 5: dyadic slopes 2^i iff i >= 5", (1, 5) in lists[5])
cone5 = [m for m in range(1, 12) if is_slope(5, 2 ** (m - 1))]
check("cone 5 + 2^m (value D + 2U, U = 2^(m-1)) supplied exactly from m = 6", cone5 == list(range(6, 12)), str(cone5))
check("grid: 7 = N(1,1,1), 21 = N(1,2,2); head of level 2 is 13, 21, 37",
      1 + 2 + 4 == 7 and 1 + 4 + 16 == 21 and [5 + 2 ** (r + 2) for r in (1, 2, 3)] == [13, 21, 37])

# ---------------------------------------------------------------- source_supply.out sections (5), (6)
txt = open(os.path.join(EXP, "source_supply.out")).read()
sec5 = txt[txt.index("(5) h:"):txt.index("(6) h:")]
sec6 = txt[txt.index("(6) h:"):]
w5 = {int(a): (int(b), int(c)) for a, b, c in re.findall(r"(\d+):(\d+)/(\d+)", sec5)}
check("section (5) of source_supply.out = windows recomputed from the slope lists", w5 == win)
m = {int(a): int(b) for a, b in re.findall(r"(\d+):(\d+)", sec6.split("Pearson")[0].split("\n", 1)[1])}
check("section (6) lists m(h) for 1 and the 34 attained 3 <= h <= 73", sorted(m) == sorted(lists), str(len(m)))


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))


for lo, n_exp, q in ((23, 26, 0.68), (15, 29, 0.66), (3, 34, 0.41)):
    hs = [h for h in sorted(lists) if lo <= h <= 73]
    c = pearson([win[h][0] / win[h][1] for h in hs], [m[h] for h in hs])
    check("Pearson %d <= h <= 73: %d offsets, %.2f" % (lo, n_exp, q), len(hs) == n_exp and abs(c - q) < 0.005,
          "%d offsets, %.4f" % (len(hs), c))
two = [h for h in sorted(m) if m[h] == 2]
check("m(h) = 2 exactly at 9, 17, 19, 35, 39, 47, 63 (h <= 73)", two == [9, 17, 19, 35, 39, 47, 63], str(two))
check("m(h) >= 6 for 63 < h <= 73, m(67) = 6", all(m[h] >= 6 for h in m if h > 63) and m[67] == 6)
rank = sorted((Fraction(*win[h]), h) for h in lists if 23 <= h <= 73)
five = [h for _, h in rank[:5]]
check("five poorest of [23, 73]: 63, 47, 67, 35, 39", five == [63, 47, 67, 35, 39], str(five))
check("the sixth share is strictly larger (27 and 55 at 0.357)", rank[5][0] > rank[4][0], str(rank[5]))
sh2 = {h: round(win[h][0] / win[h][1], 3) for h in (35, 39, 47, 63, 67, 9, 17, 19)}
check("shares 35, 39, 47, 63 = 0.278, 0.350, 0.208, 0.156; 67 = 9/34 = 0.265; 9, 17, 19 = 0.60, 0.44, 0.50",
      (sh2[35], sh2[39], sh2[47], sh2[63], sh2[67]) == (0.278, 0.35, 0.208, 0.156, 0.265)
      and win[67] == (9, 34) and (round(sh2[9], 2), round(sh2[17], 2), round(sh2[19], 2)) == (0.6, 0.44, 0.5),
      str(sh2))

# ---------------------------------------------------------------- line_density_classes.out
rows = {}
for line in open(os.path.join(EXP, "line_density_classes.out")):
    mm = re.match(r"\s*(\d+) \| ([\d.]+) \| ([\d.]+) \(([\d.]+), ([\d.]+)\) \| ([\d.]+) \| ([\d.]+)", line)
    if mm:
        rows[int(mm.group(1))] = tuple(float(x) for x in mm.groups()[1:])
check("3 | N: 0.578 against 0.557, ratio 1.038", abs(rows[3][0] - 0.578) < 6e-4 and abs(rows[3][1] - 0.557) < 6e-4
      and rows[3][4] == 1.038, str(rows[3]))
check("ratios at 5, 7, 11, 13 are 1.004, 1.001, 0.999, 1.001",
      [rows[p][4] for p in (5, 7, 11, 13)] == [1.004, 1.001, 0.999, 1.001])
check("share with p | D: 0.357, 0.201, 0.143, 0.092, 0.078",
      [rows[p][5] for p in (3, 5, 7, 11, 13)] == [0.357, 0.201, 0.143, 0.092, 0.078])
pois = math.log(1 - rows[3][0]) / math.log(1 - rows[3][1])
check("Poisson reading: -ln(1 - share) in the ratio 1.06", abs(pois - 1.06) < 0.005, "%.4f" % pois)
check("F1 ratio of class factors at 3 is 1.375", Fraction(1) - Fraction(1, 9) + Fraction(1, 3) == Fraction(11, 8) * (1 - Fraction(1, 9)))

# ---------------------------------------------------------------- line_density_order10.out
t10 = open(os.path.join(EXP, "line_density_order10.out")).read()
blocks = re.findall(r"\[2\^(\d+), 2\^\d+\): min\s+(\d+) \(N = \d+\), median\s+\d+, mean\s+[\d.]+, share ([\d.]+)", t10)
shares = {int(k): float(s) for k, _, s in blocks}
mins = {int(k): int(mn) for k, mn, _ in blocks}
check("shares k = 7..13", [shares[k] for k in range(7, 14)] == [0.296, 0.334, 0.381, 0.464, 0.518, 0.588, 0.616])
check("block minima k = 7..13", [mins[k] for k in range(7, 14)] == [5, 7, 24, 63, 146, 352, 794])
check("share of the block k = 6 is 0.271 (bound 0.27 of (A) 2)", shares[6] == 0.271)
check("shares increasing k = 6..13; k = 14 (unsaturated) is 0.602", all(shares[k] < shares[k + 1] for k in range(6, 13))
      and shares[14] == 0.602)
check("no violation of r >= N/64 on [129, 2^14)", "r(N) >= N/64 for every odd N in [129, 2^14): yes; violations 0" in t10)
check("top decile on [4096, 8192): 0.94 (odd U), 0.84 (even U)",
      "0.92 0.94   overall 0.488" in t10 and "0.84 0.84   overall 0.692" in t10)

# ---------------------------------------------------------------- head_values_check.out
hv = open(os.path.join(EXP, "head_values_check.out")).read()
check("heads at a non-sink root miss [1, 7, 21] below 9153", "not a sink: odd N < 9153 missed: [1, 7, 21]" in hv)
check("heads with D <= N/3 miss [1, 7, 11, 13, 21, 47, 67]",
      "(iii) heads with D <= N/3: odd N < 9153 missed: [1, 7, 11, 13, 21, 47, 67]" in hv)

print()
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
print("ALL VERIFIER CHECKS PASSED")
