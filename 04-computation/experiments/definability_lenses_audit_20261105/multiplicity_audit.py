#!/usr/bin/env python3
"""Independent re-reading of the P_d(1) lines of strong_ears_375.out (odd d = 3..373) and strong_ears_625.out
(odd d = 375..623), and recomputation of every multiplicity statement of THM-4884 (2) / note section 4 (corrected
text of commit 886693f53a)."""
import re, sys
base = r"C:/tmp/wt-s81/04-computation/experiments/h_spectrum_multiplicities_20261021/"
def grab(fn, lo):
    txt = open(base + fn).read()
    m = re.search(r"P_d\(1\)[^:]*:\s*\[([^\]]*)\]", txt)
    vals = [int(x) for x in m.group(1).split(',')]
    return {lo + 2 * i: v for i, v in enumerate(vals)}
m375 = grab("strong_ears_375.out", 3); m625 = grab("strong_ears_625.out", 375)
print(f"375 file: {len(m375)} values, d = {min(m375)}..{max(m375)}; 625 file: {len(m625)} values, d = {min(m625)}..{max(m625)}")
assert max(m375) == 373 and min(m625) == 375 and max(m625) == 623
m = dict(m375); m.update(m625)
# also cross-check the 375 file's own per-order polynomials for a few d (P_125 etc.) against P_d(1)
txt375 = open(base + "strong_ears_375.out").read()
polys = re.findall(r"P_(\d+) = ([^\n]*)", txt375)
bad = 0
for d, rhs in polys:
    d = int(d); coeffs = [int(c) for c in re.findall(r"(\d+)z\^", rhs)]
    if sum(coeffs) != m[d]: bad += 1; print("polynomial/P_d(1) mismatch at", d, rhs, m[d])
print(f"{len(polys)} explicit polynomials P_d in the 375 file cross-checked against P_d(1): mismatches {bad}")
def S(pred): return [h for h in sorted(m) if pred(m[h])]
print("m = 0:", S(lambda v: v == 0))
print("m = 1:", S(lambda v: v == 1))
print("m = 2:", S(lambda v: v == 2))
print("m = 3:", S(lambda v: v == 3))
print("m = 4:", S(lambda v: v == 4))
print("m = 5:", S(lambda v: v == 5))
above63 = {h: v for h, v in m.items() if 63 < h < 625 and v > 0}
print("min m(h) for attained 63 < h < 625:", min(above63.values()), "at", [h for h, v in above63.items() if v == min(above63.values())])
print("all attained h in (63,625) have m >= 6:", all(v >= 6 for v in above63.values()), "; >= 4:", all(v >= 4 for v in above63.values()))
print("largest h with m(h) <= 3:", max(h for h, v in m.items() if v <= 3), "; largest with m(h) <= 5:", max(h for h, v in m.items() if v <= 5))
r = {h: v for h, v in m.items() if 375 <= h <= 623}
print("min m(h) on [375, 623]:", min(r.values()), "at", [h for h, v in r.items() if v == min(r.values())], "; >= 73:", all(v >= 73 for v in r.values()))
vals = set(m.values())
print("multiplicity values NOT taken below 625, up to 40:", [k for k in range(0, 41) if k not in vals])
print("every value 0..38 other than 9,13,21,22,34,37 taken:", all(k in vals for k in range(0, 39) if k not in (9, 13, 21, 22, 34, 37)))
print("max multiplicity below 625:", max(m.values()), "at", [h for h, v in m.items() if v == max(m.values())])
print("values 39..80 not taken:", [k for k in range(39, 81) if k not in vals])
print("table h<=101:", " ".join(f"{h}:{m[h]}" for h in sorted(m) if h <= 101))
print("m(41) =", m[41], "; m(15) =", m[15], "; m(25) =", m[25], "; m(67) =", m[67], "; m(61) =", m[61], "; m(63) =", m[63])
# Q1 clauses
U = S(lambda v: v == 1)
print("Q1: U within [3,100]:", all(3 <= h <= 100 for h in U), "; |U| <= 12:", len(U) <= 12, "; contains 3,5,11,13:", all(x in U for x in (3, 5, 11, 13)), "; contains 41:", 41 in U,
      "; every element > 15 prime or 3*prime: vacuous" if not [h for h in U if h > 15] else "")
print("Q1: every multiplicity 1..10 occurs:", all(k in vals for k in range(1, 11)), "; missing:", [k for k in range(1, 11) if k not in vals])
