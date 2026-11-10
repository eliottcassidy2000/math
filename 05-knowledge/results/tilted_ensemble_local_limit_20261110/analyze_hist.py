#!/usr/bin/env python3
"""Local structure of the class-count function N_n(h) = #{tournament classes of order n with H = h} (n = 9, 10).
Input: out/hist_9.txt, out/hist_10.txt ("H count" lines produced by hist.c from nauty's gentourng)."""
import math, collections, os
here = os.path.dirname(os.path.abspath(__file__))
def load(n):
    d = {}
    for l in open(os.path.join(here, "out", f"hist_{n}.txt")):
        h, c = map(int, l.split()); d[h] = d.get(h, 0) + c
    return d
for n in (9, 10):
    d = load(n); tot = sum(d.values())
    mean = sum(h*c for h, c in d.items())/tot
    sd = math.sqrt(sum(h*h*c for h, c in d.items())/tot - mean**2)
    print(f"order {n}: classes {tot}, distinct values {len(d)}, max H {max(d)}, class-mean {mean:.1f}, class-sd {sd:.1f}")
    M = 1
    while (M + 2 in d) or (M + 2 in (7, 21)): M += 2
    print(f"   every odd h <= {M} occurs except 7 and 21; first missing value {M + 2}")
    lo, hi = int(mean - sd), int(mean + sd)
    zs = []; small = min(d.get(h, 0) for h in range(lo | 1, hi, 2))
    for h in range(lo | 1, hi, 2):
        nb = (d.get(h - 2, 0) + d.get(h + 2, 0))/2
        if nb > 0: zs.append((d.get(h, 0) - nb)/math.sqrt(1.5*nb))
    print(f"   bulk window [{lo},{hi}]: least class count {small}; rms normalised second difference {math.sqrt(sum(z*z for z in zs)/len(zs)):.2f} (1.00 for independent Poisson counts)")
    for q in (3, 4, 5, 7, 8):
        s = collections.defaultdict(float); k = collections.defaultdict(int)
        for h in range(lo | 1, hi, 2):
            nb = (d.get(h-2, 0) + d.get(h+2, 0) + d.get(h-4, 0) + d.get(h+4, 0))/4
            if nb > 0: s[h % q] += d.get(h, 0)/nb; k[h % q] += 1
        print(f"   relative class count by h mod {q}:", {r: round(s[r]/k[r], 3) for r in sorted(s)})
