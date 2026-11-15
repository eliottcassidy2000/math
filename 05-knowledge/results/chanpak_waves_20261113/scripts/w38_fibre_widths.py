#!/usr/bin/env python3
"""Wave 38: exact H-fibres of P2 at orders 4..10.

Input: files pairs<n>.txt with one line "H P2" per isomorphism class (from w36_pair_spectrum.c over gentourng).
For each order n this prints
  * the step g_n = gcd of all within-fibre differences of P2 (24 at even order, 12 at odd order, by the
    mod-24 law of wave 37 and its odd-order analogue),
  * the mod-g law P2 mod g as a function of H mod (lcm),
  * the widest fibre: absolute width A = max P2 - min P2, width in steps A/g, number of P2 values, fill,
    relative width A/H,
  * the median relative width A/H over the fibres whose H lies in the central half of the class-weighted
    H distribution (the "typical" fibres), and the median number of steps.
Usage: python3 w38_fibre_widths.py DIR [orders...]
"""
import sys
from collections import defaultdict
from math import gcd
from statistics import median

d = sys.argv[1]
orders = [int(x) for x in sys.argv[2:]] or list(range(4, 11))
for n in orders:
    F = defaultdict(set)
    Hs = []
    for line in open(f"{d}/pairs{n}.txt"):
        h, p = map(int, line.split())
        F[h].add(p)
        Hs.append(h)
    Hs.sort()
    g = 0
    for h, ps in F.items():
        m = min(ps)
        for p in ps:
            g = gcd(g, p - m)
    # residue law: P2 mod 24 against H mod 24
    law = defaultdict(set)
    for h, ps in F.items():
        for p in ps:
            law[h % 24].add(p % 24)
    lawtxt = " ".join(f"{k}:{sorted(v)}" for k, v in sorted(law.items()))
    lo, hi = Hs[len(Hs) // 4], Hs[3 * len(Hs) // 4]
    rows = []
    for h, ps in F.items():
        A = max(ps) - min(ps)
        rows.append((h, len(ps), A))
    w = max(rows, key=lambda t: t[2])
    gg = g if g else 1
    cen = [(h, k, A) for (h, k, A) in rows if lo <= h <= hi and k >= 2]
    medrel = median(A / h for (h, k, A) in cen) if cen else 0
    medsteps = median(A / gg for (h, k, A) in cen) if cen else 0
    medfill = median(k / (A / gg + 1) for (h, k, A) in cen) if cen else 0
    meanratio = sum(p for ps in F.values() for p in ps) / max(1, sum(h for h, ps in F.items() for p in ps))
    print(f"order {n}: classes {len(Hs)}, distinct H {len(F)}, step g = {g}")
    print(f"  P2 mod 24 by H mod 24: {lawtxt}")
    print(f"  widest fibre: H = {w[0]}, {w[1]} values, width {w[2]} = {w[2] / gg:.0f} steps, "
          f"fill {w[1] / (w[2] / gg + 1):.3f}, width/H = {w[2] / w[0]:.3f}")
    print(f"  central H in [{lo}, {hi}]: {len(cen)} fibres with >= 2 values; median width/H = {medrel:.3f}, "
          f"median steps = {medsteps:.0f}, median fill = {medfill:.3f}; mean P2/H over classes = {meanratio:.2f} "
          f"(2(n-1) = {2 * (n - 1)})")
