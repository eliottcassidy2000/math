#!/usr/bin/env python3
"""Wave 38: residue equidistribution of (H, Q) over the explicit two-apex family.

The family (chanpak_pair_line_census_20261115/TWO_APEX_COARSE_PLACEMENT.md): spine 0..m-1 transitive, apexes A, B with
A -> B, A beats spine vertex i iff a_i = 1, B beats it iff b_i = 1, a_0 = 1, a_{m-1} = 0, b != 0.  For m = 4..10 the
script computes (H, Q) for every word pair with the exact subset DP (w38_hq.c), then for primes l = 5, 7, 11, 13 the
distribution of (H mod l, Q mod l) over F_l^2: the number of empty cells, max/min cell ratio and the total-variation
distance to uniform.  A bilinear-sieve route to prime cyclic-join counts (README wave 38, GRAND_SYNTHESIS.md s.6) needs
such equidistribution for beta = (Q - H) + H sqrt(-2) modulo ideals of Z[sqrt -2].
Usage: python3 w38_two_apex_equidistribution.py HQ_BINARY
"""
import subprocess
import sys
from collections import Counter

HQ = sys.argv[1]


def two_apex_out(m, a, b):
    n = m + 2
    out = [0] * n
    for i in range(m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    A, B = m, m + 1
    for i in range(m):
        if (a >> i) & 1:
            out[A] |= 1 << i
        else:
            out[i] |= 1 << A
        if (b >> i) & 1:
            out[B] |= 1 << i
        else:
            out[i] |= 1 << B
    out[A] |= 1 << B
    return out


for m in range(4, 11):
    words = [(a, b) for a in range(1 << m) if (a & 1) and not (a >> (m - 1)) & 1 for b in range(1, 1 << m)]
    inp = "".join(f"{m + 2} " + " ".join(map(str, two_apex_out(m, a, b))) + "\n" for a, b in words)
    out = subprocess.run([HQ], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    HQs = [tuple(map(int, l.split())) for l in out if l.strip()]
    line = [f"m={m} (order {m + 2}): {len(HQs)} sources"]
    for l in (5, 7, 11, 13):
        c = Counter((h % l, q % l) for h, q in HQs)
        cells = l * l
        empty = cells - len(c)
        tv = 0.5 * sum(abs(c.get((x, y), 0) / len(HQs) - 1 / cells) for x in range(l) for y in range(l))
        mx, mn = max(c.values()), (min(c.values()) if empty == 0 else 0)
        line.append(f"l={l}: empty {empty}/{cells}, TV {tv:.3f}, max/min {mx / mn if mn else float('inf'):.2f}")
    print("; ".join(line))
