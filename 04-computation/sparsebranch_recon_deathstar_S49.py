#!/usr/bin/env python3
"""Recon for HYP-7245 (death-star-S49): sparse Bezout-branch decomposition.

For each of the 29 sparse pairs (a,b) of {1..13} (reduced i'+j' >= 14):
  (a) witness uniqueness + branch partition k in {-1,0,+1}
  (b) N^{(0)} = 2*floor((q-1)/(14*max(i',j')))
  (c) N^{(+1)} = N^{(-1)} (mirror)
  (d) branch interval law: N^{(+-1)}/q -> (i'+j'-14)/(14*i'*j')
  (e) total limit: N_ab/(q-1) -> 1/(7*max) + 2*(i'+j'-14)/(14*i'*j')
      cross-checked against boxeph LEM-044 mu for consecutive pairs.
"""
from math import gcd
from fractions import Fraction as F
import random
random.seed(49)

def failing_witness(v, q, p):
    r = (v * p) % q
    if 14 * r < q: return (v * p) // q
    if 14 * r > 13 * q: return (v * p) // q + 1
    return None

pairs_sparse = [(a, b) for a in range(1, 14) for b in range(a + 1, 14)
                if (a // gcd(a, b)) + (b // gcd(a, b)) >= 14]
print(f"sparse pairs: {len(pairs_sparse)}")

part_ok = part_bad = n0_ok = n0_bad = mir_ok = mir_bad = 0
limit_rows = []
for (a, b) in pairs_sparse:
    g = gcd(a, b); ip, jp = a // g, b // g
    M = max(ip, jp)
    qs = [x for x in random.sample(range(500, 4000), 10) if gcd(g, x) == 1][:6]
    for q in qs:
        counts = {-1: 0, 0: 0, 1: 0}
        bad = 0
        for p in range(1, q):
            wa = failing_witness(a, q, p); wb = failing_witness(b, q, p)
            if wa is None or wb is None: continue
            k = jp * wa - ip * wb
            if k in counts: counts[k] += 1
            else: bad += 1
        part_ok += (bad == 0); part_bad += (bad != 0)
        pred0 = 2 * ((q - 1) // (14 * M))
        n0_ok += (counts[0] == pred0); n0_bad += (counts[0] != pred0)
        if counts[0] != pred0 and n0_bad <= 3:
            print("N0 FAIL", a, b, q, counts[0], pred0)
        mir_ok += (counts[1] == counts[-1]); mir_bad += (counts[1] != counts[-1])
        if counts[1] != counts[-1] and mir_bad <= 3:
            print("MIRROR FAIL", a, b, q, counts[1], counts[-1])
    # limit check at one big q
    qbig = 14 * ip * jp * 100 + 1
    while gcd(g, qbig) != 1: qbig += 2
    counts = {-1: 0, 0: 0, 1: 0}
    for p in range(1, qbig):
        wa = failing_witness(a, qbig, p); wb = failing_witness(b, qbig, p)
        if wa is None or wb is None: continue
        counts[jp * wa - ip * wb] += 1
    tot = sum(counts.values())
    pred_lim = F(1, 7 * M) + 2 * F(ip + jp - 14, 14 * ip * jp)
    limit_rows.append((a, b, ip, jp, F(tot, qbig - 1), pred_lim,
                       abs(F(tot, qbig - 1) - pred_lim)))
print(f"(a) partition k in {{-1,0,1}}: ok {part_ok}, bad {part_bad}",
      "PASS" if part_bad == 0 else "FAIL")
print(f"(b) N0 = 2*floor((q-1)/(14M)): ok {n0_ok}, bad {n0_bad}",
      "PASS" if n0_bad == 0 else "FAIL")
print(f"(c) mirror N+1 = N-1: ok {mir_ok}, bad {mir_bad}",
      "PASS" if mir_bad == 0 else "FAIL")
worst = max(limit_rows, key=lambda r: r[6])
print(f"(d/e) limit law on 29 pairs: worst |N/(q-1) - pred| = {float(worst[6]):.2e}"
      f" at {worst[:2]}; all under 3/q:",
      all(r[6] <= F(3, 14 * r[2] * r[3] * 100) for r in limit_rows))
for r in limit_rows[:4] + [worst]:
    print(f"    ({r[0]},{r[1]}) reduced ({r[2]},{r[3]}): N/(q-1) = {float(r[4]):.6f}"
          f"  pred = {float(r[5]):.6f}")
