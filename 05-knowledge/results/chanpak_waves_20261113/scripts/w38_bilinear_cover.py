#!/usr/bin/env python3
"""Wave 38: the bilinear (Z[sqrt -2]) Chan-Pak test.

Cyclic join with a point: H(cyc(A, B, *)) = Re(beta_A conj beta_B) = x_A x_B + 2 y_A y_B, with
beta = x + y sqrt(-2), x = H + P2, y = H (THM-5172 (C), THM-5222).  Using every block of order <= ORD (exact class
lists), mark all odd values <= X and list the odd numbers in [1, X] that are missed.  A miss is not a proof of
non-representability (larger blocks are not enumerated); a hit is an explicit tournament of order |A| + |B| + 1.
Also: the diagonal values x^2 + 2y^2 and the Ramanujan control x^2 + y^2 + 10 z^2 (its odd exceptions).
Usage: python3 w38_bilinear_cover.py PAIRS_DIR ORD X
"""
import sys
import numpy as np

PD, ORD, X = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
pairs = {(1, 1)}                           # one vertex: H = 1, P2 = 0 -> x = 1, y = 1
pairs |= {(3, 1)}                          # order 2: H = 1, P2 = 2
pairs |= {(7, 1), (9, 3)}                  # order 3: TT3 (1, 6), C3 (3, 6)
for n in range(4, ORD + 1):
    for line in open(f"{PD}/pairs{n}.txt"):
        h, p2 = map(int, line.split())
        pairs.add((h + p2, h))
P = np.array(sorted(pairs), dtype=np.int64)
print(f"blocks of order <= {ORD}: {len(P)} distinct (x, y) = (H + P2, H)")
hit = np.zeros(X + 1, dtype=bool)
xs, ys = P[:, 0], P[:, 1]
order = np.argsort(xs)
xs, ys = xs[order], ys[order]
for i in range(len(xs)):
    xa, ya = xs[i], ys[i]
    lim = X // xa
    j = np.searchsorted(xs, lim, side="right")
    if j == 0:
        continue
    v = xa * xs[:j] + 2 * ya * ys[:j]
    v = v[v <= X]
    hit[v] = True
odd = np.arange(1, X + 1, 2)
miss = odd[~hit[odd]]
print(f"odd N <= {X} missed by cyc(A, B, *) with |A|, |B| <= {ORD}: {len(miss)}")
print(f"  first misses: {miss[:40].tolist()}")
if len(miss):
    print(f"  largest misses: {miss[-10:].tolist()}")
for lo, hi in ((1, 10 ** 3), (10 ** 3, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, X)):
    if hi > X:
        continue
    seg = odd[(odd >= lo) & (odd < hi)]
    print(f"  [{lo}, {hi}): missed {int((~hit[seg]).sum())} of {len(seg)} odd numbers")
# Ramanujan control
R = np.zeros(20001, dtype=bool)
for z in range(0, 45):
    for x in range(0, 142):
        for y in range(x, 142):
            v = x * x + y * y + 10 * z * z
            if v <= 20000:
                R[v] = True
print("odd numbers <= 20000 not of the form x^2 + y^2 + 10 z^2 (Ramanujan; Ono-Soundararajan under GRH):",
      [n for n in range(1, 20001, 2) if not R[n]])
