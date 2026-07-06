#!/usr/bin/env python3
"""
opus-2026-07-06-S95 -- HYP-4236(a): THE |B|=5 PRODUCTION CERTIFICATE TABLE.

For every 5-subset B of {1..12} (792 bases): compute the largest zone-clear component
of good(B, 1/13 + 1/500) (zones = the c<=7 block-worst radii at q=7 and degenerate
multiples, per the pinning-aware S93 closure), rationalize to the smallest denominator
Q with an integer sub-interval (A/Q, (A+LL)/Q), verify the toothMiss integer condition
for all five base runners over the widened m-window, and emit the row.

Output: 05-knowledge/results/cert_table_B5_opus_S95.tsv  (base | A | LL | Q | len | ok)
Plus: the 10 SMALLEST-component bases formatted as Lean row data for LRCClearRowsB5.
Every row is checked EXACTLY here (fractions); the Lean sample re-checks in-kernel.
"""
from fractions import Fraction as F
from itertools import combinations

RHO = F(1, 13); EPS = F(1, 500); SAFETY = F(1, 134)
L7 = F(1, 78)

zones = []
for p in range(1, 7):
    zones.append((F(p, 7), L7 / 2 + SAFETY))
for g in (2, 3):
    for p in range(1, 7 * g):
        if F(p, 7 * g).denominator == 7 * g:
            zones.append((F(p, 7 * g), L7 / (2 * g) + SAFETY))

def good_components(B, beta):
    ivs = []
    for w in B:
        for m in range(0, w + 1):
            a, b = (F(m) - beta) / w, (F(m) + beta) / w
            if b > 0 and a < 1:
                ivs.append((max(a, F(0)), min(b, F(1))))
    ivs.sort()
    merged = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    good, prev = [], F(0)
    for a, b in merged:
        if a > prev: good.append((prev, a))
        prev = max(prev, b)
    if prev < 1: good.append((prev, F(1)))
    return good

def largest_clear(B):
    best = (F(0), None, None)
    for a, b in good_components(list(B), RHO + EPS):
        segs = [(a, b)]
        for c, r in zones:
            da, db = c - r, c + r
            segs = [s for pair in segs for s in
                    ([(pair[0], min(pair[1], da))] if pair[0] < da else []) +
                    ([(max(pair[0], db), pair[1])] if pair[1] > db else [])
                    if s[0] < s[1]]
            if not segs: break
        for x, y in segs:
            if y - x > best[0]: best = (y - x, x, y)
    return best

def rationalize(lo, hi):
    for Q in range(2, 20000):
        Alo = int(lo * Q) + 1
        hiQ = hi * Q
        Ahi = int(hiQ) - (0 if hiQ != int(hiQ) else 1)
        if Alo < Ahi:
            return Q, Alo, Ahi - Alo
    return None

def toothmiss_ok(w, A, LL, Q):
    lo_m = (w * A - Q); hi_m = (w * (A + LL) + Q)
    for m in range(lo_m // Q - 1, hi_m // Q + 2):
        if w * A - Q <= Q * m <= w * (A + LL) + Q:
            if not (13 * w * (A + LL) < Q * (13 * m - 1) or Q * (13 * m + 1) < 13 * w * A):
                return False
    return True

rows, failures = [], []
for B in combinations(range(1, 13), 5):
    ln, lo, hi = largest_clear(B)
    if ln == 0:
        failures.append((B, "no clear component"))
        continue
    r = rationalize(lo, hi)
    if r is None:
        failures.append((B, "no small-Q rationalization"))
        continue
    Q, A, LL = r
    ok = all(toothmiss_ok(w, A, LL, Q) for w in B)
    rows.append((B, A, LL, Q, float(F(LL, Q)), ok))
    if not ok:
        failures.append((B, f"toothMiss FAILED at ({A},{LL},{Q})"))

with open("05-knowledge/results/cert_table_B5_opus_S95.tsv", "w") as f:
    f.write("base\tA\tLL\tQ\tlength\ttoothMiss_ok\n")
    for B, A, LL, Q, ln, ok in rows:
        f.write(f"{','.join(map(str,B))}\t{A}\t{LL}\t{Q}\t{ln:.6f}\t{ok}\n")

nok = sum(1 for r in rows if r[5])
print(f"bases: 792; rows emitted: {len(rows)}; toothMiss-verified: {nok}; failures: {len(failures)}")
for B, why in failures[:10]:
    print(f"  FAIL {B}: {why}")
rows.sort(key=lambda r: r[4])
print("\n10 smallest components (the Lean sample):")
for B, A, LL, Q, ln, ok in rows[:10]:
    print(f"  base {B}: ({A}/{Q}, {A+LL}/{Q}) len {ln:.5f} ok={ok}")
print("\nLean row data (base, A, LL, Q):")
for B, A, LL, Q, ln, ok in rows[:10]:
    print(f"  {list(B)}, {A}, {LL}, {Q}")
