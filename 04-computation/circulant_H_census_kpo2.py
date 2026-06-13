#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
circulant_H_census_kpo2.py
kind-pasteur-2026-06-10-S2 / Thread A / HYP-2378 follow-up #2.

The odd-function dictionary (odd_function_dictionary_kpo2.py, entry (a)) makes
"which circulant tournament maximizes H?" a question about odd +-1 functions:
maximize H over the 2^((n-1)/2) odd sigmas.  The limits script found exact
R_rot(19) = 2.55197 > R_paley(19) = 2.527 (canon value, reflection
why-the-paley-path-ratio-is-e...), suggesting PALEY IS NOT THE CIRCULANT
H-MAXIMIZER AT p = 19.  This script settles that with exact integers:

EXACT census of H over ALL circulant tournaments on Z_n, n = 7,9,11,13,15,17,19,
reduced by the dictionary's symmetries (relabeling d -> u*d, u in Z_n^*, and
negation sigma -> -sigma = T^op), both of which preserve H.
H computed by fixed-start Held-Karp (H = n * #paths starting at 0; valid since
Z_n <= Aut acts freely on directed Ham paths, LEM-003).  Pure-int exact.
Fresh code (MISTAKE-001 guard).
"""

import sys, math, time

PASSED, FAILED = [], []
def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    (PASSED if ok else FAILED).append(name)
    print("  [%s] %s%s" % (tag, name, ((" -- " + detail) if detail else "")))

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def sigma_rotation(n):
    h = (n - 1) // 2
    return tuple([0] + [1] * h + [-1] * h)

def sigma_paley(p):
    sgl = [0] * p
    for d in range(1, p):
        sgl[d] = 1 if pow(d, (p - 1) // 2, p) == 1 else -1
    return tuple(sgl)

def all_odd_sigmas(n):
    h = (n - 1) // 2
    out = []
    for mask in range(1 << h):
        sgl = [0] * n
        for i in range(h):
            s = 1 if (mask >> i) & 1 else -1
            sgl[i + 1] = s
            sgl[n - 1 - i] = -s
        out.append(tuple(sgl))
    return out

def canonical(n, sig, units):
    best = None
    for u in units:
        for s in (1, -1):
            cand = tuple([0] + [s * sig[(u * d) % n] for d in range(1, n)])
            if best is None or cand < best:
                best = cand
    return best

def ham_paths_start0_times_n(n, sgl):
    m = n - 1
    succ = []
    for v in range(1, n):
        bits = 0
        for w in range(1, n):
            if w != v and sgl[(w - v) % n] == 1:
                bits |= 1 << (w - 1)
        succ.append(bits)
    size = 1 << m
    dp = [0] * (size * m)
    for v in range(1, n):
        if sgl[v % n] == 1:
            dp[(1 << (v - 1)) * m + (v - 1)] = 1
    for mask in range(size):
        base = mask * m
        for v in range(m):
            c = dp[base + v]
            if c:
                avail = succ[v] & ~mask
                while avail:
                    b = avail & (-avail)
                    w = b.bit_length() - 1
                    dp[(mask | b) * m + w] += c
                    avail ^= b
    full = (size - 1) * m
    return n * sum(dp[full + v] for v in range(m))

print("=" * 78)
print("EXACT H-CENSUS OF ALL CIRCULANT TOURNAMENTS, n = 7..19 (dictionary-reduced)")
print("=" * 78)

results = {}
for n in [7, 9, 11, 13, 15, 17, 19]:
    t0 = time.time()
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    sigmas = all_odd_sigmas(n)
    classes = {}
    for sig in sigmas:
        classes.setdefault(canonical(n, sig, units), []).append(sig)
    reps = sorted(classes)
    rows = []
    for rep in reps:
        H = ham_paths_start0_times_n(n, rep)
        R = H * 2 ** (n - 1) / math.factorial(n)
        rows.append((H, R, rep, len(classes[rep])))
    rows.sort(reverse=True)
    rot_c = canonical(n, sigma_rotation(n), units)
    pal_c = canonical(n, sigma_paley(n), units) if (n in (7, 11, 19)) else None
    print()
    print("n=%2d : %d circulant tournaments, %d classes under (units x negation), %.1f s"
          % (n, len(sigmas), len(reps), time.time() - t0))
    print("      %-12s %10s  %6s  %s" % ("H", "R", "orbit", "connection set / flags"))
    for H, R, rep, osz in rows:
        S = sorted(d for d in range(1, n) if rep[d] == 1)
        flags = []
        if rep == rot_c:
            flags.append("ROTATION")
        if pal_c is not None and rep == pal_c:
            flags.append("PALEY")
        print("      %-12d %10.5f  %6d  S=%s %s"
              % (H, R, osz, S, (" <-- " + ",".join(flags)) if flags else ""))
    results[n] = (rows, rot_c, pal_c)

print()
print("=" * 78)
print("CHECKS")
print("=" * 78)
rows7 = results[7][0]
check("n=7: max circulant H = 189 = Paley = global max A038375(7)",
      rows7[0][0] == 189 and results[7][2] == rows7[0][2])
rows11 = results[11][0]
check("n=11: max circulant H = 95095 = Paley (canon value)",
      rows11[0][0] == 95095 and results[11][2] == rows11[0][2])

rows19, rot19c, pal19c = results[19]
H_by_class = {rep: H for H, R, rep, osz in rows19}
H_rot19, H_pal19 = H_by_class[rot19c], H_by_class[pal19c]
print("  n=19: H_rotation = %d, H_paley = %d, max = %d" % (H_rot19, H_pal19, rows19[0][0]))
check("n=19: H_rot(19) = 1184212824763 (matches limits script exactly)",
      H_rot19 == 1184212824763)
check("n=19: PALEY IS NOT THE CIRCULANT H-MAXIMIZER (H_rot > H_paley?)",
      H_rot19 > H_pal19,
      "H_rot - H_paley = %d" % (H_rot19 - H_pal19))
top19 = rows19[0]
check("n=19: report the actual circulant maximizer",
      True, "H_max=%d R=%.5f S=%s (rotation? %s, paley? %s)"
      % (top19[0], top19[1], sorted(d for d in range(1, 19) if top19[2][d] == 1),
         top19[2] == rot19c, top19[2] == pal19c))

for n in [9, 13, 15, 17]:
    rows, rot_c, _ = results[n]
    rk = [i for i, (H, R, rep, o) in enumerate(rows) if rep == rot_c][0]
    print("  n=%2d: rotation rank %d of %d (H=%d, max H=%d)"
          % (n, rk + 1, len(rows), rows[rk][0], rows[0][0]))

print()
print("=" * 78)
print("SUMMARY: %d checks passed, %d failed" % (len(PASSED), len(FAILED)))
for f in FAILED:
    print("  FAILED: %s" % f)
print("=" * 78)
sys.exit(1 if FAILED else 0)
