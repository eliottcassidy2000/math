#!/usr/bin/env python3
"""Audit G, THM-4609 statement 5 (Z_5 example x/5, (x+4)/5, (6x+3)/5, (11x+2)/5, (16x+1)/5).
(A) Haar pairs (y, y+e) by direct integer orbits, y uniform mod 5^(T+64): merge fraction and debt visits to M = 1.
(B) the authors' integer test: n uniform with 2000 base-5 digits, n and n+1, merge within 100/1000/4000 steps.
(C) small integers: cycles reached from n <= 20000; equal-time merges of n, n+1 (n < 20000) within 3000 steps."""
import time, math, random
from mwsim import MW, run_pairs, fmt_run
z5 = MW(5, [1, 1, 6, 11, 16], [0, 4, 3, 2, 1])
print("Lambda =", sum(math.log(x / 5) for x in z5.m) / 5)
t0 = time.time()
cps = [16, 64, 256, 1024, 4096, 16384]
res = run_pairs(z5, 1, 16384, 2000, 101, cps)
print("(A) offset e = 1, N = 2000, T up to 16384:")
for line in fmt_run(res, cps): print("   ", line)
mt = sorted(res['merge_times']); print("    merge-time quantiles (50/90/99/max):", mt[len(mt)//2], mt[int(.9*len(mt))], mt[int(.99*len(mt))], mt[-1],
                                      f"  [{time.time()-t0:.0f}s]", flush=True)
cps2 = [16, 64, 256, 1024, 4096]
for e in (2, 3, 4, 5, 6, 10, 25, 125, 1000, -1):
    res = run_pairs(z5, e, 4096, 1000, 200 + abs(e), cps2)
    q = [res['merged_by'][c] / res['N'] for c in cps2]
    vw = res['vis_win'][4096]; mt = sorted(res['merge_times'])
    print(f"(A) e = {e:5d}: merged by T = {cps2}: {[round(x, 4) for x in q]};  visits in (1024,4096] per unmerged {vw[0]/max(vw[1],1):.3f};"
          f" latest merge {mt[-1] if mt else None}  [{time.time()-t0:.0f}s]", flush=True)
# (B)
rnd = random.Random(6006)
for T_ in (100, 1000, 4000):
    merged = 0; N = 1000
    for _ in range(N):
        n = rnd.randrange(5 ** 2000, 2 * 5 ** 2000); a, b = n, n + 1
        for t in range(T_):
            a, b = z5.T(a), z5.T(b)
            if a == b: merged += 1; break
    q = merged / N
    print(f"(B) 2000-base-5-digit n, n+1 merge within {T_}: {q:.4f} +- {(q*(1-q)/N)**.5:.4f}  [{time.time()-t0:.0f}s]", flush=True)
# (C)
cycles = {}
for n in range(1, 20001):
    x = n; seen = set()
    while x not in seen:
        seen.add(x); x = z5.T(x)
    cyc = [x]; y = z5.T(x)
    while y != x: cyc.append(y); y = z5.T(y)
    key = min(cyc); cycles[key] = cycles.get(key, 0) + 1
print("(C) cycles (min element: count) from 1..20000:", dict(sorted(cycles.items())))
cyc_list = {}
for key in cycles:
    c = [key]; y = z5.T(key)
    while y != key: c.append(y); y = z5.T(y)
    cyc_list[key] = c
print("    cycles:", {k: v for k, v in sorted(cyc_list.items())})
same = 0; tot = 0
for n in range(1, 20000):
    a, b = n, n + 1; ok = False
    for t in range(3000):
        a, b = z5.T(a), z5.T(b)
        if a == b: ok = True; break
    same += ok; tot += 1
print(f"(C) consecutive pairs n < 20000 merging at equal time within 3000 steps: {same}/{tot} = {same/tot:.4f}")
same = 0; tot = 0
for n in range(1, 20001):
    a, b = n, n + 1; ok = False
    for t in range(3000):
        a, b = z5.T(a), z5.T(b)
        if a == b: ok = True; break
    same += ok; tot += 1
print(f"(C) same with n <= 20000: {same}/{tot} = {same/tot:.4f}   [{time.time()-t0:.0f}s]")
