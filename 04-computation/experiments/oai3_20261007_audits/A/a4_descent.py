#!/usr/bin/env python3
"""Audit A, item 4: from every (0, e), 0 < |e| <= 10^4, an explicit finite coin string reaches (0, 0) using ONLY the
generic chain step. Strategy: e > 0 as in the proof (depart c=0, halve at level 1 with c=0, return with c=1);
e < 0 by the mirror strategy (depart c=1 to k=-1, halve with c=0, return with c=1). Verify landing values
e -> (U(e)-1)/2 (e>0) and e -> -(U(|e|)-1)/2 (e<0), bound (U(e)-1)/2 <= (3e-1)/4, strict descent of |e|."""
from a3_drift import step, v2
def U(x):
    y = 3 * x + 1
    return y >> v2(y)
def cstep(k, N, c):
    sig = N & 1
    beta = c if k >= 0 else c ^ sig
    return step(k, N, beta)

bad = 0; maxlen = 0; maxexc = 0; worst = None
LIM = 10 ** 4
for e0 in list(range(-LIM, 0)) + list(range(1, LIM + 1)):
    k, N = 0, e0; n = 0; exc = 0
    while (k, N) != (0, 0):
        assert k == 0
        while N & 1 == 0 and N != 0:
            k, N = cstep(k, N, 0); n += 1
        if N == 0: break
        e = N
        if e > 0:
            k, N = cstep(k, N, 0); n += 1
            assert k == 1
            while N & 1 == 0:
                k, N = cstep(k, N, 0); n += 1; assert k == 1
            k, N = cstep(k, N, 1); n += 1
            if not (k == 0 and N == (U(e) - 1) // 2 and 0 <= N <= (3 * e - 1) / 4 and N < e): bad += 1
        else:
            k, N = cstep(k, N, 1); n += 1
            assert k == -1
            while N & 1 == 0:
                k, N = cstep(k, N, 0); n += 1; assert k == -1
            k, N = cstep(k, N, 1); n += 1
            if not (k == 0 and N == -((U(-e) - 1) // 2) and abs(N) < abs(e)): bad += 1
        exc += 1
        if n > 10 ** 6: bad += 1; break
    if n > maxlen: maxlen, worst = n, e0
    maxexc = max(maxexc, exc)
print(f"[descent] all 0<|e|<={LIM}: failures {bad}; longest absorbing coin string {maxlen} steps (from e={worst}); max excursions {maxexc}")
