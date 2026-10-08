#!/usr/bin/env python3
"""Audit E, task E (tail part): long-horizon pair-chain Monte Carlo (own integer-coordinate chain, validated against
direct orbits in a_table_check.py (A3)), to measure merges LATER than the author's horizon T = 400.

Each path from (0, e0) runs until absorption, until |f| > 10^CUT (recorded as escaped; the largest |f| seen before any
later merge is reported to show the cut is harmless), or until TMAX.
Usage: python3 e_chain_tail.py p NPATHS TMAX CUT SEED [e0]
"""
import sys, random, collections, math, time

def main():
    p, npaths, TMAX, CUT, seed = (int(a) for a in sys.argv[1:6])
    e0 = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    rnd = random.Random(seed)

    class Pow(list):                       # p^i, extended on demand (|k| stays ~ sqrt(t))
        def __getitem__(self, i):
            while i >= len(self): self.append(self[-1] * p if len(self) else 1)
            return list.__getitem__(self, i)
    pw = Pow([1])
    lim = 10 ** CUT
    merged = 0; escaped = 0; timeout = 0
    times = []
    maxf_before_merge = 0.0
    t0 = time.time()
    for _ in range(npaths):
        k, E = 0, e0
        maxlogf = 0.0
        t = 0
        result = None
        while t < TMAX:
            t += 1
            b = rnd.getrandbits(1)
            if E & 1 == 0:
                if b == 0: E >>= 1
                elif k >= 0: E = (p * E + 1 - pw[k]) >> 1
                else: E = (p * E + pw[-k] - 1) >> 1
            elif b == 0:
                E = (p * E + 1) >> 1 if k >= 0 else (E + pw[-k - 1]) >> 1
                k += 1
            else:
                E = (E - pw[k - 1]) >> 1 if k >= 1 else (p * E - 1) >> 1
                k -= 1
            if k == 0 and E == 0:
                result = 'm'; break
            a = abs(E); q = pw[abs(k)]
            if a > lim * q:
                result = 'e'; break
            if t % 16 == 0 and a > q:
                lf = math.log10(a) - math.log10(q)     # math.log10 accepts arbitrarily large ints
                if lf > maxlogf: maxlogf = lf
        if result == 'm':
            merged += 1; times.append(t)
            maxf_before_merge = max(maxf_before_merge, maxlogf)
        elif result == 'e': escaped += 1
        else: timeout += 1
    q = merged / npaths
    se = (q * (1 - q) / npaths) ** 0.5
    ts = sorted(times)
    gt = lambda T: sum(1 for x in ts if x > T)
    print(f"p={p} e0={e0} npaths={npaths} TMAX={TMAX} cut=1e{CUT} seed={seed}: merged {merged} -> q = {q:.6f} +- {se:.6f};"
          f" escaped {escaped}, timed out {timeout}; merges after t=100: {gt(100)}, after 200: {gt(200)}, after 400: {gt(400)},"
          f" after 1000: {gt(1000)}; latest {ts[-1] if ts else None};"
          f" max log10|f| (sampled every 16 steps) on a path that later merged: {maxf_before_merge:.1f}; runtime {time.time()-t0:.0f}s", flush=True)
    late = [x for x in ts if x > 200]
    print(f"  merge times > 200: {late[:40]}", flush=True)
    cnt = collections.Counter(ts)
    print(f"  merges with t <= 22: {sum(c for t, c in cnt.items() if t <= 22)}; after 22: {sum(c for t, c in cnt.items() if t > 22)};"
          f" per-time t<=22: {dict(sorted((t, c) for t, c in cnt.items() if t <= 22))}", flush=True)

if __name__ == '__main__':
    main()
