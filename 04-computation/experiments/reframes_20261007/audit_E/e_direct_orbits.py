#!/usr/bin/env python3
"""Audit E, task E (primary numerics): q_p = P(y and y+1 ever meet) by DIRECT big-integer orbits, no pair chain.

y is a uniform random B-bit integer (top bit forced), so its first B Terras parities are exactly i.i.d. fair; we need
N < B. For p >= 5 the integers grow, so no cycle is ever met. We record:
  * the first equal-time meeting T^t(y) = T^t(y+1), t <= N (an equal-time merge);
  * any unequal-time meeting T^m(y) = T^m'(y+1), m != m' <= N (THM-4606 statement 3 says: null set; expect none);
  * the merge-time histogram.
Usage: python3 e_direct_orbits.py p NSAMP N B SEED
"""
import sys, random, time, collections

def main():
    p, nsamp, N, B, seed = (int(a) for a in sys.argv[1:6])
    rnd = random.Random(seed)
    merged = 0
    times = collections.Counter()
    unequal = 0
    t0 = time.time()
    for s in range(nsamp):
        y = rnd.getrandbits(B) | (1 << (B - 1))
        a, b = y, y + 1
        seen_a = {a: 0}
        seen_b = {b: 0}
        hit = None
        for t in range(1, N + 1):
            a = (p * a + 1) >> 1 if a & 1 else a >> 1
            b = (p * b + 1) >> 1 if b & 1 else b >> 1
            if a == b:
                hit = t
                break
            # unequal-time meetings: a equals an earlier b, or b equals an earlier a
            if a in seen_b or b in seen_a:
                unequal += 1
            seen_a[a] = t; seen_b[b] = t
        if hit is not None:
            merged += 1; times[hit] += 1
    q = merged / nsamp
    se = (q * (1 - q) / nsamp) ** 0.5
    ts = sorted(times.elements())
    print(f"p={p} nsamp={nsamp} N={N} B={B} seed={seed}: merged {merged} -> q_p(N) = {q:.6f} +- {se:.6f};"
          f" unequal-time meetings {unequal}; latest merge {max(ts) if ts else None};"
          f" earliest {min(ts) if ts else None}; median {ts[len(ts)//2] if ts else None}; runtime {time.time()-t0:.0f}s", flush=True)
    # merge-time histogram in bins
    bins = [(1, 10), (11, 20), (21, 40), (41, 80), (81, 160), (161, 320), (321, 10**9)]
    print("  merge-time bins: " + ", ".join(f"[{lo},{hi if hi < 10**9 else 'inf'}]: {sum(c for t, c in times.items() if lo <= t <= hi)}" for lo, hi in bins), flush=True)
    print("  first merge times (count): " + ", ".join(f"{t}:{c}" for t, c in sorted(times.items())[:15]), flush=True)

if __name__ == '__main__':
    main()
