#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c93 -- HYP-8020: THE MOD-p DESCENT OF WALL A --
diagnostics for the pipeline (descend the inverse theorem to F_p, prove it
with finite-field technology, cage-lift back to Z) and its computational
twin (the SAFE-PRIME SIEVE).

(a) BUDGET: are there enough safe/unsmooth primes >= 191 for the S-T caging
    product (ln prod > ln B_13 - ln lcm(2..14) ~ 657.6)?  Tiers: safe
    (p = 2q+1, q prime) and unsmooth (largest-pf(p-1) >= (p-1)/4).
(b) SEA SAMPLES at safe primes {719, 839, 983} vs SMOOTH controls
    {631, 991, 1009} (largest-pf fractions ~0.01-0.03): greedy draw success,
    distinct minimal covers, near-shell distances -- the S6/collapse contrast
    at matched sizes.
(c) RESONANCE-CONCENTRATION (the expansion diagnostic): the S6 mechanism
    says smooth p-1 admits subgroup clustering of ratio-resonances; measure
    the distribution of r(w,w') = |A(w) & A(w')| over random pairs -- tails
    (max, p99, variance) safe vs smooth.  Heavy tails at smooth primes =
    the concentration that feeds covers; light tails at safe primes = the
    expansion that starves them (the analytic mechanism the descent's
    finite-field proof would formalize).
"""
import sys, time, random
from math import log
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build
from lrc14_I13p1_minimal_covers_kps_S128c88 import ansatz_library
from lrc14_shell_census_sample_kps_S128c92 import sample_minimal

random.seed(93)

def primes_in(lo, hi):
    sieve = bytearray([1])*(hi+1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(hi**0.5)+1):
        if sieve[i]: sieve[i*i::i] = b'\x00'*len(sieve[i*i::i])
    return [i for i in range(lo, hi+1) if sieve[i]]

def largest_pf(n):
    f = 1; d = 2
    while d*d <= n:
        while n % d == 0: f = max(f, d); n //= d
        d += 1
    return max(f, n) if n > 1 else f

def is_safe(p):
    q = (p-1)//2
    return all(q % d for d in range(2, int(q**0.5)+1)) and q > 1

def part_a():
    print("== (a) the safe-prime / unsmooth-prime caging budget ==", flush=True)
    target = 13*(12*log(91) - log(13)) - log(360360)
    for name, pred in (("SAFE (p=2q+1)", is_safe),
                      ("UNSMOOTH (lpf(p-1) >= (p-1)/4)", lambda p: largest_pf(p-1) >= (p-1)/4)):
        s = 0.0; ps = []
        for p in primes_in(191, 100000):
            if p % 7 == 0 or not pred(p): continue
            ps.append(p); s += log(p)
            if s > target: break
        if s > target:
            print(f"  {name}: {len(ps)} primes, {ps[0]}..{ps[-1]}, sum ln = {s:.1f} > {target:.1f} OK", flush=True)
        else:
            print(f"  {name}: INSUFFICIENT below 100000 (sum {s:.1f} < {target:.1f})", flush=True)

def part_b():
    print("\n== (b) sea samples: safe {719, 839, 983} vs smooth controls {631, 991, 1009} ==", flush=True)
    for group, ps in (("SAFE", (719, 839, 983)), ("SMOOTH", (631, 991, 1009))):
        for p in ps:
            t0 = time.time()
            ctx = build(p)
            h, dk, fold, inv = ctx[0], ctx[1], ctx[5], ctx[6]
            lib = ansatz_library(p, ctx); libsets = list(lib.keys())
            seen = {}; fails = 0
            for _ in range(400):
                W = sample_minimal(p, ctx, topt=random.choice((2, 3, 4)))
                if W is None: fails += 1; continue
                seen[W] = seen.get(W, 0) + 1
            k = sum(seen.values()); distinct = len(seen); coll = k - distinct
            near = {}
            for W in seen:
                V = frozenset(fold(inv[w]) for w in W)
                d = min(len(V.symmetric_difference(F))//2 for F in libsets)
                near[d] = near.get(d, 0) + 1
            est = f"N~{k*k/(2*coll):.0f}" if coll else f"N>>{k*k//2}"
            print(f"  [{group}] p={p} (lpf-frac {largest_pf(p-1)/(p-1):.3f}): ok={k}/400 distinct={distinct} "
                  f"{est} near-shell {dict(sorted(near.items()))} [{time.time()-t0:.0f}s]", flush=True)

def part_c():
    print("\n== (c) resonance-concentration diagnostic (expansion, safe vs smooth) ==", flush=True)
    import statistics
    for group, ps in (("SAFE", (719, 839, 983)), ("SMOOTH", (631, 991, 1009))):
        for p in ps:
            ctx = build(p)
            h, dk, maskA = ctx[0], ctx[1], ctx[2]
            rs = []
            for _ in range(3000):
                w1, w2 = random.sample(range(1, h+1), 2)
                rs.append((maskA[w1] & maskA[w2]).bit_count())
            rs.sort()
            mean = statistics.mean(rs); var = statistics.variance(rs)
            p99 = rs[int(0.99*len(rs))]; mx = rs[-1]
            exp_mean = dk*dk/h  # heuristic independent-overlap mean
            print(f"  [{group}] p={p}: r-mean {mean:.2f} (indep ~{exp_mean:.2f}) var {var:.2f} "
                  f"p99 {p99} max {mx}  tail-ratio max/mean {mx/mean:.1f}", flush=True)

if __name__ == "__main__":
    part_a(); part_b(); part_c()
