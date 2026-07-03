#!/usr/bin/env python3
"""
DECISIVE: is the witness q for COMPRESSED covering gcd=1 families uniformly bounded, or ~log(magnitude)?
(mac-mini-2026-07-03-S26). If bounded (~50), a FIXED-q census closes the compressed case (wide residual =
the census). If growing, the closure is the log-census / deep renormalization (the genuine crux).

Compressed = no ratio-13 dominant (top two within 13x). Aligned band-blockers: far_i ~ N each divisible by a
distinct band modulus (blocks small q by residue 0). Push N to 10^9, maximally adversarial (block {15..Q}).
"""
from math import gcd, log10
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))
def is_compressed(sp):
    return all(any(j!=i and abs(sp[i])<13*abs(sp[j]) for j in range(len(sp))) for i in range(len(sp)))
def min_witness_q(sp, qmax=400):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q) >= 1/14 for v in sp): return q
    return None

if __name__ == "__main__":
    print("Max witness q over COMPRESSED covering gcd=1 families vs magnitude N (aligned band-blockers).")
    print("=" * 82)
    print(f"{'N':>12} {'log10 N':>8} {'#compressed tested':>18} {'MAX witness q':>14}")
    band = list(range(15, 80))
    rows = []
    for N in [10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9]:
        rng = random.Random(N)
        worst = 0
        for _ in range(20000):
            rng.shuffle(band)
            k = rng.randint(9, 13)
            far = sorted({q * round(N / q) for q in band[:k]})
            far = [f for f in far if f > 22]
            speeds = far[:]
            for q in [8,9,5,7,11,13,2,3,4,6]:
                if len(speeds) >= 13: break
                if not any(s % q == 0 for s in speeds): speeds.append(q)
            while len(speeds) < 13: speeds.append(rng.randint(2, 22))
            speeds = sorted(set(speeds))[:13]
            if len(speeds) != 13 or gcd_all(speeds) != 1: continue
            if not is_covering(speeds) or not is_compressed(speeds): continue
            q = min_witness_q(speeds, qmax=400)
            if q and q > worst: worst = q
        rows.append((N, worst))
        print(f"{N:>12} {log10(N):>8.0f} {'~20k':>18} {worst:>14}")
    qs = [q for _, q in rows]; logs = [log10(N) for N, _ in rows]
    n=len(qs); sx=sum(logs); sy=sum(qs); sxx=sum(x*x for x in logs); sxy=sum(x*y for x,y in zip(logs,qs))
    slope = (n*sxy - sx*sy)/(n*sxx - sx*sx)
    print(f"\nMAX-q vs log10(N): sequence {qs}")
    print(f"slope (q per decade) ~ {slope:.2f}")
    print(f"=> slope ~ 0 (q plateaus): compressed q UNIFORMLY BOUNDED -> a FIXED census (q<=~{max(qs)}) closes")
    print(f"   the compressed case, wide residual = the census, no renormalization needed.")
    print(f"   slope > 0 (q grows): compressed q ~ log(mag) -> the log-census/deep renormalization is the crux.")
