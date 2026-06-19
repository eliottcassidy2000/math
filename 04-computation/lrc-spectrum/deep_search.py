"""
Deeper level-a search. For target level a, q=a(k+1)-1, t*=m/q.
The safe pool = {v in [1,hi] : ||v m/q|| >= a/q}. We must choose k speeds from the
pool whose EXACT M equals a/q (i.e. no other t beats it). Heuristic that worked at
k=7: a 'small block' (k - r smallest safe speeds, possibly with the run broken) + r
killers. We do a guided beam search:
  - start from the k smallest safe speeds (gives M >= a/q but usually > because the
    dense small run creates a competing t);
  - greedily REPLACE the speed that most reduces M (swap a small speed for a larger
    safe speed) until M stops decreasing or hits a/q.
This local search finds level-a sets without exponential blowup.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def safe_pool(m, q, a, hi):
    out = []
    for v in range(1, hi + 1):
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            out.append(v)
    return out


def Mval(S):
    return M_exact_fast(S)[0]


def local_min_M(k, a, m, q, hi):
    pool = safe_pool(m, q, a, hi)
    if len(pool) < k:
        return None
    poolset = pool
    S = set(pool[:k])
    curM = Mval(sorted(S))
    target = Fraction(a, q)
    improved = True
    it = 0
    while improved and it < 200:
        improved = False
        it += 1
        # try swapping each in-element for each out-element; greedy best single swap
        best_swap = None
        Slist = sorted(S)
        out_candidates = [v for v in poolset if v not in S]
        for vin in Slist:
            for vout in out_candidates:
                S2 = (S - {vin}) | {vout}
                if len(S2) != k:
                    continue
                if setgcd(S2) != 1:
                    continue
                M2 = Mval(sorted(S2))
                if M2 < curM:
                    if best_swap is None or M2 < best_swap[0]:
                        best_swap = (M2, vin, vout)
                        if M2 == target:
                            break
            if best_swap and best_swap[0] == target:
                break
        if best_swap and best_swap[0] < curM:
            _, vin, vout = best_swap
            S = (S - {vin}) | {vout}
            curM = best_swap[0]
            improved = True
    return sorted(S), curM


def best_for_k(k, amax=None, hi_mult=2):
    if amax is None:
        amax = int(k**0.5) + 4
    floor = Fraction(1, k + 1)
    best = None
    for a in range(2, amax + 1):
        q = a * (k + 1) - 1
        hi = hi_mult * q
        # try a few coprime m (m and structure matter); sample small set
        ms = [m for m in range(1, q // 2 + 1) if gcd(m, q) == 1]
        # limit to avoid blowup: try all if q small, else sample
        if len(ms) > 40:
            ms = ms[:40]
        for m in ms:
            res = local_min_M(k, a, m, q, hi)
            if res is None:
                continue
            S, M = res
            if M <= floor:
                continue
            g = M - floor
            cur = dict(S=S, M=M, g=g, gk2=g * k * k, a=level_a(M, k), why=(a, q, m))
            if best is None or cur['gk2'] < best['gk2']:
                best = cur
    return best


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [31]
    for k in ks:
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} ====", flush=True)
        b = best_for_k(k)
        if b:
            print("  S=", b['S'])
            print(f"  M={b['M']} (~{float(b['M']):.10f}) g*k^2={float(b['gk2']):.6f} "
                  f"a={float(b['a']):.4f} (a,q,m)={b['why']}", flush=True)
        else:
            print("  none", flush=True)
