"""
Lean greedy level-a search. Exact M throughout.
For target level a: q=a(k+1)-1, t*=m/q (few curated m). Start S = k smallest safe speeds.
Greedily swap one in-speed for one out-speed to reduce exact M, until target a/q reached
or no improvement. Keep pool and swap candidates small for speed.
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


def Mval(S):
    return M_exact_fast(S)[0]


def safe_pool(m, q, a, hi):
    out = []
    for v in range(1, hi + 1):
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            out.append(v)
    return out


def greedy(k, a, m, q, hi, max_iter=40):
    pool = safe_pool(m, q, a, hi)
    if len(pool) < k:
        return None
    target = Fraction(a, q)
    S = sorted(pool[:k])
    curM = Mval(S)
    if curM == target:
        return S, curM
    for _ in range(max_iter):
        Sset = set(S)
        # out-candidates: larger safe speeds not in S (these are the 'killers')
        outpool = [v for v in pool if v not in Sset]
        # prefer swapping smaller in-speeds for larger out-speeds: take first improving
        moved = False
        # try in ascending in-speed, descending no -- ascending out gives nearest killers
        for vin in S:  # ascending
            for vout in outpool:  # ascending
                if vout <= vin:
                    continue
                S2 = (Sset - {vin}) | {vout}
                if setgcd(S2) != 1:
                    continue
                M2 = Mval(sorted(S2))
                if M2 < curM:
                    S = sorted(S2)
                    curM = M2
                    moved = True
                    break
            if moved:
                break
        if not moved or curM == target:
            break
    return S, curM


def best_for_k(k, amax=None, hi_mult=2, m_samples=8):
    if amax is None:
        amax = int(k**0.5) + 3
    floor = Fraction(1, k + 1)
    best = None
    for a in range(2, amax + 1):
        q = a * (k + 1) - 1
        hi = hi_mult * q
        ms = [m for m in range(1, q) if gcd(m, q) == 1]
        # sample evenly spaced m to get variety
        if len(ms) > m_samples:
            step = len(ms) // m_samples
            ms = ms[::step][:m_samples]
        for m in ms:
            res = greedy(k, a, m, q, hi)
            if res is None:
                continue
            S, M = res
            if M <= floor:
                continue
            g = M - floor
            cur = dict(S=S, M=M, gk2=g * k * k, a=level_a(M, k), why=(a, q, m))
            if best is None or cur['gk2'] < best['gk2']:
                best = cur
    return best


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [31]
    for k in ks:
        b = best_for_k(k)
        floor = Fraction(1, k + 1)
        if b:
            print(f"k={k} sqrt(k)={k**0.5:.3f}: M={b['M']} (~{float(b['M']):.10f}) "
                  f"g*k^2={float(b['gk2']):.6f} a={float(b['a']):.4f} (a,q,m)={b['why']}")
            print(f"   S={b['S']}", flush=True)
        else:
            print(f"k={k}: none", flush=True)
