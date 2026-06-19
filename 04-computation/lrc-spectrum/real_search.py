"""
Real search for near-floor non-tight k-sets, using fast exact M.

Families:
 (DIL) dilation: S = d*{1..k} won't be gcd1; instead "stretched AP" with a unit added.
 (KILL) block + killers: dense AP {1..b} plus (k-b) killer speeds, killers chosen as
        multiples of (b+1) or of lcm to be "invisible" at the block's optimal t.
 (PRIM) primorial: t*=m/q with q=a(k+1)-1; S = the k smallest v with ||v m/q|| >= a/q
        but EXCLUDING the trivial dense prefix -- instead take v ranging so the set is
        the EXACT safe set on a full period (residues), realized by the k smallest.

The collapse we saw: dense prefix {1..k-1} always re-optimizes to t~1/k. So the
primorial set must AVOID having a long run 1,2,3,...; it should be the safe residues.
We build S = the k smallest positive integers whose residue mod q is safe, where
'safe' = residue in [a, q-a] (so ||v m/q||>=a/q with m=1, t*=1/q). With m=1 the safe
speeds are exactly integers v with (v mod q) in [a, q-a]. The k smallest such.
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


def ev(S, k):
    S = sorted(set(int(x) for x in S))
    if len(S) != k or min(S) <= 0 or setgcd(S) != 1:
        return None
    M, t = M_exact_fast(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, g=g, gk2=g * k * k, tight=(g == 0),
                a=(level_a(M, k) if g > 0 else None), maxS=max(S))


def safe_set_m1(q, a, k):
    """k smallest v>=1 with (v mod q) in [a, q-a]. (level a/q at t*=1/q, m=1)"""
    out = []
    v = 1
    while len(out) < k:
        rmod = v % q
        if a <= rmod <= q - a:
            out.append(v)
        v += 1
        if v > 50 * q:
            return None
    return out


def safe_set_m(q, m, a, k):
    out = []
    v = 1
    while len(out) < k:
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            out.append(v)
        v += 1
        if v > 50 * q:
            return None
    return out


def fmt(r):
    if r is None:
        return "INVALID/tight"
    return (f"M={r['M']} (~{float(r['M']):.10f}) g*k^2={float(r['gk2']):.6f} "
            f"a={float(r['a']):.4f} maxS={r['maxS']}")


def best_prim(k, amax=None):
    """Scan q = a(k+1)-1 ... and offsets; m=1 safe-residue construction + all coprime m."""
    if amax is None:
        amax = int(k**0.5) + 5
    floor = Fraction(1, k + 1)
    best = None
    # a is the target level; q near a(k+1)-1; try band width a around floor
    for a in range(2, amax + 1):
        for q in range(a * (k + 1) - 1 - 2, a * (k + 1) + 3):
            if q < k:
                continue
            for m in range(1, q):
                if gcd(m, q) != 1:
                    continue
                S = safe_set_m(q, m, a, k)
                if S is None:
                    continue
                r = ev(S, k)
                if r is None or r['tight'] or r['g'] <= 0:
                    continue
                if best is None or r['gk2'] < best['gk2']:
                    best = r
                    best['why'] = (a, q, m)
    return best


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [30, 31]
    for k in ks:
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} ====", flush=True)
        b = best_prim(k)
        if b:
            print("  BEST prim-style:", b['S'])
            print("   ", fmt(b), " (a,q,m)=", b['why'], flush=True)
        else:
            print("  none", flush=True)
