"""
Direct construction of near-floor sets from a chosen rational t* = m/q.

For a fixed t* = m/q in lowest terms, the set of "safe" speeds at level r is
   Safe(t*, r) = { v >= 1 : ||v * m/q|| >= r }.
If we take r = a/q (so M-level a target, q = a(k+1)-1 ideally) and collect the k
SMALLEST safe speeds, then by construction min_v ||v t*|| >= r, so M(S) >= r.
We then compute M(S) EXACTLY (it may be a bit larger, which is fine -- still close).

We want M(S) as small as possible while size = k and gcd = 1 and S non-tight.
The achievable level a = M/(M(k+1)-1).

We scan q and m (and the level r), build S, evaluate exactly.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from lrc_maxmin import M_exact
from primorial_family import level_a


def frac_norm(t):
    f = t - int(t)
    if f < 0:
        f += 1
    return f if f <= Fraction(1, 2) else 1 - f


def safe_speeds(m, q, r, k, cap):
    """k smallest v>=1 with ||v*m/q|| >= r."""
    out = []
    v = 1
    while len(out) < k and v <= cap:
        if frac_norm(Fraction(v * m, q)) >= r:
            out.append(v)
        v += 1
    if len(out) < k:
        return None
    return out


def setgcd(S):
    return reduce(gcd, S)


def evalset(S, k):
    S = sorted(set(int(x) for x in S))
    if len(S) != k or setgcd(S) != 1 or min(S) <= 0:
        return None
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, floor=floor, g=g, tight=(g == 0),
                gk2=g * k * k, a=(level_a(M, k) if g > 0 else None), maxS=max(S))


def search_k(k, q_lo, q_hi, verbose=False):
    best = None
    floor = Fraction(1, k + 1)
    for q in range(q_lo, q_hi + 1):
        # target level r slightly above floor: r = ceil(q/(k+1)+1)/q-ish.
        # We want r close to floor 1/(k+1). r = j/q for j just above q/(k+1).
        jmin = q // (k + 1) + 1  # smallest j with j/q > 1/(k+1)? check
        for j in range(jmin, jmin + 4):
            r = Fraction(j, q)
            if r <= floor:
                continue
            for m in range(1, q):
                if gcd(m, q) != 1:
                    continue
                S = safe_speeds(m, q, r, k, cap=4 * q)
                if S is None:
                    continue
                rr = evalset(S, k)
                if rr is None or rr['tight'] or rr['g'] <= 0:
                    continue
                if best is None or rr['gk2'] < best['gk2']:
                    best = rr
                    if verbose:
                        print(f"   q={q} m={m} j={j}: g*k^2={float(rr['gk2']):.5f} a={float(rr['a']):.4f} M={rr['M']}")
    return best


if __name__ == "__main__":
    targets = [29, 30, 31, 36, 40, 48, 60]
    for k in targets:
        # q ranges around a(k+1) for a up to ~ sqrt(k)+
        amax = int(k**0.5) + 3
        q_lo = k  # at least
        q_hi = amax * (k + 1)
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} scanning q in [{q_lo},{q_hi}] ====")
        best = search_k(k, q_lo, q_hi, verbose=False)
        if best:
            print(f"  BEST: S={best['S']}")
            print(f"        M={best['M']} (~{float(best['M']):.10f}) g*k^2={float(best['gk2']):.6f} "
                  f"a={float(best['a']):.4f} maxS={best['maxS']}")
        else:
            print("  none")
