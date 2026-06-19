"""
Fast screening: for chosen t*=m/q and level r=j/q, the constructed set S of the k
smallest safe speeds has min_v ||v t*|| >= r BY CONSTRUCTION, so M(S) >= r.
We use r as a CHEAP lower-proxy for M(S) (M(S) is the true value, >= r). We rank
candidates by r (smaller r => closer to floor => larger proxy-level), then run the
EXACT M_exact only on the most promising few to report rigorous M, g*k^2, a.

This decouples the (cheap) combinatorial screen from the (expensive) exact eval.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from lrc_maxmin import M_exact
from primorial_family import level_a


def frac_norm_q(x, q):
    # ||x/q|| as a Fraction, x integer
    r = x % q
    r = min(r, q - r)
    return r  # this is q*||x/q||, an integer in [0, q//2]


def safe_speeds(m, q, jr, k, cap):
    """k smallest v>=1 with q*||v m/q|| >= jr  (i.e. ||v m/q|| >= jr/q)."""
    out = []
    v = 1
    mm = m % q
    cur = 0
    while len(out) < k and v <= cap:
        cur = (cur + mm) % q
        d = min(cur, q - cur)
        if d >= jr:
            out.append(v)
        v += 1
    if len(out) < k:
        return None
    return out


def setgcd(S):
    return reduce(gcd, S)


def screen_k(k, q_lo, q_hi, top=40):
    """Return list of candidate sets ranked by proxy level (smallest r first)."""
    floor = Fraction(1, k + 1)
    cands = []
    seen = set()
    for q in range(q_lo, q_hi + 1):
        jmin = q // (k + 1) + 1
        for jr in range(jmin, jmin + 3):
            if Fraction(jr, q) <= floor:
                continue
            for m in range(1, q):
                if gcd(m, q) != 1:
                    continue
                S = safe_speeds(m, q, jr, k, cap=6 * q)
                if S is None:
                    continue
                if setgcd(S) != 1:
                    continue
                key = tuple(S)
                if key in seen:
                    continue
                seen.add(key)
                r = Fraction(jr, q)
                # proxy g*k^2 lower bound
                proxy_gk2 = float((r - floor) * k * k)
                cands.append((proxy_gk2, r, S))
    cands.sort(key=lambda x: x[0])
    return cands[:top]


def exact_eval(S, k):
    S = sorted(set(S))
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, g=g, gk2=g * k * k,
                a=(level_a(M, k) if g > 0 else None), tight=(g == 0), maxS=max(S))


def best_for_k(k, amax=None, top=100000):
    """Evaluate EXACT M on ALL screened candidates (top=large), keep min g*k^2."""
    if amax is None:
        amax = int(k**0.5) + 4
    cands = screen_k(k, max(k, k + 1), amax * (k + 1), top=top)
    best = None
    for proxy, r, S in cands:
        rr = exact_eval(S, k)
        if rr['tight'] or rr['g'] <= 0:
            continue
        if best is None or rr['gk2'] < best['gk2']:
            best = rr
    return best, cands


if __name__ == "__main__":
    import json
    ks = [int(x) for x in sys.argv[1:]] or [29, 30, 31, 36, 40, 48, 60]
    for k in ks:
        best, cands = best_for_k(k)
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} (screened {len(cands)} cands) ====")
        if best:
            print(f"  S={best['S']}")
            print(f"  M={best['M']} (~{float(best['M']):.10f}) g*k^2={float(best['gk2']):.6f} "
                  f"a={float(best['a']):.4f} maxS={best['maxS']}")
        else:
            print("  none non-tight")
        sys.stdout.flush()
