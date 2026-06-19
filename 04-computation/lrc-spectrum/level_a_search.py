"""
Targeted search for level-a family: q = a(k+1)-1, t* = m/q, M = a/q.
The optimal set at level a/q: pick t*=m/q, then S = a chosen k-subset of the safe
speeds {v : ||v m/q|| >= a/q} such that NO OTHER t beats a/q. We don't know which
subset a priori, so we search: take the safe residues mod q, and among integers
realizing them pick combinations. But exhaustive subset is huge.

Better: empirically, sigma_2 minimizers look like "small dense block with a few gaps
+ one or more killers near q". Let's parametrize:
  S = (prefix-with-gaps of size k-r) UNION (r killers), all within [1, q],
each speed having safe residue mod q (||.||>=a/q at t*=m/q). We then verify M exactly.

We brute force over: a in [2..], q=a(k+1)-1, m coprime to q with 1<=m<q (m and q-m
equivalent so m in [1,q//2]). For each (a,q,m): compute safe residues; the safe speeds
in [1,q] = those integers; if exactly... we need k of them. Take the k SMALLEST safe
speeds in [1, q] (period). Evaluate. Also take k smallest safe in [1, 2q], etc.
KEY DIFFERENCE from before: we now also require the construction to be checked at the
RIGHT q (=a(k+1)-1), and we look at whether a>1 is actually realized (not collapsing).
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


def safe_speeds_in_range(m, q, a, hi):
    """all v in [1,hi] with min((v*m)%q, q-(v*m)%q) >= a, in increasing order."""
    out = []
    for v in range(1, hi + 1):
        r = (v * m) % q
        gg = r if r <= q - r else q - r
        if gg >= a:
            out.append(v)
    return out


def best_level_a(k, amax=None, hi_mult=3, verbose=False):
    if amax is None:
        amax = int(k**1.0)  # allow large a target; we'll see what materializes
    floor = Fraction(1, k + 1)
    best = None
    for a in range(2, amax + 1):
        q = a * (k + 1) - 1
        hi = hi_mult * q
        for m in range(1, q // 2 + 1):
            if gcd(m, q) != 1:
                continue
            safe = safe_speeds_in_range(m, q, a, hi)
            if len(safe) < k:
                continue
            # take k smallest safe speeds
            S = safe[:k]
            r = ev(S, k)
            if r and not r['tight'] and r['g'] > 0:
                if best is None or r['gk2'] < best['gk2']:
                    best = r
                    best['why'] = (a, q, m)
                    if verbose:
                        print(f"   a={a} q={q} m={m}: gk2={float(r['gk2']):.5f} a_act={float(r['a']):.3f} M={r['M']} S={S}", flush=True)
    return best


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [7, 8, 12, 20, 30, 31]
    for k in ks:
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} ====", flush=True)
        b = best_level_a(k, verbose=True)
        if b:
            print("  BEST:", b['S'])
            print(f"   M={b['M']} (~{float(b['M']):.10f}) g*k^2={float(b['gk2']):.6f} "
                  f"a={float(b['a']):.4f} (a,q,m)={b['why']} maxS={b['maxS']}", flush=True)
        else:
            print("  none", flush=True)
