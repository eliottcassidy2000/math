"""
Faster targeted constructions, exact M via crossing candidates.

Focus families:
 (S) "skip" sets: {1,...,k+1} \ {one element}  (k speeds, near-AP, NON-tight)
 (D) dilation: take AP {1..k} -> not gcd1 after scaling, so use known dilation trick
 (K) multi-killer: small dense AP block {1..b} plus killers = j*(b+1) for j=2.. to fill k
 (P) primorial-style: speeds related to k-1 primorial

We optimize for smallest g*k^2.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from lrc_maxmin import M_exact
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def evalset(S, k=None):
    S = sorted(set(int(x) for x in S))
    if k is None:
        k = len(S)
    if len(S) != k or setgcd(S) != 1 or min(S) <= 0:
        return None
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, floor=floor, g=g, tight=(g == 0),
                gk2=g * k * k, a=(level_a(M, k) if g > 0 else None))


def fmt(r):
    if r is None:
        return "INVALID"
    a = r['a']
    return (f"S={r['S']}\n      M={r['M']} (~{float(r['M']):.10f}) g*k^2={float(r['gk2']):.6f} "
            f"a={(f'{float(a):.4f}' if a else None)} tight={r['tight']} maxS={max(r['S'])}")


if __name__ == "__main__":
    for k in [29, 30, 31, 35, 40, 47, 59]:
        print(f"==== k={k}  sqrt(k)={k**0.5:.3f} ====")
        results = []
        # (S) skip sets {1..k+1}\{j}
        for j in range(1, k + 2):
            S = [x for x in range(1, k + 2) if x != j]
            if len(S) != k:
                continue
            r = evalset(S, k)
            if r and not r['tight'] and r['g'] > 0:
                results.append(("skip%d" % j, r))
        results.sort(key=lambda x: x[1]['gk2'])
        if results:
            tag, r = results[0]
            print(f"  best SKIP ({tag}):", fmt(r))
        else:
            print("  SKIP: none non-tight")
