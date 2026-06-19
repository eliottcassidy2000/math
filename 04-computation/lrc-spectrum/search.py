"""
Search for k-speed gcd-1 sets with M(S) close to floor 1/(k+1), beating primorial.

Strategy ideas from the prompt:
 - multiple large "killer" speeds (several multiples of (k-1) or lcm(1..m))
 - speeds = lcm(1..m)*j patterns
 - dilations
 - perturbations of the AP {1..k}

Background fact: The AP {1,...,k} achieves M = 1/(k+1) exactly at t=1/(k+1):
  v*t = v/(k+1), ||.|| min over v=1..k. v=1 gives 1/(k+1); v=k gives k/(k+1)->1/(k+1).
So the AP is tight. To get a NON-tight set close to floor, we perturb a few speeds.

Key construction (known "near-AP" extremals): take {1,2,...,k} and replace some
speeds by larger ones that are "almost invisible" at t ~ 1/(k+1) because they are
multiples of (k+1) shifted, keeping the min near 1/(k+1).

We'll search systematically.
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
    S = sorted(set(S))
    if k is None:
        k = len(S)
    if len(S) != k:
        return None
    if setgcd(S) != 1:
        return None
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    if M <= floor:
        # tight or below floor (below shouldn't happen for valid sets per LRC)
        is_tight = (M == floor)
        return dict(S=S, M=M, t=t, floor=floor, g=M - floor, tight=is_tight,
                    gk2=(M - floor) * k * k, a=level_a(M, k))
    return dict(S=S, M=M, t=t, floor=floor, g=M - floor, tight=False,
                gk2=(M - floor) * k * k, a=level_a(M, k))


def fmt(r):
    if r is None:
        return "INVALID"
    return (f"M={r['M']} (~{float(r['M']):.9f}) g={r['g']} g*k^2={float(r['gk2']):.6f} "
            f"a={float(r['a']) if r['a'] else None} tight={r['tight']}")


if __name__ == "__main__":
    # 1) AP perturbations at moderate k: replace top element k by k+something
    for k in [29, 30, 31, 37, 41, 43, 47, 53, 59]:
        base = list(range(1, k + 1))
        best = None
        # replace one element with a "killer" multiple of (k+1) +/- small, or large value
        cand_results = []
        # try replacing the largest few entries with multiples of (k+1)
        for repl_idx in range(k - 1, max(k - 6, 0) - 1, -1):
            for mult in range(1, 6):
                for off in range(-3, 4):
                    newv = mult * (k + 1) + off
                    if newv <= 0:
                        continue
                    S = base[:repl_idx] + base[repl_idx + 1:] + [newv]
                    S = sorted(set(S))
                    if len(S) != k:
                        continue
                    r = evalset(S, k)
                    if r and not r['tight'] and r['M'] > r['floor']:
                        cand_results.append(r)
        if cand_results:
            cand_results.sort(key=lambda r: r['M'])
            print(f"k={k}: best non-tight via single-killer replacement:")
            print("   S=", cand_results[0]['S'])
            print("   ", fmt(cand_results[0]))
        else:
            print(f"k={k}: no non-tight found in this family")
