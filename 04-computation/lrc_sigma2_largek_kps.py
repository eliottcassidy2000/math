#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sigma_2(k) at large k via low-perturbation families off the AP, and the
lower-bound test for the spectral gap g(k)=sigma_2(k)-1/(k+1).
kind-pasteur-2026-06-19-S9.

KEY LEMMA (denominator bound): M(S)=max_t min_v||vt|| is attained at t*=m/(v_i+-v_j)
for a binding pair, so its reduced denominator q divides (v_i+-v_j) => q <= 2 max(S).
Hence  g(k) = (p(k+1)-q)/(q(k+1)) >= 1/(q(k+1)) >= 1/(2 max(S) (k+1)).
=> g(k) dips below Theta(1/k^2)  <=>  the sigma_2 extremizer needs max(S)/k -> infinity.

So we track, for the minimal-M low-perturbation config at each k:
   M, (p,q), a (if M=a/(a(k+1)-1)), max(S), max(S)/k, q/k, g*k^2.
If max(S)/k stays bounded => g=Theta(1/k^2) (gap does NOT dip).  If it grows => dips.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def frac_norm(x):
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)

def maxmin(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(i, len(S)):
            for d in (S[i] + S[j], S[i] - S[j]):
                if d == 0: continue
                d = abs(d)
                for m in range(1, d): cand.add(Fraction(m, d))
    best = Fraction(0); bestt = Fraction(0)
    for t in cand:
        mv = min(frac_norm(v * t) for v in S)
        if mv > best: best, bestt = mv, t
    return best, bestt

def primitive(S): return reduce(gcd, S) == 1

def search_sigma2(k, W=None, depth=2):
    """Minimal M over: AP {1..k} with up to `depth` of the top speeds replaced by
       free values in (k, W].  Returns (M, S, tau)."""
    if W is None: W = 4 * k + 6
    floor = Fraction(1, k + 1)
    best = (Fraction(1, 1), None, None)
    base = list(range(1, k + 1))
    # depth-1: replace position i (any) with w
    for i in range(k):
        for w in range(k + 1, W + 1):
            S = base[:i] + base[i+1:] + [w]
            S = tuple(sorted(set(S)))
            if len(S) != k or not primitive(S): continue
            M, t = maxmin(S)
            if floor < M < best[0]:
                best = (M, S, t)
    if depth >= 2:
        # replace the top two AP speeds k-1,k with free w1<w2 (both > some threshold)
        for i in range(k):
            for j in range(i+1, k):
                rest = [base[x] for x in range(k) if x != i and x != j]
                for w1 in range(k + 1, W + 1):
                    for w2 in range(w1 + 1, W + 1):
                        S = tuple(sorted(rest + [w1, w2]))
                        if len(S) != k or not primitive(S): continue
                        M, t = maxmin(S)
                        if floor < M < best[0]:
                            best = (M, S, t)
        # this double loop is O(k^2 W^2) maxmin calls -> only feasible small; gate it
    return best

def search_sigma2_fast(k, W=None):
    """Faster: depth-1 over ALL replace positions, plus depth-2 only replacing the
       single top speed by two free speeds (apex-split), bounded."""
    if W is None: W = 3 * k + 6
    floor = Fraction(1, k + 1)
    best = (Fraction(1, 1), None, None)
    base = list(range(1, k + 1))
    # depth-1
    for i in range(k):
        for w in range(k + 1, W + 1):
            S = tuple(sorted(set(base[:i] + base[i+1:] + [w])))
            if len(S) != k or not primitive(S): continue
            M, t = maxmin(S)
            if floor < M < best[0]: best = (M, S, t)
    # depth-2 apex-split: drop top two {k-1,k}, add two free w1<w2
    rest = base[:-2]
    for w1 in range(k, W + 1):
        for w2 in range(w1 + 1, W + 1):
            S = tuple(sorted(rest + [w1, w2]))
            if len(S) != k or not primitive(S): continue
            M, t = maxmin(S)
            if floor < M < best[0]: best = (M, S, t)
    # depth-2 drop {one middle, apex}: only a few middle positions j (cheap, catches
    # the (1,2,3,4,5,7,18)-type configs that skip one interior speed)
    for j in (1, 2, k // 2):
        if not (0 <= j < k - 1): continue
        rest2 = [base[x] for x in range(k) if x != j and x != k-1]
        for w1 in range(k, W + 1):
            for w2 in range(w1 + 1, W + 1):
                S = tuple(sorted(rest2 + [w1, w2]))
                if len(S) != k or not primitive(S): continue
                M, t = maxmin(S)
                if floor < M < best[0]: best = (M, S, t)
    return best

def express(M, k):
    p, q = M.numerator, M.denominator
    e = p * (k + 1) - q
    a = p if e == 1 else None
    return a, e, p, q

if __name__ == "__main__":
    print("=== sigma_2(k) lower-bound test via low-perturbation search ===")
    print(f"  {'k':>2} {'min M':>10} {'g*k^2':>8} {'a':>4} {'e':>3} {'q':>5} {'q/k':>6} {'maxS':>5} {'maxS/k':>7}  witness")
    data = []
    for k in range(5, 17):
        M, S, t = search_sigma2_fast(k)
        if S is None:
            print(f"  {k:>2}  (none found below mediant)")
            continue
        floor = Fraction(1, k + 1)
        g = M - floor
        a, e, p, q = express(M, k)
        maxS = max(S)
        data.append((k, M, g, a, e, q, maxS))
        print(f"  {k:>2} {str(M):>10} {float(g)*k*k:>8.4f} {str(a):>4} {e:>3} {q:>5} {q/k:>6.2f} {maxS:>5} {maxS/k:>7.3f}  {S}")
    print("\n  Trend check: if maxS/k and q/k stay bounded => g(k)=Theta(1/k^2) (gap does NOT dip).")
    print("  a increasing slowly, e usually 1 => sigma_2 on the Stern-Brocot mediant path a/(a(k+1)-1).")
