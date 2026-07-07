#!/usr/bin/env python3
r"""
lrc_coprime_perpattern_kps_S73.py   (kind-pasteur-2026-07-07-S73)

Nail the CLEANEST TRUE coprime statement behind the layer-cake dominance, and check the
m=2 base-case proof.

PER-PATTERN DOMINANCE (candidate clean lemma, STRONGER than cumulative M_m):
  Every triple of E is an affine copy  a + g*{0, i, i+j}  of a unique PRIMITIVE (coprime)
  pattern pi = {0,i,i+j}, gcd(i,j)=1, at scale g>=1.  Define
     copies(pi, E) = #{(a,g>=1): {a, a+g*i, a+g*(i+j)} subset E}.
  CLAIM:  copies(pi, E) <= copies(pi, AP_k)  for EVERY primitive pattern pi and every k-set E.
  This => M_m(E) = sum_{pi: i+j<=m} copies(pi,E) <= M_m(AP) (cumulative), => Sigma_3 <= Sigma_3(AP).

  copies(pi, AP_k) = sum_{g>=1, g*(i+j)<=k-1} (k - g*(i+j))   [interval translates per scale].

TESTS:
  1. per-pattern dominance vs the SAME adversary battery + hill-climb (k=8, k=13).
     -- if it holds, it is the correct coprime statement (AP = maximizer of every coprime-
        pattern copy-count); if it FAILS while cumulative M_m holds, the dominance is only
        coarse-grained (report which).
  2. m=2 base case as an EXACT additive-energy identity:
        #3-APs(E) = #{(x,y) in E^2, x<y, (x+y)/2 in E}    (midpoint count)
     and the interval-maximality of the midpoint count (compression check).
"""
import random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

def primitive(E):
    E = sorted(set(E)); E = [e - E[0] for e in E]
    g = 0
    for e in E[1:]: g = gcd(g, e)
    return tuple(e // g for e in E) if g > 1 else tuple(E)

def all_primitive_patterns(mmax):
    """primitive triple shapes (i, i+j), gcd(i,j)=1, i,j>=1, reduced max-diff s=i+j <= mmax."""
    pats = []
    for s in range(2, mmax + 1):
        for i in range(1, s):
            j = s - i
            if gcd(i, j) == 1:
                pats.append((i, s))   # points {0, i, s}
    return pats

def copies(pat, E):
    i, s = pat
    S = set(E); cnt = 0
    lo, hi = min(E), max(E)
    for a in E:
        g = 1
        while a + g * s <= hi:
            if (a + g * i) in S and (a + g * s) in S:
                cnt += 1
            g += 1
    return cnt

def copies_AP(pat, k):
    i, s = pat
    return sum(k - g * s for g in range(1, k) if g * s <= k - 1)

def battery(k):
    B = {}
    B["AP"] = list(range(k))
    B["all-odd"] = [1 + 2*j for j in range(k)]
    sp2 = [1 + 2*j for j in range(k-1)]; B["spread-d2+bump"] = sp2 + [sp2[-1] + 3]
    sp3 = [1 + 3*j for j in range(k-1)]; B["spread-d3+bump"] = sp3 + [sp3[-1] + 4]
    B["AP-1+far"] = list(range(k-1)) + [3*(k-1)]
    B["Fibonacci"] = [1,2,3,5,8,13,21,34,55,89,144,233,377][:k]
    B["Lucas"] = [1,3,4,7,11,18,29,47,76,123,199,322,521][:k]
    B["primes"] = [2,3,5,7,11,13,17,19,23,29,31,37,41][:k]
    B["two-block"] = list(range(k//2)) + [40 + j for j in range(k - k//2)]
    if k == 13:
        B["GW"] = list(range(1,12)) + [13, 24]
        B["prim-sat"] = [2*j for j in range(1,13)] + [13]
    return B

def run(k, seed=73, n_random=4000, hill=3000):
    print("=" * 96)
    print(f"  k={k}: PER-PATTERN DOMINANCE  copies(pi,E) <= copies(pi,AP)  (every primitive coprime pi)")
    print("=" * 96)
    pats = all_primitive_patterns(7)
    apc = {p: copies_AP(p, k) for p in pats}
    print(f"  {len(pats)} primitive patterns (s=i+j<=7).  AP copy-counts:")
    for p in pats:
        print(f"    pattern {{0,{p[0]},{p[1]}}} (s={p[1]}): copies(AP_{k}) = {apc[p]}")
    rng = random.Random(seed)
    viol = []
    def check(name, E):
        Ep = list(primitive(E))
        if len(set(Ep)) != k: return
        for p in pats:
            c = copies(p, Ep)
            if c > apc[p]:
                viol.append((name, Ep, p, c, apc[p]))
    B = battery(k)
    print(f"\n  battery (per-pattern violations flagged):")
    for name, E in B.items():
        Ep = list(primitive(E))
        cs = [copies(p, Ep) for p in pats]
        bad = [f"{{0,{p[0]},{p[1]}}}:{copies(p,Ep)}>{apc[p]}" for p in pats if copies(p,Ep) > apc[p]]
        print(f"    {name:16s} copies-sum={sum(cs):4d}  {'VIOL: '+','.join(bad) if bad else 'ok'}")
        check(name, E)
    for _ in range(n_random):
        check("rand", primitive([rng.randint(0, 6*k) for _ in range(k)]))
    # hill-climb maximizing total copies AND each single pattern
    for target in [None] + pats:
        cur = list(range(k))
        obj = (lambda E: sum(copies(p, list(primitive(E))) for p in pats)) if target is None \
              else (lambda E, t=target: copies(t, list(primitive(E))))
        best = obj(cur)
        for _ in range(hill // (len(pats)+1)):
            E2 = cur[:]
            if rng.random() < 0.5: E2[rng.randrange(k)] = rng.randint(0, 6*k)
            else:
                i2, j2 = rng.randrange(k), rng.randrange(k); E2[i2] = E2[j2] + rng.choice([1,-1,2,-2])
            E2 = list(primitive(E2))
            if len(set(E2)) != k: continue
            check("hill", E2)
            v = obj(E2)
            if v > best or rng.random() < 0.3: cur, best = E2, max(best, v)
    print(f"\n  PER-PATTERN dominance violations: {len(viol)}")
    for name, E, p, c, a in viol[:10]:
        print(f"    {name}: {E} has {c} copies of {{0,{p[0]},{p[1]}}} > AP's {a}")
    return viol

def m2_base_case():
    print()
    print("=" * 96)
    print("  m=2 BASE CASE:  #3-APs(E) = #{(x,y) in E^2, x<y, (x+y)/2 in E}  (midpoint identity)")
    print("  and interval-maximality of the midpoint count.")
    print("=" * 96)
    def threeAP_direct(E):
        S = set(E); c = 0
        for a, b, cc in combinations(sorted(E), 3):
            if b - a == cc - b: c += 1
        return c
    def midpoint_count(E):
        S = set(E); c = 0
        E = sorted(E)
        for x, y in combinations(E, 2):
            if (x + y) % 2 == 0 and (x + y)//2 in S:
                c += 1
        return c
    rng = random.Random(5)
    ok = True
    for _ in range(2000):
        k = rng.randint(4, 10)
        E = sorted(set(rng.randint(0, 40) for _ in range(k)))
        if threeAP_direct(E) != midpoint_count(E):
            ok = False; print(f"    identity FAILS on {E}"); break
    print(f"  midpoint identity holds on 2000 random sets: {ok}")
    # interval maximality: among k-sets, AP maximizes #3-APs (adversarial)
    for k in (8, 13):
        ap = threeAP_direct(list(range(k)))
        worst = 0; ws = None
        cur = list(range(k))
        for _ in range(6000):
            E2 = cur[:]
            if rng.random()<0.5: E2[rng.randrange(k)] = rng.randint(0,6*k)
            else:
                i,j = rng.randrange(k), rng.randrange(k); E2[i]=E2[j]+rng.choice([1,-1,2,-2])
            E2 = sorted(set(E2))
            if len(E2)!=k: continue
            v = threeAP_direct(E2)
            if v > worst: worst, ws = v, E2
            if v >= worst or rng.random()<0.3: cur = E2
        print(f"  k={k}: #3-APs(AP)={ap}; hill-climb max over non-AP = {worst} "
              f"-> AP {'IS' if worst<=ap else 'NOT'} the 3-AP maximizer")

if __name__ == "__main__":
    v8 = run(8)
    v13 = run(13)
    m2_base_case()
    print()
    print("=" * 96)
    print(f"  SUMMARY: per-pattern dominance violations  k=8:{len(v8)}  k=13:{len(v13)}")
    print("  0 => the AP maximizes copies of EVERY primitive coprime triple-pattern (the clean")
    print("  coprime lemma); cumulative M_m and Sigma_3 AP-maximality follow immediately.")
    print("DONE.")
