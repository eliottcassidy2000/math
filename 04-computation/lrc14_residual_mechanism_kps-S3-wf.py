#!/usr/bin/env python3
"""
LRC(14) S3 residual -- finding the PROVABLE mechanism behind C(S) via v=max.

We test several candidate provable lower bounds on W(S\\{max}):

(A) "small-part dominant": is the widest arc of S\\{max} essentially the widest arc
    of the SMALL PART alone, minus cluster-tooth bites? I.e. does the small part
    (subset of {1..13}) carve out a wide arc that the (few, tightly-clustered) large
    teeth fail to fully cover?

(B) "cluster-collapse window": near a tau where Vmin*tau mod 1 = 1/2 (cluster center
    parked at the antipode), all cluster runners u=V+d have u*tau approx 1/2 + d*tau,
    spread d*tau<= s*tau. If s*tau < 6/7 (band fits a single gap of width 6/7 in the
    safe middle), the whole cluster + small part is safe. Where can this hold AND the
    small part is safe? This is the cluster-as-one-runner reduction. Test its width.

(C) The decisive reframing: M(S) >= 1/14 directly. Replace the cluster L by a SINGLE
    surrogate runner. Claim: M(small ∪ {Vmin}) and M(small ∪ {Vmax}) bracket; if the
    small part ∪ one cluster runner already has a wide arc that the rest of the cluster
    cannot bite below 1/(7max), done.

We compute, per S3 set:
  - W_small  = Wwidth(small part only)
  - W_one    = Wwidth(small ∪ {one cluster runner})   (min over cluster choices)
  - W_rest   = Wwidth(S\\{max})  (actual)
  - thr      = 1/(7*max)
and check which simple quantity provably lower-bounds W_rest > thr.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
from collections import Counter

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def gen_S3(seed=0, target=2000, vmax_cap=1500):
    rng = random.Random(seed); out = []; tries = 0
    while len(out) < target and tries < target * 400:
        tries += 1
        c = rng.choice([2, 2, 3, 3, 4, 5, 6]); nsmall = 13 - c
        if nsmall < 1: continue
        V = rng.choice([rng.randint(40, 80), rng.randint(80, 200),
                        rng.randint(200, 600), rng.randint(600, vmax_cap)])
        spread = rng.choice([14, 20, 28, 35, 42, 56, 70])
        cluster = set()
        while len(cluster) < c: cluster.add(V + rng.randint(0, spread))
        small = set(); pool = list(range(1, 14)); rng.shuffle(pool)
        for x in pool:
            if len(small) >= nsmall: break
            small.add(x)
        if len(small) < nsmall: continue
        S = sorted(small | cluster)
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    S3 = gen_S3(seed=3, target=400)
    print(f"n S3 = {len(S3)}")
    thr = F(1, 14)
    # Per-set quantities
    Wsmall_vs = []; Wone_vs = []; Wrest_vs = []
    one_fires = 0; small_only_fires = 0; cluster_count = Counter()
    fail_one = []
    for S in S3:
        Vmax = S[-1]; rest = [u for u in S if u != Vmax]
        small = [u for u in S if u <= 13]
        cluster = [u for u in S if u > 13]
        cluster_count[len(cluster)] += 1
        t = F(1, 7*Vmax)
        Wrest = Wwidth(rest)
        # (A) small part alone
        Wsmall = Wwidth(small)
        # (B) small ∪ one cluster runner (min over choices = worst)
        Wone = min(Wwidth(small + [u]) for u in cluster)
        Wsmall_vs.append(float(Wsmall * 7 * Vmax))
        Wone_vs.append(float(Wone * 7 * Vmax))
        Wrest_vs.append(float(Wrest * 7 * Vmax))
        if Wone > t: one_fires += 1
        else: fail_one.append((S, float(Wone*7*Vmax), float(Wrest*7*Vmax)))
        if Wsmall > t: small_only_fires += 1
    n = len(S3)
    def stats(a): return f"min={min(a):.3f} mean={sum(a)/len(a):.3f} max={max(a):.3f}"
    print(f"cluster sizes among S\\max: {dict(sorted(cluster_count.items()))}")
    print(f"W(small)*7*max         : {stats(Wsmall_vs)}   (>1 in {small_only_fires}/{n})")
    print(f"W(small+1cluster)*7*max: {stats(Wone_vs)}   (>1 in {one_fires}/{n})")
    print(f"W(S\\max)*7*max  ACTUAL : {stats(Wrest_vs)}")
    print(f"\n(B) 'small + worst single cluster runner' fires C on {one_fires}/{n}")
    if fail_one:
        print(f"  {len(fail_one)} sets where worst single-cluster surrogate < threshold:")
        for S, wone, wrest in fail_one[:8]:
            print(f"     S={S}  Wone*7max={wone:.3f}  Wrest*7max={wrest:.3f}")
