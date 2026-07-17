# opus-2026-07-16-S334 -- HYP-7170 / THM-936: THE GAP-FREE CORE TAXONOMY.
# The last unstructured regime: 13-speed primitive packets with NO internal
# scale gap (all consecutive ratios <= 15 -- no cascade step, no gluing
# junction applies).  TRUE LRC(14) normalization: lam = 1/14.
#
# Per packet (all EXACT):
#   * uncovered set W = circle minus union of 13 combs; mu, kappa;
#   * TIGHT detection: mu = 0 => check the 14-grid t = p/14 for closed
#     loneliness (all ||x t|| >= 1/14) -- the mu6 locus;
#   * WITNESS: simplest rational (Stern-Brocot) in the largest component ->
#     witness denominator q_w (the Lean-certifiable object: 13 rational
#     distance checks);
#   * resonance profile: Theta* = min_{q,p<=13} |qb-pa|(q+p-1) over pairs;
#     #3-term APs; gcd structure.
# CENSUS: strata sizes, witness-denominator law vs family and vs Theta*,
# and the headline: is EVERY sampled gap-free packet lonely (mu > 0 or
# tight-attained)?  [LRC(14) predicts yes; each instance gets a kernel-
# checkable certificate.]
from fractions import Fraction
from math import floor, gcd
import random, itertools
from functools import reduce

F = Fraction
LAM = F(1, 14)

def subtract_comb(V, x):
    w = F(1, 14 * x)
    out = []
    for (a, b) in V:
        cur = a
        j0 = floor((a - w) * x); j1 = floor((b + w) * x) + 1
        for j in range(j0, j1 + 1):
            lo, hi = F(j, x) - w, F(j, x) + w
            if hi <= cur: continue
            if lo >= b: break
            if lo > cur: out.append((cur, lo))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return out

def uncovered(speeds):
    V = [(F(0), F(1))]
    for x in sorted(speeds): V = subtract_comb(V, x)
    return V

def mu(V): return sum(b - a for a, b in V)

def simplest_in(a, b):
    """simplest rational strictly inside (a, b), 0 <= a < b: classic
    continued-fraction descent, O(log denominator)."""
    fa = a.numerator // a.denominator          # floor(a); a < fa+1 always
    if F(fa + 1) < b:                          # integer strictly inside
        return F(fa + 1)
    if a == fa:                                # left endpoint an integer
        inv = F(1) / (b - fa)                  # want smallest q > 1/(b-fa)
        q = inv.numerator // inv.denominator + 1
        return fa + F(1, q)
    sub = simplest_in(F(1) / (b - fa), F(1) / (a - fa))
    return fa + F(1) / sub

def lonely_at(speeds, t):
    return all(min((x * t) % 1, 1 - (x * t) % 1) >= LAM for x in speeds)

def theta_star(speeds):
    best = None
    for a, b in itertools.combinations(speeds, 2):
        lo, hi = min(a, b), max(a, b)
        for q in range(1, 14):
            for p in range(1, 14):
                v = abs(q * hi - p * lo) * (q + p - 1)
                if best is None or v < best: best = v
    return best

def count_3ap(speeds):
    s = set(speeds)
    return sum(1 for a, b in itertools.combinations(speeds, 2)
               if a != b and (a + b) % 2 == 0 and (a + b) // 2 in s
               and (a + b) // 2 not in (a, b))

def gapfree(speeds):
    s = sorted(speeds)
    return all(s[i+1] <= 15 * s[i] for i in range(12)) and len(set(s)) == 13

def primitive(speeds):
    return reduce(gcd, speeds) == 1

# ---------------- sample generator
random.seed(334)
families = {}
families['tight'] = [list(range(1, 14))]
aps = []
for d in (1, 2, 3, 5, 7, 11):
    for a in (1, 2, 3, 5, 8, 14, 20, 50):
        P = [a + k * d for k in range(13)]
        if gapfree(P) and primitive(P): aps.append(P)
families['AP'] = aps[:20]
naps = []
for base_d in (2, 3, 5, 7):
    for a in (3, 5, 9, 15, 25):
        for pert in (1, -1, 2):
            P = [a + k * base_d for k in range(13)]
            P[6] += pert
            if len(set(P)) == 13 and gapfree(P) and primitive(P):
                naps.append(sorted(P))
families['nearAP'] = naps[:20]
geos = []
for qn, qd in ((4, 3), (3, 2), (5, 4), (7, 5)):
    for c in (7, 11, 20, 33, 50):
        P, v = [], c
        for _ in range(13):
            P.append(round(v)); v = v * qn / qd
        P = sorted(set(int(x) for x in P))
        if len(P) == 13 and gapfree(P) and primitive(P): geos.append(P)
families['geometric'] = geos[:16]
pool = []
tries = 0
while len(pool) < 200 and tries < 20000:
    tries += 1
    base = random.randint(50, 600)
    P = sorted(random.sample(range(base, min(base * 10, 5000)), 13))
    if gapfree(P) and primitive(P): pool.append(P)
pool_th = sorted(pool, key=theta_star)
families['random'] = pool_th[:120]
families['random_cleanest'] = pool_th[-80:]   # top-Theta* quartile
print("=" * 74)
print("THE GAP-FREE CORE TAXONOMY (lam = 1/14 TRUE LRC; all exact)")
print("=" * 74)
from collections import defaultdict
strata = defaultdict(int)
rows = []
not_lonely = []
for fam, packs in families.items():
    qws, mus = [], []
    for P in packs:
        V = uncovered(P)
        m = mu(V)
        th = theta_star(P)
        n3 = count_3ap(P)
        if m == 0:
            # tight or dead: check the 14-grid closed attainment
            att = [p for p in range(1, 14) if lonely_at(P, F(p, 14))]
            if att:
                strata['TIGHT (mu=0, 14-grid attained)'] += 1
                rows.append((fam, P if len(str(P)) < 55 else 'big', 0, 'p/14', th, n3))
            else:
                # search other rationals q <= 2002 quickly? honest: record
                not_lonely.append((fam, P))
                strata['NO-WITNESS-FOUND (mu=0, grid fails)'] += 1
        else:
            big = max(V, key=lambda iv: iv[1] - iv[0])
            wt = simplest_in(big[0], big[1])
            ok = lonely_at(P, wt)
            assert ok, (P, wt)
            strata['GENERIC (mu>0, rational witness)'] += 1
            qws.append(wt.denominator); mus.append(float(m))
    if qws:
        qws.sort()
        print(f"  {fam:15s}: n={len(packs):3d}  mu range [{min(mus):.4f}, "
              f"{max(mus):.4f}]  witness q: med {qws[len(qws)//2]:4d} "
              f"max {qws[-1]:5d}")
    else:
        print(f"  {fam:15s}: n={len(packs):3d}  (no mu>0 members)")
print()
print("STRATA CENSUS:")
for k, v in sorted(strata.items()): print(f"  {k}: {v}")
print()
if not_lonely:
    print("  !! packets with mu=0 and no 14-grid witness (INVESTIGATE):")
    for fam, P in not_lonely[:6]: print(f"     [{fam}] {P}")
else:
    print("  EVERY sampled gap-free packet is LONELY (mu>0 witness or tight")
    print("  14-grid attainment) -- each carries a kernel-checkable rational")
    print("  certificate: 13 exact distance checks. The accessible gap-free")
    print("  core is WITNESS-CERTIFIABLE end to end.")
print()
print("WITNESS-DENOMINATOR LAW (all mu>0 packets, by Theta* bins):")
bins = defaultdict(list)
for fam, packs in families.items():
    for P in packs:
        V = uncovered(P)
        if mu(V) == 0: continue
        th = theta_star(P)
        big = max(V, key=lambda iv: iv[1] - iv[0])
        wt = simplest_in(big[0], big[1])
        key = 0 if th < 30 else (1 if th < 60 else (2 if th < 120 else 3))
        bins[key].append(wt.denominator)
lbl = {0: 'Theta*<30', 1: '30-59', 2: '60-119', 3: '>=120'}
for k in sorted(bins):
    v = sorted(bins[k])
    print(f"  {lbl[k]:10s}: n={len(v):3d}  witness q med {v[len(v)//2]:4d}  "
          f"p90 {v[int(len(v)*0.9)]:5d}  max {v[-1]:5d}")
