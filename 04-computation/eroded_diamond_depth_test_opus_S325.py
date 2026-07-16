# opus-2026-07-16-S325 -- HYP-7025: (A) THE DEFINITIVE ERODED-DIAMOND DEPTH
# TEST -- exact THM-789 folded objects H_{a,b} = {t : ||at|| + ||bt|| >= 11/13}
# for the twins (a,b) = (11,2) [(x,y)=(13,9)] vs (15,2) [(17,13)], plus the
# U0-eroded versions E = H cap G_deep(U0); scalar calibration; alpha to depth 25.
# (B) the relation-second-moment probe: Kendall-Wei/Perron lambda per class
# vs H and x at n = 5, 6.
from fractions import Fraction
from math import pi, gcd
import cmath, sys
sys.path.insert(0, '04-computation')

def norm_dist(fr):   # ||fr|| for Fraction in [0,1)
    fr = fr % 1
    return min(fr, 1 - fr)

def folded_set(a, b, thresh):
    # exact arcs of {t in [0,1): ||a t|| + ||b t|| >= thresh}
    # kinks of ||a t||: t = j/(2a); of ||b t||: j/(2b)
    ks = sorted(set([Fraction(j, 2*a) for j in range(2*a + 1)] +
                    [Fraction(j, 2*b) for j in range(2*b + 1)]))
    arcs = []
    for i in range(len(ks) - 1):
        lo, hi = ks[i], ks[i+1]
        # Q is linear on [lo, hi]: Q(t) = ||at|| + ||bt||
        qlo = norm_dist(a*lo) + norm_dist(b*lo)
        qhi = norm_dist(a*hi) + norm_dist(b*hi)
        if qlo >= thresh and qhi >= thresh:
            arcs.append((lo, hi))
        elif qlo < thresh and qhi < thresh:
            continue
        else:
            # crossing: solve linear
            tt = lo + (thresh - qlo)/(qhi - qlo) * (hi - lo)
            if qlo >= thresh: arcs.append((lo, tt))
            else: arcs.append((tt, hi))
    # merge adjacent
    merged = []
    for (lo, hi) in arcs:
        if merged and lo == merged[-1][1]: merged[-1] = (merged[-1][0], hi)
        else: merged.append((lo, hi))
    return merged

def intersect(A, B):
    out = []
    for (a1, b1) in A:
        for (a2, b2) in B:
            lo, hi = max(a1, a2), min(b1, b2)
            if lo < hi: out.append((lo, hi))
    return sorted(out)

def deep_set(U, depth):
    # {t : min_u ||u t|| >= depth}: the core's deep times
    ivs = [(Fraction(0), Fraction(1))]
    for v in U:
        bands = []
        for k in range(v):
            lo = (Fraction(k) + depth)/v
            hi = (Fraction(k+1) - depth)/v
            if lo < hi: bands.append((lo, hi))
        ivs = intersect(ivs, bands)
    return ivs

def moments_arcs(arcs, K):
    tot = float(sum(b - a for a, b in arcs))
    ms = {0: 1.0+0j}
    for k in range(1, K + 1):
        z = 0j
        for (a, b) in arcs:
            z += (cmath.exp(-2j*pi*k*float(b)) - cmath.exp(-2j*pi*k*float(a)))/(-2j*pi*k)
        ms[k] = z / tot
    return ms

def opuc_float(ms, N):
    polys, alphas = [], []
    def c(k): return ms[abs(k)] if k >= 0 else ms[abs(k)].conjugate()
    for n in range(N + 1):
        coef = [0j]*(n+1); coef[n] = 1.0+0j
        for m in range(n):
            pm = polys[m]
            num = sum(pm[j].conjugate() * c(n - j) for j in range(len(pm)))
            den = sum(pm[i].conjugate() * pm[j] * c(i - j)
                      for i in range(len(pm)) for j in range(len(pm)))
            if abs(den) < 1e-13: return alphas
            r = num / den
            for j in range(len(pm)): coef[j] -= r * pm[j]
        polys.append(coef)
        if n >= 1: alphas.append(-coef[0].conjugate())
    return alphas

TH = Fraction(11, 13)
U0 = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12]
print("(A) THE DEFINITIVE ERODED-DIAMOND DEPTH TEST")
objs = {}
for name, (a, b) in (('TWIN-A (11,2)', (11, 2)), ('TWIN-B (15,2)', (15, 2))):
    H = folded_set(a, b, TH)
    muH = sum(hi - lo for lo, hi in H)
    print(f"   {name}: H has {len(H)} components, mu(H) = {muH} = {float(muH):.6f}")
    objs[name] = H
eq = objs['TWIN-A (11,2)'] and (sum(b-a for a,b in objs['TWIN-A (11,2)']) ==
                                 sum(b-a for a,b in objs['TWIN-B (15,2)']))
print(f"   scalar measures EQUAL (canon calibration): {eq}")

# eroded versions: intersect with U0's 1/11-deep set (per THM-789's depth)
D = deep_set(U0, Fraction(1, 11))
print(f"   U0 deep set (depth 1/11): {len(D)} components, mu = "
      f"{float(sum(b-a for a,b in D)):.6f}")
for name in list(objs):
    E = intersect(objs[name], D)
    objs['ERODED ' + name] = E
    print(f"   ERODED {name}: {len(E)} components, mu = "
          f"{float(sum(b-a for a,b in E)):.6f}")

print("\n   alpha-vectors (depth 25) and per-depth twin separations:")
als = {}
for name, arcs in objs.items():
    if not arcs: continue
    ms = moments_arcs(arcs, 26)
    als[name] = opuc_float(ms, 25)
for pair, (k1, k2) in (('RAW H', ('TWIN-A (11,2)', 'TWIN-B (15,2)')),
                        ('ERODED', ('ERODED TWIN-A (11,2)', 'ERODED TWIN-B (15,2)'))):
    if k1 not in als or k2 not in als: continue
    A, B = als[k1], als[k2]
    L = min(len(A), len(B))
    d = [abs(A[i] - B[i]) for i in range(L)]
    print(f"   {pair}: max sep depth 0-9: {max(d[:10]):.4f}; "
          f"10-17: {max(d[10:18]) if L > 10 else 0:.4f}; "
          f"18-{L-1}: {max(d[18:]) if L > 18 else 0:.4f}")

print("\n(B) THE RELATION SECOND MOMENT: Kendall-Wei/Perron per class vs H, x")
from smith_diagram_of_the_metagraph_opus_S307 import build
import numpy as np
from collections import defaultdict
for n in (5, 6):
    B5 = build(n)
    cls_of, H_of, x_of = B5['cls_of'], B5['H_of'], B5['x_of']
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    rep = {}
    for t in range(1 << len(tiles)):
        c = cls_of[t]
        if c not in rep: rep[c] = t
    rows = []
    for c, t in rep.items():
        A = np.zeros((n, n))
        for k in range(2, n+1): A[k-1][k-2] = 1
        for i, (xx, yy) in enumerate(tiles):
            if (t >> i) & 1: A[xx-1][yy-1] = 1
            else: A[yy-1][xx-1] = 1
        lam = max(np.linalg.eigvals(A).real)
        rows.append((c, lam, H_of[c], x_of[c]))
    lams = np.array([r[1] for r in rows]); Hs = np.array([float(r[2]) for r in rows])
    xs = np.array([float(r[3]) for r in rows])
    print(f"   n={n}: corr(lambda, H) = {np.corrcoef(lams, Hs)[0,1]:+.4f}; "
          f"corr(lambda, -x) = {np.corrcoef(lams, -xs)[0,1]:+.4f}; "
          f"lambda range [{lams.min():.4f}, {lams.max():.4f}]")
    # is lambda an x-refinement? per-level spread
    lev = defaultdict(list)
    for c, lam, H, x in rows: lev[x].append(lam)
    sp = {x: (min(v), max(v)) for x, v in lev.items() if len(v) > 1}
    print(f"      per-level lambda spreads (x: min..max): "
          f"{ {k: (round(a,3), round(b,3)) for k,(a,b) in sorted(sp.items())} }")
