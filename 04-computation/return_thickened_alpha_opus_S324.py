# opus-2026-07-16-S324 -- HYP-7000: THE RETURN-THICKENED ALPHA.
# Thicken each row's safe set (delta = 1/14) by the THM-789 return interval
# R = (-1/858, 1/858) (Minkowski: expand arcs by 1/858, merge overlaps).
# The merge pattern = the return-thickened component incidence (where THM-789
# says the twins differ). Tests:
#   (1) merge census: #arcs before/after, merged-gap positions, mu growth;
#   (2) alpha to depth 12 on the thickened measures: twin separation?
#   (3) depth-20 tails on the RAW measures (codex's growth question).
from fractions import Fraction
from math import gcd, pi
import cmath

def safe_arcs(S, delta):
    ivs = [(Fraction(0), Fraction(1))]
    for v in S:
        bands = []
        for k in range(v):
            lo = (Fraction(k) + delta)/v
            hi = (Fraction(k+1) - delta)/v
            if lo < hi: bands.append((lo, hi))
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def thicken(arcs, r):
    # Minkowski sum with (-r, r) on the circle, then merge
    ex = sorted([((a - r) % 1, (a - r) % 1 + (b - a) + 2*r) for (a, b) in arcs])
    # unroll circle: treat as linear with wrap merge at the end
    merged = []
    for (a, b) in ex:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    # wrap: if last extends past 1 and overlaps first
    if len(merged) > 1 and merged[-1][1] % 1 >= merged[0][0] and merged[-1][1] > 1:
        merged[0] = (merged[-1][0] - 1, max(merged[0][1], merged[-1][1] - 1))
        merged.pop()
    return merged

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

U0 = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12]
R = Fraction(1, 858)
TW = [('TWIN-A (13,9)', sorted([2*u for u in U0] + [13, 9])),
      ('TWIN-B (17,13)', sorted([2*u for u in U0] + [17, 13]))]

print("(1) MERGE CENSUS (thicken by R = (-1/858, 1/858)):")
data = {}
for name, S in TW:
    raw = safe_arcs(S, Fraction(1, 14))
    th = thicken(raw, R)
    mu_r = float(sum(b - a for a, b in raw))
    mu_t = float(sum(float(b) - float(a) for a, b in th))
    merges = len(raw) - len(th)
    # gap spectrum: raw inter-arc gaps < 2/858 (these merged)
    gaps = sorted(float(raw[i+1][0] - raw[i][1]) for i in range(len(raw)-1))
    small = [g for g in gaps if g < float(2*R)]
    print(f"   {name}: raw arcs = {len(raw)}, thickened = {len(th)} "
          f"(merges = {merges}); mu {mu_r:.5f} -> {mu_t:.5f}; "
          f"merged-gap count = {len(small)}")
    print(f"      smallest 6 gaps: {[f'{g:.6f}' for g in gaps[:6]]} "
          f"(threshold 2/858 = {float(2*R):.6f})")
    data[name] = (raw, th)

print("\n(2) THE RETURN-THICKENED ALPHA (depth 12):")
als = {}
for name, (raw, th) in data.items():
    ms = moments_arcs(th, 13)
    al = opuc_float(ms, 12)
    als[name] = al
    print(f"   {name}: {' '.join(f'{a.real:+.3f}' for a in al)}")
a, b = als['TWIN-A (13,9)'], als['TWIN-B (17,13)']
L = min(len(a), len(b))
dist = max(abs(a[i] - b[i]) for i in range(L))
print(f"   MAX |alpha_A - alpha_B| (thickened, depth {L}): {dist:.4f}")

print("\n(3) RAW depth-20 tails (the growth question):")
for name, (raw, th) in data.items():
    ms = moments_arcs(raw, 21)
    al = opuc_float(ms, 20)
    als['raw ' + name] = al
    print(f"   {name}: tail 10..19: {' '.join(f'{x.real:+.3f}' for x in al[10:20])}")
ra, rb = als['raw TWIN-A (13,9)'], als['raw TWIN-B (17,13)']
L2 = min(len(ra), len(rb))
per_depth = [abs(ra[i] - rb[i]) for i in range(L2)]
print(f"   raw |diff| by depth: {' '.join(f'{d:.3f}' for d in per_depth)}")
print(f"   raw max separation depth 0..9: {max(per_depth[:10]):.4f}; "
      f"depth 10..{L2-1}: {max(per_depth[10:]):.4f}")
