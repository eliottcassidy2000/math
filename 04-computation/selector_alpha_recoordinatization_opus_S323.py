# opus-2026-07-16-S323 -- HYP-6995 / THM-885: THE SELECTOR RE-COORDINATIZATION.
# Two-sheet rows S = 2U u {x, y} through the THM-883 alpha-map.
# Named rows (canon): U0 = {1,2,3,5,7,8,9,10,11,12}:
#   (13,9)  THM-789's trapped liar (return set trapped at t0 = 4/17)
#   (17,13) the fingerprint twin (identical raw fingerprints, diff incidence)
#   (13,5)  THM-824's diagnostic pair
# Controls: perturbed cores and generic odd pairs (gate-failing/loose-ish).
# Statistics per row: exact M; #arcs; mu(G); alpha_0..alpha_9; HUG = max|alpha_n|
# for n >= 2; ALT = period-2 alternation score; verdict table.
from fractions import Fraction
from math import gcd, pi
import cmath, itertools

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
            if abs(den) < 1e-14: return alphas
            r = num / den
            for j in range(len(pm)): coef[j] -= r * pm[j]
        polys.append(coef)
        if n >= 1: alphas.append(-coef[0].conjugate())
    return alphas

def exact_M(S, qmax=400):
    best = (Fraction(0), None)
    for q in range(2, qmax+1):
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            t = Fraction(p, q)
            m = min(min((v*t) % 1, 1 - (v*t) % 1) for v in S)
            if m > best[0]: best = (m, t)
    return best

U0 = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12]
ROWS = [
    ('U0+(13,9) TRAPPED-LIAR', sorted([2*u for u in U0] + [13, 9])),
    ('U0+(17,13) FP-TWIN',     sorted([2*u for u in U0] + [17, 13])),
    ('U0+(13,5) THM824-PAIR',  sorted([2*u for u in U0] + [13, 5])),
    ('U0+(15,9) control-odd',  sorted([2*u for u in U0] + [15, 9])),
    ('Upert+(13,9) core-ctrl', sorted([2*u for u in [1,2,3,4,7,8,9,10,11,12]] + [13, 9])),
    ('Uap+(13,9) AP-core',     sorted([2*u for u in range(1, 11)] + [13, 9])),
]
print(f"{'row':26s} {'M':10s} {'arcs':4s} {'mu(G)':8s} {'HUG':6s} {'ALT':6s}  alpha (|.| shown, sign of Re)")
for name, S in ROWS:
    M, targ = exact_M(S)
    arcs = safe_arcs(S, Fraction(1, 14))
    if not arcs:
        print(f"{name:26s} {str(M):10s}   -- safe set EMPTY at 1/14 (covering row)")
        continue
    ms = moments_arcs(arcs, 11)
    al = opuc_float(ms, 10)
    mu = float(sum(b - a for a, b in arcs))
    hug = max(abs(a) for a in al[2:]) if len(al) > 2 else 0
    # period-2 alternation: correlation of sign(Re alpha_n) with (-1)^n over tail
    signs = [1 if a.real >= 0 else -1 for a in al[1:]]
    alt = abs(sum(s*(-1)**i for i, s in enumerate(signs)))/len(signs) if signs else 0
    astr = ' '.join(f"{'+' if a.real>=0 else '-'}{abs(a):.2f}" for a in al)
    print(f"{name:26s} {str(M):10s} {len(arcs):4d} {mu:8.5f} {hug:6.3f} {alt:6.2f}  {astr}")
print()
print("hypothesis: deep/gate-surviving rows -> HUG ~ 1 with high ALT (period-2");
print("boundary-hugging, the two-band signature); loose/escaping rows -> spread.")
