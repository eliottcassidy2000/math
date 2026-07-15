# opus-2026-07-15-S314 -- HYP-6925 / THM-864: beat-localization lemma battery.
# Planted relations: choose (q, p, y), set B = (p*A + s*y)/q (s = +-1), require
# integrality + coprimality + Y*-minimality of the planted relation.
# Measure err = |mu(D_A cap D_B cap E) - mu(E) rho(A,B)| against the candidate
# bound  err <= C * kappa_E * (p+q) / (169 y)   [the sub-orbit AP argument's
# shape; C to be pinned and then DOMINATED by the proven constant].
from fractions import Fraction
from math import gcd
import math

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def mu(ivs): return sum(b - a for (a, b) in ivs)

def comb_teeth_in(x, a, b):
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def pair_measure_on(E, x1, x2):
    tot = Fraction(0)
    for (a, b) in E:
        for (lo, hi) in comb_teeth_in(x1, a, b):
            for (l2, h2) in comb_teeth_in(x2, lo, hi): tot += h2 - l2
    return tot

def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        T += (min(2*a, s - 13*m)) * (1 if m == 0 else 2)
        m += 1
    return T

def rho(x1, x2):
    g = gcd(x1, x2)
    a, b = x1//g, x2//g
    if a > b: a, b = b, a
    return Fraction(T_of(a, b), 13*a*b)

def Ystar(x1, x2, Q=13):
    best = None
    for q in range(1, Q+1):
        p0 = round(q*x2/x1)
        for p in (p0-1, p0, p0+1):
            if p < 1: continue
            v = abs(q*x2 - p*x1)
            if best is None or v < best[0]: best = (v, q, p)
    return best

PREFIXES = [[1,2,3,4,5], [2,3,5,6,8], [8,9,10,11,12]]
print("planted-relation battery: err vs kappa*(p+q)/(169 y)")
print("prefix kap | q  p    y    A     B    | Y*(q,p)     err        ratio err*169y/(kap(p+q))")
worst = 0
for P in PREFIXES:
    E = safe_set(P)
    muE = mu(E)
    kap = len(E)
    for (q, p) in [(1,1), (1,2), (1,3), (2,3), (1,12), (3,4), (2,5)]:
        for y in (5, 11, 23, 47, 95):
            for A0 in (611, 1200):
                # find A near A0 with q | (p*A + y) and gcd conditions
                found = None
                for A in range(A0, A0 + 6*q + 1):
                    if (p*A + y) % q: continue
                    B = (p*A + y)//q
                    if B <= A: continue
                    if gcd(A, B) != 1: continue
                    v, qq, pp = Ystar(A, B)
                    if v != y or qq != q: continue   # planted relation must be minimal
                    found = (A, B); break
                if not found: continue
                A, B = found
                er = abs(pair_measure_on(E, A, B) - muE*rho(A, B))
                ratio = float(er) * 169 * y / (kap * (p + q))
                worst = max(worst, ratio)
                print(f"{str(P):16s} {kap:2d} | {q:2d} {p:2d} {y:4d} {A:5d} {B:6d} | "
                      f"({q},{p})@{y}   {float(er):.6f}   {ratio:7.4f}")
print(f"\nWORST ratio err*169y/(kappa(p+q)) = {worst:.4f}")
print(f"=> empirical C in err <= C kappa (p+q)/(169 y): C ~ {worst:.3f}")
