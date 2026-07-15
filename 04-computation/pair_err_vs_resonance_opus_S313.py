# opus-2026-07-15-S313 -- the E-restriction error vs resonance distance.
# err(x1,x2) = |mu(D_x1 cap D_x2 cap E) - muE * rho(x1,x2)|
# Y(x1,x2)  = min over 1 <= q,p <= 60 of |q*x2 - p*x1| (resonance distance)
# Conjectured sharp law: err <= c / Y  (ET gives ~ log Y / Y with Q0 ~ 2500).
from fractions import Fraction
from math import gcd
import math, random

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

def Ydist(x1, x2, Q=60):
    best = None
    for q in range(1, Q+1):
        # p nearest to q*x2/x1
        p0 = round(q*x2/x1)
        for p in (p0-1, p0, p0+1):
            if p < 1: continue
            v = abs(q*x2 - p*x1)
            if best is None or v < best: best = v
    return best

P = [1, 2, 3, 4, 5]
E = safe_set(P)
muE = mu(E)
random.seed(3)
rows = []
# spread of pair types: random, near-dilate (planted), consecutive
pairs = [(random.randint(300, 900), random.randint(300, 2000)) for _ in range(14)]
pairs += [(500, 6001), (500, 6013), (700, 8397), (611, 613), (400, 407)]
print(" x1    x2     g   Y(res)   err          err*Y")
for (x1, x2) in pairs:
    if x1 == x2: continue
    er = abs(pair_measure_on(E, x1, x2) - muE * rho(x1, x2))
    Y = Ydist(x1, x2)
    rows.append((x1, x2, Y, er))
    print(f"{x1:5d} {x2:6d} {gcd(x1,x2):3d} {Y:6d}   {float(er):.6f}   {float(er)*Y:8.3f}")
mx = max(float(er)*Y for (_, _, Y, er) in rows)
print(f"\nmax err*Y over sample: {mx:.3f}  (conjectured sharp law err <= c/Y, c ~ {mx:.2f})")
print(f"budget per tree edge for phi* = 17/546 with muE >= 3/13: "
      f"{float(Fraction(3,13)*Fraction(17,546)/6/2):.6f} => need Y >= {mx/float(Fraction(3,13)*Fraction(17,546)/6/2):.0f}")
