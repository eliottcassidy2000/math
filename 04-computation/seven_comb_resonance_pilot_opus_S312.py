# opus-2026-07-15-S312 -- the seven-comb resonance pilot (HYP-6895 LRC leg).
# Setting: scale-one radius-7 chart. Prefix P (5 unreplaced speeds of [12]),
# strict-safe set E = {t : min_{q in P} ||qt|| > 1/13} (exact rational intervals).
# Combs D_x, x = r + 13h for the 7 replaced residues r.
# CLAIMS TESTED:
#  (L1: periodicity lemma) x * (|E cap D_x| - (2/13) mu(E)) is EXACTLY periodic
#      in h with period dividing Lam = lcm(denominators of E's endpoints)/13-ish
#      -- verify by exact computation across h and h + Lam.
#  (P2: pair defects) mu(D_i cap D_j cap E) ~ (2/13)^2 mu(E) with defect
#      shrinking in min(x_i, x_j) EXCEPT near resonances (|x_i - x_j| small,
#      or x_i +- x_j hitting E's denominators).
#  (P3: quasi-independence) uncovered mass mu(E \ union of 7 combs) close to
#      mu(E) * (11/13)^7 for generic lifts -- coverage (tightness) requires
#      large positive dependence.
from fractions import Fraction
from math import gcd, lcm
import itertools, random

DELTA = Fraction(1, 13)

def safe_set(P):
    # components of {t in (0,1): ||qt|| > 1/13 for all q}, exact
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
    # danger teeth of speed x intersecting (a,b): [j/x - 1/(13x), j/x + 1/(13x)]
    w = Fraction(1, 13*x)
    j0 = (a - w) * x
    j1 = (b + w) * x
    import math
    out = []
    for j in range(math.floor(j0), math.floor(j1) + 2):
        lo, hi = Fraction(j, x) - w, Fraction(j, x) + w
        lo2, hi2 = max(lo, a), min(hi, b)
        if lo2 < hi2: out.append((lo2, hi2))
    return out

def inter_measure(E, x):
    tot = Fraction(0)
    for (a, b) in E:
        for (lo, hi) in comb_teeth_in(x, a, b): tot += hi - lo
    return tot

def pair_measure(E, x1, x2):
    tot = Fraction(0)
    for (a, b) in E:
        for (lo, hi) in comb_teeth_in(x1, a, b):
            for (l2, h2) in comb_teeth_in(x2, lo, hi):
                tot += h2 - l2
    return tot

def subtract_comb(ivs, x):
    out = []
    for (a, b) in ivs:
        cuts = comb_teeth_in(x, a, b)
        cur = a
        for (lo, hi) in sorted(cuts):
            if lo > cur: out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return [iv for iv in out if iv[0] < iv[1]]

P = [1, 2, 3, 4, 5]
R = [6, 7, 8, 9, 10, 11, 12]
E = safe_set(P)
muE = mu(E)
dens = [d for (a, b) in E for d in (a.denominator, b.denominator)]
Lam = 1
for d in dens: Lam = lcm(Lam, d)
print(f"prefix P={P}: E has {len(E)} components, mu(E) = {muE} ~ {float(muE):.5f}")
print(f"endpoint denominator lcm = {Lam}")

# L1: periodicity of the anomaly
r = 6
print("\nL1 anomaly x*(mass - (2/13)muE) for x = r+13h, r=6:")
base = {}
for h in range(1, 9):
    x = r + 13*h
    an = (inter_measure(E, x) - Fraction(2,13)*muE) * x
    base[h] = an
    print(f"   h={h:3d} x={x:4d}: anomaly*x = {an} ~ {float(an):+.6f}")
# periodic check: h vs h + Lam/13? test candidate periods
for T in (5, 10, 20, 60):
    ok = all((inter_measure(E, r+13*(h+T)) - Fraction(2,13)*muE)*(r+13*(h+T)) == base[h]
             for h in range(1, 5))
    print(f"   period T={T}: {'EXACT' if ok else 'no'}")

# P2: pair defects (fix x1 = 6+13*4 = 58; scan x2 = 7+13h)
x1 = 6 + 13*4
print(f"\nP2 pair defect mu(D1 D2 E) - (4/169)muE, x1={x1}, x2=7+13h:")
for h in list(range(1, 8)) + [4]:
    x2 = 7 + 13*h
    d = pair_measure(E, x1, x2) - Fraction(4, 169)*muE
    near = abs(x2 - x1)
    print(f"   x2={x2:4d} (|x2-x1|={near:3d}): defect = {float(d):+.7f} "
          f"(x{float(d)*x2:+.4f} per 1/x2)")

# P3: uncovered mass for random 7-lifts vs the independence prediction
pred = muE * Fraction(11, 13)**7
print(f"\nP3 uncovered mass vs independence prediction {float(pred):.6f}:")
random.seed(42)
for trial in range(6):
    hs = [random.randint(2, 40) for _ in R]
    xs = [r + 13*h for r, h in zip(R, hs)]
    U = list(E)
    for x in xs: U = subtract_comb(U, x)
    got = mu(U)
    print(f"   xs={xs}: uncovered = {float(got):.6f} ratio = {float(got/pred):.4f}")
# adversarial: all near-equal x (max resonance)
xs = [499 + i for i in range(7)]  # 7 consecutive integers ~ maximally correlated?
U = list(E)
for x in xs: U = subtract_comb(U, x)
print(f"   consecutive xs={xs}: uncovered = {float(mu(U)):.6f} "
      f"ratio = {float(mu(U)/pred):.4f}")
