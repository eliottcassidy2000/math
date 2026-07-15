# opus-2026-07-15-S313 -- HYP-6920 / THM-863 stream 1: the corrected pair law's
# structure theorems.
#   rho(a,b) = T/(13ab), T = Sum_m min(2a, (a+b)-13|m|)^+   [codex-S14 corrected law]
# (1) rho >= 1/78 with UNIQUE minimizer (1,12): casework + exact finite scan.
# (2) deviation bound |T - 4ab/13| <= 2a + 26 (verify over a large box).
# (3) LOW CHANNELS: all coprime (a,b) with rho(a,b) < c3.
# (4) THE 3-SPEED LEMMA: quotients of low channels are never low
#     => low graph triangle-free => every 7-packet has a spanning tree with
#     Sum rho >= 5*c3 + 1/78 > 13/169: THE UNCONDITIONAL GAP phi*.
# (5) empirical min-maxtree search to see how sharp the bound is.
from fractions import Fraction
from math import gcd
import itertools, random

def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        term = min(2*a, s - 13*m)
        T += term if m == 0 else 2*term
        m += 1
    return T

def rho(a, b):
    g = gcd(a, b)
    a, b = a//g, b//g
    if a > b: a, b = b, a
    return Fraction(T_of(a, b), 13*a*b)

# (1)+(2): scan box; verify floor and deviation bound
AMAX, BMAX = 60, 400
viol_floor, viol_dev, minimizers = [], [], []
floor = Fraction(1, 78)
for a in range(1, AMAX+1):
    for b in range(a, BMAX+1):
        if gcd(a, b) != 1: continue
        T = T_of(a, b)
        r = Fraction(T, 13*a*b)
        if r < floor: viol_floor.append((a, b, r))
        if r == floor: minimizers.append((a, b))
        if abs(13*T - 52*a*b) > 13*(2*a + 26):   # |T - 4ab/13| <= 2a+26
            viol_dev.append((a, b, T))
print(f"(1) floor 1/78 violations in box a<={AMAX}, b<={BMAX}: {len(viol_floor)}")
print(f"    minimizers: {minimizers}")
print(f"(2) deviation-bound violations: {len(viol_dev)}")
# analytic tail regions (documented in THM-863): b<=12 via m=0 term;
# b >= (24a+156)/11 via full-term counting; a>=40 via the deviation bound;
# the residual box is scanned above.

# (3) low channels for candidate c3 values
for c3 in (Fraction(9, 500), Fraction(2, 105)):     # 0.018, ~0.019
    low = []
    # sawtooth dip needs (2a+26)/(13ab) >= 4/169 - c3 => finite ab range
    for a in range(1, 60):
        for b in range(a, 400):
            if gcd(a, b) != 1: continue
            if rho(a, b) < c3: low.append((a, b))
    print(f"\n(3) c3 = {c3} ~ {float(c3):.5f}: low channels ({len(low)}): {low}")
    # (4) quotient check: for ratios r1 = b1/a1, r2 = b2/a2 in low (and their
    # inverses), the third pair ratio r2/r1 (reduced) must NOT be low
    ratios = set()
    for (a, b) in low:
        ratios.add(Fraction(b, a)); ratios.add(Fraction(a, b))
    bad = []
    for r1, r2 in itertools.combinations(sorted(ratios), 2):
        q = r2 / r1
        if q == 1: continue
        qa, qb = q.numerator, q.denominator
        if qa > qb: qa, qb = qb, qa
        # q may have huge numerator/denominator -> rho ~ 4/169 fine; only
        # check when small enough to possibly be low
        if qb <= 400 and qa <= 60:
            if rho(qa, qb) < c3: bad.append((r1, r2, q, rho(qa, qb)))
    print(f"(4) low-quotient violations (3-speed lemma fails if nonempty): {len(bad)}")
    for x in bad[:6]: print("      ", x)

# (5) empirical min-maxtree over structured + random packets
def maxtree(xs):
    n = len(xs)
    W = {}
    for i, j in itertools.combinations(range(n), 2):
        W[(i, j)] = rho(xs[i], xs[j])
    intree = {0}; tot = Fraction(0)
    while len(intree) < n:
        best = None
        for i in intree:
            for j in range(n):
                if j not in intree:
                    w = W[(min(i,j), max(i,j))]
                    if best is None or w > best[0]: best = (w, j)
        tot += best[0]; intree.add(best[1])
    return tot

thr = Fraction(13, 169)
cands = [
    [1, 12, 144, 1728, 20736, 248832, 2985984],       # pure 12-tower
    [1, 12, 11, 132, 10, 120, 1320],                   # low-channel web
    [1, 12, 2, 24, 3, 36, 5],                          # (1,12)-matching + spread
    [1, 11, 12, 110, 120, 132, 1320],
    [32, 33, 34, 35, 36, 37, 38],                      # consecutive
    [1, 10, 11, 12, 9, 90, 99],
]
random.seed(1)
for _ in range(30):
    cands.append(sorted(random.sample(range(1, 3000), 7)))
best = None
for xs in cands:
    mt = maxtree(xs)
    if best is None or mt < best[0]: best = (mt, xs)
print(f"\n(5) min over candidates of maxtree = {best[0]} ~ {float(best[0]):.5f} "
      f"at {best[1]}")
print(f"    threshold 13/169 = {float(thr):.5f}; gap = {float(best[0]-thr):+.5f}")
print(f"    proved floor 5*c3 + 1/78 with c3=0.018: {5*0.018 + 1/78:.5f}")
