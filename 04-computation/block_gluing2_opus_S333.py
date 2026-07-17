# opus-2026-07-16-S333 -- HYP-7160 / THM-932 part 2: THE SHARP-INPUT GLUING
# + coercive composition demos. Revision from part 1's negative demo:
#   the theorem's loss terms use kappa(V_{i-1}) EXACT (computable data),
#   not the crude accumulated block bound (30x overcount); and the junction
#   gaps must satisfy the RELATIVE-LOSS LAW loss_{i+1} ~ c * b_i / G_{i+1}
#   (scale-free), so G ~ 300 with c ~ 1..2 is the coercive regime.
# Demos:
#   (D2) TWO-BLOCK 7+6: B1 = first 7 of X; B2 = 300 * (last-6 base);
#        exact mu + exact kappa(V_1) + eta profile => bound; verify sign +
#        exact >= bound.
#   (D3) THREE-BLOCK 5+4+4 at G = 300, 300: all block data exact; exact
#        total mu (background-grade compute); verify.
from fractions import Fraction
from math import floor
import bisect

F = Fraction
lam = F(1, 13)

def subtract_comb(V, x, lam):
    w = F(lam.numerator, lam.denominator * x)
    out = []
    for (a, b) in V:
        cur = a
        j0 = floor((a - w) * x)
        j1 = floor((b + w) * x) + 1
        for j in range(j0, j1 + 1):
            lo, hi = F(j, x) - w, F(j, x) + w
            if hi <= cur: continue
            if lo >= b: break
            if lo > cur: out.append((cur, lo))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return out

def W_block(B):
    V = [(F(0), F(1))]
    for x in sorted(B): V = subtract_comb(V, x, lam)
    return V

def mu(V): return sum(b - a for a, b in V)

def circ_kappa(V):
    if len(V) >= 2 and V[0][0] == 0 and V[-1][1] == 1: return len(V) - 1
    return len(V) if V else 0

class Prefix:
    def __init__(self, W):
        self.s = [a for a, b in W]; self.e = [b for a, b in W]
        self.c = [F(0)]
        for a, b in W: self.c.append(self.c[-1] + (b - a))
        self.tot = self.c[-1]
    def cdf(self, t):
        k = bisect.bisect_right(self.s, t) - 1
        if k < 0: return F(0)
        return self.c[k] + min(max(t - self.s[k], F(0)), self.e[k] - self.s[k])
    def arc(self, a, l):
        b = a + l
        if b <= 1: return self.cdf(b) - self.cdf(a)
        return (self.tot - self.cdf(a)) + self.cdf(b - 1)

def eta(W, l):
    P = Prefix(W)
    cands = set()
    for a, b in W:
        for p in (a, b):
            cands.add(p % 1); cands.add((p - l) % 1)
    return min(P.arc(a, l) for a in cands) / l

print("=" * 72)
print("ETA PROFILES of the base blocks (exact; scale unit 1/300):")
BASE7 = [300, 406, 511, 652, 862, 963, 1074]
BASE6 = [862, 963, 1074, 1357, 1459, 1571]   # used dilated as B2 in D2
BASE5 = [300, 406, 511, 652, 862]
BASE4 = [300, 406, 511, 652]
profiles = {}
for name, B in [('B7', BASE7), ('B6', BASE6), ('B5', BASE5), ('B4', BASE4)]:
    W = W_block(B)
    profiles[name] = (W, mu(W), circ_kappa(W))
    es = {}
    for num, den in [(1, 600), (1, 300), (2, 300), (4, 300), (8, 300)]:
        es[(num, den)] = eta(W, F(num, den))
    print(f"  {name} ({len(B)} combs): mu = {float(mu(W)):.6f}  kappa = "
          f"{circ_kappa(W)}  eta(1/600, 1/300, 2/300, 4/300, 8/300) = "
          + ", ".join(f"{float(es[k]):.4f}" for k in
                      [(1,600),(1,300),(2,300),(4,300),(8,300)]))
    profiles[name] = (W, mu(W), circ_kappa(W), es)

print()
print("=" * 72)
print("(D2) TWO-BLOCK COMPOSITION 7+6, junction gap G = 300:")
G = 300
B1 = BASE7
B2 = [G * x for x in BASE6]
print(f"  B1 = {B1}")
print(f"  B2 = 300 * {BASE6}  (gap min(B2)/max(B1) = {B2[0]/B1[-1]:.1f})")
W1, m1, k1, es1 = profiles['B7']
# choose the scale l2 = c/min(B2) optimizing bound = m1*eta2(c) - k1*l2
WB6, mB6, kB6, esB6 = profiles['B6']
best = None
for num, den in [(1, 600), (1, 300), (2, 300), (4, 300), (8, 300)]:
    # eta_{B2}(l) = eta_{base6}(G l); pick l = (num/den)/G * (862/300 adj):
    # dilation: B2 = G*base6 => eta_B2(l) = eta_base6(G*l); base6 scale unit
    # is 1/862.. but we computed eta_base6 at k/300 grids -- use them directly
    l2 = F(num, den * G)
    e2 = esB6[(num, den)]
    bnd = m1 * e2 - k1 * l2
    if best is None or bnd > best[0]: best = (bnd, l2, e2, (num, den))
bnd2, l2, e2, sc2 = best
V = [(F(0), F(1))]
for x in sorted(B1 + B2): V = subtract_comb(V, x, lam)
mu2 = mu(V)
print(f"  block data: m1 = {float(m1):.6f}, kappa(V1) = {k1} (exact); "
      f"eta2 = {float(e2):.6f} at l2 = ({sc2[0]}/{sc2[1]})/300")
print(f"  GLUED BOUND = m1*eta2 - kappa1*l2 = {float(bnd2):+.6f}"
      f"   {'COERCIVE' if bnd2 > 0 else 'not coercive'}")
print(f"  EXACT mu(V) = {float(mu2):.6f}   bound holds: {mu2 >= bnd2}")

print()
print("=" * 72)
print("(D3) THREE-BLOCK COMPOSITION 5+4+4, gaps G2 = G3 = 300:")
B1 = BASE5
B2 = [300 * x for x in BASE4]
B3 = [90000 * x for x in BASE4]
print(f"  B1 = {B1}")
print(f"  B2 = {B2}  (gap {B2[0]/B1[-1]:.1f})")
print(f"  B3 = 90000*{BASE4}  (gap {B3[0]/B2[-1]:.1f}; max speed {B3[-1]:,})")
W1, m1, k1, es1 = profiles['B5']
WB4, mB4, kB4, esB4 = profiles['B4']
# stage 1+2 exact
V2 = [(F(0), F(1))]
for x in sorted(B1 + B2): V2 = subtract_comb(V2, x, lam)
m12, k12 = mu(V2), circ_kappa(V2)
# choose scales
best = None
for n2, d2 in [(1, 600), (1, 300), (2, 300), (4, 300)]:
    for n3, d3 in [(1, 600), (1, 300), (2, 300), (4, 300)]:
        l2 = F(n2, d2 * 300); e2v = esB4[(n2, d2)]
        l3 = F(n3, d3 * 90000); e3v = esB4[(n3, d3)]
        b = m1 * e2v * e3v - k1 * l2 * e3v - k12 * l3
        if best is None or b > best[0]:
            best = (b, l2, e2v, l3, e3v, (n2, d2), (n3, d3))
bnd3, l2, e2v, l3, e3v, sc2, sc3 = best
print(f"  block data: m1 = {float(m1):.6f} kappa(V1) = {k1}; after B2: "
      f"mu(V2) = {float(m12):.6f} kappa(V2) = {k12} (exact)")
print(f"  scales: eta2 = {float(e2v):.4f} at ({sc2[0]}/{sc2[1]})/300; "
      f"eta3 = {float(e3v):.4f} at ({sc3[0]}/{sc3[1]})/90000")
print(f"  GLUED BOUND = m1 e2 e3 - k1 l2 e3 - kappa(V2) l3 = "
      f"{float(bnd3):+.6f}   {'COERCIVE' if bnd3 > 0 else 'not coercive'}")
print(f"  intermediate check: mu(V2) >= m1*eta2 - k1*l2 = "
      f"{float(m1 * e2v - k1 * l2):+.6f}: {m12 >= m1 * e2v - k1 * l2}")
V3 = V2
for x in sorted(B3): V3 = subtract_comb(V3, x, lam)
mu3 = mu(V3)
print(f"  EXACT mu(V3) = {float(mu3):.6f}   bound holds: {mu3 >= bnd3}   "
      f"NONEMPTY: {mu3 > 0}")
print()
print("  13 speeds, three bands, max ratio "
      f"{(B3[-1]/300):,.0f}; within-band ratios < 3 (NOT lacunary),")
print("  gaps 300 (NOT single-band): certified by COMPOSITION alone.")
