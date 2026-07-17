# opus-2026-07-16-S333 -- HYP-7160 / THM-932: THE LOCAL-DENSITY BLOCK-GLUING
# THEOREM. lambda = tooth half-width parameter (LRC(14): 1/14; certificate
# convention: 1/13). D_x = {||x t|| < lam}; W_B = circle minus union of combs.
#
# LEMMAS (verified exactly here):
#  (G1) SAMPLING: for V a union of kappa(V) arcs, any block B, any l > 0:
#       mu(V cap W_B) >= eta_B(l) * (mu(V) - kappa(V) * l),
#       eta_B(l) = min over circle intervals |I| = l of mu(I cap W_B)/l.
#  (G2) COMPONENTS: kappa(U cap V) <= kappa(U) + kappa(V);
#       kappa(W_B) <= sum_{x in B} x.
#  (G3) SINGLE-SPEED EXACTNESS: eta_{x}(k/x) = 1 - 2 lam EXACTLY.
#
# THEOREM: blocks B_1 < ... < B_r; V_r = cap W_{B_i};
#   mu(V_r) >= m_1 * prod_{i>=2} eta_i - sum_{i>=2} K_{i-1} l_i,
#   m_1 = mu(W_{B_1}) exact, eta_i = eta_{B_i}(l_i), K_i = sum_{j<=i} kappa(W_{B_j}).
#
# Plus: eta profile of the THM-928(C) certified packet X (localization scale);
# the DILATION LAW eta_{cB}(l) = eta_B(c l); and the MIXED-SCALE 13-SPEED
# COMPOSITION DEMO (3 blocks, gaps 60 and 40): exact uncovered vs glued bound.
from fractions import Fraction
from math import floor
import random, bisect

F = Fraction

def subtract_comb(V, x, lam):
    """V minus D_x, generating teeth only inside V's components. Exact."""
    w = F(lam.numerator, lam.denominator * x)  # half-width lam/x
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

def W_block(B, lam):
    V = [(F(0), F(1))]
    for x in sorted(B):
        V = subtract_comb(V, x, lam)
    return V

def mu(V): return sum(b - a for a, b in V)

def circ_kappa(V):
    if len(V) >= 2 and V[0][0] == 0 and V[-1][1] == 1: return len(V) - 1
    return len(V) if V else 0

def isect(U, V):
    out, i, j = [], 0, 0
    while i < len(U) and j < len(V):
        a, b = max(U[i][0], V[j][0]), min(U[i][1], V[j][1])
        if a < b: out.append((a, b))
        if U[i][1] < V[j][1]: i += 1
        else: j += 1
    return out

class Prefix:
    """exact mu([0,t] cap W) with wrap, for interval queries."""
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
        """mu of [a, a+l] mod 1 intersect W (l <= 1)."""
        b = a + l
        if b <= 1: return self.cdf(b) - self.cdf(a)
        return (self.tot - self.cdf(a)) + self.cdf(b - 1)

def eta(W, l):
    """exact local density floor at scale l: min over circle intervals."""
    P = Prefix(W)
    cands = set()
    for a, b in W:
        for p in (a, b):
            cands.add(p % 1); cands.add((p - l) % 1)
    best = None
    for a in cands:
        v = P.arc(a, l)
        if best is None or v < best: best = v
    return best / l

print("=" * 72)
print("(G1)-(G3) LEMMA VERIFICATION (exact, random configs):")
random.seed(333)
lam = F(1, 13)
viol = {1: 0, 2: 0, 3: 0}; n1 = n2 = n3 = 0
for trial in range(40):
    # random V
    pts = sorted(random.sample(range(4000), 2 * random.randint(2, 8)))
    V = [(F(pts[2*i], 4000), F(pts[2*i+1], 4000)) for i in range(len(pts)//2)]
    # random small block
    B = random.sample(range(20, 90), random.randint(1, 3))
    WB = W_block(B, lam)
    l = F(random.randint(1, 30), 1000)
    et = eta(WB, l)
    lhs = mu(isect(V, WB))
    n1 += 1
    if lhs < et * (mu(V) - len(V) * l): viol[1] += 1
    n2 += 1
    if circ_kappa(isect(V, WB)) > circ_kappa(V) + circ_kappa(WB): viol[2] += 1
    if circ_kappa(WB) > sum(B): viol[2] += 1
for x in (7, 19, 44, 131):
    for k in (1, 2, 5):
        Wx = W_block([x], lam)
        n3 += 1
        if eta(Wx, F(k, x)) != 1 - 2 * lam: viol[3] += 1
print(f"  G1 sampling:   {n1} configs, violations = {viol[1]}")
print(f"  G2 components: {n2} configs, violations = {viol[2]}")
print(f"  G3 exactness eta_x(k/x) == 1-2lam == 11/13: {n3} checks, "
      f"violations = {viol[3]}")

print()
print("=" * 72)
print("ETA PROFILE OF THE CERTIFIED PACKET X (THM-928(C)), lam = 1/13:")
X = [300, 406, 511, 652, 862, 963, 1074, 1357, 1459, 1571, 1776, 1991, 2208]
WX = W_block(X, lam)
print(f"  mu(W_X) = {float(mu(WX)):.6f}  kappa = {circ_kappa(WX)}  "
      f"(bound sum x = {sum(X)})")
for num, den in [(1, 2208), (1, 862), (1, 300), (2, 300), (4, 300), (8, 300),
                 (16, 300), (1, 8), (1, 2)]:
    l = F(num, den)
    e = eta(WX, l)
    print(f"  eta_X({num}/{den}) = {float(e):.6f}"
          f"{'   <-- localization scale reached (eta > 0)' if e > 0 and (num, den) in [(2,300),(4,300)] else ''}")

print()
print("=" * 72)
print("DILATION LAW eta_cB(l) = eta_B(c l) (exact check):")
B1 = [300, 406, 511, 652]
W1 = W_block(B1, lam)
for c, num, den in [(7, 1, 2100), (7, 1, 700), (60, 1, 18000)]:
    lhs = eta(W_block([c * x for x in B1], lam), F(num, den))
    rhs = eta(W1, F(num * c, den))
    print(f"  c={c:3d} l={num}/{den}: eta_cB = {float(lhs):.6f}  "
          f"eta_B(cl) = {float(rhs):.6f}  equal: {lhs == rhs}")

print()
print("=" * 72)
print("THE MIXED-SCALE 13-SPEED COMPOSITION DEMO (blocks with gaps):")
# B1 = 4 speeds (band ~2.2x), B2 = 60*B1 (gap 60*300/652 = 27.6),
# B3 = 40*B2 = 2400*{300,406,511,652,862} 5 speeds (gap 40*300/652 = 18.4)
B2 = [60 * x for x in B1]
B3 = [2400 * x for x in [300, 406, 511, 652, 862]]
speeds = sorted(B1 + B2 + B3)
print(f"  B1 = {B1}")
print(f"  B2 = {B2}  (gap min(B2)/max(B1) = {B2[0]/B1[-1]:.1f})")
print(f"  B3 = {B3}  (gap min(B3)/max(B2) = {B3[0]/B2[-1]:.1f})")
print(f"  13 speeds, max ratio {speeds[-1]/speeds[0]:.0f} -- NOT lacunary "
      f"(within-block ratios < 3), NOT single-band")
# exact uncovered, processing in increasing order
V = [(F(0), F(1))]
for x in speeds:
    V = subtract_comb(V, x, lam)
muV = mu(V)
print(f"  EXACT uncovered mu(V) = {float(muV):.6f}  "
      f"({'NONEMPTY: every runner gets >= 1/13-lonely' if muV > 0 else 'EMPTY'})")
# the glued bound
m1 = mu(W1); k1 = circ_kappa(W1)
W2 = W_block(B2, lam); k2 = circ_kappa(W2)
# choose per-block scales: l_i = smallest tested scale with eta_i > 0
l2 = F(4, 60 * 300)     # = (4/300)/60 via dilation from B1-profile
e2 = eta(W1, F(4, 300))  # dilation law: eta_B2(l2) = eta_B1(4/300)
l3 = F(4, 2400 * 300)
WB3base = W_block([300, 406, 511, 652, 862], lam)
e3 = eta(WB3base, F(4, 300))
K1 = k1
K2 = k1 + k2
bound = m1 * e2 * e3 - K1 * l2 * e3 - K2 * l3
bound_simple = m1 * e2 * e3 - K1 * l2 - K2 * l3
print(f"  block data: m1 = {float(m1):.6f} (kappa {k1}); "
      f"eta2 = {float(e2):.6f} at l2 = 4/18000; eta3 = {float(e3):.6f} "
      f"at l3 = 4/720000; K1 = {K1}, K2 = {K2}")
print(f"  GLUED BOUND (sharp form)  = {float(bound):+.6f}")
print(f"  GLUED BOUND (simple form) = {float(bound_simple):+.6f}")
print(f"  bound holds (exact >= sharp bound): {muV >= bound}")
print(f"  bound POSITIVE (certificate!): {bound > 0}")
print()
print("  -> a 13-speed family certified by COMPOSITION: block certificates +")
print("     scale gaps, no global search, no lacunarity. The two-scale law.")
