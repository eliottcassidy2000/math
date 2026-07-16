# opus-2026-07-16-S328 -- HYP-7045 / THM-897:
# (A) THE LEVEL-5 WALL: W5(m) = sum_{k<=5} (-1)^k C(m,k) 2^k 13^(5-k) exact;
#     + the exact order-5 Bonferroni referee on a FULL 13-comb packet (no
#     prefix: the whole scale-one problem) -- S1..S5 exact.
# (B) THE TRIPLE BEAT LEMMA T1 (scale-separated; proof = THM-815 discrepancy
#     per component): |mu(W cap D3) - (2/13) mu(W)| <= 22 kappa_W / (169 x3).
#     Referee: exact triples, scale-ordered; check the bound dominates.
# (C) T2 empirics (same-scale incl. resonant): the enhancement vs the pair
#     cascade bound mu3 <= (2/13 + eps) mu2.
from fractions import Fraction
from math import comb, gcd
import math, itertools, random

# ---------- (A) the wall
print("(A) THE LEVEL-5 WALL W5(m)/13^5:")
for m in range(11, 18):
    W = sum((-1)**k * comb(m, k) * 2**k * 13**(5-k) for k in range(6))
    print(f"   m' = {m:2d}: W5 = {W:7d} {'COERCIVE' if W > 0 else 'dead'}"
          + ("   <-- THE FULL PROBLEM" if m == 13 else ""))
W13 = sum((-1)**k * comb(13, k) * 2**k * 13**(5-k) for k in range(6))
print(f"   m'=13 margin: {W13}/371293 = {W13/371293:.4f}; factor: 5x7^2x11x13 = {5*7*7*11*13}")

# exact referee: a full 13-comb packet, no prefix
def teeth(x, a, b):
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def inter(arcs, x):
    out = []
    for (a, b) in arcs: out.extend(teeth(x, a, b))
    return out

def mu(arcs): return sum(b - a for a, b in arcs)

random.seed(28)
xs = sorted(r + 13*random.randint(4, 25) for r in range(1, 14))
print(f"\n   13-comb packet (residues 1..13 lifted): S1..S5 exact "
      f"(this is heavy: C(13,5) = 1287 quintuple intersections)")
FULL = [(Fraction(0), Fraction(1))]
lvl = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
for i in range(13): lvl[1][(i,)] = inter(FULL, xs[i])
S = {1: sum(mu(v) for v in lvl[1].values())}
for k in range(2, 6):
    tot = Fraction(0)
    for idx in itertools.combinations(range(13), k):
        prev = lvl[k-1][idx[:-1]]
        cur = inter(prev, xs[idx[-1]])
        if k < 5: lvl[k][idx] = cur
        tot += mu(cur)
    S[k] = tot
    print(f"   S{k} = {float(tot):.6f}", flush=True)
bound = 1 - S[1] + S[2] - S[3] + S[4] - S[5]
# actual uncovered
U = FULL
for x in xs:
    out = []
    for (a, b) in U:
        cur = a
        for (lo, hi) in sorted(teeth(x, a, b)):
            if lo > cur: out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    U = [iv for iv in out if iv[0] < iv[1]]
print(f"   BONF5 >= {float(bound):+.6f}   actual uncovered = {float(mu(U)):.6f}"
      f"   {'COERCIVE' if bound > 0 else 'not coercive'}")

# ---------- (B) T1 referee
print("\n(B) T1 (scale-separated triple beat): |mu(W cap D3) - (2/13)muW| <= 22 kappa_W/(169 x3)")
E = FULL
for (x1, x2, x3) in [(40, 67, 1201), (33, 90, 2003), (51, 52, 1500), (29, 41, 601)]:
    W = inter(inter(E, x1), x2)
    m3 = mu(inter(W, x3))
    err = abs(m3 - Fraction(2, 13)*mu(W))
    cap = Fraction(22*len(W), 169*x3)
    print(f"   ({x1},{x2},{x3}): err = {float(err):.7f}  cap = {float(cap):.7f}  "
          f"dominates: {err <= cap}")

# ---------- (C) T2 empirics: same-scale, the cascade bound
print("\n(C) same-scale triples: mu3 vs the cascade (2/13)*mu2 (enhancement ratio):")
for (x1, x2, x3) in [(150, 151, 152), (99, 199, 300), (77, 143, 169),
                     (200, 201, 403), (101, 137, 211)]:
    W = inter(inter(E, x1), x2)
    m2 = mu(W); m3 = mu(inter(W, x3))
    print(f"   ({x1},{x2},{x3}): mu2 = {float(m2):.6f}  mu3 = {float(m3):.6f}  "
          f"mu3/((2/13)mu2) = {float(m3/(Fraction(2,13)*m2)):.3f}")
