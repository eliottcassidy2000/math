# opus-2026-07-17-S336 -- HYP-7190 part 2: THE CLUSTER GAP LEMMA (exact).
# Whenever the displayed numerator is positive, any k <= 6 positive speeds B
# leave, inside the window [a,b], a safe subinterval
# (all ||w t|| >= 1/14) of width
#   delta >= [ (1 - k/7)(b-a) - k/(7 m) ] / (1 + k + (b-a) * Sum_B x)
# (union bound on covered measure + pigeonhole on <= 1 + k + (b-a)Sum x
# components).  This powers the BLOCK nested cascade: junctions of ratio
# >= G feed the next block's window condition min(B') * delta >= 2.
from fractions import Fraction
from math import floor
import random

F = Fraction

def subtract_comb(V, x):
    w = F(1, 14 * x)
    out = []
    for (a, b) in V:
        cur = a
        for j in range(floor((a - w) * x), floor((b + w) * x) + 2):
            lo, hi = F(j, x) - w, F(j, x) + w
            if hi <= cur: continue
            if lo >= b: break
            if lo > cur: out.append((cur, lo))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return out

random.seed(3360)
print("(1) THE CLUSTER GAP LEMMA, exact battery:")
worst = None
viol = 0
positive = 0
positive_by_k = {k: 0 for k in range(1, 7)}
for trial in range(400):
    k = random.randint(1, 6)
    m = random.randint(2, 400)
    rho = random.choice([2, 5, 13, 50, 200])
    B = sorted(random.sample(range(m, max(m + k, m * rho)), k))
    a = F(random.randint(0, 1000), 1009)
    L = F(random.randint(2, 40), 7 * B[0])   # window length ~ scale of min
    b = a + L
    V = [(a, b)]
    for x in B: V = subtract_comb(V, x)
    widest = max((hi - lo for lo, hi in V), default=F(0))
    Sx = sum(B)
    bound = (F(7 - k, 7) * L - F(k, 7 * B[0])) / (1 + k + L * Sx)
    # A nonpositive bound makes no existence claim; audit only theorem rows.
    if bound > 0:
        positive += 1
        positive_by_k[k] += 1
        if widest < bound:
            viol += 1
    r = float(widest / bound) if bound > 0 else float('inf')
    if bound > 0 and (worst is None or r < worst[0]):
        worst = (r, k, B[:3], float(L))
print(f"   positive-bound rows = {positive}/400 by k = {positive_by_k}; violations = {viol}")
print(f"   worst widest/bound ratio = {worst[0]:.3f} (k={worst[1]})")

# The original length range has no positive k=6 row.  Force B[0]=m and
# mL>6, the exact positivity threshold for k=6, to exercise that boundary.
rng6 = random.Random(3361)
k6_viol = 0
k6_worst = None
for _ in range(100):
    m = rng6.randint(2, 400)
    rho = rng6.choice([2, 5, 13, 50, 200])
    B = [m] + sorted(rng6.sample(range(m + 1, max(m + 6, m * rho)), 5))
    a = F(rng6.randint(0, 1000), 1009)
    L = F(rng6.randint(43, 84), 7 * m)
    V = [(a, a + L)]
    for x in B:
        V = subtract_comb(V, x)
    widest = max((hi - lo for lo, hi in V), default=F(0))
    bound = (F(1, 7) * L - F(6, 7 * m)) / (7 + L * sum(B))
    assert bound > 0
    if widest < bound:
        k6_viol += 1
    ratio = float(widest / bound)
    if k6_worst is None or ratio < k6_worst:
        k6_worst = ratio
print(f"   targeted positive k=6 rows = 100; violations = {k6_viol}; "
      f"worst ratio = {k6_worst:.3f}")

print()
print("(2) THE BLOCK NESTED CASCADE 6+6+1 (end-to-end exact):")
# blocks: B1 = 6 arbitrary comparable-ish, B2 = G*scale 6 more, B3 = singleton
for G in (60, 200):
    B1 = sorted(random.sample(range(10, 60), 6))
    m2 = B1[-1] * G
    B2 = sorted(random.sample(range(m2, 3 * m2), 6))
    B3 = [B2[-1] * G]
    speeds = B1 + B2 + B3
    # nested construction: window [0,1] -> B1 widest gap -> B2 widest gap in it -> B3
    V = [(F(0), F(1))]
    ok_chain = True
    for Bi in (B1, B2, B3):
        for x in Bi: V = subtract_comb(V, x)
        if not V: ok_chain = False; break
        big = max(V, key=lambda iv: iv[1] - iv[0])
        V = [big]   # nest into the widest component
    if ok_chain:
        lo, hi = V[0]
        t = (lo + hi) / 2
        dists = sorted(min((x * t) % 1, 1 - (x * t) % 1) for x in speeds)
        print(f"   G={G:4d}: 13 speeds (max {speeds[-1]:,}); nested witness "
              f"t = {t.numerator}/{t.denominator}; min dist = {dists[0]} "
              f">= 1/14: {dists[0] >= F(1, 14)}")
    else:
        print(f"   G={G:4d}: chain DIED (window emptied) -- G too small")
print()
print("(3) junction-constant probe: smallest G with 20/20 success (6+6+1):")
for G in (4, 8, 15, 30, 60):
    succ = 0
    for _ in range(20):
        B1 = sorted(random.sample(range(10, 60), 6))
        m2 = B1[-1] * G
        B2 = sorted(random.sample(range(m2, 3 * m2), 6))
        B3 = [B2[-1] * G]
        V = [(F(0), F(1))]
        alive = True
        for Bi in (B1, B2, B3):
            for x in Bi: V = subtract_comb(V, x)
            if not V: alive = False; break
            V = [max(V, key=lambda iv: iv[1] - iv[0])]
        if alive: succ += 1
    print(f"   G = {G:3d}: {succ}/20 nested chains survive")
