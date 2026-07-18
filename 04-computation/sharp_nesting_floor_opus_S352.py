from fractions import Fraction as F
from math import gcd, floor
import random
LAM = F(1,14)

def pair_overlap_exact(a, b):
    g = gcd(a, b); tot = F(0); Mb = LAM*(a+b); m = 0
    while True:
        if F(m) >= Mb*g and m > 0: break
        for mm in ([0] if m == 0 else [m,-m]):
            if mm % g: continue
            if F(abs(mm)) >= Mb: continue
            w = min(2*LAM*a, 2*LAM*b, LAM*(a+b)-abs(mm))
            if w > 0: tot += F(w, a*b)*g
        m += g
        if m > int(Mb)+2: break
    return tot

print("TEST: mu(D_a n D_b) >= 4*lam^2 - 2*lam*(a/b) ?   [lam=1/14: 1/49 - (1/7)(a/b)]")
random.seed(352)
viol = 0; tot = 0; worst = None
for _ in range(4000):
    a = random.randint(1, 200); b = random.randint(a, 40*a)
    mu = pair_overlap_exact(a, b)
    bound = 4*LAM**2 - 2*LAM*F(a, b)
    tot += 1
    slack = mu - bound
    if slack < 0:
        viol += 1
        if worst is None or slack < worst[0]: worst = (slack, a, b, mu, bound)
print(f"  {tot} random pairs (ratio up to 40): violations = {viol}")
if worst: print(f"  worst violation: a={worst[1]} b={worst[2]} mu={worst[3]} bound={worst[4]} slack={float(worst[0]):.6f}")

# separated regime specifically
print("\n  SEPARATED regime b >= 7a (where the bound is positive):")
viol2 = 0; n2 = 0; minslack = None
for _ in range(3000):
    a = random.randint(1, 150); b = random.randint(7*a, 60*a)
    mu = pair_overlap_exact(a, b); bound = 4*LAM**2 - 2*LAM*F(a,b)
    n2 += 1
    s = mu - bound
    if s < 0: viol2 += 1
    if minslack is None or s < minslack[0]: minslack = (s, a, b)
print(f"  {n2} pairs with b >= 7a: violations = {viol2}; min slack = {float(minslack[0]):.6f} at a={minslack[1]} b={minslack[2]}")
print(f"\n  sanity: 4*lam^2 = {4*LAM**2} = 1/49; bound positive iff a/b < 2*lam = 1/7")
