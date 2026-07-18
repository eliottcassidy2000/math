# opus-2026-07-17-S360 -- HYP-7430: (A) THE S4 QUADRUPLE FLOOR by iterated
# containment; (B) the divisor-count reduction for the sieve/existence step.
from fractions import Fraction as F
from math import gcd, floor
import random, itertools
LAM = F(1,14)

def teeth(x, lo, hi):
    w = LAM / x; out = []
    for j in range(floor((lo - w)*x), floor((hi + w)*x) + 2):
        a, b = max(F(j,x) - w, lo), min(F(j,x) + w, hi)
        if a < b: out.append((a,b))
    return out

def inter(u, v):
    out, i, j = [], 0, 0
    while i < len(u) and j < len(v):
        a, b = max(u[i][0], v[j][0]), min(u[i][1], v[j][1])
        if a < b: out.append((a,b))
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out

def mu_k(S):
    """exact measure of the k-fold comb intersection over [0,1]."""
    cur = teeth(S[0], F(0), F(1))
    for x in S[1:]: cur = inter(cur, teeth(x, F(0), F(1)))
    return sum(b - a for a, b in cur)

print("(A) THE S4 QUADRUPLE FLOOR.  independence value (2*lam)^4 =", (2*LAM)**4,
      "=", float((2*LAM)**4))
print("    THM-1012 generalises by ITERATED counting: inside each a-arc count")
print("    b-cells, inside each b-arc count c-cells, etc.  Each level loses at")
print("    most one cell, so the natural floor is")
print("      mu_4 >= (2lam)^4 * prod_i (1 - a_i/(2*lam*a_{i+1}))   [separated]")
random.seed(360)
def floor4(S):
    p = (2*LAM)**4
    for i in range(3):
        f = 1 - F(S[i], 1) / (2*LAM*S[i+1])
        if f <= 0: return F(0)
        p *= f
    return p
viol = n = 0; sep_ok = sep_n = 0
for _ in range(150):
    a = random.randint(1, 30)
    S = [a]
    for _ in range(3): S.append(S[-1] * random.randint(8, 25))   # separated
    ex = mu_k(S); fl = floor4(S)
    n += 1
    if ex < fl: viol += 1
    if fl > 0:
        sep_n += 1
        if ex >= fl: sep_ok += 1
print(f"    separated quadruples: {n} tested, floor violated {viol}/{n};"
      f" floor positive in {sep_n}")
print()
print("    COMPARABLE quadruples (the dense-core regime): does the floor survive?")
viol2 = n2 = 0; pos = 0
for _ in range(60):
    m0 = random.randint(20, 60)
    S = sorted(random.sample(range(m0, 13*m0), 4))
    ex = mu_k(S); fl = floor4(S)
    n2 += 1
    if ex < fl: viol2 += 1
    if fl > 0: pos += 1
print(f"    {n2} comparable quadruples: violations {viol2}, floor positive in {pos}")
print("    => as with pairs, the INDEPENDENCE-type floor is positive only when")
print("       separated; comparable quadruples need the sawtooth analogue.")
print()
print("    exact S4 behaviour on comparable quadruples (what a floor must beat):")
vals = sorted(float(mu_k(sorted(random.sample(range(50, 650), 4)))) for _ in range(40))
print(f"      exact mu_4 range [{vals[0]:.6f}, {vals[-1]:.6f}], median {vals[len(vals)//2]:.6f}")
print(f"      independence (2lam)^4 = {float((2*LAM)**4):.6f}")
