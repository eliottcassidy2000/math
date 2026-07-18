# opus-2026-07-17-S363 -- HYP-7460: IS THE LAST OPEN REGIME BOUNDED?
# S362 found BONF5 failures concentrate at small speeds.  Before chasing a
# speed threshold, settle the right FORM of the question: every S_k is
# dilation-invariant (mu(cap D_{k v_i}) = mu(cap D_{v_i})), so BONF5 is too --
# hence a failing family dilates to arbitrarily large speeds while still
# failing, and the residual is UNBOUNDED as a set but possibly bounded
# UP TO DILATION (i.e. finitely many PRIMITIVE residual families).
# That is the same shape THM-1025 used: reduce to primitive, finite table.
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random, itertools
LAM = F(1,14)

def teeth(x, lo, hi):
    w = LAM / x; out = []
    for j in range(floor((lo-w)*x), floor((hi+w)*x)+2):
        a,b = max(F(j,x)-w, lo), min(F(j,x)+w, hi)
        if a < b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j = [],0,0
    while i < len(u) and j < len(v):
        a,b = max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a < b: out.append((a,b))
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out
def mu(S):
    cur = teeth(S[0], F(0), F(1))
    for x in S[1:]: cur = inter(cur, teeth(x, F(0), F(1)))
    return sum(b-a for a,b in cur)

print("(1) DILATION INVARIANCE OF THE S_k (hence of BONF5):")
random.seed(363)
bad = n = 0
for _ in range(120):
    k = random.randint(2, 7)
    S = sorted(random.sample(range(3, 60), random.randint(2, 4)))
    n += 1
    if mu(S) != mu([k*x for x in S]): bad += 1
print(f"    {n} (subset, dilation) pairs, sizes 2-4: failures = {bad}")
print("    => every S_k is dilation-invariant, so BONF5(k*V) = BONF5(V).")
print("       A BONF5-failing family therefore dilates to ARBITRARILY LARGE")
print("       speeds while still failing: the residual is UNBOUNDED as a set.")
print("       The correct question is whether it is bounded UP TO DILATION,")
print("       i.e. whether only finitely many PRIMITIVE families fail.")

print()
print("(2) SO: DO THE S362 FAILURES SURVIVE PRIMITIVITY?")
print("    (if the small-speed failures were merely dilates, primitivity kills")
print("     them; if they are genuinely primitive, the residual is real)")
fails = [[32,60,95,144],[49,55,59,72],[63,88,107,186],[71,72,136,158],
         [84,91,167,200],[89,118,219,322]]
for V in fails:
    g = reduce(gcd, V)
    print(f"    {V}: gcd = {g}  {'PRIMITIVE' if g == 1 else 'dilate'}")
print()
print("    => the S362 failures are primitive, so they are genuine residual")
print("       members, not artefacts of a common factor.")
