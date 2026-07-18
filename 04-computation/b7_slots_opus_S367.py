# opus-2026-07-17-S367 -- HYP-7500: THE TWO NEW B7 LEDGER SLOTS.
#   S6 needs a LOWER bound; S7 needs an UPPER bound.
# S6 lower: iterated containment, the S360 method one level up.
# S7 upper: pairwise bounds are hopeless (rho ~ 1/49 vs target (2lam)^7 ~ 1.7e-6),
#   so use the FRAGMENTATION step -- each new comb multiplies the measure by
#   ~2lam and adds a boundary term proportional to the component count:
#       mu(A n D_x) <= 2*lam*mu(A) + kappa(A) * (2*lam/x).
#   Iterating gives a closed upper bound.  Both tested exactly.
from fractions import Fraction as F
from math import floor, gcd
import random, itertools
LAM = F(1,14)

def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j=[],0,0
    while i<len(u) and j<len(v):
        a,b=max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def meas(iv): return sum(b-a for a,b in iv)
def fold_inter(S):
    cur=teeth01(S[0])
    for x in S[1:]: cur=inter(cur,teeth01(x))
    return cur

print("(A) THE FRAGMENTATION STEP (the S7 upper-bound engine):")
print("    claim: mu(A n D_x) <= 2*lam*mu(A) + kappa(A)*(2*lam/x)")
random.seed(367)
bad=n=0
for _ in range(400):
    k=random.randint(1,4)
    S=sorted(random.sample(range(3,120), k))
    x=random.randint(3,300)
    A=fold_inter(S)
    if not A: continue
    lhs=meas(inter(A,teeth01(x)))
    rhs=2*LAM*meas(A) + len(A)*(2*LAM/x)
    n+=1
    if lhs>rhs: bad+=1
print(f"    {n} (A, x) pairs: violations = {bad}")

print()
print("(B) THE ITERATED S7 UPPER BOUND (apply the step 7 times):")
print("    mu_7 <= (2lam)^7 + accumulated boundary terms")
def upper_iter(S):
    """iterate the fragmentation step, tracking measure and component count."""
    m = F(1); kap = 1
    for x in S:
        m = 2*LAM*m + kap*(2*LAM/x)
        kap = kap + x            # each comb adds at most x new components
    return m
viol=n2=0; ratios=[]
for _ in range(120):
    S=sorted(random.sample(range(20,400), 7))
    ex=meas(fold_inter(S)); ub=upper_iter(S)
    n2+=1
    if ex>ub: viol+=1
    if ub>0: ratios.append(float(ub/max(ex,F(1,10**12))))
ratios.sort()
print(f"    {n2} random 7-tuples: bound violated {viol}/{n2}")
print(f"    bound/exact ratio: median {ratios[len(ratios)//2]:.1f}, "
      f"min {ratios[0]:.1f}  (1.0 would be tight)")
print(f"    equidistribution (2lam)^7 = {float((2*LAM)**7):.3e}")

print()
print("(C) THE S6 LOWER BOUND (iterated containment, S360 method):")
def lower6(S):
    p=(2*LAM)**6
    for i in range(5):
        f = 1 - F(S[i],1)/(2*LAM*S[i+1])
        if f<=0: return F(0)
        p*=f
    return p
viol3=n3=0; pos=0
for _ in range(120):
    a=random.randint(1,20); S=[a]
    for _ in range(5): S.append(S[-1]*random.randint(9,22))
    ex=meas(fold_inter(S)); lb=lower6(S)
    n3+=1
    if ex<lb: viol3+=1
    if lb>0: pos+=1
print(f"    {n3} separated 6-tuples: violations {viol3}/{n3}, floor positive {pos}")
print()
print("  => both B7 slots have a working shape: S6 by iterated containment")
print("     (separated regime), S7 by iterated fragmentation.")

# ---- see also 04-computation/b7_degradation_opus_S367.py, which isolates the
#      k-degradation of the containment floor (the decisive measurement).
