# opus-2026-07-17-S376 -- IS THERE A MEASURE GAP ABOVE THE TIGHT FAMILY?
# The trade-off hypothesis is refuted (uncovered is independent of q0).  The
# remaining question worth settling: {1,...,13} has uncovered EXACTLY 0, and my
# adversarial minimiser reached 0.063.  If uncovered were either 0 or >= c > 0,
# that GAP would be a strong theorem -- it would say the tight family is
# isolated and everything else is robustly lonely.  Test: perturb {1,...,13}.
from fractions import Fraction as F
from functools import reduce
from math import gcd
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))

print("(4) PERTURBING THE TIGHT FAMILY {1,...,13}")
print("    family                                    uncovered")
base=list(range(1,14))
print(f"    {str(base):42s} {float(uncovered(base)):.8f}   <- TIGHT")
for repl,new in [(13,14),(13,15),(13,16),(13,20),(13,26),(12,14),(1,14),(1,15),(7,14)]:
    V=sorted([x for x in base if x!=repl]+[new])
    if len(set(V))!=13: continue
    print(f"    swap {repl:2d}->{new:2d}: {str(V):30s} {float(uncovered(V)):.8f}")

print()
print("(5) SCALED PERTURBATIONS -- k*{1..13} with one speed nudged")
for k in [2,3]:
    base=[k*i for i in range(1,14)]
    print(f"    {k}*{{1..13}} itself: {float(uncovered(base)):.8f}")
    for d in [1,-1,2]:
        V=sorted(base[:-1]+[base[-1]+d])
        if len(set(V))!=13: continue
        print(f"      last speed {base[-1]}->{base[-1]+d}: uncovered = {float(uncovered(V)):.8f}")

print()
print("(6) VERDICT ON THE GAP")
vals=[]
base=list(range(1,14))
for repl in base:
    for new in range(14,30):
        V=sorted([x for x in base if x!=repl]+[new])
        if len(set(V))!=13: continue
        vals.append((float(uncovered(V)), repl, new))
vals.sort()
print(f"    smallest uncovered among single-swap perturbations: {vals[0][0]:.8f}")
print(f"      (swap {vals[0][1]} -> {vals[0][2]})")
print(f"    next few: {[('%.6f'%v) for v,_,_ in vals[1:5]]}")
print()
print("    if these descend continuously toward 0, there is NO gap and the")
print("    tight family is a limit point, not an isolated one.")
