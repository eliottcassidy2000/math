from fractions import Fraction as F
from itertools import combinations
V=[1,2,3,4,5,6,7,8,9,10,11,13,36]
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
print(f"V = {V}   distinct: {len(set(V))==13}")
# exhaustive over ALL critical denominators (pinch lemma set)
D=set()
for a,b in combinations(V,2):
    D.add(a+b)
    if a!=b: D.add(abs(a-b))
for v in V: D.add(2*v)
best=F(0); arg=None
for q in sorted(D):
    if q<2: continue
    for p in range(1,q):
        g=gap_at(V,F(p,q))
        if g>best: best,arg=g,F(p,q)
print(f"M = {best} = {float(best):.8f}   at t* = {arg}")
print(f"  1/14 = {float(F(1,14)):.8f}   2/27 = {float(F(2,27)):.8f}")
print(f"  strictly inside the gap: {F(1,14) < best < F(2,27)}")
# active pair and determinant
up=[];dn=[]
for v in V:
    x=v*arg; a=int(x); r=x-a
    if r<0: a-=1; r+=1
    if r==best: up.append((v,a))
    if 1-r==best: dn.append((v,a+1))
print(f"  straddling: up={up}  down={dn}")
if up and dn:
    vi,ai=up[0]; vj,aj=dn[0]; Dt=vi*aj-vj*ai
    print(f"  active pair ({vi},{vj})  D={Dt}  s={vi+vj}  D/s = {F(Dt,vi+vj)}"
          f"   matches M: {F(Dt,vi+vj)==best}")
    print(f"  s <= 14D ?  {vi+vj} <= {14*Dt}  -> {vi+vj<=14*Dt}  (LRC(14) holds here)")
# independent cross-check: brute force over a fine rational grid
print()
print("  independent check: scan ALL p/q with q <= 200")
b2=F(0); a2=None
for q in range(1,201):
    for p in range(1,q):
        g=gap_at(V,F(p,q))
        if g>b2: b2,a2=g,F(p,q)
print(f"    best over q<=200: {b2} at {a2}   agrees with pinch-set answer: {b2==best}")
print()
print("  relation to the second extremal {1..11,13,24}:")
print("    same shape {1,...,11,13,x}: x=24 -> M=1/14 (extremal); x=36 -> M=3/41 (in gap)")
print(f"    36 = 3*12, 24 = 2*12; denominators 41 = 3*14-1, 14 = 1*14")
