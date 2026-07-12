from math import gcd
from fractions import Fraction as F
from itertools import combinations
def lcm(a,b): return a//gcd(a,b)*b
def bunch(E):
    nz=[abs(e) for e in E if e]; has0=len(nz)<len(E)
    L=1
    for e in nz: L=lcm(L,e)
    D=7*L; pts=set([0,D])
    for e in nz:
        st=L//e; pts.update(range(0,D+1,st))
    pts=sorted(pts); pn=[0]*8; b0=1 if has0 else 0
    for t1,t2 in zip(pts,pts[1:]):
        s=t1+t2; hit=b0
        for e in nz: hit|=1<<((7*e*s//(2*D))%7)
        pn[7-bin(hit).count("1")]+=t2-t1
    p=[F(x,D) for x in pn]; T=[sum(p[j:]) for j in range(8)]
    return 2*T[5]+4*T[6]
def norm(E):
    g=0
    for e in E: g=gcd(g,e)
    return tuple(e//g for e in E) if g>1 else tuple(E)
# focused: mod-7 offset r, internal set = all 9-subsets of {0..M} for small M -- exhaustive small
maxb=None;maxb_s=None
for r in range(1,7):
    for M in range(9,13):
        for js in combinations(range(0,M),9):
            E=norm(tuple(r+7*j for j in js))
            b=bunch(E)
            if maxb is None or b>maxb: maxb=b;maxb_s=E
print(f"focused mod-7 exhaustive (offset r=1..6, internal 9-of-[0..{12}]):")
print(f"  MAX BUNCH = {maxb} ~ {float(maxb):.6f} at {maxb_s}")
# is {1,8..57} it?
print(f"  bunch({{1,8..57}}) = {bunch(tuple(1+7*i for i in range(9)))} = {float(bunch(tuple(1+7*i for i in range(9)))):.6f}")
# corrected separation
minPOS=F(4717,882); FL=F(432,91)
print(f"\nCORRECTED SEPARATION:")
print(f"  minPOS (consec) = 4717/882 = {float(minPOS):.5f}")
print(f"  maxBUNCH (mod-7 pole) = {maxb} = {float(maxb):.5f}")
sep=minPOS-maxb
print(f"  minPOS - maxBUNCH = {sep} = {float(sep):.5f}")
print(f"  vs floor 432/91 = {float(FL):.5f}: margin {float(sep-FL):+.5f}  {'CLOSES' if sep>=FL else 'FAILS'}")
