from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd=0,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd=qb,q
    return Fraction(bn,bd)
def saturated(v): return all(any(x%q==0 for x in v) for q in range(2,15))
# extend exhaustive SINGLE-SCALE saturated census to larger max, find min M and any below 1/12
for hi in (24,26):
    mins=None; minv=None; cnt=0; below12=[]
    for combo in combinations(range(1,hi+1),13):
        if not saturated(combo): continue
        cnt+=1; M=Mw(list(combo))
        if M<Fraction(1,12): below12.append((combo,M))
        if mins is None or M<mins: mins=M; minv=combo
    print('{1..%d}: %d saturated, min M=%s=%.5f, #below 1/12=%d'%(hi,cnt,str(mins),float(mins),len(below12)))
    for c,M in below12[:5]:
        print('   BELOW 1/12:',c,'M=%s=%.5f ratio=%.1f'%(str(M),float(M),c[-1]/c[0]))
    if hi==26: print('   extremal:',minv)
