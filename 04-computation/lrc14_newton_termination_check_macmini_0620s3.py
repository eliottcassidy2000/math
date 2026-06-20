import itertools
from fractions import Fraction as F
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def DeltaS(B,S):
    S=list(S); tot=F(0)
    for r in range(len(S)+1):
        for T in itertools.combinations(S,r):
            tot+=(-1)**(len(S)-len(T))*measS7(list(B)+list(T))
    return tot
# my functional, 0 in B. Test |S|=7 and |S|=8 far runners.
B=(0,1,2,3)
for S in [(15,16,17,18,19,20,21),(15,17,19,23,29,31,37)]:
    print(f"|S|={len(S)}: Delta_S(B={B}, S={S}) = {DeltaS(B,S)} = {float(DeltaS(B,S)):+.6f}")
# also check the SUM over all S subset F equals p0 (Newton identity, finite 2^r terms)
F7=(17,19,23,29,31,37,41)
direct=measS7(list(B)+list(F7))
newton=sum(DeltaS(B,S) for r in range(len(F7)+1) for S in itertools.combinations(F7,r))
print(f"\nNewton identity (r=7): p0 direct={float(direct):.6f}  sum_S Delta_S={float(newton):.6f}  match={direct==newton}")
# magnitude of |S|=7 term relative to |S|=2
d7=DeltaS(B,F7); d2=max(abs(DeltaS(B,p)) for p in itertools.combinations(F7,2))
print(f"|Delta_F (|S|=7)|={float(abs(d7)):.6f}  max|Delta (|S|=2)|={float(d2):.6f}  ratio={float(abs(d7)/d2) if d2 else 0:.4f}")
