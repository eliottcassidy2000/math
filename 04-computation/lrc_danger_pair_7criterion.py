from fractions import Fraction as F
import math
delta=F(1,14)
def danger_measure(T):
    bps=set([F(0),F(1)])
    for v in T:
        for k in range(0,v+1):
            for off in (delta,1-delta,F(0)):
                t=(k+off)/v
                if 0<=t<=1: bps.add(t)
    bps=sorted(bps);tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0:continue
        mid=(x0+x1)/2
        if all(((v*mid)%1<delta or (v*mid)%1>1-delta) for v in T):tot+=x1-x0
    return tot
bad=0
for a in range(1,20):
    for b in range(a+1,22):
        if math.gcd(a,b)!=1:continue
        m=danger_measure([a,b]); amp=m*49
        crit=(a%7==0 or b%7==0)
        if (amp==1)!=crit: bad+=1; print("MISMATCH",a,b,float(amp))
print("criterion '7|a or 7|b => independence' mismatches:",bad)
