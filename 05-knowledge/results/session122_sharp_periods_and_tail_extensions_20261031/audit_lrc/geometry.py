"""Independent full-comb geometry; no producer or copied geometry imports."""
from fractions import Fraction as F
from functools import lru_cache
R=F(1,14)
def need(ok,label):
    if not ok: raise RuntimeError(label)
@lru_cache(None)
def comb(v,r=R):
    need(isinstance(v,int) and v>0,'positive integer speed')
    return tuple(((j+r)/v,(j+1-r)/v) for j in range(v))
def intersection(a,b):
    i=j=0; out=[]
    while i<len(a) and j<len(b):
        l=max(a[i][0],b[j][0]); h=min(a[i][1],b[j][1])
        if l<=h:
            if out and l<=out[-1][1]:out[-1]=(out[-1][0],max(h,out[-1][1]))
            else:out.append((l,h))
        if a[i][1]<b[j][1]:i+=1
        elif a[i][1]>b[j][1]:j+=1
        else:i+=1;j+=1
    return tuple(out)
def safe(speeds,r=R):
    out=((F(0),F(1)),)
    for v in sorted(set(speeds)):out=intersection(out,comb(v,r))
    return out
def width(intervals):return max((h-l for l,h in intervals),default=F(0))
def check_witness(t,speeds,r=R):
    need(all(r <= (v*t)%1 <= 1-r for v in speeds),'original-speed witness')
