import itertools, random, sys
from fractions import Fraction as F
from math import gcd
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass
def bps(E):
    E=sorted(set(E)); b={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): b.add(F(a,7*e))
    return sorted(x for x in b if 0<=x<=1)
def p0(E):
    E=sorted(set(E)); B=bps(E); s=F(0)
    for lo,hi in zip(B,B[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if not (set(range(1,7))-set(int((e*mid)%1*7) for e in E)): s+=hi-lo
    return s
def prim(E):
    E=sorted(set(e-min(E) for e in E)); g=0
    for e in E: g=gcd(g,e)
    if g>1: E=[e//g for e in E]
    return sorted(set(E))
CAPS={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
QB={8:0.19660,9:0.36210,10:0.44789,11:0.53125,12:0.60224}
MAXE=80
best={k:(F(-1),None) for k in CAPS}; overcap=[]; overqb=[]; nck={k:0 for k in CAPS}
def consider(E):
    E=prim(E); k=len(E)
    if k not in CAPS: return
    if E[0]!=0 or sorted(E)[-2]<=14 or max(E)>MAXE: return
    nck[k]+=1; v=p0(E)
    if v>best[k][0]: best[k]=(v,E)
    if v>CAPS[k]: overcap.append((k,E,float(v)))
    if float(v)>QB[k]+1e-9: overqb.append((k,E,float(v)-QB[k]))
random.seed(7)
for k in CAPS:
    for jump in (16,18,20,24,30,40,55,70): consider(list(range(k-1))+[k-2+jump])
    for s1 in range(1,k):
        s2=k-s1
        for M in (16,18,22,28,38,55): consider(list(range(s1))+[M+i for i in range(s2)])
    for parts in [(3,3,k-6),(2,2,k-4),(4,2,k-6),(1,1,k-2),(k-4,2,2),(2,k-4,2)]:
        if min(parts)<1 or sum(parts)!=k: continue
        a,b,c=parts
        for (M1,M2) in [(18,40),(16,35),(22,55)]:
            consider(list(range(a))+[M1+i for i in range(b)]+[M1+M2+i for i in range(c)])
    for M in (15,21,28,49):
        if k>=7: consider(list(range(7))+[M+7*j for j in range(k-7)])
        consider(list(range(k-3))+[M,M+1,M+3])
    s1=k//2; s2=k-s1
    for M in (16,20,28,40,55): consider(list(range(s1))+[M+x for x in range(s2)])
for k in CAPS:
    for _ in range(1200):
        nclu=random.randint(2,min(4,k-1))
        cuts=sorted(random.sample(range(1,k),nclu-1))
        sizes=[b-a for a,b in zip([0]+cuts,cuts+[k])]
        base=0; E=[]
        for s in sizes:
            M=random.choice([15,16,18,22,28,40,60]); start=base+M
            if s==1: clu=[start]
            else:
                span=random.randint(s-1,min(13,3*s))
                pts=sorted(random.sample(range(1,span+1),s-1)) if span>=s-1 else []
                clu=[start]+[start+p for p in pts] if len(pts)==s-1 else list(range(start,start+s))
            E+=clu; base=max(E)+1
        consider(E)
print("k  #wide   best_p0    cap      QB     >cap >QB")
for k in CAPS:
    bv,be=best[k]
    print(f"{k} {nck[k]:6d} {float(bv):.5f} {float(CAPS[k]):.5f} {QB[k]:.5f} {float(bv)>float(CAPS[k])} {float(bv)>QB[k]+1e-9}")
    print("    witness",be)
print("violations p0>cap:",len(overcap))
for r in overcap[:8]: print("   !!!",r)
overqb.sort(key=lambda t:-t[2])
print("rows over QB wide-bound:",len(overqb))
for r in overqb[:10]: print("   over-QB k=%d +%.5f"%(r[0],r[2]),r[1])
