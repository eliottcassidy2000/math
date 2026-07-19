from fractions import Fraction as F
from itertools import combinations
LO=F(1,14); HI=F(3,41)
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def Mfull(V):
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0); arg=None
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def act(V,t,g):
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))
        if 1-r==g: dn.append((v,a+1))
    return up,dn
print(f"SEARCHING (1/14, 3/41) = ({float(LO):.8f}, {float(HI):.8f})")
print("  shapes: {1..11,13,x}, {1..12,x}, {1..10,12,13,x}, {1..11,14,x}   x <= 400")
hits=[]; n=0
shapes=[("{1..11,13,x}", list(range(1,12))+[13]),
        ("{1..12,x}",    list(range(1,13))),
        ("{1..10,12,13,x}", list(range(1,11))+[12,13]),
        ("{1..11,14,x}", list(range(1,12))+[14])]
for nm,base in shapes:
    for x in range(2,401):
        V=sorted(set(base+[x]))
        if len(V)!=13: continue
        n+=1
        M,t=Mfull(V)
        if LO<M<HI:
            up,dn=act(V,t,M)
            if up and dn:
                vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai; s=vi+vj
            else: D=s=None
            hits.append((M,nm,x,D,s,14*D-s if D else None))
print(f"  {n} families scanned; hits strictly inside: {len(hits)}")
seen=set()
for M,nm,x,D,s,sl in sorted(hits)[:12]:
    if M in seen: continue
    seen.add(M)
    print(f"    M={str(M):9s} ({float(M):.8f})  {nm} x={x}   D={D} s={s} slack={sl}")
if hits:
    mn=min(h[0] for h in hits)
    print(f"  SMALLEST value found strictly above 1/14: {mn} = {float(mn):.8f}")
else:
    print("  none -- (1/14, 3/41) appears empty for these shapes")
