"""Independent confirmation that the deep low-mu champions are real (exact vs high-res sample)."""
from fractions import Fraction as F
def mu_exact(E):
    E=sorted(set(int(e) for e in E)); k=len(E)
    if k<=1: return F(1) if k==1 else F(0)
    TH=F(2,7); bp=set([F(0),F(1)]); diffs=set()
    for i in range(k):
        for j in range(i+1,k): diffs.add(E[j]-E[i])
    for d in diffs:
        for m in range(0,d+1): bp.add(F(m,d))
    obp=sorted(b for b in bp if F(0)<=b<=F(1)); refined=set(obp)
    for a,b in zip(obp,obp[1:]):
        mid=(a+b)/2; floors={e:(e*mid).__floor__() for e in E}
        order=sorted(E,key=lambda e:e*mid-floors[e])
        for t in range(k):
            if t==k-1:
                ef,el=order[0],order[-1]; slope=ef-el; const=F(1)-floors[ef]+floors[el]
            else:
                eh,elo=order[t+1],order[t]; slope=eh-elo; const=-floors[eh]+floors[elo]
            if slope!=0:
                xb=(TH-const)/slope
                if a<xb<b: refined.add(xb)
    refined=sorted(refined); tot=F(0)
    for a,b in zip(refined,refined[1:]):
        mid=(a+b)/2; pts=sorted(set((e*mid)%1 for e in E))
        if len(pts)==1: mg=F(1)
        else:
            g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1]); mg=max(g)
        if mg>TH: tot+=(b-a)
    return tot
def mu_sample(E,N):
    E=sorted(set(E)); th=2.0/7.0; cnt=0
    for s in range(N):
        x=(s+0.5)/N; pts=sorted((e*x)%1.0 for e in E)
        mg=0.0; prev=pts[0]
        for t in range(1,len(pts)):
            g=pts[t]-prev
            if g>mg: mg=g
            prev=pts[t]
        w=pts[0]+1-pts[-1]
        if w>mg: mg=w
        if mg>th: cnt+=1
    return cnt/N
champions=[
 [0,7,13,14,15,16,17,19,21,24,27,29,40],
 [0,7,11,20,21,23,28,33,35,36,39,42,45],
 [0,11,21,28,33,35,36,37,39,42,45,49,62],
 [0,12,17,20,24,26,27,28,32,36,37,47,60],
]
for E in champions:
    ex=mu_exact(E); sa=mu_sample(E,2000000)
    print(f"E={E}\n   exact={ex}={float(ex):.6f}  sample(2e6)={sa:.6f}  diff={abs(float(ex)-sa):.2e}  <1/14? {ex<F(1,14)}",flush=True)
