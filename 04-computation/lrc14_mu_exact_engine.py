from fractions import Fraction as F

def _frac(q):
    return q - q.__floor__()

def mu_exact(E):
    E = sorted(set(E))
    k = len(E)
    diffs = set()
    for i in range(k):
        for j in range(k):
            d = E[i]-E[j]
            if d>0: diffs.add(d)
    bps = set([F(0),F(1)])
    for d in diffs:
        for t in range(0, d+1):
            x = F(t,d)
            if 0<=x<=1: bps.add(x)
    bps = sorted(bps)
    total = F(0)
    g0 = F(2,7)
    for a,b in zip(bps, bps[1:]):
        if a==b: continue
        mid = (a+b)/2
        fr = [_frac(F(E[i])*mid) for i in range(k)]
        order = sorted(range(k), key=lambda i: fr[i])
        n = [(F(E[i])*mid).__floor__() for i in range(k)]
        cross = set([a,b]); m=k
        for r in range(m):
            i1=order[r]; i2=order[(r+1)%m]
            wrap = 1 if r==m-1 else 0
            slope=E[i2]-E[i1]; const=-n[i2]+n[i1]+wrap
            if slope!=0:
                xc=(g0-const)/slope
                if a<xc<b: cross.add(xc)
        cross=sorted(cross)
        for u,v in zip(cross,cross[1:]):
            if u==v: continue
            mm=(u+v)/2
            P=[F(E[i])*mm-n[i] for i in range(k)]
            Ps=sorted(P)
            gaps=[Ps[r+1]-Ps[r] for r in range(m-1)]+[Ps[0]+1-Ps[-1]]
            if max(gaps)>g0: total+=(v-u)
    return total
