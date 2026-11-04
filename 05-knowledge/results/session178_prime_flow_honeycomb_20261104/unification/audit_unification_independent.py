"""Independent structural audit using random loop/multiedge graphs and exact paths."""
import random
from fractions import Fraction

def need(ok):
    if not ok: raise RuntimeError('independent audit failed')
rng=random.Random(178)
for case in range(1000):
    n=rng.randrange(1,15)
    # First a spanning tree; add arbitrary oriented parallel edges and loops.
    es=[(rng.randrange(v),v) for v in range(1,n)]
    es=[e if rng.randrange(2) else e[::-1] for e in es]
    es += [(rng.randrange(n),rng.randrange(n)) for _ in range(rng.randrange(30))]
    obs=[rng.randrange(n) for _ in range(n)]
    anchors={}
    for v,y in enumerate(obs): anchors.setdefault(y,v)
    for v,y in enumerate(obs):
        if v==anchors[y]:continue
        root=anchors[y];prev={root:None};queue=[root]
        for u in queue:
            for k,(a,b) in enumerate(es):
                if a==u and b not in prev:prev[b]=(u,k,1);queue.append(b)
                if b==u and a not in prev:prev[a]=(u,k,-1);queue.append(a)
        z=[0]*len(es);u=v
        while prev[u] is not None:
            u,k,sgn=prev[u];z[k]+=sgn
        debt=[0]*n
        for (a,b),c in zip(es,z):debt[a]-=c;debt[b]+=c
        need(debt==[int(i==v)-int(i==root) for i in range(n)])
        collapsed={}
        for i,d in enumerate(debt):collapsed[obs[i]]=collapsed.get(obs[i],0)+d
        need(not any(collapsed.values()))
        # Independent potential-telescoping check for signed random chains.
    z=[rng.randrange(-4,5) for _ in es];h=[rng.randrange(-10,11) for _ in range(n)]
    debt=[0]*n
    for (a,b),c in zip(es,z):debt[a]-=c;debt[b]+=c
    need(sum(c*(h[b]-h[a]) for (a,b),c in zip(es,z))==sum(x*y for x,y in zip(h,debt)))
    ratios=Fraction(1)
    for (a,b),c in zip(es,z): ratios*=Fraction(b+1,a+1)**c
    need(ratios==__import__('functools').reduce(lambda a,t:a*Fraction(t[0]+1)**t[1], enumerate(debt),Fraction(1)))
print('PASS: 1000 independent random connected multigraphs; all fibre-anchor path witnesses; signed additive and rational multiplicative identities')
