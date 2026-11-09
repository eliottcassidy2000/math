# exhaustive small-order data with flags for the seed/good analysis.
# per connected even class: n d U pend_any pendNU(count, capped 2) iso r3
import sys
from collections import defaultdict
NB=int(sys.argv[2])
comps=defaultdict(set)
for line in open(sys.argv[1]):
    p=line.split(); n=int(p[0]); d=int(p[2])
    if p[3]!='E' or n>NB: continue
    edges=p[5] if len(p)>5 else ''
    adj=[set() for _ in range(n)]
    for i in range(0,len(edges),2):
        u=int(edges[i],16); v=int(edges[i+1],16); adj[u].add(v); adj[v].add(u)
    deg=[len(a) for a in adj]
    U = (n==1) or any(x==n-1 for x in deg)
    # pendant vertices and whether neighbour universal
    pend=[v for v in range(n) if deg[v]==1]
    pnu=sum(1 for v in pend if deg[next(iter(adj[v]))] < n-1)
    pany=len(pend)
    iso=sum(1 for v in range(n) if deg[v]==0)
    r3=any(x==2 for x in deg)
    comps[n].add((d,int(U),min(pany,2),min(pnu,2),min(iso,2),int(r3)))
# unions: flags combine; a pendant's neighbour is non-universal in a union of >=2 components
res=defaultdict(set)
for n,L in comps.items(): res[n]|=L
changed=True
while changed:
    changed=False
    new=defaultdict(set)
    for n1,L1 in list(res.items()):
        for n2,L2 in comps.items():
            if n1+n2>NB: continue
            for a in L1:
                for b in L2:
                    pany=min(a[2]+b[2],2)
                    new[n1+n2].add((a[0]*b[0],0,pany,pany,min(a[4]+b[4],2),max(a[5],b[5])))
    for n,L in new.items():
        before=len(res[n]); res[n]|=L
        if len(res[n])!=before: changed=True
for n in sorted(res):
    for t in sorted(res[n]): print(n,*t)
