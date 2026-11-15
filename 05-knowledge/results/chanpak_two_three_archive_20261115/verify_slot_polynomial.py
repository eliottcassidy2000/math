from itertools import combinations,permutations
from math import log,isclose
n=5
pairs=list(combinations(range(n),2))
perms=list(permutations(range(n)))
def need(q):
    if not q:raise RuntimeError("slot polynomial mismatch")
for code in range(1<<len(pairs)):
    edges={ab if code>>b&1 else ab[::-1] for b,ab in enumerate(pairs)}
    h=[1]+[0]*31
    for s in range(1,32):
        h[s]=sum(all(e in edges for e in zip(p,p[1:]))
                 for p in permutations(i for i in range(n) if s>>i&1))
    P=[0]*32
    for s in range(32):
        P[s]=sum(h[a]*h[a^s] for a in range(32) if not a & ~s)
    q2=P[31]
    q3=sum(h[a]*P[31^a] for a in range(32))
    c3=sum((a,b) in edges and (b,c) in edges and (c,a) in edges or
           (a,c) in edges and (c,b) in edges and (b,a) in edges
           for a,b,c in combinations(range(n),3))
    c5=sum(all(e in edges for e in zip(p, p[1:]+(p[0],))) for p in perms if p[0]==0)
    need(q2==2**5+2*c3*2**3+2*c5*2)
    need(q3==3**5+2*c3*3**3+2*c5*3)
    need((q2+q3)%5==0)
    L=log(3,2)
    # This float control is supplemental to the exact termwise proof.
    need(q3<=q2**L*(1+c3/2)**(1-L)+1e-8)
print("1024 independent ordered-slot and cycle-packing controls")
