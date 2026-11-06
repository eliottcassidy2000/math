"""Literal cover sets for the two union bounds and the outgoing-cover bijection."""
from itertools import permutations
import json

def need(q,m):
    if not q:raise RuntimeError(m)
def good(rows,p):return all(rows[x]>>y&1 for x,y in zip(p,p[1:]))
def paths(rows,vertices):return {p for p in permutations(vertices) if good(rows,p)}
def covers(rows,vertices):
    out=set()
    for p in permutations(vertices):
        for k in range(len(p)+1):
            a,b=p[:k],p[k:]
            if good(rows,a) and good(rows,b):out.add((a,b))
    return out
inputs=cover_checks=0
for n in range(1,5):
    for bits in range(1<<(n*(n-1)//2)):
        rows=[0]*n;k=0
        for i in range(n):
            for j in range(i+1,n):
                if bits>>k&1:rows[i]|=1<<j
                else:rows[j]|=1<<i
                k+=1
        H=paths(rows,tuple(range(n)));C=covers(rows,tuple(range(n)))
        for r in range(n):
            body=tuple(v for v in range(n) if v!=r);A=covers(rows,body)
            L={(a,b) for a,b in C if a and a[-1]==r}
            R={(a,b) for a,b in C if a and a[0]==r}
            singleton={((r,),p) for p in paths(rows,body)}
            after={(p[:p.index(r)+1],p[p.index(r)+1:]) for p in H}
            before={(p[p.index(r):],p[:p.index(r)]) for p in H}
            need(after|singleton <= L and before|singleton <= R,'literal union inclusion')
            need(len(after&singleton)==sum(p[0]==r for p in H),'after intersection is starts')
            need(len(before&singleton)==sum(p[-1]==r for p in H),'before intersection is ends')
            extra=R-(before|singleton)
            need(extra=={(a,b) for a,b in R if len(a)>=2 and b and rows[r]>>b[-1]&1},'exact f complement')
            deletion={(a[:-1],b) for a,b in L}
            need(len(deletion)==len(L) and deletion<=A,'L deletion is injection')
            complement={(a,b) for a,b in A if a and rows[r]>>a[-1]&1}
            need(A-deletion==complement,'exact outgoing L complement')
            outdegree=sum(rows[r]>>v&1 for v in range(n))
            if outdegree==1:
                need(not extra,'unique outneighbor forces f zero')
            cover_checks+=len(C)+len(A);inputs+=1
print(json.dumps(dict(status='PASS',marked_inputs=inputs,literal_cover_visits=cover_checks,
    checks=['both union injections and exact intersections','exact f complement','L deletion/complement bijection','unique-outneighbor f=0']),sort_keys=True,indent=2))
