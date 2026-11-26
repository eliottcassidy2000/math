from collections import Counter
from itertools import combinations, product
from pathlib import Path
import json, sys
from sympy import Poly, Symbol, oo, discriminant

def read(path):
    it=iter(map(int,Path(path).read_text().split()))
    n=next(it); rows=[]
    for _ in range(n):
        y,k=next(it),next(it)
        rows.append((y,[tuple(next(it) for _ in range(4)) for _ in range(k)]))
    marker=rows[0][1][0][1:]; assert n==53 and next(it,None) is None
    return [(y,[r for r in q if all(s!=t for s,t in zip(r[1:],marker))]) for y,q in rows[1:]]

def main(path):
    owners=[{} for _ in range(3)]; holes=[]; refs=[]
    for k,(y,q) in enumerate(read(path)):
        if y<0: holes.append(q);continue
        refs.append(k)
        for i,s in enumerate(next(r[1:] for r in q if r[0]==y)):
            assert s not in owners[i];owners[i][s]=k
    assert len(refs)==50 and list(map(len,holes))==[32,33]
    pool=[]
    for a,b in product(*holes):
        m=sum(x==y for x,y in zip(a[1:],b[1:]))
        hit=Counter(owners[i][s] for x in (a,b) for i,s in enumerate(x[1:]) if s in owners[i])
        pool.append((sum(hit.values())+m,m,hit))
    Z=Symbol("z")
    def polys(m):
        p=Counter(q for q,s,_ in pool if m is None or m==s)
        f=sum(v*Z**(k-min(p)) for k,v in p.items())
        return [str(f),int(Poly(f,Z).count_roots(-oo,oo)),int(discriminant(f,Z))]
    U=[polys(None),polys(0)]
    def result(n):
        lc=0; mini=100; defects=Counter(); floors=Counter()
        for D in combinations(refs,n):
            q=Counter(s-sum(h[i] for i in D) for s,_,h in pool)
            floors[min(q)]+=1
            delta=min((q[k]**2-q[k-1]*q[k+1],k) for k in range(min(q)+1,max(q)))
            if delta[0]<0:
                mini=min(mini,delta[0]);defects[delta[1]]+=1;lc+=1
        return {"releases":n,"nonlogconcave":lc,"defect_centers":dict(defects),"min_raw_defect":mini,"relaxed_floors":dict(floors)}
    return {"polynomials":U,"released_tests":[result(k) for k in range(3)]}
if __name__=="__main__":
    print(json.dumps(main(sys.argv[1]),indent=2))
