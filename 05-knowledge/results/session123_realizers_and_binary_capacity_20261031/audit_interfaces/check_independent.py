from itertools import combinations, product
from fractions import Fraction as Q
import json

def require(b,msg):
    if not b: raise RuntimeError(msg)

def matching(edges):
    if not edges:return 0
    a,b=edges[0]
    return max(matching(edges[1:]),1+matching([e for e in edges[1:] if a not in e and b not in e]))

def forest(edges,n):
    p=list(range(n))
    def root(a):
        while p[a]!=a:a=p[a]
        return a
    for a,b in edges:
        a,b=root(a),root(b)
        if a==b:return False
        p[a]=b
    return True

counts=dict(forests=0,protections=0,forest_assignments=0,quadratics=0,compiler_instances=0,compiler_vectors=0)
# Independent constructive leaf elimination on every labelled forest on <=5 vertices.
for n in range(2,6):
    all_e=list(combinations(range(n),2))
    for mask in range(1,1<<len(all_e)):
        edges=[e for i,e in enumerate(all_e) if mask>>i&1]
        if not forest(edges,n):continue
        counts['forests']+=1
        for pm in range(1<<len(edges)):
            prot=[e for i,e in enumerate(edges) if pm>>i&1]
            un=[e for i,e in enumerate(edges) if not pm>>i&1]
            remaining=un[:]; factors=[]; xorcost=0
            while remaining:
                deg=[sum(v in e for e in remaining) for v in range(n)]
                x=next(i for i,d in enumerate(deg) if d==1)
                y=next(b if a==x else a for a,b in remaining if x in (a,b))
                neigh=[b if a==y else a for a,b in remaining if y in (a,b)]
                factors.append((y,neigh));xorcost+=len(neigh)-1
                remaining=[e for e in remaining if x not in e and y not in e]
            nu=matching(un)
            require(len(factors)==nu,'leaf matching')
            require(xorcost==len(un)-nu,'inner xor count')
            count=len(prot)+xorcost+len(factors)+len(prot)+len(factors)-1
            require(2*len(edges)-1-count==len(un)-nu,'gate count')
            cover=all((pm>>i&1) or (pm>>j&1) for i,j in combinations(range(len(edges)),2) if set(edges[i])&set(edges[j]))
            require(cover==(len(un)==nu),'line graph cover')
            counts['protections']+=1
            for bits in product(range(2),repeat=n):
                original=sum(bits[a]*bits[b] for a,b in edges)%2
                repaired=(sum(bits[a]*bits[b] for a,b in prot)+sum(bits[y]*(sum(bits[x] for x in ns)%2) for y,ns in factors))%2
                require(original==repaired,'forest root')
                # All protected products are literally retained, not reconstructed.
                counts['forest_assignments']+=1
# All scalar quadratics on 3 variables, all boundary splits; selected pairs.
n=3;pairs=list(combinations(range(n),2))
polys=[]
for coefficients in product(range(2),repeat=1+n+len(pairs)):
    def eval_at(x,c=coefficients):
        return (c[0]+sum(c[1+i]*((x>>i)&1) for i in range(n))+sum(c[1+n+j]*((x>>a)&1)*((x>>b)&1) for j,(a,b) in enumerate(pairs)))%2
    vals=[eval_at(x) for x in range(1<<n)]
    polys.append(vals)
for p in range(n+1):
    m=n-p
    for idx,f in enumerate(polys):
        for g in [None,polys[(37*idx+11)%len(polys)]]:
            fs=[f] if g is None else [f,g]
            truth={h for h in range(1<<m) if all(v[x^(h<<p)]==v[x] for v in fs for x in range(1<<n))}
            # Polar zero is tested on basis vectors; then constant derivatives removed.
            rad={h for h in range(1<<m) if all(v[(h<<p)^(1<<i)]^v[h<<p]^v[1<<i]^v[0]==0 for v in fs for i in range(n))}
            kernel={h for h in rad if all(v[h<<p]==v[0] for v in fs)}
            require(truth==kernel,'polar-kernel formula')
            require(all((a^b) in kernel for a in kernel for b in kernel),'kernel linear')
            counts['quadratics']+=1
# Deliberately incorrect heuristics must fail.
triangle=[sum(((x>>a)&1)*((x>>b)&1) for a,b in pairs)%2 for x in range(8)]
require(all(triangle[7^(1<<i)]^triangle[7]^triangle[1<<i]^triangle[0]==0 for i in range(3)) and triangle[7]!=triangle[0],'derivative-one hostile')
# Independent endpoint-cell scan vs matrix feasible floor cells, including exact 1/2.
for n in range(1,4):
    for speeds in product(range(1,6),repeat=n):
        for radius in [Q(1,2),Q(1,3),Q(1,4),Q(2,5)]:
            p,q=radius.numerator,radius.denominator
            def safe(t):
                return all(radius<=v*t-(v*t).__floor__()<=1-radius for v in speeds)
            endpoints=sorted({Q(0),Q(1)}|{(Q(z)+r)/v for v in speeds for z in range(v) for r in [radius,1-radius]})
            direct_points={t for t in endpoints if safe(t)}
            cells=[]
            for z in product(*(range(v) for v in speeds)):
                matrix=all(q*speeds[j]*z[i]-q*speeds[i]*z[j]<=(q-p)*speeds[i]-p*speeds[j] for i in range(n) for j in range(n))
                lo=max((Q(z[i])+radius)/speeds[i] for i in range(n))
                hi=min((Q(z[i])+1-radius)/speeds[i] for i in range(n))
                require(matrix==(lo<=hi),'compiler equivalence')
                if matrix:
                    cells.append((lo,hi));require(safe(lo) and safe(hi) and safe((lo+hi)/2),'cell safety')
                counts['compiler_vectors']+=1
            require(bool(cells)==bool(direct_points),'endpoint scan feasibility')
            require(all(any(lo<=t<=hi for lo,hi in cells) for t in direct_points),'all boundary witnesses retained')
            counts['compiler_instances']+=1
print(json.dumps({'status':'PASS','counts':counts,'arithmetic':'exact rational; no producer imports'},sort_keys=True))
