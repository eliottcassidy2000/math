#!/usr/bin/env python3
"""Independent exact forest quotient and pullback audit; no source imports."""
import itertools as it,json,math,sys
from collections import defaultdict
from pathlib import Path

CHECKS=0
def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def classify(n,mask):
    edges=list(it.combinations(range(n),2)); adj=[[] for _ in range(n)]
    for j,(a,b) in enumerate(edges):
        if mask>>j&1:
            adj[a].append(b);adj[b].append(a)
            if len(adj[a])>2 or len(adj[b])>2:return None
    seen=set(); old=[];new=[];sgn=1;isol=[]
    for v in range(n):
        if not adj[v]:isol.append(v);continue
        if v in seen or len(adj[v])!=1:continue
        path=[v];prev=-1;cur=v
        while True:
            nxt=[w for w in adj[cur] if w!=prev]
            if not nxt:break
            prev,cur=cur,nxt[0];path.append(cur)
        seen.update(path)
        if path[0]>path[-1]:path.reverse()
        for x,y in zip(path,path[1:]):sgn*=1 if x<y else -1
        old.append((tuple(sorted(path)),path[0],path[-1]))
        new.append((path[0],path[-1],(len(path)-2)%2))
    if len(seen)+len(isol)!=n:return None
    return (tuple(sorted(old)),(tuple(isol),tuple(sorted(new))),sgn)

def groups(n):
    old=defaultdict(list);new=defaultdict(list)
    for mask in range(1<<(n*(n-1)//2)):
        ans=classify(n,mask)
        if ans is not None:
            a,b,s=ans;old[a].append((mask,s));new[b].append((mask,s))
    return old,new

def endpoint_count(n):
    b=[1]
    for j in range(1,n+1):b.append(sum(math.comb(j-1,k-1)*(1 if k==1 else math.comb(k,2))*b[j-k] for k in range(1,j+1)))
    return b[n]

def reduced_count(n):
    total=1
    for k in range(1,n//2+1):
        matching=math.factorial(2*k)//(2**k*math.factorial(k))
        for r in range(n-2*k+1):
            parity=sum(math.comb(k,j) for j in range(min(k,r)+1) if (r-j)%2==0)
            total+=math.comb(n,2*k)*matching*math.comb(n-2*k,r)*parity
    return total

def closure(n,gn,gn1):
    """Substitute every one-vertex cut directly in enlarged R polynomials."""
    eo=list(it.combinations(range(n),2));en=list(it.combinations(range(n+1),2))
    old_index={en.index(e):j for j,e in enumerate(eo)}
    new_index={en.index((v,n)):v for v in range(n)}
    old_lookup={m:(key,s) for key,ms in gn.items() for m,s in ms}
    verified=0
    for cut in range(1<<n):
        for large,monomials in gn1.items():
            poly=defaultdict(int)
            for mask,s in monomials:
                om=0;co=s
                for j in range(len(en)):
                    if not (mask>>j&1):continue
                    if j in old_index:om|=1<<old_index[j]
                    else:co*=1 if cut>>new_index[j]&1 else -1
                poly[om]+=co
            values={}
            for m,co in poly.items():
                if not co:continue
                require(m in old_lookup,'pullback support forest')
                key,s=old_lookup[m];v=co*s
                require(key not in values or values[key]==v,'pullback signed quotient consistency')
                values[key]=v
            for key,v in values.items():
                require(all(poly.get(m,0)*s==v for m,s in gn[key]),'pullback missing coefficients')
            verified+=1
    return verified

def main():
    gs={};census=[]
    for n in range(7):
        old,new=groups(n);gs[n]=new
        require(len(old)==endpoint_count(n),'endpoint census')
        require(len(new)==reduced_count(n),'parity boundary census')
        require(sum(len(x) for x in old.values())==sum(len(x) for x in new.values()),'forest partition')
        histogram=defaultdict(int)
        for key,g in old.items():histogram[classify(n,g[0][0])[1]]+=1
        census.append({'n':n,'edge_subsets':1<<(n*(n-1)//2),'forests':sum(map(len,new.values())),'endpoint_vertex_blocks':len(old),'boundary_parity_blocks':len(new),'old_coordinates_per_new':{str(k):sum(v==k for v in histogram.values()) for k in sorted(set(histogram.values()))}})
    cl=[{'old_n':n,'cuts':1<<n,'basis_pullbacks':closure(n,gs[n],gs[n+1])} for n in (2,3,4,5)]
    # The first strict loss, two assigned interior vertices.
    edges=list(it.combinations(range(6),2))
    def fm(es):return sum(1<<edges.index(tuple(sorted(e))) for e in es)
    f=fm([(0,4),(4,1),(2,5),(5,3)])
    g=fm([(0,5),(5,1),(2,4),(4,3)])
    cf,cg=classify(6,f),classify(6,g)
    require(cf[0]!=cg[0] and cf[1:]==cg[1:],'first hidden interior reassignment')
    out={'source_commit':'7b76044ca13e6684d1746923285b83dbf83828d2','checks':CHECKS,'census':census,'closures':cl,'b_sequence':[endpoint_count(n) for n in range(13)],'c_sequence':[reduced_count(n) for n in range(13)],'hostile':{'n':6,'edge_order':edges,'masks':[f,g],'old_signatures':[cf[0],cg[0]],'new_signature':cf[1],'canonical_signs':[cf[2],cg[2]]}}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
