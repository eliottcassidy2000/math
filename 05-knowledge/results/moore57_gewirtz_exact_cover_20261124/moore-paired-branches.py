"""Check and prepare a two-anchor exact-cover test for the Moore fiber gauge.

Run with no inputs for a complete Hoffman–Singleton witness replay.
When a Gewirtz sharp-set candidate becomes available, call paired_columns
with the candidate row, a selected g, and the complete proposed group.
"""
from itertools import combinations, permutations

def inverse(p):
    a=[0]*len(p)
    for i,j in enumerate(p):
        a[j]=i
    return tuple(a)

def mul(p,q):
    return tuple(p[t] for t in q)

def sharp(S):
    n=len(S)
    return all(len({p[x] for p in S})==n for x in range(n))

def paired_columns(S,g,K):
    """Return exact-cover options j -> (q, bitmask). One per slot, masks disjoint.
    Bits index the two separated Latin tables. The fixed maps occupy masks not
    present in any option; max coverage of the remaining masks gives sharpness.
    """
    n=len(S)
    Id=tuple(range(n))
    assert len(set(S))==n and Id in S and g in S and g!=Id and sharp(S)
    fixedT=(Id,inverse(g))
    fixedR=(Id,g)
    banT=[{p[x] for p in fixedT} for x in range(n)]
    banR=[{p[x] for p in fixedR} for x in range(n)]
    def pack(p,q):
        return sum(1<<(n*x+y) for x,y in enumerate(p)) | sum(1<<(n*n+n*x+y) for x,y in enumerate(q))
    pp=[p for p in S if p not in fixedR]
    pool=[(q,inverse(q)) for q in K if all(y not in banT[x] for x,y in enumerate(q))]
    options=[]
    for p in pp:
        col=[]
        for q,qi in pool:
            r=mul(qi,p)
            if all(y not in banR[x] for x,y in enumerate(r)):
                col.append((q,pack(q,r)))
        options.append(col)
    return pp,options

def hs():
    V=[(t,i,j) for t in range(2) for i in range(5) for j in range(5)]
    E={v:set() for v in V}
    for v,w in combinations(V,2):
        t,i,j=v; s,k,l=w
        if (t==s and i==k and (j-l)%5 in ({1,4} if t==0 else {2,3})) or (t==0 and s==1 and l==(i*k+j)%5):
            E[v].add(w);E[w].add(v)
    assert all(len(a)==7 for a in E.values())
    assert all(len(E[v]&E[w])==(w not in E[v]) for v,w in combinations(V,2))
    return E

def gauge(E,r):
    ns=sorted(E[r]); d=len(ns)
    L=[sorted(E[v]-{r}) for v in ns]
    for i in range(1,d):
        L[i]=[next(y for y in L[i] if y in E[x]) for x in L[0]]
    F={}
    for i in range(d):
        for j in range(d):
            if i!=j:
                F[i,j]=tuple(next(k for k,y in enumerate(L[j]) if y in E[x]) for x in L[i])
    return L,F

def main():
    E=hs()
    r=min(E)
    L,F=gauge(E,r)
    n=6
    I=tuple(range(n))
    cells=[(i,x) for i in range(1,n+1) for x in range(n)]
    adj={(i,x):{(j,F[i,j][x]) for j in range(1,n+1) if i!=j} for i,x in cells}
    for (i,x),(j,y) in combinations(cells,2):
        expected=(i!=j and x!=y and (j,y) not in adj[i,x])
        assert len(adj[i,x]&adj[j,y])==expected
    a,b=1,2
    S=[F[a,j] for j in range(n+1) if a!=j]
    g=F[a,b]
    assert F[a,0]==I and F[b,0]==I
    p,opts=paired_columns(S,g,list(permutations(range(n))))
    actual={F[a,j]:F[b,j] for j in range(3,n+1)}
    taken=0
    for t,choices in zip(p,opts):
        mask=next(m for q,m in choices if q==actual[t])
        assert not taken&mask
        taken|=mask
    assert taken.bit_count()==2*n*(n-2)
    print('HS(7):', [len(x) for x in opts], 'choices per second-anchor slot, verified known cover and grid identity')
if __name__=='__main__':
    main()
