from itertools import combinations,product
from collections import Counter,defaultdict

def graph(V,test):
    return [set(j for j,w in enumerate(V) if i!=j and test(v,w)) for i,v in enumerate(V)]

def petersen():
    V=list(map(frozenset,combinations(range(5),2)))
    return graph(V,lambda x,y:not(x&y))

def hs():
    V=[(t,i,j) for t in range(2) for i in range(5) for j in range(5)]
    def test(v,w):
        t,i,j=v;s,k,l=w
        if t>s:return test(w,v)
        if t==s:return i==k and (j-l)%5 in ((1,4) if t==0 else (2,3))
        return l==(i*k+j)%5
    return graph(V,test)

def automorphisms(G,pairs):
    """Exact simultaneous individualize/refine isomorphism search."""
    n=len(G); V=list(range(n))
    def extend(F):
        a=[0]*n;b=[0]*n
        for c,(x,y) in enumerate(F,1):a[x]=c;b[y]=c
        while True:
            c={}; keys=[]
            for z in (a,b):
                keys.append([c.setdefault((z[i],tuple(sorted(z[j] for j in G[i]))),len(c)) for i in V])
            u,v=keys
            sizes=Counter(u)
            if sizes!=Counter(v):return
            if len(sizes)==len(set(a)):break
            a,b=u,v
        if len(sizes)==n:
            q={k:i for i,k in enumerate(v)}
            p=tuple(q[k] for k in u)
            assert all((p[i] in G[p[j]])==(i in G[j]) for i,j in combinations(V,2))
            yield p
            return
        cls=min((k for k,v in sizes.items() if v>1),key=lambda k:(sizes[k],k))
        i=u.index(cls)
        for j in V:
            if v[j]==cls:
                yield from extend(F+[(i,j)])
    return extend(list(pairs))

def part(p):
    s=set(range(len(p)));cc=[]
    while s:
        a=min(s);b=a;k=0
        while b in s:
            k+=1;s.remove(b);b=p[b]
        cc.append(k)
    return tuple(sorted(cc))

def data(G):
    r=0;v=min(G[r])
    rows=sorted(G[r]-{v});cols=sorted(G[v]-{r})
    cell=[]
    for x in rows:
        cells=[]
        for y in cols:
            z,=(G[x]&G[y])
            cells.append(z)
        cell.append(cells)
    assert len({z for C in cell for z in C})==len(rows)**2
    P=list(automorphisms(G,[(r,r),(v,v)]))
    flip=list(automorphisms(G,[(r,v),(v,r)]))
    K={(tuple(rows.index(p[x]) for x in rows),tuple(cols.index(p[y]) for y in cols)) for p in P}
    for p in P:
        for i,x in enumerate(rows):
            for j,y in enumerate(cols):
                assert p[cell[i][j]]==cell[rows.index(p[x])][cols.index(p[y])]
    for p in flip:
        for i,x in enumerate(rows):
            for j,y in enumerate(cols):
                assert p[cell[i][j]]==cell[rows.index(p[y])][cols.index(p[x])]
    corr=Counter((part(s),part(t)) for s,t in K)
    fixstat=Counter((sum(i==x for i,x in enumerate(s)),sum(i==x for i,x in enumerate(t))) for s,t in K)
    print('degree',len(G[r]),'edge',r,v,'ordered',len(P),'swapping',len(flip),
       'projections',len({s for s,t in K}),len({t for s,t in K}),flush=True)
    kinds=Counter((part(tuple(cols.index(p[x]) for x in rows)),part(tuple(rows.index(p[y]) for y in cols))) for p in flip)
    print('fixed-axis counts',sorted(fixstat.items()),flush=True)
    print('paired cycle types',sorted(corr.items()),flush=True)
    fliporders=Counter(part(p) for p in flip)
    print('swap involutions',sum(all(p[p[i]]==i for i in range(len(G))) for p in flip), 'hop full cycle shapes',sorted(fliporders.items()),flush=True)
    return K,corr
K,corr=data(petersen()); assert len(K)==4
K,corr=data(hs()); assert len(K)==720
assert corr[(1,1,1,1,2),(2,2,2)]==15

G=hs();r=0;v=min(G[r]);R=sorted(G[r]-{v})
swap=next(p for p in automorphisms(G,[(r,v),(v,r)])
          if all(p[p[i]]==i for i in range(len(G))))
C=[swap[x] for x in R]
F={tuple(R.index(p[x]) for x in R):tuple(C.index(p[y]) for y in C)
   for p in automorphisms(G,[(r,r),(v,v)])}
assert all(F[F[p]]==p for p in F)
assert (R,C)==([4,25,30,35,40,45],[2,26,31,36,41,46])
expected=[(1,0,5,4,3,2),(4,5,3,2,0,1),(5,2,1,4,3,0),
          (1,0,3,2,5,4),(2,5,0,4,3,1)]
for i,w in enumerate(expected):
    p=list(range(6));p[i],p[i+1]=p[i+1],p[i]
    assert F[tuple(p)]==w
print("HS aligned axes and involutory outer generators verified:",R,C,expected)
for G in (petersen(),hs()):
    r=0;v=min(G[r])
    assert all(next(automorphisms(G,[(r,x),(v,y)]),None) is not None
               for x in range(len(G)) for y in G[x])
    print('confirmed directed-edge transitivity',len(G),len(G)*len(G[0]))

def recover_pairs(G):
    from itertools import permutations
    r=0;v=min(G[r]);X=sorted(G[r]-{v});Y=sorted(G[v]-{r})
    b=len(X);I=tuple(range(b))
    C=[[next(iter(G[x]&G[y])) for y in Y] for x in X]
    F={}
    for i in I:
        F[i,i]=I
        for j in I:
            if i!=j:
                F[i,j]=tuple(next(y for y,t in enumerate(C[j]) if t in G[s])
                               for s in C[i])
    S=set()
    for p in permutations(range(b)):
        for target0 in I:
            q=[-1]*b
            for j in I:q[F[0,j][0]]=F[p[0],p[j]][target0]
            if all(tuple(q[x] for x in F[i,j])==tuple(F[p[i],p[j]][x] for x in q)
                   for i,j in combinations(range(b),2)):
                S.add((p,tuple(q)))
    T={(tuple(X.index(f[x]) for x in X),tuple(Y.index(f[y]) for y in Y))
       for f in automorphisms(G,[(r,r),(v,v)])}
    assert S==T
    print('second algorithm: sharp first row determined opposite action:',len(S))
recover_pairs(petersen());recover_pairs(hs())
