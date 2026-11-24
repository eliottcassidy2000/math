from itertools import combinations, permutations
from math import isqrt

def build(V,test):
    return {v:{w for w in V if v!=w and test(v,w)} for v in V}

def cyc5():
    return build(range(5),lambda a,b:(a-b)%5 in (1,4))

def petersen():
    return build(tuple(combinations(range(5),2)),lambda a,b:not(set(a)&set(b)))

def hs():
    V=[(t,i,j) for t in range(2) for i in range(5) for j in range(5)]
    def test(v,w):
        t,i,j=v;s,k,l=w
        if t>s:return test(w,v)
        if t==s:return i==k and (j-l)%5 in ((1,4) if t==0 else (2,3))
        return l==(i*k+j)%5
    return build(V,test)

def order(G,d):
    r=min(G)
    N={v:G[v]|{v} for v in G}
    if d==2:return [r,*sorted(G[r]),*sorted(set(G)-N[r])]
    if d==3:
        s=min(G[r]);F={r,s};M=N[r]|N[s]
        return [*sorted(set(G)-M),*sorted(M-F),*sorted(F)]
    F=N[r]
    S=sorted(set(G)-F)+sorted(F)
    return S

def rows(G,d):
    m=len(G);S=order(G,d);N=[G[v]|{v} for v in S]
    r=[0]*(2*m)
    for i,j in combinations(range(m),2):
        r[i]|=1<<j
        r[m+j]|=1<<(m+i)
    for i,j in ((i,j) for i in range(m) for j in range(m)):
        if S[j] in N[i]:r[i]|=1<<(m+j)
        else:r[m+j]|=1<<i
    return S,r

def verify(G,d):
    m=len(G)
    assert m==d*d+1 and all(len(x)==d for x in G.values())
    S,r=rows(G,d)
    assert all(len(G[v]&G[w])==(w not in G[v]) for v,w in combinations(G,2))
    code=[sum(1<<j for j,w in enumerate(S) if w==v or w in G[v]) for v in S]
    for i,j in combinations(range(m),2):
        D=(code[i]^code[j]).bit_count()
        assert D==(2*d-2 if S[j] in G[S[i]] else 2*d)
    for i,j in combinations(range(2*m),2):
        assert (r[i]>>j&1)!=(r[j]>>i&1)
        a,b=(i+m)%(2*m),(j+m)%(2*m)
        assert (r[i]>>j&1)!=(r[a]>>b&1)
    spine=list(range(m))+list(range(2*m-1,m-1,-1))
    assert all(r[i]>>j&1 for i,j in zip(spine,spine[1:]))
    reach=[1<<i for i in range(2*m)]
    changed=True
    while changed:
        changed=False
        for i in range(2*m):
            new=reach[i]
            for j in range(2*m):
                if reach[i]>>j&1:new|=r[j]
            if new!=reach[i]:changed=True;reach[i]=new
    assert len(set(reach))==1 and reach[0]==(1<<(2*m))-1
    score=[x.bit_count() for x in r]
    if d==2:
        assert score.count(7)==score.count(2)==1
        for i,j in combinations(range(2*m),2):
            if score[i]==score[j]:
                f=lambda x:(r[0]>>x&1,r[x]>>m&1)
                assert f(i)!=f(j)
    else:
        F=list(range(2*d+2,m))
        for z in F:
            assert score.count(score[z])==score.count(score[m+z])==1
        f=lambda x:(all(r[x]>>z&1 for z in F),all(r[m+z]>>x&1 for z in F))
        for i,j in combinations(range(2*m),2):
            if score[i]==score[j]:assert f(i)!=f(j)
    print(d,m,2*m,"verified Moore code, strong blue observer, rigidity")

for G,d in ((cyc5(),2),(petersen(),3),(hs(),7)):verify(G,d)

a,b=((16,1,3,7,14),(20,1,2,7,14))
paths=list(permutations(range(5)))
for r in (a,b):
    assert sum(all(r[x]>>y&1 for x,y in zip(p,p[1:])) for p in paths)==9
    p=(0,4,3,2,1)
    assert all(r[x]>>y&1 for x,y in zip(p,p[1:]))
assert sum((x^y).bit_count() for x,y in zip(a,b))==2
assert not any(all((a[i]>>j&1)==(b[p[i]]>>p[j]&1)
                   for i in range(5) for j in range(5)) for p in paths)
for d,n in ((2,5),(3,6),(7,12),(57,83)):
    m=d*d+1
    assert (n-2)*(n-3)//2 < m <= (n-1)*(n-2)//2 and isqrt(m)**2!=m
