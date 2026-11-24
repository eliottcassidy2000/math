from itertools import combinations
from functools import cache

def graph(V,adj):
    return [{j for j,w in enumerate(V) if i!=j and adj(v,w)}
            for i,v in enumerate(V)]

def pentagon():
    return graph(list(range(5)),lambda v,w:(v-w)%5 in (1,4))

def petersen():
    V=list(map(frozenset,combinations(range(5),2)))
    return graph(V,lambda v,w:not(v&w))

def HS():
    V=[(t,i,j) for t in range(2) for i in range(5) for j in range(5)]
    def adj(v,w):
        t,i,j=v;s,k,l=w
        if t>s:return adj(w,v)
        if t==s:return i==k and (j-l)%5 in ((1,4) if t==0 else (2,3))
        return l==(i*k+j)%5
    return V,graph(V,adj)

def data(m,d,tri=0,quad=0):
    k=m-1-d
    ana=m*d*(d-1)-6*tri
    nan=m*d*(m-2*d)+6*tri
    an2=m*m*d*d-2*d*m*(d+d*d)+2*m*d*d+12*tri+8*quad
    return m*d*k*k-nan,m*d*d*k*k-k*ana-d*nan-an2

def motifs(G):
    m=len(G);V=set(range(m))
    C=[V-{v}-G[v] for v in V]
    F4=sum(len(C[z]-{x,y}) for x in V for y in C[x] for z in G[y] if z!=x)
    F5=sum(len(C[w]-{x,y,z}) for x in V for y in G[x] for z in C[y]
           if z!=x for w in G[z] if w not in (x,y,z))
    return F4,F5

def small_hamilton(G):
    n=len(G)
    def rooted(start,closed):
        @cache
        def rec(mask,j):
            if mask==(1<<n)-1:return int(not closed or start in G[j])
            return sum(rec(mask|1<<k,k) for k in G[j] if not(mask>>k&1))
        return rec(1<<start,start)
    return sum(rooted(i,False) for i in range(n)),rooted(0,True)

V,H=HS()
for G,d in [(pentagon(),2),(petersen(),3),(H,7)]:
    m=len(G)
    assert all(len(a)==d for a in G)
    assert all(len(G[x]&G[y])==(y not in G[x]) for x,y in combinations(range(m),2))
    assert motifs(G)==data(m,d)
    print('case',d,'F4,F5',motifs(G))
    if m<20:
        C=[set(range(m))-{x}-G[x] for x in range(m)]
        print('graph/complement oriented Hamilton paths and cycles',small_hamilton(G),small_hamilton(C))
print('degree57 conditional targets',data(3250,57))

p=[];j=0
for i in range(5):
    p.extend((0,i,(j+s)%5) for s in range(5))
    u=(i*i+j-1)%5
    p.extend((1,i,(u+2*s)%5) for s in range(5))
    j=(j-3-i)%5
assert j==0 and len(set(p))==50
assert all(V.index(b) in H[V.index(a)] for a,b in zip(p,p[1:]+p[:1]))
print('explicit Hamilton cycle in degree7 verified')

# Test non-Moore d-regular examples and count obstructions independently.
for m,steps in [(6,(1,2,4,5)),(8,(1,3,5,7)),(9,(1,2,7,8))]:
    G=graph(list(range(m)),lambda v,w:(v-w)%m in steps)
    t=sum(b in G[a] and c in G[a] and c in G[b] for a,b,c in combinations(range(m),3))
    q=sum(len(G[a]&G[b])*(len(G[a]&G[b])-1)//2
          for a,b in combinations(range(m),2))//2
    assert motifs(G)==data(m,len(steps),t,q)
