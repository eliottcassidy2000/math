from itertools import permutations,combinations,product
from collections import Counter, defaultdict
V=tuple(range(5))
TRI=tuple(combinations(V,3)); ED=tuple(combinations(V,2))
def parity(p):return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2
def tclass(bits):
    # z_ij=1 for i->j; normalize one star
    z={e:False for e in ED}; z.update(zip(combinations(V[1:],2),bits))
    def f(a,b):return z[a,b] if a<b else not z[b,a]
    return tuple(f(a,b)^f(b,c)^f(c,a) for a,b,c in TRI)
def perm(x,p):
    d=dict(zip(TRI,x))
    return tuple(d[tuple(sorted((p[a],p[b],p[c])))]^parity((p[a],p[b],p[c])) for a,b,c in TRI)
G=tuple(p for p in permutations(V) if not parity(p))
classes={tclass(b) for b in product((False,True),repeat=6)}
left=classes.copy(); orbits=[]
while left:
    x=min(left); o={perm(x,p) for p in G}
    left-=o; orbits.append(o)
print('A5 orbit sizes',sorted(map(len,orbits)))
for o in orbits:
    print(len(o), 'converse orbit index',next(j for j,jj in enumerate(orbits) if tuple(not t for t in min(o)) in jj))
assert sorted(map(len,orbits))==[12,12,20,20]

def mask(edges):
    return sum(1<<ED.index(tuple(sorted(e))) for e in edges)
def cut(X):
    return mask((a,b) for a,b in ED if (a in X)!=(b in X))
def euler(x):
    d=set()
    for i,(a,b) in enumerate(ED):
        if x>>i&1:d.symmetric_difference_update((a,b))
    return x^cut(d)
K=(1<<10)-1
mov=Counter(euler(K^mask(zip(p,p[1:]))) for p in permutations(V) if p[0]<p[-1])
assert set(mov.values())=={4} and len(mov)==15 and all(x.bit_count()==4 for x in mov)
ev={euler(x) for x in range(1024) if x.bit_count()%2==0}
assert len(ev)==32 and (set(mov)|{0})=={x for x in ev if x.bit_count()%4==0}
assert {x^K for x in mov}.isdisjoint(mov) and K not in mov
print('32-vertex double-cover sheet: 15 Petersen moves, 4 unlabelled-direction paths per move')

P=tuple(combinations(V,2))
PE=tuple((a,b) for a,b in combinations(P,2) if not set(a)&set(b))
M=[S for S in combinations(PE,5) if len(set(sum(S,())))==10]
assert len(M)==6
assert all(len(set(a)&set(b))==1 for a,b in combinations(M,2))
for x in mov:
    active={v for i,e in enumerate(ED) if x>>i&1 for v in e}
    diag=tuple(e for i,e in enumerate(ED) if set(e)<=active and not x>>i&1)
    assert diag in PE
print('Six matchings, pairwise intersections give 15 K6 edges')

def count(a):
    return sum(all(a[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(V))
for i,o in enumerate(orbits):
    c=min(o)
    b=next(b for b in product((False,True),repeat=6) if tclass(b)==c)
    rb=dict(zip(combinations(V[1:],2),b))
    hs=Counter()
    for sx in range(16):
        a=[0]*5
        for u,v in ED:
            f=rb.get((u,v),False) ^ (sx>>u&1) ^ (sx>>v&1)
            if f:a[u]|=1<<v
            else:a[v]|=1<<u
        hs[count(a)]+=1
    print('orbit index',i,'count histogram',sorted(hs.items()))
    assert sorted(hs.items())==([(1,5),(9,5),(11,5),(15,1)] if len(o)==12 else [(3,3),(5,6),(9,3),(13,3),(15,1)])

def reflection(e):
    p=list(V)
    for a,b in e:p[a],p[b]=b,a
    return tuple(p)
def ord_product(p,q):
    t=tuple(p[q[i]] for i in V); u=tuple(t)
    k=1
    while u!=V:
        u=tuple(t[u[i]] for i in V);k+=1
    return k
def suppress(edges):
    adj={u:set() for u in P}
    for u,v in edges:adj[u].add(v);adj[v].add(u)
    while any(len(z)==2 for z in adj.values()):
        u=next(u for u in adj if len(adj[u])==2)
        a,b=adj.pop(u)
        assert b not in adj[a]
        adj[a].remove(u);adj[b].remove(u);adj[a].add(b);adj[b].add(a)
    return adj
kinds=Counter()
for e,f in combinations(PE,2):
    o=ord_product(reflection(e),reflection(f))
    if set(e)&set(f):
        assert o==3; kinds[3]+=1; continue
    a=suppress(set(PE)-{e,f})
    triangles=sum(v in a[u] and w in a[u] and w in a[v] for u,v,w in combinations(a,3))
    assert len(a)==6 and set(map(len,a.values()))=={3}
    assert (o,triangles) in {(2,0),(5,2)}
    kinds[o]+=1
assert kinds=={2:15,3:30,5:60}
print('Petersen edge-pair product orders',sorted(kinds.items()))
