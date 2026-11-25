from itertools import product,combinations
from collections import Counter
from sys import argv

def mul(a,b):
    c=(a if b&1 else 0)^((a<<1) if b&2 else 0)
    return c^(7 if c&4 else 0)
M=[[mul(a,b) for b in range(4)] for a in range(4)]
def norm(x):
    a=next(v for v in x if v)
    return tuple(M[a][M[a][v]] for v in x)
P=sorted({norm(x) for x in product(range(4),repeat=3) if any(x)})
idx={v:i for i,v in enumerate(P)}
N=[[not (M[a][d]^M[b][e]^M[c][f]) for d,e,f in P] for a,b,c in P]
arcs=[]; line_bits=[sum((1<<i) for i in range(21) if row[i]) for row in N]
for C in combinations(range(21),6):
    mask=sum(1<<a for a in C)
    if all((mask&m).bit_count()<=2 for m in line_bits):
        arcs.append(mask)
assert len(arcs)==168

def trans(i,j,t):
    out=[]
    for x in P:
        y=list(x);y[i]^=M[t][y[j]]
        out.append(idx[norm(y)])
    return tuple(out)
def apply(p,mask):
    return sum(1<<p[j] for j in range(21) if mask>>j&1)
def pl_perm(p):
    out=[apply(p,h) for h in orbit]
    return bytes(ind[o] for o in out)
gens=[trans(i,j,t) for i,j in [(0,1),(1,2),(2,0)] for t in [1,2]]
orb={arcs[0]}; stack=[arcs[0]]
for s in stack:
    for g in gens:
        h=apply(g,s)
        if h not in orb: orb.add(h);stack.append(h)
assert len(orb)==56
orbit=sorted(orb);ind={h:i for i,h in enumerate(orbit)}
frob=tuple(idx[tuple(M[a][a] for a in x)] for x in P)
diag=tuple(idx[norm((M[2][a],b,c))] for a,b,c in P)
def dual_mask(h):
    return sum(1<<j for j,l in enumerate(line_bits) if h&l==0)
ds=[dual_mask(h) for h in orbit]
for _ in range(3):
    if set(ds)==orb:break
    ds=[apply(diag,s) for s in ds]
else:raise AssertionError('duality orbit')
du=bytes(ind[d] for d in ds)
assert sorted(du)==list(range(56))
Adj=[sum(1<<j for j,h in enumerate(orbit) if (g&h)==0) for g in orbit]
for i in range(56):
    assert Adj[i].bit_count()==10
    for j in range(i):
        assert (Adj[i]&Adj[j]).bit_count()==(0 if Adj[i]>>j&1 else 2)
for p in [du,pl_perm(frob)]:
    assert all(((Adj[i]>>j)&1)==((Adj[p[i]]>>p[j])&1) for i in range(56) for j in range(56))
G=[pl_perm(g) for g in gens]

def gbfs(gen):
    I=bytes(range(56));S=[I];cache={I}
    for s in S:
        for t in gen:
            u=s.translate(t+bytes(range(56,256)))
            if u not in cache:cache.add(u);S.append(u)
    return S
G1=gbfs(G+[pl_perm(frob)])
G2=gbfs(G+[pl_perm(frob),du])
assert (len(G1),len(G2))==(40320,80640)
assert len({frozenset(apply(diag,h) for h in orb),frozenset(apply(diag,apply(diag,h)) for h in orb),frozenset(orb)})==3
print('group sizes',len(G1),len(G2),flush=True)

def meet_families():
    point=[{h for h in orbit if h>>i&1} for i in range(21)]
    line=[{h for h in orbit if h&m==0} for m in line_bits]
    assert all(len(s)==16 for s in point+line)
    for i in range(21):
        for j in range(21):
            assert len(point[i]&line[j])==(0 if N[i][j] else 6)
            assert len(point[i]&point[j])==(16 if i==j else 4)
    return point,line

def mod2_shadow(G):
    basis={}
    def reduce(v):
        while v:
            k=v.bit_length()
            if k not in basis:break
            v ^= basis[k]
        return v
    J=(1<<(56**2))-1
    for i,p in enumerate(G):
        v=reduce(sum(1<<(56*j+k) for j,k in enumerate(p)))
        if v: basis[v.bit_length()]=v
        if v and reduce(J)==0:return True,len(basis),i+1
    return False,len(basis),len(G)

def mod3_shadow(G):
    # Pack coefficient +1 and -1 at disjoint places.
    basis={}
    def add(a,b,c,d):
        z=~(a|b);w=~(c|d)
        return (a&w)|(z&c)|(b&d),(b&w)|(z&d)|(a&c)
    def reduce(a,b):
        while a|b:
            k=(a|b).bit_length()
            if k not in basis:break
            c,d=basis[k]
            if (a>>(k-1))&1: a,b=add(a,b,d,c)
            else: a,b=add(a,b,c,d)
        return a,b
    J=(1<<(56**2))-1
    for i,p in enumerate(G):
        a,b=reduce(sum(1<<(56*j+k) for j,k in enumerate(p)),0)
        if a|b:
            k=(a|b).bit_length()
            if not (a>>(k-1))&1:a,b=b,a
            basis[k]=a,b
            if reduce(J,0)==(0,0):return True,len(basis),i+1
    return False,len(basis),len(G)


def cycles(p):
    left=set(range(len(p)));ors=[]
    while left:
        x=left.pop();c=[x];y=p[x]
        while y!=x:c.append(y);left.remove(y);y=p[y]
        ors.append(c)
    return ors

def check_derangements():
    pointbits=[sum(1<<j for j,h in enumerate(orbit) if h>>i&1) for i in range(21)]
    paths=[[j for j,h in enumerate(orbit) if h>>i&1] for i in range(21)]
    hist=Counter()
    for p in G1:
        if any(i==j for i,j in enumerate(p)):continue
        k=sum(sum(1<<p[j] for j in b)==pointbits[i] for i,b in enumerate(paths))
        hist[k]+=1
    assert hist=={0:11520,1:7560}
    print('fixed projective points of hyperoval derangements:',hist)

def check_cosets():
    p=next(p for p in G1 if sorted(map(len,cycles(p)))==[7]*8)
    ors=cycles(p)
    def images(t):
        return tuple(sum(1<<t[i] for i in C) for C in ors)
    C=images(bytes(range(56)))
    candidates=set()
    for g in G2:
        w=images(g)
        if all(not a&b for a,b in zip(C,w)):
            candidates.add(w)
    S=sorted(candidates)
    pairs=[set(j for j,w in enumerate(S) if all(not a&b for a,b in zip(v,w))) for v in S]
    edges=sum(map(len,pairs))//2
    triangles=sum(len(pairs[i]&pairs[j]) for i,R in enumerate(pairs) for j in R)//6
    assert (len(S),edges,triangles)==(201,252,0)
    print('C7 normalized coset graph:',len(S),edges,triangles)

if __name__=='__main__':
    meet_families()
    check_derangements()
    check_cosets()
    if '--spans' in argv:
        result=[(mod2_shadow(G),mod3_shadow(G)) for G in (G1,G2)]
        assert [[a[0],b[0]] for a,b in result]==[[True,False],[True,True]]
        print('integer-shadow modulus checks:',result)
