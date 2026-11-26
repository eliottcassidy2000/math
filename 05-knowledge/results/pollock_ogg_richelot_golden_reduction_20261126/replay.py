from itertools import combinations, permutations, product

OGG=(2,3,5,7,11,13,17,19,23,29,31,41,47,59,71)
V=tuple(range(5))
pairs=list(combinations(V,2))
E=[frozenset((u,v)) for u,v in combinations(pairs,2) if not set(u)&set(v)]
cycles=set()
for tail in permutations(V[1:]):
    seq=(0,)+tail
    cycles.add(frozenset(tuple(sorted((seq[i],seq[(i+1)%5]))) for i in V))
assert len(cycles)==12
M={frozenset(i for i,e in enumerate(E) if len(e&c)==1) for c in cycles}
assert len(M)==6
assert all(len(m)==5 and len(set().union(*(E[i] for i in m)))==10 for m in M)
assert all(len(m&n)==1 for m,n in combinations(M,2))
def omit(e):
    return next(v for v in V if all(v not in pair for pair in e))
labels=[frozenset(j for j,m in enumerate(sorted(M,key=sorted)) if i in m) for i in range(15)]
assert len(set(labels))==15 and all(len(l)==2 for l in labels)
for i,j in combinations(range(15),2):
    assert (not labels[i]&labels[j])==(bool(E[i]&E[j]) or omit(E[i])==omit(E[j]))
print('15 chart labels on six Petersen matchings; 105 intersection checks')

vals=(-3,-1,0,1,3)
def delta(e):
    a,b=sorted(e)
    d=vals[omit(e)]
    return (vals[a[0]]-d)*(vals[a[1]]-d)-(vals[b[0]]-d)*(vals[b[1]]-d)
print('absolute determinant multiset:',sorted(abs(delta(e)) for e in E))
assert sum(delta(e)==0 for e in E)==4
nondeg5=frozenset(i for i,e in enumerate(E) if delta(e)%5)
gold=frozenset(tuple(sorted((i,j))) for i,j in combinations(V,2) if (vals[i]-vals[j])%5 in (1,4))
assert nondeg5==frozenset(i for i,e in enumerate(E) if len(e&gold)==1)
print('Richelot decompositions at good primes 4, except 10 at 5; remaining five = marked golden matching')

def points(p):
    n=0
    for i in range(4):
        for y in product(range(p),repeat=3-i):
            x=(0,)*i+(1,)+y
            x+=((-sum(x))%p,)
            if p==3:
                v=sum(x[a]*x[b]*x[c] for a,b,c in combinations(V,3))
            else:
                v=sum(a**3 for a in x)
            n+=v%p==0
    return n
for p in sorted(set(OGG)|{37,61,73,101}):
    chi=0 if p==5 else 1 if p%5 in (1,4) else -1
    n=points(p)
    assert n==p*p+(6+chi)*p+1
    print(p, 'Ogg' if p in OGG else 'non-Ogg', 'geometric lines over Fp:',{0:21,1:27,-1:15}[chi], 'surface points:',n)
