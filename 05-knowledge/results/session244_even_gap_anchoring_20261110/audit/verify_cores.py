#!/usr/bin/env python3
"""Independent forward-generation census, not producer deletion recursion."""
import itertools,json

def need(x,msg):
    if not x:raise RuntimeError(msg)

def reaches(a,r):
    seen={r}
    while True:
        t=seen|{j for i in seen for j in range(len(a)) if a[i]>>j&1}
        if t==seen:return len(t)==len(a)
        seen=t

def insert(a,root,pos):
    n=len(a);b=list(a)+[0 if root is None else 1<<root]
    for i in range(n):
        if i!=root:b[i]|=1<<n
    # Relabel newborn to pos and old vertices in their original order.
    perm=list(range(n));perm.insert(pos,n)
    return tuple(sum(1<<j for j,w in enumerate(perm) if b[v]>>w&1) for v in perm)

def advance(bank,guard):
    result=set()
    for a in bank:
        roots=[None]+[r for r in range(len(a)) if not guard or reaches(a,r)]
        for r in roots:
            for pos in range(len(a)+1):result.add(insert(a,r,pos))
    return result

unrestricted=guarded={(0,)};counts=[]
for n in range(1,7):
    if n>1:unrestricted=advance(unrestricted,False);guarded=advance(guarded,True)
    counts.append([n,len(unrestricted),len(guarded)])
need([x[1] for x in counts]==[1,2,8,64,1000,26720],'unrestricted forward census')
need([x[2] for x in counts]==[1,2,8,56,760,18960],'guarded forward census')
cores=[(6,12,24,17,3),(6,12,48,5,11,27),(6,12,24,33,11,23),(6,12,24,33,43,7),(6,24,26,33,41,7),(6,24,42,33,13,19)]
known={5:set(),6:set()};sizes=[]
for a in cores:
    orbit={tuple(sum(1<<j for j,w in enumerate(p) if a[v]>>w&1) for v in p) for p in itertools.permutations(range(len(a)))}
    need(not known[len(a)]&orbit,'orbit overlap');known[len(a)]|=orbit;sizes.append(len(orbit))
for n in (5,6):
    found=set();edges=list(itertools.combinations(range(n),2))
    for code in range(1<<len(edges)):
        a=[0]*n
        for k,(i,j) in enumerate(edges):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
        if min(x.bit_count() for x in a)<2:continue
        # Test every eligible proper subset directly, with no recursive cache.
        if any(all((a[v]&mask).bit_count()>=2 for v in range(n) if mask>>v&1) for mask in range(1,(1<<n)-1) if mask.bit_count()>=5):continue
        found.add(tuple(a))
    need(found==known[n],('minimal cores',n))
print(json.dumps({'status':'PASS','method':'forward graph insertion plus direct induced-subset obstruction checks','counts_order_unrestricted_guarded':counts,'core_orbit_sizes':sizes},indent=2))
