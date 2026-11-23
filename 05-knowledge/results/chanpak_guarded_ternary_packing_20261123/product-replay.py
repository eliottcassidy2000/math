#!/usr/bin/env python3
"""Independent small Held–Karp replay for shared-apex identities."""
from itertools import combinations, product
from functools import lru_cache

@lru_cache(None)
def subsets(rows):
    n = len(rows)
    end = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        end[1<<v][v] = 1
    for mask in range(1,1<<n):
        for v,count in enumerate(end[mask]):
            if count:
                nxt = rows[v] & ~mask
                while nxt:
                    x = nxt & -nxt
                    nxt -= x
                    end[mask|x][x.bit_length()-1] += count
    return end, [1]+[sum(x) for x in end[1:]]

def tournaments(n):
    edges = list(combinations(range(n),2))
    for bits in product((0,1),repeat=len(edges)):
        rows=[0]*n
        for (i,j),s in zip(edges,bits):
            if s: i,j=j,i
            rows[i] |= 1<<j
        yield tuple(rows)

def converse(rows):
    return tuple(sum(1<<i for i,r in enumerate(rows) if r>>j&1) for j in range(len(rows)))

def h(rows):
    return subsets(rows)[1][-1]

def join(w,v):
    a,b=len(w)-1,len(v)-1
    return tuple((r&((1<<a)-1)) | (((1<<b)-1)<<a) | ((r>>a)<<(a+b)) for r in w[:-1]) + tuple(
        ((r&((1<<b)-1))<<a) | ((r>>b)<<(a+b)) for r in v[:-1]) + (w[-1]|(v[-1]<<a),)

def sig(w):
    t=tuple(r&((1<<(len(w)-1))-1) for r in w[:-1])
    o=w[-1]
    end,cnt=subsets(t)
    st,_=subsets(converse(t))
    n=len(t); full=(1<<n)-1
    i=full^o
    ao=bi=0
    for s in range(1,1<<n):
        if s.bit_count()%2==0:
            # Remove a chosen start and run a DP for paths that start there.
            for start in range(n):
                if not (s>>start&1): continue
                # Tiny replay (order <= 3 in join census).
                from itertools import permutations
                for tail in permutations([p for p in range(n) if p!=start and s>>p&1]):
                    p=(start,)+tail
                    if all(t[a]>>b&1 for a,b in zip(p,p[1:])):
                        if o>>p[0]&1 and o>>p[-1]&1: ao+=cnt[full^s]
                        if i>>p[0]&1 and i>>p[-1]&1: bi+=cnt[full^s]
    ro=sum(end[full][v] for v in range(n) if o>>v&1)
    ri=sum(st[full][v] for v in range(n) if i>>v&1)
    return cnt[full],2*ao+ro,2*bi+ri,ro,ri

root_checks=0
for n in range(1,6):
    for rows in tournaments(n):
        incoming,_=subsets(rows)
        starts,_=subsets(converse(rows))
        cnt=subsets(rows)[1]; full=(1<<n)-1
        for v in range(n):
            value=sum((-1)**s.bit_count()*starts[s][v]*cnt[full^s] for s in range(1,full+1))
            assert value==-incoming[full][v]
            root_checks+=1

all4=[(w,sig(w),h(w)) for n in range(2,5) for w in tournaments(n)]
join_checks=0
for w,(d,x,_,ro,_),hw in all4:
    for v,(e,_,y,_,ri),hv in all4:
        assert h(join(w,v))==hw*hv+x*(hv-e)+y*(hw-d)+x*y+ro*ri
        join_checks+=1
c=(2,4,1)
for n in range(3,10):
    a=tuple(sum(1<<j for j in range(n) if j==i+1 or j<=i-2) for i in range(n))
    w=tuple(r|((1<<n) if i<n-1 else 0) for i,r in enumerate(a))+(1<<(n-1),)
    assert subsets(a)[0][-1][-1]==1
    assert h(join(w,c))==3*h(w)+2
    print(n,h(a),h(w),h(join(w,c)))
print(root_checks,join_checks,'identity checks passed')
