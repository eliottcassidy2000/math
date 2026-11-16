"""Tiny weighted cycle-pack certificate; no SAT solver or isomorphism atlas."""

from functools import lru_cache
from itertools import combinations, permutations

@lru_cache(None)
def packing_masks(n):
    edge = {e:i for i,e in enumerate(combinations(range(n),2))}
    def mcycle(path):
        p=w=0
        for u,v in zip(path,path[1:]+path[:1]):
            bit=1<<edge[min(u,v),max(u,v)]
            p |= bit
            if u<v: w |= bit
        return p,w
    cycles=[]
    for k in range(3,n+1,2):
        for vertices in combinations(range(n),k):
            v,*others=vertices
            for tail in permutations(others):
                p,w=mcycle((v,)+tail)
                cycles.append((sum(1<<u for u in vertices),p,w,k))
    ans=[(p,w,2) for _,p,w,_ in cycles]
    for (v,p,w,j),(vv,pp,ww,k) in combinations(cycles,2):
        if not v&vv:
            ans.append((p|pp,w|ww,4))
    assert n<=8  # at most two odd directed cycles
    return list(edge),ans

def cycle_root(n):
    edge,_=packing_masks(n)
    z=sum(1<<edge.index((v,v+1)) for v in range(n-1))
    return z|(1<<edge.index((0,n-1))), z

def bounds(packs,assigned,positive):
    lo=hi=1
    for e,expected,value in packs:
        if (expected^positive)&(e&assigned):continue
        hi+=value
        if e&assigned==e:lo+=value
    return lo,hi

def solve(n,N):
    edge,packs=packing_masks(n)
    allbits=(1<<len(edge))-1
    root=cycle_root(n)
    leaves={}
    def recur(known,positive,proof):
        lo,hi=bounds(packs,known,positive)
        if not lo<=N<=hi:
            why="high" if lo>N else "low"
            leaves[proof]=(lo,hi,why)
            return
        free=allbits^known
        if not free:
            raise AssertionError((N,n,known,positive))
        # Pick the edge which immediately forces the most branches,
        # then shrink total residual interval as secondary criterion.
        options=[]
        while free:
            p=free&-free;free-=p
            ab=[bounds(packs,known|p,k) for k in (positive,positive|p)]
            options.append((sum(not a<=N<=b for a,b in ab),-sum(b-a for a,b in ab),p))
        _,_,p=max(options)
        i=p.bit_length()-1
        recur(known|p,positive,proof+((edge[i],0),))
        recur(known|p,positive|p,proof+((edge[i],1),))
    recur(*root,())
    assert sum(2**(len(edge)-root[0].bit_count()-len(p)) for p in leaves)==2**(len(edge)-root[0].bit_count())
    return leaves

if __name__=="__main__":
    for N,m in [(7,4),(21,6)]:
        for n in range(3,m+1):
            data=solve(n,N)
            print(N,n,len(data), 'leaves',sum(v[2]=="high" for v in data.values()),'high',max(map(len,data)),'depth')
