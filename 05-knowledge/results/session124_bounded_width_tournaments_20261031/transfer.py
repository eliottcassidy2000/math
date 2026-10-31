"""Exact path-cover transfer for forward tournaments of backedge bandwidth k."""
from collections import Counter
from functools import lru_cache

@lru_cache(None)
def step(state, r, k, mask):
    # Frontier labels 0,...,r-1, newest new vertex r; bit=1 means new->old.
    out=Counter()
    incoming=[None]+[i for i,(_,e) in enumerate(state) if e==-1 or not(mask>>e&1)]
    outgoing=[None]+[i for i,(s,_) in enumerate(state) if s>=0 and mask>>s&1]
    for i in incoming:
        for j in outgoing:
            if i is not None and i==j: continue
            parts=[p for a,p in enumerate(state) if a!=i and a!=j]
            parts.append((r if i is None else state[i][0], r if j is None else state[j][1]))
            if r==k:
                parts=[(s-1 if s>=0 else -1,e-1 if e>=0 else -1) for s,e in parts]
            if sum(s==-1 for s,e in parts)>1: continue
            out[tuple(sorted(parts))]+=1
    return tuple(out.items())

def advance(vec,r,k,mask):
    out=Counter()
    for state,n in vec.items():
        for dest,m in step(state,r,k,mask): out[dest]+=n*m
    return out

def count(word,k):
    vec={():1}
    for n,mask in enumerate(word): vec=advance(vec,min(n,k),k,mask)
    return sum(v for s,v in vec.items() if len(s)==1)

def adjacency(word,k):
    n=len(word); a=[[False]*n for _ in range(n)]
    for j,mask in enumerate(word):
        r=min(j,k)
        for i in range(j):
            backward=i>=j-r and (mask>>(i-j+r)&1)
            a[j][i]=bool(backward);a[i][j]=not backward
    return a

def dp(word,k):
    a=adjacency(word,k);n=len(a)
    if not n:return 0
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if not d[s][i]:continue
            for j in range(n):
                if not s>>j&1 and a[i][j]:d[s|1<<j][j]+=d[s][i]
    return sum(d[-1])

def carrier(k):
    states={()};r=0
    for _ in range(k):
        states={t for s in states for a in range(1<<r) for t,m in step(s,r,k,a)};r+=1
    frontier=set(states)
    while frontier:
        new={t for s in frontier for a in range(1<<k) for t,m in step(s,k,k,a)}-states
        states|=new;frontier=new
    return tuple(sorted(states))

if __name__=='__main__':
    import itertools,json,random
    checks=0
    for k,maxn in [(0,8),(1,9),(2,7),(3,5)]:
        for n in range(1,maxn+1):
            for word in itertools.product(*(range(1<<min(i,k)) for i in range(n))):
                x=count(word,k);y=dp(word,k)
                if x!=y:raise ValueError((k,word,x,y))
                checks+=1
        print('k',k,'states',len(carrier(k)),flush=True)
    rng=random.Random(124)
    for k in [2,3,4]:
        for _ in range(50):
            n=rng.randrange(6,11);w=tuple(rng.randrange(1<<min(i,k)) for i in range(n))
            if count(w,k)!=dp(w,k):raise ValueError(w)
            checks+=1
    print(json.dumps({'checks':checks,'k1_backward':[count((0,)+(1,)*(n-1),1) for n in range(1,14)]}))
