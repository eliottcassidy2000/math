"""Independent literal tournament kernel; no producer/repository imports.
Mask gauge: bit k=1 means i->j in lexicographic unordered pairs i<j.
"""
from itertools import permutations
from math import comb, factorial
from functools import lru_cache

def tournament(n, mask):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if mask & (1<<k) else (j,i)
            a[u]|=1<<v;k+=1
    return tuple(a)

def validate(a):
    n=len(a)
    if any(row<0 or row >= 1<<n or row>>i&1 for i,row in enumerate(a)):
        raise ValueError('Invalid outgoing row')
    if any(((a[i]>>j&1)+(a[j]>>i&1))!=1 for i in range(n) for j in range(i+1,n)):
        raise ValueError('Pair is not uniquely oriented')
    return a

def induced(a, vs):
    return tuple(sum(((a[v]>>w)&1)<<j for j,w in enumerate(vs)) for v in vs)

def literal_paths(a):
    return [p for p in permutations(range(len(a))) if all((a[u]>>v)&1 for u,v in zip(p,p[1:]))]

def h_literal(a):return len(literal_paths(a))

def h_recursive(a):
    """Backwards deletion recursion independent of forward producer DP."""
    n=len(a)
    @lru_cache(None)
    def ending(s,v):
        s^=1<<v
        if not s:return 1
        return sum(ending(s,u) for u in range(n) if s>>u&1 and a[u]>>v&1)
    if not n:return 1
    return sum(ending((1<<n)-1,v) for v in range(n))

def path_cover_literal(a):
    """Permute vertices, force cuts at bad arcs, choose optional other cuts."""
    n=len(a)
    if not n:return (1,)
    counts=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum(not(a[u]>>v&1) for u,v in zip(p,p[1:]))
        for k in range(1,n+1):
            if bad<=k-1:counts[k]+=comb(n-1-bad,k-1-bad)
    for k in range(1,n+1):
        if counts[k]%factorial(k):raise ArithmeticError('not integral')
        counts[k]//=factorial(k)
    return tuple(counts)

def scc(a):
    n=len(a);reach=[a[v]|1<<v for v in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    left=set(range(n));out=[]
    while left:
        v=min(left);block=tuple(u for u in left if reach[v]>>u&1 and reach[u]>>v&1)
        left.difference_update(block);out.append(block)
    return tuple(sorted(out,key=lambda b:-sum(bool(a[b[0]]>>c[0]&1) for c in out if b!=c)))

def core(a):return induced(a,sorted(v for c in scc(a) if len(c)>1 for v in c))

def odd_cycles(a):
    n=len(a);cycles=[]
    for k in range(3,n+1,2):
        for p in permutations(range(n),k):
            if p[0]!=min(p):continue
            if all(a[u]>>v&1 for u,v in zip(p,p[1:]+p[:1])):
                cycles.append((p,sum(1<<v for v in p)))
    return tuple(cycles)

def ocf(a):
    cycles=odd_cycles(a)
    def collect(i,used):
        value=1
        for j in range(i,len(cycles)):
            if not used & cycles[j][1]:value+=2*collect(j+1,used|cycles[j][1])
        return value
    return collect(0,0)

def ordinal(*blocks):
    n=sum(map(len,blocks));rows=[];off=0
    for b in blocks:
        top=((1<<n)-1)^((1<<(off+len(b)))-1)
        rows.extend((r<<off)|top for r in b);off+=len(b)
    return tuple(rows)

def substitute(q,blocks):
    offsets=[];off=0
    for b in blocks:offsets.append(off);off+=len(b)
    out=[]
    for i,b in enumerate(blocks):
        external=sum(((1<<len(c))-1)<<offsets[j] for j,c in enumerate(blocks) if q[i]>>j&1)
        out.extend((r<<offsets[i])|external for r in b)
    return validate(tuple(out))

def extend(a,outgoing):
    n=len(a)
    return tuple(r|((0 if outgoing>>i&1 else 1)<<n) for i,r in enumerate(a))+(outgoing,)

def require(ok,message):
    if not ok:raise AssertionError(message)
