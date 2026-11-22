from itertools import combinations,permutations
from functools import cache
@cache
def chart(n):
    E=tuple(combinations(range(n),2)); return E
def rows(n,m):
    A=[0]*n
    for i,(u,v) in enumerate(chart(n)):
        if m>>i&1:A[u]|=1<<v
        else:A[v]|=1<<u
    return tuple(A)
def h(A):
    n=len(A); d=[[0]*n for _ in range(1<<n)]
    for j in range(n):d[1<<j][j]=1
    for s in range(1,1<<n):
        for j,c in enumerate(d[s]):
            if not c:continue
            b=A[j]&~s
            while b:
                e=b&-b;b-=e
                d[s|e][e.bit_length()-1]+=c
    return sum(d[-1])
def strongly(A):
    n=len(A)
    for j in range(n):
        s=1<<j
        for _ in range(n):
            t=s
            for a in range(n):
                if s>>a&1:t|=A[a]
            if t==s:break
            s=t
        if s!=(1<<n)-1:return False
    return True
def marked_canonical(A,u,v):
    return min(sum(((A[P[a]]>>P[b])&1)<<j for j,(a,b) in enumerate(chart(len(A))))
               for W in permutations(set(range(len(A)))-{u,v}) for P in ((u,v)+W,(v,u)+W))
def scan(n):
    A=[rows(n,m) for m in range(1<<len(chart(n)))]
    H=[h(a) for a in A]
    keep={(5,9):set(),(19,23):set()}
    if n <= 5: assert max(H) < 19
    for m,z in enumerate(H):
        if not any(z==a for a,b in keep):continue
        if z != 5 and not strongly(A[m]):continue
        for i,(u,v) in enumerate(chart(n)):
            other=m^(1<<i); k=(z,H[other])
            if k not in keep or not strongly(A[other]):continue
            keep[k].add(marked_canonical(A[m],u,v))
    return {k: sorted(x) for k,x in keep.items()}
expected = {3:([],[]),4:([],[]),5:([9,111,259,319],[]),
            6:([],[57,248,2085,2099,2255,2294])}
for n, result in expected.items():
    got=scan(n)
    answer=(got[5,9],got[19,23])
    assert answer==result, (n,answer)
    print(n, *answer)
    for S in answer:
        for m in S:
            A=rows(n,m)
            C=tuple(((1<<n)-1) ^ (1<<j) ^ row for j,row in enumerate(A))
            assert marked_canonical(C,0,1) in S
        print("  converse orbits:",sorted({tuple(sorted((m,marked_canonical(tuple(((1<<n)-1)^(1<<j)^row for j,row in enumerate(rows(n,m))),0,1)))) for m in S}))
