from shell import *
from collections import defaultdict,deque
from math import ceil

@lru_cache(None)
def canon(R):
    return canonical(R)
def stat(n,blue_only=False):
    F=tuple((i,j) for i in range(n) for j in range(i+2,n))
    idx={e:k for k,e in enumerate(F)}
    mirror=tuple(idx[n-1-j,n-1-i] for i,j in F)
    base=[(1<<i)-1 for i in range(n)]
    for i in range(n-1):base[i]^=1<<(i+1);base[i+1]^=1<<i
    def read(bits):
        rows=base[:]
        for k,(i,j) in enumerate(F):
            if bits>>k&1:rows[i]^=1<<j;rows[j]^=1<<i
        return min(canon(tuple(rows)),canon(tuple(((1<<n)-1)^(1<<i)^r for i,r in enumerate(rows))))
    G=defaultdict(set)
    root=read((1<<len(F))-1)
    for a in range(1<<(len(F)-1)):
        if blue_only and any((a>>i&1)!=(a>>j&1) for i,j in enumerate(mirror)):continue
        A,B=read(a),read(((1<<len(F))-1)^a)
        G[A];G[B]
        if A!=B:G[A].add(B);G[B].add(A)
    C={root};todo=[root]
    for a in todo:
        for b in G[a]:
            if b not in C:todo.append(b);C.add(b)
    e=sum(len(G[a]) for a in C)//2;v=len(C)
    return v,e,(max(0,ceil((e-3*v+6)/6)) if v>=3 else 0)
if __name__=='__main__':
    import sys
    for n in range(3,1+(int(sys.argv[1]) if len(sys.argv)>1 else 6)):
        print(n,'full',stat(n),'blue',stat(n,True),flush=True)
