"""Check the triangle packing tradeoff and reconstruct the sparse witness code."""
from itertools import combinations
from math import comb

def arcs(n,s):
    A=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        if s>>b&1:A[j]|=1<<i
        else:A[i]|=1<<j
    return tuple(A)

def height(A):
    n=len(A);M=(1<<n)-1
    D=[[0]*n for _ in range(M+1)]
    for v in range(n):D[1<<v][v]=1
    for s in range(1,M+1):
        for i,d in enumerate(D[s]):
            if not d:continue
            t=A[i]&~s
            while t:
                bit=t&-t;t^=bit
                D[s|bit][bit.bit_length()-1]+=d
    return sum(D[M])

def pivot(A,S):
    if not S:return ()
    v=S[0]
    return pivot(A,tuple(x for x in S[1:] if A[x]>>v&1))+(v,)+pivot(A,tuple(x for x in S[1:] if A[v]>>x&1))

def check(A):
    n=len(A)
    cy=[]
    for x,y,z in combinations(range(n),3):
        if all(sum(A[i]>>j&1 for j in (x,y,z))==1 for i in (x,y,z)):
            cy.append((1<<x)|(1<<y)|(1<<z))
    F=k=0
    for c in cy:
        if not c&F:F|=c;k+=1
    h=height(A)
    assert h>=3**k+2*(len(cy)-k)*3**max(k-3,0)
    assert len(cy)*2**k<=15*h
    order=pivot(A,tuple(range(n)))
    M={i for i,x in enumerate(order) if F>>x&1}
    normal=[i for i in range(n) if i not in M]
    support={(i,j) for i,j in combinations(range(n),2)
             if (i in M or j in M) and A[order[j]]>>order[i]&1}
    W=sorted(normal, key=lambda i: -sum(A[order[i]]>>order[j]&1 for j in normal))
    ranks=[W.index(i) for i in normal]
    inv=[sum(q>r for q in ranks[:i]) for i,r in enumerate(ranks)]
    assert sum(inv)+len(support)<=len(cy)
    # Decode the ranks by inserting sequentially into the canonical transitive order.
    decoded=[]
    for j,i in enumerate(normal):decoded.insert(j-inv[j],i)
    assert decoded==W
    for i,j in combinations(range(n),2):
        back= ((i,j) in support) if (i in M or j in M) else (decoded.index(i)>decoded.index(j))
        assert back==bool(A[order[j]]>>order[i]&1)

if __name__=='__main__':
    total=0
    for n in range(1,7):
        for s in range(1<<comb(n,2)):
            check(arcs(n,s));total+=1
    print(total,'labeled tournaments: pack budget, inversion code and exact decode passed')
