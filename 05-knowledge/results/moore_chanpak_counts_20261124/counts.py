"""Exact small Moore counts, a capped 50-vertex count, and port signatures."""
from collections import Counter
from itertools import combinations, permutations
from math import comb

def tournament(n, forward):
    rows=[0]*n
    for i,j in combinations(range(n),2):
        a,b=(i,j) if forward(i,j) else (j,i)
        rows[a]|=1<<b
    return tuple(rows)

def masks(n): return [(i,j) for i,j in combinations(range(n),2) if j>i+1]

def from_mask(n, t):
    back={e for k,e in enumerate(masks(n)) if (t>>k)&1}
    return tournament(n,lambda i,j:(i,j) not in back)

def block(T, x):
    return tuple(T[i]>>j&1 for i,j in zip(x,x[1:]))

def flipped(T, i, j):
    R=list(T); R[i]^=1<<j; R[j]^=1<<i; return tuple(R)

def dp(T):
    n=len(T);full=(1<<n)-1
    F=[[0]*n for _ in range(full+1)]
    for i in range(n):F[1<<i][i]=1
    for S in range(1,full+1):
        for i,p in enumerate(F[S]):
            if not p:continue
            A=T[i]&(full^S)
            while A:
                j=A&-A;A^=j
                F[S|j][j.bit_length()-1]+=p
    return F

def path_cycle(T):
    n=len(T)
    F=dp(T); H=sum(F[-1])
    A=[sum(F[U]) if U else 1 for U in range(1<<n)]
    Z2=sum(A[S]*A[((1<<n)-1)^S] for S in range(1<<n))
    # Appending the rest to distinguished vertex 0 computes cycles up to rotation.
    B=[[0]*n for _ in range(1<<n)]; B[1][0]=1
    for S in range(1,1<<n,2):
        for i,p in enumerate(B[S]):
            if not p:continue
            R=T[i]&~S
            while R:
                j=R&-R;R^=j;B[S|j][j.bit_length()-1]+=p
    C=sum(p for j,p in enumerate(B[-1]) if T[j]&1)
    return H,C,Z2

def cycle5():
    return tournament(5,lambda i,j:(j-i)%5 in (1,2))

def compact(q):
    m=5*q
    def f(a,b):
        i,x=divmod(a%m,5);j,y=divmod(b%m,5)
        if a//m==b//m:
            return (1 <= (y-x)%5<=2) if i==j else (a<m)
        return (y-x-i*j)%5 !=0
    return tournament(2*m,f)

def observer(n, adjacent, order):
    def f(i,j):
        if j<n:return True
        if i>=n:return False
        a,b=order[i],order[j-n]
        return a==b or adjacent(a,b)
    return tournament(2*n,f)

def hs_edge(a,b):
    if a==b:return False
    if a>=25 and b<25:a,b=b,a
    i,x=divmod(a%25,5); j,y=divmod(b%25,5)
    if a//25==b//25:return i==j and (y-x)%5 in ((1,4) if a<25 else (2,3))
    return (y-x-i*j)%5==0

def hs_spine():
    P=[];a=0
    for i in range(5):
        P.extend(5*i+(a+j)%5 for j in range(5))
        P.extend(25+5*i+(a+4+i*i+2*j)%5 for j in range(5))
        a=(a-3-i)%5
    assert a==0 and all(hs_edge(a,b) for a,b in zip(P,P[1:]+P[:1]))
    return P

def fixed_block(T, blocks):
    """Count paths in which every displayed block is visited contiguously."""
    n=len(T); m=len(blocks)
    local={}
    for B in blocks:
        for i,j in permutations(B,2):
            rest=tuple(x for x in B if x!=i and x!=j)
            local[i,j]=sum(all(block(T,(i,*P,j))) for P in permutations(rest))
    F=[[0]*n for _ in range(1<<m)]
    for b,B in enumerate(blocks):
        for j in B:F[1<<b][j]=sum(local[i,j] for i in B if i!=j)
    for S in range(1,1<<m):
        for b,B in enumerate(blocks):
            if (S>>b)&1:continue
            G=F[S|(1<<b)]
            for i in B:
                pre=sum(p for j,p in enumerate(F[S]) if (T[j]>>i)&1)
                if pre:
                    for j in B:
                        if i!=j:G[j]+=pre*local[i,j]
    return sum(F[-1])

def chord_subsets(mask, n):
    E=[e for k,e in enumerate(masks(n)) if (mask>>k)&1]
    base=from_mask(n,0);H=[]
    for s in range(1<<len(E)):
        T=base
        for k,e in enumerate(E):
            if (s>>k)&1:T=flipped(T,*e)
        H.append(sum(dp(T)[-1]))
    hits=[(i,i^(1<<k),h,H[i^(1<<k)]) for i,h in enumerate(H)
          for k in range(len(E)) if h+2==H[i^(1<<k)]]
    return H,hits

if __name__=='__main__':
    eye={5:(1<<len(masks(5)))-1,10:(1<<len(masks(10)))-1}
    values={
        'C5 circulant':cycle5(),
        'C5 graph-spine sparse':from_mask(5,1<<masks(5).index((0,4))),
        'C5 graph-spine dense':from_mask(5,eye[5]^(1<<masks(5).index((0,4)))),
        'C5 observer 01243':observer(5,lambda i,j:(i-j)%5 in (1,4),(0,1,2,4,3)),
        'C5 observer 01432':observer(5,lambda i,j:(i-j)%5 in (1,4),(0,1,4,3,2)),
        'Petersen compact':compact(1),
        'Petersen graph-spine sparse':from_mask(10,4304406596),
        'Petersen graph-spine dense':from_mask(10,4304406596^eye[10]),
    }
    for name,T in values.items():print(name,path_cycle(T))
    for m in (4304406596,4311892004):
        H,hits=chord_subsets(m,10)
        print('Petersen chord subcube:',m,'set',sorted(set(H)),'exact +2 witnesses',hits)
    T=compact(5)
    print('HS compact, contiguous-block paths',fixed_block(T,[list(range(i,i+5)) for i in range(0,50,5)]))
    Q=hs_spine()
    for dense in (False,True):
        T=tournament(50,lambda i,j:j==i+1 or (hs_edge(Q[i],Q[j])==dense))
        print('HS graph-spine',dense,'contiguous-block paths',fixed_block(T,[list(range(i,i+5)) for i in range(0,50,5)]))
