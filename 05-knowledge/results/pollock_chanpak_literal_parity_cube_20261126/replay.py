from itertools import combinations, permutations, product
from functools import cache
from math import comb

P=tuple(combinations(range(7),2))
MARK=((0,1),(2,3),(4,5))
BASE=42232
MASKS=tuple(BASE|sum(1<<P.index(MARK[j]) for j in range(3) if e>>j&1) for e in range(8))
G=((0,1,2,3,4,5,6),(0,1,3,2,5,4,6),(1,0,4,5,3,2,6),(1,0,5,4,2,3,6))

def dig(n,mask):
    A=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=not(mask>>k&1)
        A[j][i]=not A[i][j]
    return A

def runvec(A,v):
    out=[z for z in range(len(A)) if z!=v]
    C=[0]*(len(A)+1)
    for path in permutations(out):
        L=(None,)+path; R=path+(None,)
        necessary=[i for i,(a,b) in enumerate(zip(L,R)) if a is not None and b is not None and not A[a][b]]
        allowed={i for i,(a,b) in enumerate(zip(L,R)) if (a is None or A[a][v]) and (b is None or A[v][b])}
        if not set(necessary)<=allowed: continue
        z=len(necessary); f=len(allowed)-z
        for j in range(f+1):
            C[z+j]+=comb(f,j)
    return tuple(C)

def hp(A):
    n=len(A)
    @cache
    def dp(mask,j):
        if mask==1<<j:return 1
        return sum(dp(mask^(1<<j),i) for i in range(n) if (mask>>i)&1 and A[i][j])
    return sum(dp((1<<n)-1,j) for j in range(n))

def strong(A):
    n=len(A)
    for B in (A,tuple(zip(*A))):
        U={0}
        for i in range(n): U|={v for u in U for v in range(n) if B[u][v]}
        if len(U)!=n:return False
    return True

def inflate(A,B):
    k=len(B)
    return [[A[i][j] if i<6 and j<6 else A[min(i,6)][min(j,6)] if min(i,j)<6 else B[i-6][j-6]
             for j in range(6+k)] for i in range(6+k)]

D=tuple(dig(7,m) for m in MASKS)
for m,A in enumerate(D):
    parity=m.bit_count()%2
    assert strong(A)
    C=runvec(A,6)
    assert C==(43,143-2*parity,137-parity,38,0,0,0,0)
    targets=set()
    for g in G:
        B=[[False]*7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                B[g[i]][g[j]]=A[i][j]
        targets.add(D.index(B))
    assert targets=={i for i in range(8) if i.bit_count()%2==parity}
print('Klein group: two regular four-point orbits; exact run vectors match')

for k in range(1,6):
    l=38*3**k+22*2**k-17
    assert tuple(hp(inflate(A,dig(k,0))) for A in D)==tuple(l+(0 if i.bit_count()%2 else 2**k) for i in range(8))
print('five transitive modules verified independently by subset path DP')

def F(x):return (x**3-x)//6
for k,mask,t,a,b,d,z,l,delta in [(4,4,11,2,2,4,9,5065,32),(5,40,1,1,6,9,38,20247,108)]:
    T=(t,t+a+b,t+a+d,t+b+d,z)
    U=(t+a,t+b,t+d,t+a+b+d,z)
    n0,n1=sum(map(F,T)),sum(map(F,U))
    assert n1-n0==a*b*d
    assert sum(T)==sum(U) and sum(x*x for x in T)==sum(x*x for x in U)
    assert (2*n0+1,2*n1+1)==(l,l+delta)
    J=[inflate(A,dig(k,mask)) for A in D]
    assert all(map(strong,J))
    assert tuple(map(hp,J))==tuple(l+(0 if i.bit_count()%2 else delta) for i in range(8))
    print('common Pollock endpoints',n0,n1,'strong tournament endpoints',l,l+delta,'indices',T,U)

N=(141*225-1)//2
assert N==sum(map(F,(5,11,33,33,28)))
assert N+225==sum(map(F,(8,8,30,36,28)))
assert hp(dig(3,2))==3 and hp(dig(4,4))==5
print('ordered-join control:',N,N+225)
