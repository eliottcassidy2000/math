from itertools import permutations, product

def path(A,P):
    return all(A[u]>>v&1 for u,v in zip(P,P[1:]))

def hc(A):
    return sum(path(A,P) for P in permutations(range(len(A))))

def cover(A):
    n=len(A)
    W=[[0]*n for _ in A]
    for P in permutations(range(n)):
        for k in range(1,n):
            if path(A,P[:k]) and path(A,P[k:]):
                x,y=P[k-1:k+1]
                W[x][y]+=1
                W[y][x]+=1
    return W

def extend(A,S):
    n=len(A)
    return tuple(x|(0 if S>>i&1 else 1<<n) for i,x in enumerate(A))+(S,)

def z(A,k):
    n=len(A); c=0
    for f in product(range(k),repeat=n):
        m=1
        for i in range(k):
            m*=sum(path(A,P) for P in permutations(j for j in range(n) if f[j]==i))
        c+=m
    return c

def wrap(A,k):
    a=tuple(A)
    for _ in range(k):a=extend(a,0)
    n=len(a)
    return tuple(x|(1<<n) for x in a)+(1<<(n+1),(1<<n)-1)

for A in ((2,4,1),(2,8,3,5),(2,4,17,7,11)):
    h=hc(A);W=cover(A)
    assert all(W[i][j]>=2 and W[i][j]%2==0 for i in range(len(A)) for j in range(i+1,len(A)))
    assert sum(map(sum,W))==2*(z(A,2)-2*h)
    for S in range(1<<len(A)):
        D=sum(W[i][j] for i in range(len(A)) for j in range(len(A)) if A[i]>>j&1 and S>>i&1 and not S>>j&1)
        assert hc(extend(A,S))==h+D
    for k in range(2):
        for i in range(1,4):
            assert z(wrap(A,k),i)==i*((i+1)**k*z(A,i+1)+i**(k+1)*z(A,i)-(i-1)**k*z(A,i-1))

r={1:1};expect=[[1],[1,1],[-1,5,2],[-11,13,28,6],[-37,-43,254,180,24]]
for e in expect:
    assert list(v for k,v in sorted(r.items()))==e
    nxt={}
    for j,v in r.items():
        for i,c in ((j+1,j),(j,j*j),(j-1,-j)):
            if i:nxt[i]=nxt.get(i,0)+v*c
    r=nxt
for signs,d5 in (([2,-2,2,-2,2],199),([2,2,-2,2,-2],-27)):
    ans=[]
    for row,target in zip(expect,signs):
        t=target-sum(x*y for x,y in zip(row,ans))
        assert t%row[-1]==0
        ans.append(t//row[-1])
    assert ans[-1]==d5 and d5%5
print('literal mixed cuts and cyclic transfers pass')
