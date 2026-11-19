from functools import lru_cache
from itertools import combinations

def profiles(A,x):
    n=len(A)-1
    p=[i for i in range(n+1) if i!=x]
    B=[sum((1<<j) for j,v in enumerate(p) if A[u]>>v&1) for u in p]
    F=(1<<n)-1
    O=sum(1<<j for j,v in enumerate(p) if A[x]>>v&1)
    I=F^O
    def path(st):
        dp=[[0]*n for _ in range(F+1)]
        for i in range(n):
            if st>>i&1:dp[1<<i][i]=1
        for s,row in enumerate(dp):
            for i,v in enumerate(row):
                if not v:continue
                e=B[i]&~s
                while e:
                    b=e&-e;e-=b
                    dp[s|b][b.bit_length()-1]+=v
        return dp
    calls=path(O)
    C=[0]*(F+1)
    for i in range(n):
        D=path(1<<i)
        for s in range(1<<i,F+1,1<<(i+1)):
            if s.bit_count()>=3 and s.bit_count()%2:
                C[s]=sum(D[s][j] for j in range(n) if B[j]>>i&1)
    d=2+n//3
    @lru_cache(None)
    def Z(s):
        if not s:return [1]+[0]*d
        b=s&-s
        a=Z(s^b)[:]
        t=s
        while t:
            if t&b and C[t]:
                v=Z(s^t)
                for k in range(d):
                    a[k+1]+=C[t]*v[k]
            t=(t-1)&s
        return a
    U=Z(F)[:];V=[0]*(d+1);W=[0]*(d+1)
    for s in range(1,F+1):
        c=sum(calls[s][i] for i in range(n) if I>>i&1)
        if not c:continue
        X=V if s.bit_count()%2 else W
        for j,v in enumerate(Z(F^s)[:-1]):X[j+1]+=c*v
    return U,V,W

def at2(P):return sum(c<<i for i,c in enumerate(P))

def require(x,detail):
    if not x:
        raise AssertionError(detail)

def blow(A,x,m):
    pts=[i for i in range(len(A)) for k in range(m if i==x else 1)]
    return [sum(1<<j for j,v in enumerate(pts)
                if (i<j if u==v else A[u]>>v&1)) for i,u in enumerate(pts)]

def literal_polynomial(A):
    return profiles(A+[0],len(A))[0]

def ham(A):
    n=len(A);F=(1<<n)-1
    d=[[0]*n for _ in range(1<<n)]
    require(all(not(A[i]>>i&1) and all((A[i]>>j&1)+(A[j]>>i&1)==1
                for j in range(i)) for i in range(n)),A)
    for i in range(n):d[1<<i][i]=1
    for s,row in enumerate(d):
        for i,x in enumerate(row):
            if not x:continue
            J=A[i]&~s
            while J:
                j=J&-J;J-=j
                d[s|j][j.bit_length()-1]+=x
    return sum(d[-1])

BANK=[
    ([4,5,48,7,43,75,159,63],0,(1,0,0),(0,36,4)),
    ([0,21,73,67,45,15,179,63],7,(1,2,-1),(0,34,5)),
]
for A,x,expectedC,expectedD in BANK:
    U,V,W=profiles(A,x)
    C=tuple(a-b for a,b in zip(U,V))
    D=tuple(a+b for a,b in zip(V,W))
    require(C==expectedC+(0,)*(len(C)-3) and D==expectedD+(0,)*(len(D)-3),(C,D))
    for m in range(1,5):
        T=blow(A,x,m)
        Z=literal_polynomial(T)
        rhs=tuple(a+2**(m-1)*b for a,b in zip(C,D))
        require(tuple(Z)==rhs+(0,)*(len(Z)-len(rhs)),(m,Z,rhs))
        require(at2(Z)==ham(T)==1+44*(1<<m),(m,Z))
# Minimum exterior at matching number two; verify the pole coefficients
# from independently computed cycle packings at m=1,2,3,4.
from fractions import Fraction as Q
A=[6,28,24,17,1]
seq=[literal_polynomial(blow(A,0,m)) for m in range(1,5)]
bases=[-1,1,2,3]
for z in (-1,0,1,2,3,4):
    mat=[[Q(b**m) for b in bases]+[Q(sum(t*z**i for i,t in enumerate(seq[m-1])))]
         for m in range(1,5)]
    for i in range(4):
        j=next(j for j in range(i,4) if mat[j][i])
        mat[i],mat[j]=mat[j],mat[i]
        mat[i]=[x/mat[i][i] for x in mat[i]]
        for j in range(4):
            if j!=i:
                mat[j]=[x-mat[j][i]*y for x,y in zip(mat[j],mat[i])]
    require(mat[0][-1]==Q(z*(z-2),2) and mat[3][-1]==Q(z*(z+2),2),(z,mat))
print('PASS: two identical Hamilton rays with distinct cycle profiles, and tight-cut poles')
