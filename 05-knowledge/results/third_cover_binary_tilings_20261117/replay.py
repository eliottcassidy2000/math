"""Independent subset path/cover convolution and literal matrix checks."""
from functools import cache
from itertools import combinations,product,permutations
from math import prod,gcd
from pathlib import Path
import json

def need(z):
    if not z:raise ArithmeticError('certificate failed')

def cover(T,k=4):
    n=len(T);f=(1<<n)-1
    need(all(0<=r<=f and not r>>i&1 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i)))
    dp=[[0]*n for _ in range(f+1)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,f+1):
        for j,v in enumerate(dp[S]):
            if not v:continue
            a=T[j]&~S
            while a:
                x=a&-a; a^=x
                dp[S|x][x.bit_length()-1]+=v
    H=[1]+[sum(a) for a in dp[1:]]
    @cache
    def cnt(S,j):
        if j==1:return H[S]
        z=0;A=S
        while True:
            z+=H[A]*cnt(S^A,j-1)
            if not A:return z
            A=(A-1)&S
    return tuple(cnt(f,j) for j in range(1,k+1))

def framed(n,t,u):
    M=[[int(i==j) for j in range(n)] for i in range(n)]
    for i in range(n-1):M[i][i+1]=1;M[i+1][i]=-1
    E=[(i,j) for i,j in combinations(range(n),2) if j-i>1]
    for (i,j),a,b in zip(E,t,u):M[i][j]=a;M[j][i]=-b
    return M

def had(M):
    return all(sum(a*b for a,b in zip(M[i],M[j]))==0 for i,j in combinations(range(len(M)),2))

def sylvester_path(m):
    n=1<<m; E=[i^(i>>1) for i in range(n)];M=[]
    for x in E:
        M.append([(-1)**((x&y).bit_count()+y.bit_count()) for y in E])
    d=[1]
    for i in range(n-1):d.append(d[-1]*M[i][i+1])
    return [[d[i]*d[j]*v for j,v in enumerate(row)] for i,row in enumerate(M)]

def check_fano():
    S=(1,2,4); T=tuple(sum(1<<((i+d)%7) for d in S) for i in range(7))
    F=[{frozenset((a+t)%7 for a in B) for t in range(7)} for B in ((0,1,3),(0,2,3))]
    C={frozenset((i,j,k)) for i,j,k in combinations(range(7),3)
       if (T[i]>>j&1)==(T[j]>>k&1)==(T[k]>>i&1)}
    need(F[0].isdisjoint(F[1]) and F[0]|F[1]==C and len(C)==14)
    for D in F:
        need(len(D)==7)
        need(all(sum(i in L and j in L for L in D)==1 for i,j in combinations(range(7),2)))
        need(all(a&b for a,b in combinations(D,2)))
    for X,Y in ((F[0],F[1]),(F[1],F[0])):
        for L in X:need(sum(L.isdisjoint(K) for K in Y)==1)
    need(cover(T,3)==(189,2016,12159))

@cache
def canon(rows):
    n=len(rows)
    cells=[tuple(i for i,r in enumerate(rows) if r.bit_count()==k) for k in range(n)]
    return min(tuple(sum((rows[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        for choice in product(*(permutations(c) for c in cells if c))
        for p in [sum(choice,())])

def converse(a):
    f=(1<<len(a))-1
    return tuple(f^(1<<i)^r for i,r in enumerate(a))

def merged(a):return min(canon(a),canon(converse(a)))

def edge_types(T,U):
    n=len(T);f=(1<<n)-1; goal=merged(U);types=set()
    for A in (T,converse(T)):
        for p in permutations(range(n)):
            if not all(A[i]>>j&1 for i,j in zip(p,p[1:])):continue
            X=tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(p))for i in p)
            Y=tuple(f^(1<<i)^x^((1<<(i-1) if i else 0)|(1<<(i+1) if i<n-1 else 0)) for i,x in enumerate(X))
            if merged(Y)==goal:
                mirrored=all((X[i]>>j&1)+(X[n-1-i]>>(n-1-j)&1)==1 for i,j in combinations(range(n),2))
                types.add('blue' if mirrored else 'black')
    return types

def metagraph_cases():
    examples = [('black only', (78, 20, 56, 114, 225, 67, 134, 47), (234, 100, 9, 50, 167, 68, 156, 46), 527, 441), ('no edge', (254, 84, 72, 18, 36, 78, 152, 62), (2, 84, 105, 243, 165, 195, 145, 7), 103, 253), ('same node', (26, 124, 217, 176, 160, 197, 153, 3), (202, 12, 105, 80, 167, 203, 146, 14), 405, 405), ('blue only', (142, 108, 88, 48, 227, 197, 137, 14), (194, 164, 105, 19, 39, 73, 154, 60), 565, 637)]
    for label,T,U,h,k in examples:
        n=len(T); E=[(i,j)for i,j in combinations(range(n),2) if j-i>1]
        need(all((T[i]>> (i+1)&1) and (U[i]>>(i+1)&1) for i in range(n-1)))
        need(had(framed(n,[2*(T[i]>>j&1)-1 for i,j in E],[2*(U[i]>>j&1)-1 for i,j in E])))
        need(cover(T,1)==(h,) and cover(U,1)==(k,))
        if label=='same node':need(merged(T)==merged(U))
        else:
            need(merged(T)!=merged(U))
            expectation={'blue only':{'blue'},'black only':{'black'},'no edge':set()}[label]
            need(edge_types(T,U)==expectation)
        print(label,h,k)

def main():
    Z=json.loads(Path(__file__).with_name('witness.json').read_text())
    P=[]
    for L in Z:
        n=0; R=[]
        for T,w in L:
            p=cover(T);need(tuple(w)==p[:3]);n+=len(T);R.append(p)
        need(n==137);P.append(tuple(prod(x) for x in zip(*R)))
    a,b=P
    need(tuple(7*x for x in a[:3])==b[:3] and 7*a[3]!=b[3])
    need(a[0]%3==a[1]%3==0 and gcd(a[0],a[1])>1)
    print("Exact ×7 at q=1,2,3; q=4 differs, and lower H and P2 are multiples of 3")
    for n in (2,4):
        M=list(product((-1,1),repeat=(n-1)*(n-2)//2))
        S=[(t,u) for t,u in product(M,repeat=2) if had(framed(n,t,u))]
        need(len(S)==(1 if n==2 else 3))
        if n==4:
            need(sum(t==u for t,u in S)==2)
            need(not any(tuple(-a for a in t)==u for t,u in S))
            need(not any(t[::-1]==u for t,u in S))
    for m in range(1,6):
        M=sylvester_path(m); n=len(M)
        need(had(M) and all(M[i][i]==1 for i in range(n)))
        need(all(M[i][i+1]==1 and M[i+1][i]==-1 for i in range(n-1)))
    check_fano()
    metagraph_cases()
    print('Two-sided framing, Sylvester controls and both Paley Fano designs: PASS')

if __name__=='__main__':main()
