#!/usr/bin/env python3
"""Three independent constructions: literal tournament DP, positive walk,
and branched first-return equations. Explicit checks survive python -O."""

def require(v,msg):
    if not v: raise RuntimeError(msg)

def literal(m):
    # each iteration adds a,b with old -> a -> b -> old
    n=1+2*m
    a=[[False]*n for _ in range(n)]
    for t in range(m):
        i=2*t+1; j=i+1
        for v in range(i):
            a[v][i]=True; a[j][v]=True
        a[i][j]=True
    dp=[[0]*n for _ in range(1<<n)]
    dp[1][0]=1
    # all possible starts, independently, not just vertex 0
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not dp[s][v]: continue
            for w in range(n):
                if not s>>w&1 and a[v][w]: dp[s|1<<w][w]+=dp[s][v]
    return sum(dp[-1])

def walks(M):
    p={1:1}; out=[]
    for m in range(M+1):
        out.append(p.get(1,0)); q={}
        for j,c in p.items():
            for k,w in ((j-1,j*(j-1)),(j,j*(j+2)),(j+1,2*j+1),(j+2,1)):
                if w: q[k]=q.get(k,0)+c*w
        p=q
    return out

def mul(a,b,M):
    c=[0]*(M+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:M+1-i]): c[i+j]+=x*y
    return c

def invert_one_minus(f,M):
    g=[1]+[0]*M
    for n in range(1,M+1): g[n]=sum(f[k]*g[n-k] for k in range(1,n+1))
    return g

def fraction(M):
    e={j:[1]+[0]*M for j in range(2*M+3,2*M+5)}
    for j in range(2*M+2,0,-1):
        a=j*(j+2); b=j*(j+1)*(2*j+1); c=j*(j+1)**2*(j+2)
        prod=mul(e[j+1],e[j+2],M)
        f=[0]*(M+1)
        for n in range(1,M+1):
            f[n]=(a if n==1 else 0)+(b*e[j+1][n-2] if n>=2 else 0)+(c*prod[n-3] if n>=3 else 0)
        e[j]=invert_one_minus(f,M)
    return e[1]

if __name__=='__main__':
    x=walks(15); y=fraction(15)
    require(x==y, f'walk/fraction mismatch {x} != {y}')
    for m in range(7): require(literal(m)==x[m],f'literal mismatch at {m}')
    require(x[:8]==[1,3,15,123,1479,23955,496479,12725547],f'initial terms: {x[:8]}')
    print('PASS: first-return fraction equals forward weighted walks through z^15')
    print('PASS: independent literal tournament endpoint DP for m=0,...,6')
    print('Hamilton counts:',x[:8])

# Independent exhaustive skeleton walk enumeration and Lagrange coefficients.
def skeletons(M):
    # d maps (level,flat,one-rise,two-rise) to count of prefixes.
    from math import factorial, comb
    d={(1,0,0,0):1}
    for m in range(M+1):
        for i in range(m+1):
            for j in range(m+1):
                for k in range(m+1):
                    if i+2*j+3*k!=m: continue
                    r=i+j+k
                    formula=(factorial(r)//(factorial(i)*factorial(j)*factorial(k))*comb(m,r-1)//r) if r else 1
                    require(d.get((1,i,j,k),0)==formula, f'skeleton mismatch {(i,j,k)}')
        q={}
        def add(key,v): q[key]=q.get(key,0)+v
        for (level,i,j,k),v in d.items():
            add((level,i+1,j,k),v)
            add((level+1,i,j+1,k),v)
            add((level+2,i,j,k+1),v)
            if level>1: add((level-1,i,j,k),v)
        d=q
    print(f'PASS: Lagrange skeleton formula versus walks through z^{M}')

if __name__=='__main__': skeletons(12)

def literal_transitive(n,m):
    N=n+2*m
    a=[[False]*N for _ in range(N)]
    for i in range(n):
        for j in range(i+1,n): a[i][j]=True
    for t in range(m):
        i=n+2*t; j=i+1
        for v in range(i): a[v][i]=True; a[j][v]=True
        a[i][j]=True
    dp=[[0]*N for _ in range(1<<N)]
    for i in range(N): dp[1<<i][i]=1
    for s in range(1,1<<N):
        for i in range(N):
            if not dp[s][i]: continue
            for j in range(N):
                if not s>>j&1 and a[i][j]: dp[s|1<<j][j]+=dp[s][i]
    return sum(dp[-1])

if __name__=='__main__':
    for n in range(1,5):
        forms=[1+2**n,2*3**n+5*2**n-1,6*4**n+28*3**n+13*2**n-11]
        for m,want in enumerate(forms,1):
            require(literal_transitive(n,m)==want, f'transitive gadget mismatch {n=} {m=}')
    print('PASS: three transitive-carrier gadget formulas, n=1,...,4')

if __name__=='__main__':
    c={1:1}
    for m in range(6):
        for n in range(1,5):
            require(sum(v*q**n for q,v in c.items())==literal_transitive(n,m),f'exponential gadget mismatch {n=} {m=}')
        require(c.get(m+1)==__import__('math').factorial(m),f'leading coefficient mismatch {m=}')
        c={q:(q-1)*c.get(q-1,0)+q*q*c.get(q,0)-(q+1)*c.get(q+1,0) for q in range(1,m+3)}
    print('PASS: finite exponential formula and leading coefficient for m=0,...,5 and n=1,...,4')
