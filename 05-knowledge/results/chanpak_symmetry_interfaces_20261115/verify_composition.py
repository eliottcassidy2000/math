from itertools import permutations, product
from math import comb
def cover(rows):
    b=len(rows); f=[0]*(b+1)
    for w in permutations(range(b)):
        j=sum(not(rows[u]>>v&1) for u,v in zip(w,w[1:]))
        for k in range(j+1,b+1): f[k]+=comb(b-1-j,k-1-j)
    return f
def h(rows):
    n=len(rows); m=1<<n; d=[[0]*n for _ in range(m)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,m):
        for i in range(n):
            if d[s][i]:
                for j in range(n):
                    if s>>j&1==0 and rows[i]>>j&1:
                        d[s|1<<j][j]+=d[s][i]
    return sum(d[-1])
def cyclic(R):
    b=[len(r) for r in R]; x=[0,b[0],b[0]+b[1]]; n=sum(b); out=[0]*n
    for i,r in enumerate(R):
        j=(i+1)%3
        for u,v in enumerate(r):
            out[x[i]+u]=(v<<x[i]) | (((1<<b[j])-1)<<x[j])
    return out
def ch(F):
    z=0
    for p in product(*(range(1,len(f)) for f in F)):
        if max(p)-min(p)>1:continue
        a,b,c=p
        z+=(3 if a==b==c else 1)*F[0][a]*F[1][b]*F[2][c]
    return z
def ordinal_profile(f,t):
    N=len(f)+t-1
    E=[r**t*sum(comb(r,j)*f[j] for j in range(1,min(len(f)-1,r)+1)) for r in range(N+1)]
    return [sum((-1)**(k-r)*comb(k,r)*E[r] for r in range(k+1)) for k in range(N+1)]
if __name__=='__main__':
    R=(6,12,8,1); f=cover(R); assert f==[0,5,22,36,24]
    assert sum(all((R[i]>>j&1)==(R[p[i]]>>p[j]&1) for i in range(4) for j in range(4)) for p in permutations(range(4)))==1
    A=(24,1,3,6,14); B=(18,8,3,5,14)
    assert cover(A)==[0,11,78,198,240,120] and cover(B)==[0,13,78,198,240,120]
    K=(0,); U=(2,0); T=(6,4,0); C=(2,4,1)
    assert cover(T)==[0,1,6,6] and cover(C)==[0,3,6,6]
    Z=[K,U,T,C,R]
    for a,b,c in product(Z,repeat=3):
        if len(a)+len(b)+len(c)<=9:
            assert h(cyclic((a,b,c)))==ch([cover(a),cover(b),cover(c)])
    def optional(f,r):
        return sum(comb(r,k)*f[k] for k in range(1,min(r,len(f)-1)+1))
    for core in Z:
        f=cover(core)
        for B in (K,U):
            g=cover(cyclic((core,B,K)))
            for r in range(1,7):
                L,M,H=[optional(f,s) for s in (r-1,r,r+1)]
                target=r*(H+r*M-L) if B==K else r*((2*r+1)*H+(r*r-2)*M-(2*r-1)*L)
                assert optional(g,r)==target
    for core in Z:
        for B in Z:
            if len(core)+len(B)+1>7:continue
            f,b,g=cover(core),cover(B),cover(cyclic((core,B,K)))
            for r in range(1,7):
                x=[optional(f,j) for j in (r-1,r,r+1)]
                z=[optional(b,j) for j in (r-1,r,r+1)]
                assert optional(g,r)==r*((z[2]-z[1])*x[2]+(3*z[1]-z[0]-z[2])*x[1]-(z[1]-z[0])*x[0])
    X=cyclic((cyclic((K,U,K)),K,K))
    B=R+(0,); B=tuple((v|(1<<4)) if j<4 else v for j,v in enumerate(B))
    assert cover(X)[1]==37 and optional(cover(X),2)==352
    assert cover(B)[1]==5 and optional(cover(B),2)==64
    assert h(cyclic((X,B,K)))==18955
    P=ordinal_profile(f,3)
    assert P==[0,5,246,2406,9360,17280,15120,5040]
    assert ch([P]*3)==72465879938757 and ch([P]*3)%9==6
    print('Literal profile, rigidity and cyclic schedule controls passed.')
