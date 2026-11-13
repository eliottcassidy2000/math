from functools import cache
def join(A,B):
    n=len(A);m=len(B)
    return [a|(((1<<m)-1)<<n) for a in A]+[b<<n for b in B]
def tt(n):return [((1<<n)-1)^((1<<(j+1))-1) for j in range(n)]
def T(r,i):
    A=[]
    for _ in range(i): A=join(A,[2,4,1])
    return join(A,tt(2*(r-i)))
def Y(r,i,L):
    A=T(r,i); B=join(T(r,r-i),tt(L)); n=len(A);m=len(B)
    out=join(A,B)
    for j in range(n,n+m):out[j]|=1<<(n+m)
    return out+[(1<<n)-1]
def counts(A):
    n=len(A);N=(1<<n)-1
    D=[[0]*n for _ in range(N+1)]
    for j in range(n):D[1<<j][j]=1
    h=[1]
    for S in range(1,N+1):
        h.append(sum(D[S]))
        for j,e in enumerate(D[S]):
            if e:
                R=A[j]&~S
                while R:
                    b=R&-R;R^=b;D[S|b][b.bit_length()-1]+=e
    return h[-1],sum(h[S]*h[N^S] for S in range(N+1))
def qformula(r,i,L):
    a,b=3**i,3**(r-i); v=4**r
    p=(a,a*v,a*11**i*9**(r-i))
    q=(b,b*v*2**L,b*11**(r-i)*9**i*3**L)
    return 2*(p[1]*q[1]+sum((p[k]-p[1])*(q[k]-q[1]) for k in (0,2)))
for r,L in [(1,1),(1,3),(2,2),(2,3),(3,3)]:
    h=3**r*((4**r-1)*(4**r*2**L-1)+2)
    vals=[]
    for i in range(r+1):
        H,Q=counts(Y(r,i,L))
        if H!=h or Q!=qformula(r,i,L):raise ValueError((r,i,L,H,Q))
        vals.append(Q)
    if vals!=sorted(set(vals)):raise ValueError(vals)
    print(r,L,h,vals)
