from math import isqrt
from pathlib import Path

def recurrent(u):
    f=1
    while u:
        q,p=divmod(u,2);f+=(q+(q&1))*(p==0);u=q
    return f

def closed(u):
    bits=bin(u)[2:] if u else ''
    return 1+sum(2**(j-i-1)+(j==i+1) for i in range(len(bits)) for j in range(i+1,len(bits))
                 if bits[i]=='1' and bits[j]=='0')

def core(u):
    n=u.bit_length()
    return tuple(((1<<n)-(1<<(i+1)))+((u>>i&1)<<n) for i in range(n))+(((1<<n)-1)^u,)
def build(b,e):
    L=core(b-1);R=core(e-1);a=len(L);n=len(R)
    return L[:-1]+(L[-1]|1<<(a+n-1),)+tuple((((1<<a)-1)^((i==n-1)<<(a-1)))|
        (sum((row>>i&1)<<j for j,row in enumerate(R))<<a) for i in range(n))
def strong(T):
    n=len(T)
    for j in range(n):
        S=1<<j
        for _ in range(n):
            O=S
            for i in range(n):
                if S>>i&1:S|=T[i]
            if O==S:break
        if S!=(1<<n)-1:return False
    return True
def basic(N,b,e):
    u,v=b-1,e-1
    f,g=closed(u),closed(v)
    assert (f,g)==(recurrent(u),recurrent(v))
    assert N==b*e+f*g-(b&1)*(e&1)
    T=build(b,e)
    assert all(row>>i&1==0 and row>>len(T)==0 for i,row in enumerate(T))
    assert all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(len(T)) for j in range(i))
    assert strong(T)
    return len(T),f,g

def columns(T):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,len(D)):
        for i,w in enumerate(D[S]):
            m=T[i]&~S
            while w and m:
                x=m&-m;m-=x;D[S|x][x.bit_length()-1]+=w
    return D
def rooted(rows,r):
    D=columns(rows);P=columns(tuple(sum((row>>i&1)<<j for j,row in enumerate(rows)) for i in range(len(rows))))
    H=[1]+[sum(q) for q in D[1:]];M=len(H)-1
    return H[M],D[M][r],sum(H[S]*P[M^S][r] for S in range(M+1) if not(S>>r&1))
def F2(N,v,h,e,R,r,G):
    assert rooted(G,r)==(h,e,R)
    f=closed(v);assert f==recurrent(v)
    assert N==f*h-int(v%2==0)*e+(v+1)*R
    L=core(v);a=len(L)
    T=L[:-1]+(L[-1]|1<<(a+r),)+tuple((((1<<a)-1)^((i==r)<<(a-1)))|(g<<a) for i,g in enumerate(G))
    assert all(row>>i&1==0 and row>>len(T)==0 for i,row in enumerate(T))
    assert all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(len(T)) for j in range(i))
    assert strong(T)
    return len(T)
def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))
for k in (32,35):
    base=1<<k
    pp=next(v for v in range(base+1,base+500) if prime(v))
    q=next(v for v in range((base-1)//7+1,(base-1)//7+500) if prime(v))
    path=Path(__file__).with_name(f'band{k}.tsv')
    V=[tuple(map(int,L.split())) for L in path.read_text().splitlines()]
    extra=[]
    if k==35:
        for line in Path(__file__).with_name('band35_f2.tsv').read_text().splitlines():
            *p,z=line.split();N,*vals=map(int,p);rows=tuple(map(int,z.split(',')))
            print('F2 gap',N,'prime',prime(N),'order',F2(N,*vals,rows));extra.append(N)
    K=sorted([N for N,b,e in V]+extra)
    assert K==list(range(base+1,base+2*len(K),2))
    for N,b,e in V:
        row=basic(N,b,e)
        if N in (pp,7*q):
            print('target',N, 'N prime' if N==pp else f'7 times prime {q}', 'seam',b,e,'order, F values',row)
    print('PASS physical witnesses in complete detached band',base+1,base+2*len(K)-1)
