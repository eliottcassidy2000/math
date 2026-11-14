"""Two distinct physical three-point windows, with independent outer wings."""
from math import isqrt
from itertools import product
import sys

H9 = (
    ((14,28,24,16,1),4,2,2,1,0),
    ((14,28,24,16,1),22,1,3,0,1),
    ((14,28,24,16,1),13,3,1,2,0),
)
H45 = (
    ((22,44,24,49,34,5),5,0,0,0,0),
    ((14,12,56,48,35,3),3,0,0,0,0),
    ((14,12,56,48,35,3),9,0,0,0,0),
)

def rc(u, n):
    result=0
    for _ in range(n):
        result=2*result+1-(u&1)
        u>>=1
    return result

def attach(rows,U,a,b,x,y):
    n=len(rows)
    N=a+n+b
    M=(1<<N)-1
    Y=rc(y,b)
    contact=0
    for i in range(a):
        if x>>(a-1-i)&1:contact|=1<<i
    contact|=U<<a
    for i in range(b):
        if Y>>(b-1-i)&1:contact|=1<<(a+n+i)
    R=[]
    for i in range(N):
        if a<=i<a+n:
            r=(rows[i-a]<<a)|(M^((1<<(a+n))-1))
        else:
            r=M^((1<<(i+1))-1)
        R.append(r|((1-((contact>>i)&1))<<N))
    return tuple(R)+(contact,)

def paths(rows):
    n=len(rows)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        D[1<<i][i]=1
    for S in range(1,1<<n):
        for i in range(n):
            r=D[S][i]
            if not r:continue
            t=rows[i]&~S
            while t:
                b=t&-t
                D[S|b][b.bit_length()-1]+=r
                t^=b
    return D

def profile(rows):
    n=len(rows)-1
    assert all(not (r>>i&1) for i,r in enumerate(rows))
    assert all(((r>>j&1)+(rows[j]>>i&1)==1) for i,r in enumerate(rows) for j in range(i+1,len(rows)))
    F=(1<<n)-1
    R=tuple(r&F for r in rows[:-1]);U=rows[-1]
    end=paths(R)
    start=paths(tuple(F^r^(1<<i) for i,r in enumerate(R)))
    H=[sum(v) for v in end];H[0]=1
    h=H[-1]
    f=sum(H[A]*H[A^F] for A in range(F+1))
    P=sum(start[A][i]*H[A^F] for i in range(n) if U>>i&1 for A in range(F+1))
    Q=sum(end[A][i]*H[A^F] for i in range(n) if not (U>>i&1) for A in range(F+1))
    S=sum(end[F][i] for i in range(n) if U>>i&1)
    T=sum(start[F][i] for i in range(n) if not (U>>i&1))
    K=sum(paths(rows)[-1])
    return h,f,P,Q,S,T,K

def g(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:h+=x+(x&1)
    return h

def value(x,y,Z):
    h,f,P,Q,S,T,K=Z
    return K+h*(g(x)+g(y)-2)+Q*x+P*y+f*x*y+T*(x&1)+S*(y&1)

def contexts(kind,x,y,u):
    inner=attach(*(kind[u]))
    F=(1<<(len(inner)-1))-1
    return attach(tuple(r&F for r in inner[:-1]),inner[-1],x.bit_length(),y.bit_length(),x,y)

def solve(N,profile_center):
    """Every exact wing solution in a branch, not just a congruence hit."""
    if N<1:return
    h,f,P,Q,S,T,K=profile_center
    for rev in (False,True):
        if rev: P,Q,S,T=Q,P,T,S
        for x in range(1,isqrt(N//f)+1):
            for d in (-2,0,2):
                D=K+d+h*(g(x)-1)+Q*x+T*(x&1)
                lam=f*x+P
                Y=(N-D)//lam
                min_y=x+int(rev)
                if Y<min_y:continue
                M=2
                while lam*M<=h*(Y//M+1):M*=2
                for v in range(min(M,Y+1)):
                    lo=max(0,(min_y-v+M-1)//M)
                    hi=(Y-v)//M
                    while lo<=hi:
                        z=(hi+lo)//2
                        y=M*z+v
                        w=D+lam*y+h*(g(y)-1)+S*(y&1)
                        if w==N:
                            yield (y,x,1+d//2) if rev else (x,y,1+d//2)
                            break
                        if w<N: lo=z+1
                        else: hi=z-1

def join(A,B):
    """Ordinal join of the old cores, retaining one common apex."""
    n,m=len(A)-1,len(B)-1
    F,G=(1<<n)-1,(1<<m)-1
    return (tuple((r&F)|(G<<n)|((r>>n)<<(n+m)) for r in A[:-1])
            +tuple(((r&G)<<n)|((r>>m)<<(n+m)) for r in B[:-1])
            +(A[-1]|(B[-1]<<n),))

def merge(A,B):
    h,f,P,Q,S,T,K=A
    j,g,p,q,s,t,L=B
    J=(K-(h-S))*j+h*(L-(j-t))+P*q+(h-S)*(j-t)
    return h*j,f*g,P*g+h*p,f*q+Q*j,h*s,T*j,J

if __name__=='__main__':
    for kind,Z in [(H9,(9,1344,427,245,9,9,91)),(H45,(45,432,129,258,15,30,137))]:
        for i,spec in enumerate(kind):
            W=profile(attach(*spec))
            print(W,attach(*spec))
            assert W==(*Z[:-1],Z[-1]+2*(i-1))
        for x,y in [(1,1),(2,3),(3,2),(4,5),(5,4)]:
            for i in range(3):
                N=value(x,y,Z)+2*(i-1)
                assert sum(paths(contexts(kind,x,y,i))[-1])==N
                assert (x,y,i) in set(solve(N,Z))
        if len(sys.argv)>1:
            N=int(sys.argv[1])
            for x,y,i in solve(N,Z):print('WITNESS',N,x,y,i,contexts(kind,x,y,i))
    print('PASS: literal profiles, movable windows, and exact height inverse.')
    for a,b in product(H9,H45):
        A,B=attach(*a),attach(*b)
        assert profile(join(A,B))==merge(profile(A),profile(B))
    print('PASS: all nine literal mixed shared-apex profile products.')
