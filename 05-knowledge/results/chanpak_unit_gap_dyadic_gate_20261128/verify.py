"""Independent literal tournament checks of the two-arc rectangle and fixed-root block law."""
from itertools import product
from random import Random

def dp(A):
    n=len(A);lim=1<<n
    I=[sum(1<<i for i in range(n) if A[i]>>j&1) for j in range(n)]
    H=[1]+[0]*(lim-1);E=[[0]*n for _ in range(lim)]
    for V in range(1,lim):
        z=V
        while z:
            b=z&-z;z-=b;j=b.bit_length()-1;T=V^b
            if not T:E[V][j]=1
            else:
                w=T&I[j]
                while w:
                    q=w&-w;w-=q;E[V][j]+=E[T][q.bit_length()-1]
            H[V]+=E[V][j]
    return H,E,I

def profile(A,a):
    H,E,I=dp(A);_,S,_=dp(I);all_=len(H)-1
    return H[-1],E[-1][a],sum(row[a]*H[all_^T] for T,row in enumerate(S)),H[all_^(1<<a)]

def append(A,r,w):
    for c in w:
        b=len(A);A=[v|((1<<b) if c=='P' or a!=r else 0) for a,v in enumerate(A)]
        A.append((1<<r) if c!='P' else 0)
        if c=='M':r=b
    return A,r

def flip(A,i,j):
    A=list(A);A[i]^=1<<j;A[j]^=1<<i;return A

def extend(A,r,C):
    m=len(A);n=len(C)-1
    out=[a|(((1<<n)-1)<<m) if i!=r else a|((C[0]>>1)<<m) for i,a in enumerate(A)]
    for row in C[1:]:out.append(((row>>1)<<m)|(row&1)<<r)
    return out

def predicted(A,r,C):
    h,e,u,g=profile(A,r)
    j,b,w,d=profile(C,0);_,_,R=dp(C);_,s,v,_=profile(R,0)
    B=[r>>1 for r in C[1:]]; Q,_E,_I=dp(B);m=len(Q)-1
    q=sum(y*Q[m^i] for i,y in enumerate(Q))
    return (h-e)*d+e*s+(u-g)*(v-d)+g*(j-s),g*b,(u-g)*q+g*w,g*d

def random_tour(n,rng):
    A=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.randrange(2):A[i]|=1<<j
            else:A[j]|=1<<i
    return A

if __name__=='__main__':
    A,r=append([0],0,"PKKKMP")
    assert A==[66,124,121,113,97,65,0] and r==5
    B=flip(A,1,3)
    assert profile(B,r)==(15,0,36,9)
    assert profile(A,r)==(17,0,36,9)
    print('Flat single-arc seed',profile(B,r),profile(A,r))
    for old,gate in product((0,1),"PK"):
        got,root=append((B,A)[old],r,gate+"PPK")
        assert root==r
        assert profile(got,root)==(303+2*old,9,567 if gate=="P" else 495,9)
        print("Stored",old,gate,profile(got,root))
    for k in (1,2):
        low=303+495*2**k; all_=[]
        for old,gate in product((0,1),"KP"):
            T,j=append((B,A)[old],r,gate+"PPK"+"P"*k+"K")
            value=dp(T)[0][-1]
            assert value==low+2*old+(gate=='P')*72*2**k
            all_.append(value)
            ref,_=append(B,r,"KPPK"+"P"*k+"K")
            diff={(a,b) for a in range(len(T)) for b in range(a+1,len(T))
                  if (T[a]^ref[a])>>b&1}
            assert diff==({(1,3)} if old else set())|({(5,7)} if gate=='P' else set())
        print('Revealed',k,sorted(all_))
    rng=Random(20261128)
    for _ in range(180):
        X=random_tour(rng.randrange(1,6),rng);r=rng.randrange(len(X));Y=random_tour(rng.randrange(2,6),rng)
        value=profile(extend(X,r,Y),r)
        assert value==predicted(X,r,Y),(X,r,Y,value,predicted(X,r,Y))
    print("180 arbitrary rooted suffix-block tests")
