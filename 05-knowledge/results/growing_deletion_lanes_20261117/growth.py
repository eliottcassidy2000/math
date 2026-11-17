from math import isqrt
from functools import cache

@cache
def profile(D,U):
    n=len(D);N=1<<n
    F=[1]+[0]*(N-1);G=[0]*N
    dp=[[0]*n for _ in range(N)]
    sp=[[0]*n for _ in range(N)]
    for i in range(n):
        dp[1<<i][i]=1
        if U>>i&1:sp[1<<i][i]=1
    for s in range(1,N):
        for i in range(n):
            h,p=dp[s][i],sp[s][i]
            F[s]+=h;G[s]+=p
            if not h:continue
            o=D[i]&~s
            while o:
                b=o&-o;o-=b;j=b.bit_length()-1
                dp[s|b][j]+=h;sp[s|b][j]+=p
    se=so=te=to=J=0
    for s in range(1,N):
        h=F[(N-1)^s]
        # T is the weighted end menu on the incontacts, complement of U.
        end=sum(dp[s][i] for i in range(n) if not U>>i&1)*h
        if s.bit_count()%2:so+=G[s]*h;to+=end
        else:
            se+=G[s]*h;te+=end
            J+=sum(sp[s][i] for i in range(n) if not U>>i&1)*h
    return F[-1],sum(F[s]*F[(N-1)^s] for s in range(N)),se,so,te,to,J

B0=(62,124,120,112,32,0,49)
B1=(126,60,120,112,96,0,34)

def host(a,t,B,U):
    w=([1,1] if a==3 else [1,0,1])+[0]*t
    n=len(w);k=n+1;N=k+len(B)+1
    A=[sum(1<<j for j in range(i+1,n))+((not v)<<n) for i,v in enumerate(w)]
    A+=[sum(v<<i for i,v in enumerate(w))]
    suffix=(1<<(N-1))-(1<<k)
    return tuple(
      [x|suffix|((i!=n)<<(N-1)) for i,x in enumerate(A)]
      +[(v<<k)|((not (U>>j&1))<<(N-1))for j,v in enumerate(B)]
      +[(1<<n)|(U<<k)])

def add(T,r,z):
    A=list(T)
    for j in range(z.bit_length()):
        i=len(A);b=z>>j&1
        A=[v+((k!=r or not b)<<i)for k,v in enumerate(A)]+[(1<<r)if b else 0]
    return tuple(A)

def f(z):
    a=1
    while z:
        if not z&1:
            u=z>>1;a+=u+(u&1)
        z>>=1
    return a

def predict(a,t,z,delta):
    h=(a<<t)-1
    return h*(536*z+477+9*f(z))+672*z+581+delta

def inverse(P):
    # The proper target inverse splits the binary tape and retains both weights.
    ans=[]
    for a in (3,5):
        for t in range(1,P.bit_length()):
            h=(a<<t)-1; A=536*h+672; g=9*h
            if 486*h+581>P:break
            for d in (0,2):
                Y=P+59*h+91-d
                k=max(1,isqrt(max(1,Y//A)).bit_length()-1);q=1<<k
                for r in range(q):
                    rev=0
                    for i in range(k):rev=2*rev+(r>>i&1)
                    K=A*q+g*(q-1-rev);C=A*(r+1)+g*(f(r)-1)
                    lo=max(0,(Y-C-2*g+K+g-1)//(K+g))
                    hi=(Y-C-g)//K
                    for u in range(lo,hi+1):
                        z=q*u+r
                        if predict(a,t,z,d)==P: ans.append((a,t,z,d))
    return ans

if __name__=='__main__':
    for B,U in ((B0,70),(B1,9)):
        print('tile',profile(B,U))
    for a,t in ((3,1),(5,1),(5,3)):
        for B,U,d in ((B1,9,0),(B0,70,2)):
            T=host(a,t,B,U);r=len(T)-1
            h=(a<<t)-1
            # gap-independent literal count via full path DP in profile
            C,R,E,O,_,_,_=profile(T,0)
            assert C==predict(a,t,0,d)
            D=tuple(x&((1<<r)-1) for x in T[:-1]);U=T[-1]
            rprof=profile(D,U)
            assert rprof[0]==9*h and rprof[3]==268*h+336
            assert rprof[2]+9*h==268*h+336+9*h
    assert predict(5,3,0,2)==19537==7*2791 and all(2791%i for i in range(2,isqrt(2791)+1))
    for x in (19537,7*10007,7*10009,7*30011,7*30013,1000000007,7*1000000007):
        print('inverse',x,inverse(x),flush=True)
    C=(2,4,9,3)
    assert profile(C,0)[:2]==(5,32)
    assert profile(tuple((1<<5)-(1<<(i+1)) for i in range(5)),0)[:2]==(1,32)
    for B,U,d in ((B1,9,0),(B0,70,2)):
        T=host(3,1,B,U);r=len(T)-1
        for z in (1,2,3):
            assert profile(add(T,r,z),0)[0]==predict(3,1,z,d)
