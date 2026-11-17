from math import gcd
from random import Random

def ck(x):
    if not x:raise ArithmeticError('control failed')
def F(x):
    c=1
    while x:
        y,b=divmod(x,2)
        if not b:c+=y+y%2
        x=y
    return c
def rev(x,t):
    s=0
    for _ in range(t):s=(s<<1)+(x&1);x>>=1
    return s
def profiles(m):
    q=2**m;h=q+1;P=4*q+2*(3**m-1);L=2*q-1;d=2*q-2
    return q,h,P,L,d
def coefficients(m):
    q,h,P,L,d=profiles(m)
    return 5*h,256*P,147*P-10*h,142*h+256*L
def H(m,u,w):
    q,h,P,L,d=profiles(m)
    lx=P*u+L
    dx=d+h*(F(u)-1)+(P-2*h)*u+q*(u%2)
    return 91*h+5*(dx+h*(F(w)-1))+142*(lx+h*w)+256*lx*w
def displacement(m,G,r,z,y):
    s=G-1-r;x=y+z
    u=G*(16*x+5)+r
    w=G*(16*y+4)+s
    return H(m,u+4*G,w-4*G)-H(m,u,w)
def B(m,G,r):
    h,P,A,C=coefficients(m)
    t=G.bit_length()-1
    return 16*P*G*G+4*P*G+8*P*G*r+6*h-4*h*G+8*h*rev(r,t)-4*G*(A-C)
def select(m,G,z,D,j=0):
    ck(m%2==0 and m>=2 and G&(G-1)==0 and z>0 and D%2==0)
    h,P,A,C=coefficients(m)
    J=P*G*G
    g=gcd(h,P);mod=3*h
    if (D-3072*G)%g:return None
    ck(G>=2*mod//g and abs(D)<=J)
    t=G.bit_length()-1
    E=8*(P*G+h*(-1)**(t-1));M=mod//g
    ck(gcd(E,mod)==g)
    # r=2a+1 and modular B(r)=B(0)+E*r.
    rem=-D-(64*J-h)*z-B(m,G,0)-E
    ck(rem%g==0)
    a=(rem//g)*pow(2*E//g,-1,M)%M
    r=2*(a+j*M)+1
    if not(1<=r<G):return None
    num=D+(64*J-h)*z+B(m,G,r)
    ck(num%mod==0)
    y=num//mod;x=y+z
    ck(y>0 and not x%2)
    ck(displacement(m,G,r,z,y)==D)
    return r,x,y
rng=Random(341)
tot=0
for m in (2,4,6,8,12,16,20,24):
    h,P,*_=coefficients(m);g=gcd(h,P)
    G=1
    while G<6*h//g:G*=2
    for gg in (G,2*G,8*G):
        J=P*gg*gg
        for z in (1,2,7):
            for j in range(20):
                D=2*rng.randrange(-J//2,J//2)
                D+=2*(((3072*gg-D)*pow(2,-1,g))%g) if g>1 else 0
                if abs(D)>J:continue
                ans=select(m,gg,z,D)
                ck(ans is not None)
                if select(m,gg,z,D,1) is not None:
                    ck(select(m,gg,z,D,1)[2]>ans[2])
                tot+=1
print('PASS',tot,'graded physical arithmetic repairs')
for D in (0,2,-2,4,100,-100):
    m,G=2,256
    r,x,y=select(m,G,1,D)
    s=G-1-r;u=G*(16*x+5)+r;w=G*(16*y+4)+s
    print(m,G,D,r,x,y,H(m,u,w),H(m,u+4*G,w-4*G))

def tt(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def joined(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def lit(m,u,w,reflect):
    C=(54,32,2,119,6,20,55);ct=26
    n=m+2;R=list(tt(n));R[0]^=1<<(n-1);R[-1]=1;R=tuple(R)
    def word(x):
        s=bin(x)[2:] if x else '';return tt(len(s)),sum(int(ch)<<i for i,ch in enumerate(s))
    def dual(A,s):
        return tuple(sum(((A[j]>>i)&1)<<j for j in range(len(A))) for i in range(len(A))),((1<<len(A))-1)^s
    T,p=word(u);X,a=joined(T,R),p|(1<<len(T));Y,b=word(w)
    if reflect:X,a,Y,b=Y,b,X,a
    Y,b=dual(Y,b);c=a|(ct<<len(X))|(b<<(len(X)+len(C)))
    Z=joined(joined(X,C),Y);N=len(Z)
    G=tuple(x|((1-(c>>i&1))<<N) for i,x in enumerate(Z))+(c,)
    n=len(G);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for t in range(1,1<<n):
        for i,v in enumerate(dp[t]):
            if not v:continue
            b=G[i]&~t
            while b:
                s=b&-b;b-=s;dp[t|s][s.bit_length()-1]+=v
    return sum(dp[-1])
for m,u,w in [(1,0,2),(1,1,2),(2,0,4),(2,1,2),(2,3,2),(3,0,2)]:
    v=H(m,u,w)
    ck(lit(m,u,w,False)==v and lit(m,u,w,True)==v+2)
print('PASS 12 independent literal tournament path counts')
