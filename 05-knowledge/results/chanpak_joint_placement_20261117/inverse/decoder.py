from fractions import Fraction as Q
from math import ceil,isqrt

def F(x):
    v,p,s=0,0,1
    for b in bin(x)[2:]:
        d=int(b)
        if d==0:s+=v+p
        v,p=2*v+d,d
    return s

def check(a,b):
    if a!=b:raise AssertionError((a,b))

def rev(r,G):
    return int(format(r,f'0{G.bit_length()-1}b')[::-1],2)

def src(u,w):
    return 81+142*(u+w)+256*u*w+5*(F(u)+F(w))

def assemble(G,r,x,y):
    u=G*(16*x+5)+r;w=G*(16*y+4)+(G-1-r)
    return src(u,w),src(u+4*G,w-4*G)-src(u,w)

def inv_one(N,delta,G):
    ep=Q(1,25)
    a=1024*G+Q(5,32*G);b=Q(15,32*G)
    c=-delta-4096*G*G+20*G-30+b*(1-10*G)
    k=b/a;l=c/a
    lo=max(1,ceil((l-ep)/(1-k)),ceil((-l-ep)/(1+k)))
    hi=N
    if lo>hi:return None
    def p(S):return k*S+l
    def H0(S):return 64*(S*S-p(S)**2)+142*S+81
    if H0(hi)+21*hi<N:return None
    while hi>lo:
        md=(hi+lo)//2
        if H0(md)+21*md<N:lo=md+1
        else:hi=md
    S=lo
    if H0(S)-21*S>N:return None
    d=round(p(S))
    if abs(d-p(S))>ep:return None
    z,r2=divmod(d,16*G)
    if not(3<=r2<=2*G-1 and r2%4==3):return None
    r=(r2-1)//2
    a0,m0=divmod(S+1,G)
    if m0:return None
    y,ex=divmod(a0-16*z-10,32)
    if ex:return None
    x=y+z
    if min(x,y)<1 or x&1:return None
    return (G,r,x,y) if assemble(G,r,x,y)==(N,delta) else None

def inverse(N,delta):
    if N<1 or delta%2:return []
    out=[];G=4
    while 256*G*G<=N:
        q=inv_one(N,delta,G)
        if q:out.append(q)
        G*=2
    return out

def window_one(N,E,G,allow_zero=False):
    A=16384*G*G-5
    Bmax=6144*G*G-1004*G-10
    ymax=isqrt(N//(65536*G*G))
    Z=max((E+Bmax)//A,(15*ymax+E)//A)
    dmax=(16*Z+2)*G
    def quad(S):return 64*S*S+142*S+81
    def lower(pred):
        L,R=1,N+1
        while L<R:
            S=(L+R)//2
            if pred(S):R=S
            else:L=S+1
        return L
    lo=lower(lambda S:quad(S)+5*S+10>=N)
    hi=lower(lambda S:quad(S)-64*dmax*dmax+10>N)
    result=[]
    for S in range(lo,hi):
        T,rem=divmod(S+1,G)
        if rem or (T-10)%16:continue
        modz=(-(T-10)//16)%4
        for z in range(-Z+((modz+Z)%4),Z+1,4):
            if z==0 and not allow_zero:continue
            y,rem=divmod(T-16*z-10,32)
            x=y+z
            if rem or min(x,y)<1 or x&1:continue
            a0,b0=16*x+5,16*y+4
            base=F(a0)+F(b0)+(G-1)*a0-2
            radius=(abs(16*z+1)+1)*G+2
            sq0=quad(S)+5*base-N
            sqmin=max(0,(sq0-5*radius+63)//64)
            sqmax=(sq0+5*radius)//64
            if sqmax<sqmin:continue
            d0=isqrt(sqmin)
            if d0*d0<sqmin:d0+=1
            d1=isqrt(sqmax)
            if z<0:d0,d1=-d1,-d0
            r0=max(1,(d0-16*G*z)//2)
            r1=min(G-1,(d1-16*G*z-1)//2)
            for r in range(r0+(r0%2==0),r1+1,2):
                NN,delta=assemble(G,r,x,y)
                if NN==N and abs(delta)<=E:result.append((G,r,x,y))
    return result

def replay():
    checks=0
    for G in (4,8,16,32,64,128):
        for r in range(1,G,2):
            for x,y in ((2,1),(2,2),(2,11),(12,6),(14,19),(2**72,2**72),(2**72,2**70),(2**70,2**72)):
                N,D=assemble(G,r,x,y)
                A=16384*G*G-5
                B=4096*G*G+1004*G+30+2048*G*r+40*rev(r,G)
                check(D,15*y-A*(x-y)-B)
                check(inv_one(N,D,G),(G,r,x,y))
                check(inv_one(N+2,D,G),None)
                checks+=1
    for D in (-505,0,4,58,693,4000):
        G=32;r=5
        if D%2:continue
        A=16384*G*G-5;B=4096*G*G+1004*G+30+2048*G*r+40*rev(r,G)
        z=next(z for z in range(1,16) if (D+A*z+B)%15==0);y=(D+A*z+B)//15;x=y+z
        N,D0=assemble(G,r,x,y);check(D,D0);check((G,r,x,y) in inverse(N,D),True)
    tested=0
    for G in (4,8,32,256,2**10,2**16,2**20):
        for r in (1,G//2-1,G-1):
            for z in (0,1,3,15):
                y=(G*G*600+z)%((1<<100)+1);x=y+z
                N,D=assemble(G,r,x,y);check((G,r,x,y) in window_one(N,abs(D),G,True),True)
                tested+=1
    print(tested,'free-shift native window controls pass.')
    print(checks,'direct forward/inverse and adjacent ghost checks; full-height external controls pass.')
if __name__=='__main__':replay()
