from math import isqrt

def ck(a,b):
    if a!=b:raise AssertionError((a,b))

def TT(n):return tuple((1<<n)-(1<<(j+1)) for j in range(n))
def reverse(a):return tuple(sum((x>>i&1)<<j for j,x in enumerate(a)) for i in range(len(a)))
def mark(a,q):return tuple(x|(((q>>i&1)^1)<<len(a)) for i,x in enumerate(a))+(q,)
def tag(u):return int(bin(u)[2:][::-1],2)
def U(u):
    b=[i for i in range(u.bit_length()) if u>>i&1]
    return sum((1<<(i-j-1))+(i-j==1) for i in b for j in b if i>j)
def build(I,O,B,arcs):
    n,m,k=len(I),len(O),len(B)
    r=[x|(((1<<(m+k))-1)<<n) for x in I]+[x<<n for x in O]+[((x<<(n+m))|(((1<<m)-1)<<n)) for x in B]
    for j,i in arcs:
        r[i]^=1<<(n+j);r[n+j]^=1<<i
    return tuple(r)
def strong(a):
    for i in range(len(a)):
        r=1<<i
        for _ in a:
            for j,x in enumerate(a):
                if r>>j&1:r|=x
        if r!=(1<<len(a))-1:return False
    return True
def tables(a):
    n=len(a);N=1<<n;D=[[0]*n for _ in range(N)];h=[1]+[0]*(N-1)
    for j in range(n):D[1<<j][j]=1
    for m in range(1,N):
        h[m]=sum(D[m])
        for j,v in enumerate(D[m]):
            if not v:continue
            z=a[j]&~m
            while z:
                b=z&-z;z-=b;D[m|b][b.bit_length()-1]+=v
    return h,D
def split(h):
    P=[]
    for m in range(len(h)):
        j=m;c=0
        while True:
            c+=h[j]*h[m^j]
            if not j:break
            j=(j-1)&m
        P.append(c)
    return P
def hpq(a,third=True):
    h=tables(a)[0];p=split(h)
    if not third:return h[-1],p[-1]
    f=len(h)-1
    return h[-1],p[-1],sum(h[j]*p[f^j] for j in range(f+1))
def pos(u,B):
    I=mark(TT(u.bit_length()),tag(u))
    return build(I,(0,),B,[(0,len(I)-1)])
def neg(u,B):
    return build(TT(u.bit_length()),(0,),B,[(0,i) for i in range(u.bit_length()) if tag(u)>>i&1])
def readp(u,h,p):
    return (u+1)*p+(u+1-(u&1)-U(u))*h
def readn(u,h,p):
    return u*p+(1-U(u))*h
P=(156,65,162,54,198,211,13,74);Q=(30,28,208,36,200,23,43,107)
A=(138,136,11,112,135,151,183,12);B=(100,237,48,37,107,192,12,93)
C=(216,29,81,52,32,135,186,30);D=(86,128,42,67,46,11,182,61)
def P2pos(u,Z):
    I=mark(TT(u.bit_length()),tag(u));_,q=hpq(I,False);h,p,r=hpq(Z)
    i=len(I)-1;H,E=tables(reverse(I));IP=split(H);f=len(H)-1
    L=sum(E[m][i]*IP[f^m] for m in range(f+1))
    Lm=sum((-1)**m.bit_count()*E[m][i]*IP[f^m] for m in range(f+1))
    return 2*q*p+2*L*r+2*Lm*h
def maxruns(a,slot):
    n=len(a);N=1<<n;D=[[-1]*n for _ in range(N)]
    for j in range(n):D[1<<j][j]=(slot>>j)&1
    for m in range(1,N):
        for j,v in enumerate(D[m]):
            if v<0:continue
            z=a[j]&~m
            while z:
                b=z&-z;z-=b;k=b.bit_length()-1
                D[m|b][k]=max(D[m|b][k],v+bool(slot&b and not(slot>>j&1)))
    return max(D[-1])

def main():
    for u,Z in [(1,(0,)),(2,TT(2)),(3,(2,4,1)),(5,TT(2)),(7,TT(3)),(11,(2,4,9,3))]:
        for g in (pos,neg):ck(strong(g(u,Z)),True)
        h,p=hpq(Z,False)
        ck(hpq(neg(u,Z),False)[0],readn(u,h,p))
        ck(hpq(pos(u,Z),False),(readp(u,h,p),P2pos(u,Z)))
    for z,v in [(P,(527,7448)),(Q,(535,7456)),(A,(199,4120)),(B,(247,4144)),(C,(387,6024)),(D,(419,6008))]:
        ck(hpq(z,False),v);ck(strong(z),True)
    ck(hpq(C)[2],46899);ck(hpq(D)[2],45603)
    for t in range(2,25):
        u=(1<<t)-1
        g=(1<<t)*6921-5867
        ck((readn(u,527,7448),readn(u,535,7456)),(g,g+16))
        u=(1<<t)+1
        g=8041*(1<<(t-1))+4319
        ck((readn(u,199,4120),readn(u,247,4144)),(g,g+72))
        ck(readp(u,387,6024),readp(u,419,6008))
        ck(P2pos(u,C)>P2pos(u,D),True) if t<7 else None
    for Z in (C,D):
        ck(hpq(pos(5,Z),False),(37305,P2pos(5,Z)))
    for p in (437077,437093,261631):
        ck(all(p%d for d in range(2,isqrt(p)+1)),True)
    for depth in range(1,5):
        for size in (1,depth,depth+1):
            Z=TT(size);mask=(1<<size)-1
            for _ in range(depth):Z=neg(1,Z);mask<<=2
            ck(maxruns(Z,mask),min(depth+1,size))
    print('12 nested-read run-rank controls pass.')
    print('Six raw two-reader comparisons, six strong source profiles, 23 pumped levels, third-cover slot response and three finite primality checks pass.')
if __name__=='__main__':main()
