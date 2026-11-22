from math import gcd
from itertools import product
def z(t):
    D,U,V=t
    return D+2*U,V,D+U
def iz(t):
    D,U,V=t
    return 2*V-D,D-V,U
def vp(z,p):
    if not z:return 10**9
    n=0
    while z%p==0:n+=1;z//=p
    return n
x,y,w=(1,0,0),(0,0,1),(1,0,-1)
T={0:0,1:1,2:1}
for k in range(-1,-60,-1):T[k]=T[k+3]-T[k+2]-T[k+1]
for j in range(30):
    a,b,c,d=sum(x),sum(y),x[1],y[1]
    assert a*d-b*c==T[-j-4]==w[1]
    for D,G,p in product(range(2,20,2),range(1,12,2),(3,5,7,11,13)):
        if gcd(D,G)==1:
            A,B=a*D+b*G,c*D+d*G
            assert vp(gcd(A,B),p)==min(vp(w[1],p),vp(D*w[2]-G*w[0],p))
    for h,Q in ((1,2),(1,32),(5,32),(37,352)):
        P=h*a+Q*(a+c);R=h*b+Q*(b+d);d0=gcd(P,R)
        assert R%2==1
        for N in (5,13,81,123,341,777,8125,2*P+R,6*P+3*R):
            B=N//(2*P+R)
            actual={G for G in range(1,B+1,2) if (N-R*G)%(2*P)==0}
            if N%d0: predicted=set()
            else:
                M=2*P//d0; G0=pow(R//d0,-1,M)*(N//d0)%M
                predicted=set(range(G0,B+1,M))
            assert predicted==actual
            if 2*P*(2*P+R)>d0*N: assert len(actual)<=1
    x,y,w=z(x),z(y),iz(w)
x=(2,0,1)
for _ in range(6):x=z(x)
assert x==(48,20,37) and (2*(sum(x)+x[1]),0,sum(x))==(250,0,105)
print('rank indexing, content and last-gap congruence controls pass')
