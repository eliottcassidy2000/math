from fractions import Fraction as Q
from contextlib import redirect_stdout
from io import StringIO
with redirect_stdout(StringIO()):
    from verify import coefficients,H,B
from random import Random

def decode(N,Delta,m,G):
    h,P,A,C=coefficients(m)
    a=4*P*G+Q(h,32*G);b=Q(3*h,32*G)
    c=b*(1-10*G)-Delta-16*P*G*G-6*h+4*h*G+4*G*(A-C)
    k,l=b/a,c/a
    K=96*(1<<m)+71+142*(2**(m+1)-1)
    p=lambda s:k*s+l
    phi=lambda s:P*(s*s-p(s)**2)/4+Q(A+C,2)*s+Q(A-C,2)*p(s)+K
    ceil=lambda x:-((-x.numerator)//x.denominator)
    u=max(1,ceil((l-Q(1,100))/(1-k)),ceil((-l-Q(1,100))/(1+k)))
    v=N
    if u>v:return None
    while u<v:
        s=(u+v)//2
        if phi(s)+Q(P*s,30)>=N:v=s
        else:u=s+1
    S=u
    if abs(phi(S)-N)>Q(P*S,30):return None
    z,b=divmod(round(p(S)),16*G)
    r=(b-1)//2
    if not(1<=r<G and r%2 and b==2*r+1 and z>=1):return None
    y,rem=divmod(S+1-G*(16*z+10),32*G);x=y+z
    if rem or min(y,x)<=0 or x%2:return None
    U=G*(16*x+5)+r;V=G*(16*y+4)+G-1-r
    if H(m,U,V)!=N or H(m,U+4*G,V-4*G)-N!=Delta:return None
    return x,y,r

if __name__=='__main__':
    rng=Random(320)
    for m in range(2,23,2):
        for t in (2,6,20,40):
            G=1<<t
            for _ in range(5):
                r=2*rng.randrange(G//2)+1
                z=rng.randrange(1,300)
                y=(rng.getrandbits(80)<<1)+(z%2)+2;x=y+z
                U=G*(16*x+5)+r;V=G*(16*y+4)+G-1-r
                N=H(m,U,V);Delta=H(m,U+4*G,V-4*G)-N
                if decode(N,Delta,m,G)!=(x,y,r):raise ArithmeticError('exact joint failure')
    print('220 exact changed-flank joint inversions pass.')
