from math import isqrt
from random import Random

def F(n):
    x,z=0,1
    for a in bin(n)[2:]:
        b=int(a)
        if not b:z+=x+x%2
        x=2*x+b
    return z

def H(u,v):
    return 81+142*(u+v)+256*u*v+5*(F(u)+F(v))

def B(G,r):
    rev=int(format(r,'0'+str(G.bit_length()-1)+'b')[::-1],2)
    return 4096*G*G+1004*G+30+2048*G*r+40*rev

def brackets(lo,hi,yes):
    while lo<hi:
        a=(lo+hi)//2
        if yes(a):hi=a
        else:lo=a+1
    return lo

def decode(N,G,Z,dest=False,upper=False):
    if G<(16*Z+12)**2:raise ValueError('use bounded small-scale enumeration')
    N-=2*upper
    T=lambda S:64*S*S+142*S+81
    S=brackets(1,N,lambda S:T(S)+10*S>=N)
    if abs(T(S)-N)>10*S or (S+1)%G:return None
    k,t=divmod((S+1)//G-10,16)
    if t or k<3:return None
    A=16384*G*G-5
    den=2*A+15
    # 15 k - den*z = 2(Delta+B_r).
    z0=max(1,-(-(15*k-2*(6145*G*G+1044*G+30))//den))
    z1=min(Z,(15*k-2*4095*G*G)//den)
    if z0>z1:return None
    if z0!=z1:raise ArithmeticError('multiple coarse imbalances')
    z=z0
    if (k-z)%2:return None
    y=(k-z)//2;x=y+z
    if y<1 or x%2:return None
    def addresses(i):
        r=2*i+1
        return G*(16*x+5)+r+dest*4*G,G*(16*y+4)+G-1-r-dest*4*G
    i=brackets(0,G//2-1,lambda i:H(*addresses(i))<=N)
    if H(*addresses(i))!=N:return None
    r=2*i+1
    d=15*y-A*z-B(G,r)
    if abs(d)>G*G:return None
    return x,y,r,d

def check():
    R=Random(815)
    Z=20
    for t in range(20,58):
        G=1<<t
        for _ in range(2):
            z=R.randint(1,Z)
            base=32*R.randrange(G//32)
            d=2*R.randrange(-G*G//2,G*G//2+1)
            A=16384*G*G-5
            r=next(a for a in range(base+1,base+32,2) if (d+A*z+B(G,a))%15==0)
            y=(d+A*z+B(G,r))//15;x=y+z
            for dest in (False,True):
                a,b=G*(16*x+5)+r+4*dest*G,G*(16*y+4)+G-1-r-4*dest*G
                for upper in (False,True):
                    if decode(H(a,b)+2*upper,G,Z,dest,upper)!=(x,y,r,d):
                        raise ArithmeticError('joint target failed')
    print('304 target-only exact decodes pass (both endpoints and both reflections).')
if __name__=='__main__':check()
