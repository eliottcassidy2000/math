from math import isqrt
from random import Random

def F(n):
    x,z=0,1
    for c in bin(n)[2:]:
        b=int(c)
        if not b:z+=x+x%2
        x=2*x+b
    return z

def H(u,w):
    return 81+142*(u+w)+256*u*w+5*(F(u)+F(w))+2*(w%2)

def require(ok):
    if not ok:raise ArithmeticError('control failed')

def template(G,r,z,d):
    v=int(format(r,'0'+str(G.bit_length()-1)+'b')[::-1],2)
    return d+(16384*G*G-5)*z+4096*G*G+1004*G+30+2048*G*r+40*v

def check(G,base,z,d):
    opt=[r for r in range(base+1,base+32,2) if template(G,r,z,d)%15==0]
    require(len(opt)>=1)
    r=opt[0];y=template(G,r,z,d)//15;x=y+z
    require(x>0 and y>0 and x%2==0)
    u,w=G*(16*x+5)+r,G*(16*y+4)+(G-1-r)
    require(H(u+4*G,w-4*G)-H(u,w)==d and H(w,u)-H(u,w)==2)
    # Multiply radius by G for exact integer estimate.
    top=15*isqrt(H(u,w))*2048
    denom=524288*2048
    expect=((8*z+2)*2048*G*G+2048*r*G+d)*524288*G
    require(abs(top-expect)<=200*denom*G*G)

rng=Random(1514)
for t in range(5,75):
    G=1<<t
    for _ in range(20):
        b=32*rng.randrange(G//32)
        z=rng.randrange(100)
        d=2*rng.randrange(-G*G//2,G*G//2+1)
        check(G,b,z,d)
# Check a rational octal point against the finitely many scaled band phases.
for j in range(80):
    q=192*8**j # (27+3/7)*7
    require(q//7%8==3 and q%7==3)
print('1400 exact phase/arc/radius controls and octal controls pass.')
