"""Independent arithmetic replay of the signed entry rules and mod-24 shell."""
from math import isqrt
from random import Random

PAIRS = ((17,6), (19,4), (21,2), (23,0))

def v2(n):
    return (n & -n).bit_length() - 1

def p(n):
    w = 0
    while n:
        h = n // 2
        if n % 2 == 0:
            w += h + h % 2
        n = h
    return w

def H(x,y,j):
    return 432*x*y + 129*x + 258*y + 15*(x%2) + 30*(y%2) + 45*(p(x)+p(y)) + 135 + 2*j

def choose(u,U,N):
    v = U-u
    j = 2*N%3
    for b,c in PAIRS:
        x,y=32*u+b,32*v+c
        if H(x,y,j)%24 == N%24:
            return x,y,j
    raise RuntimeError(("suffix",u,v,N))

def band(N):
    L = isqrt(N//110592)
    g = 1 << ((L.bit_length()-1)//3)
    U = ((L+2+g-1)//g)*g-1
    begin = (U+1)//2
    end = begin + 4*isqrt(U*g)
    if end >= U:
        raise RuntimeError(("test target too small",N,U,g))
    def pair(u):
        x,y,j = choose(u,U,N)
        z=H(x,y,j)
        if z%24!=N%24 or x+y!=32*U+23:
            raise RuntimeError((N,U,u,x,y,j))
        return H(x,y,0),z,H(x,y,2)
    distances=[]
    for idx in (0,2):
        a,b=begin,end
        if not pair(a)[idx]>N>pair(b)[idx]:
            raise RuntimeError(("initial sign",N,a,b))
        while b-a>1:
            c=(a+b)//2
            if pair(c)[idx]>N:
                a=c
            else:
                b=c
        distances.extend((abs(pair(a)[1]-N),abs(pair(b)[1]-N)))
        if idx==0 and pair(a)[0]<=N or idx==2 and pair(b)[2]>N:
            raise RuntimeError(("whole windows",N,idx,a,b))
    for dist in distances:
        if dist**3>10**21*N:
            raise RuntimeError(("scale control",N,dist))

for m in range(3,10001,2):
    e = 1 if m%4==1 else -1
    x = (m-e)//4
    y = m-x
    a,b=(m-1)//2,(m+1)//2
    if sorted((a*a%m,a*b%m,b*a%m,b*b%m))!=[x,x,y,y] or y!=3*x+e:
        raise RuntimeError(("table",m,x,y))
    if x%2==0:
        left=2*x+e
        if 3*left-e!=2*y or 3*m+e!=4*y or y%2==0:
            raise RuntimeError(("bridge",m,x,y,left))
    elif y%2!=0:
        raise RuntimeError(("raw",m,x,y))

for u in range(1,81):
    for v in range(1,81):
        if (u+v)%2:
            for j in range(3):
                got={H(32*u+b,32*v+c,j)%8 for b,c in PAIRS}
                if got!={1,3,5,7}:
                    raise RuntimeError(("clock",u,v,j,got))
rng=Random(20261115)
for bits in (40,64,96,128,256,512):
    for _ in range(32):
        band((1<<(bits-1))|rng.getrandbits(bits-1)|1)
print("4999 modular tables, 9600 fixed-j clocks, 192 mod-24 brackets through 512 bits passed")
