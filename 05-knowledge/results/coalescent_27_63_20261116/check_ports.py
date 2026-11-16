from itertools import product
from math import isqrt

def square(n):
    return n>=0 and isqrt(n)**2==n

def shortcut(n):
    return (3*n+1)//2 if n%2 else n//2

def trajectory(n,steps):
    ans=[n]
    for _ in range(steps):
        n=shortcut(n);ans.append(n)
    return ans

def meeting(start):
    ans=[start]
    while ans[-1]!=182:
        ans.append(shortcut(ans[-1]))
    return ans

bases=[meeting(27),meeting(63)]
assert [len(x)-1 for x in bases]==[10,8]
shifts=[1024,2304]
for base,shift in zip(bases,shifts):
    # Check linear parity preservation at every pre-meeting step. Since
    # the slope is even, the base parity itinerary persists for every j.
    for u,v in zip(base,base[1:]):
        assert shift%2==0
        shift=(3 if u%2 else 1)*shift//2
        assert v==shortcut(u)
    assert shift==6561
assert set(bases[0][:-1]).isdisjoint(bases[1][:-1])

for j in range(80):
    x,y,D=27+1024*j,63+2304*j,9+320*j
    z=37+4096*j*(1+25*j)
    assert x+z==(D-1)**2 and y+z==(D+1)**2
for h in range(20):
    t=1+160*h
    assert 9*(t*t-1)%320==0 and (176*t*t+9)%5==0
    j=9*(t*t-1)//320
    x,y,z=27+1024*j,63+2304*j,(176*t*t+9)//5
    assert x<z<y and x+z==(8*t)**2 and y+z==(10*t)**2

def nb(x,N):
    return {z for z in range(1,N+1) if z!=x and square(x+z)}
assert sorted((a,b) for a,b in product(nb(8,41),nb(9,41)) if square(a+b))==[(41,40)]
assert 8+9==17 and 8*8+9*9==145 and 2*8*9==144
assert 1729==1+12**3==9**3+10**3
assert 1472==2*32**2-24**2
assert 1729-1472==257==1+48**2-2*32**2
assert all(257%k for k in range(2,isqrt(257)+1))
assert not square(1472+1729)
assert 14912+1472==128**2 and 14912+1729==129**2
assert (1472+49,49+120,120+1729)==(39**2,13**2,43**2)

w,u=2,1
for _ in range(8):
    w,u=3*w+4*u,2*w+3*u
    v=w-1
    a,b,c,d=v-1,(v*v+1)//2,(v*v-1)//2,v
    assert w*w-2*u*u==2 and c*c+v*v==b*b
    assert (a+b,b+c,c+d)==(u*u,v*v,u*u)

F=[0,1]
for _ in range(40):
    F.append(F[-1]+F[-2])
for k in range(5,36,3):
    u,v=F[k+1],F[k]
    p=(v*v-1)//2;q=p+1;x=u*u-q;y=u*u-p
    assert len({p,q,x,y})==4 and 0<p<q<x<y
    assert (x+q,q+p,p+y)==(u*u,v*v,u*u)
    assert y-x==1 and (x+y)%8==7
    assert x+y==v*F[k+3]+2*(-1)**k
print("Affinely preserved Collatz itineraries, quadratic/zero overhead routes, 8/9 onset, taxicab ports, and golden routing replay passed.")
