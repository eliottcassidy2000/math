"""Exact hex-spiral labels and finite-state Collatz edge detector; no dependencies."""
from math import isqrt
def require(ok, detail="check failed"):
    if not ok: raise AssertionError(detail)

D=((0,-1),(-1,0),(-1,1),(0,1),(1,0),(1,-1))
def xy(n):
    if not isinstance(n,int) or n<1: raise ValueError("positive integer required")
    if n==1: return (0,0)
    r=(isqrt(12*n-3)-3)//6
    while 1+3*r*(r+1)<n:r+=1
    while r and n<=1+3*(r-1)*r:r-=1
    k=n-(2+3*r*(r-1)); side,j=divmod(k,r)
    q,s=r,0
    for dq,ds in D[:side]:q+=r*dq;s+=r*ds
    dq,ds=D[side]
    return q+j*dq,s+j*ds

def label(q,s):
    if not isinstance(q,int) or not isinstance(s,int): raise ValueError("integer coordinates required")
    r=max(abs(q),abs(s),abs(q+s))
    if not r:return 1
    starts=((r,0),(r,-r),(0,-r),(-r,0),(-r,r),(0,r))
    for side,((a,b),(dq,ds)) in enumerate(zip(starts,D)):
        j=(q-a)//dq if dq else (s-b)//ds
        if 0<=j<r and (q,s)==(a+j*dq,b+j*ds):return 2+3*r*(r-1)+side*r+j
    raise AssertionError((q,s))

def C(n):
    if not isinstance(n,int) or n<1: raise ValueError("positive integer required")
    return n//2 if n%2==0 else 3*n+1

def directions(a,b):
    return tuple((x,y) for x,y in ((a,b),(b,a)) if C(x)==y)

def certify(n,m,a,b):
    """Return LSD-first carries for m=a*n+b, or None. Positive inputs."""
    if not all(isinstance(x,int) for x in (n,m,a,b)) or n<1 or m<1 or a not in (2,3) or b not in (0,1): raise ValueError("invalid affine checker inputs")
    L=max(n.bit_length(),m.bit_length())
    c=b;cs=[c]
    for i in range(L):
        t=a*((n>>i)&1)+c
        if t%2 != ((m>>i)&1): return None
        c=t//2;cs.append(c)
    return cs if c==0 else None

def machine_directions(n,m):
    out=[]
    for x,y in ((n,m),(m,n)):
        # x -> y either x=2*y, or x odd and y=3*x+1
        if certify(y,x,2,0) is not None or (x%2 and certify(x,y,3,1) is not None):out.append((x,y))
    return tuple(out)

def neighboring_edges(R):
    out=set()
    for n in range(1,2+3*R*(R+1)):
        q,s=xy(n)
        for dq,ds in D:
            m=label(q+dq,s+ds)
            out.update(directions(n,m))
    return sorted(out)

def test():
    for n in range(1,100001):require(label(*xy(n))==n)
    for q in range(-100,101):
        for s in range(-100,101):require(xy(label(q,s))==(q,s))
    for n in range(1,501):
        for m in range(1,501):require(directions(n,m)==machine_directions(n,m),(n,m))
    for n in range(1,100001):
        require((n,C(n)) in machine_directions(n,C(n)))
        require(all(c in (0,1) for c in certify(n,2*n,2,0)))
        require(all(c in (0,1,2) for c in certify(n,3*n+1,3,1)))
    es=neighboring_edges(30)
    require(es==neighboring_edges(6)==neighboring_edges(60)==[(1,4),(2,1),(3,10),(28,14),(30,15),(32,16),(34,17),(36,18)])
    print('PASS: 100000 label round trips; 40401 coordinate round trips; 250000 pair checks; 100000 true transitions and carry ranges')
    print('Adjacent spiral Collatz edges:',es)
    print('13 -> 40 certificate:',certify(13,40,3,1))
if __name__=='__main__':test()
