from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod

def rows(w):
    m=len(w)
    return tuple(sum(1<<(j+1) for j in range(m) if w[j]) if i==0 else
                 (int(not w[i-1]) | sum(1<<j for j in range(i+1,m+1)))
                 for i in range(m+1))

def profile(r):
    @lru_cache(None)
    def start(s,i):
        t=s^(1<<i)
        return 1 if not t else sum(start(t,j) for j in range(len(r)) if t&r[i]&(1<<j))
    @lru_cache(None)
    def h(s):
        return 1 if not s else sum(start(s,i) for i in range(len(r)) if s&(1<<i))
    full=(1<<len(r))-1
    return h(full),sum(h(s)*h(full^s) for s in range(full+1))

def closed(w):
    n=len(w)+1
    h,q=1,1<<n
    for i,bit in enumerate(w):
        if not bit:continue
        for j in range(i+1,len(w)):
            if w[j]:continue
            r=j-i-1
            h+=(1<<r) if r else 2
            q+=(1<<(n-2-r))*(3**r+1)
    return h,q

def glue(r,w,u):
    n,m,k=len(r),len(w),len(u)
    wm=((1<<m)-1)<<n
    um=((1<<k)-1)<<(n+m)
    return (tuple(a|wm|(um if i else 0) for i,a in enumerate(r))
            +tuple((a<<n)|um for a in w)
            +tuple((a<<(n+m))|1 for a in u))

def main():
    counts=[]
    seen=set()
    for n in range(3,16):
        cur=set()
        for mid in product(range(2),repeat=n-3):
            w=(1,*mid,0)
            value=closed(w)
            if not n>8 and value!=profile(rows(w)):
                raise AssertionError((n,w,value))
            cur.add(value)
        counts.append(len(cur))
        seen.update(cur)
    Z15=prod(Fraction(Q,Q-1) for H,Q in seen)
    tail=Fraction(65,128)*Fraction(81,16)*Fraction(2,3)**16
    E=sum(Fraction(1,2**n-1) for n in range(1,13))+Fraction(2,8191)
    if not (Z15<Fraction(3041,2500) and tail<Fraction(4,1000)
            and Fraction(3041,2500)*Fraction(1000,996)<Fraction(611,500)
            and E<Fraction(1607,1000)
            and Fraction(111,500)*(Fraction(1607,1000)+Fraction(8,3)*Fraction(111,500))<Fraction(49,100)):
        raise AssertionError((float(Z15),tail,E))
    print('Profiles by order n=3..15:',counts)
    print('Union',len(seen),'and Z15 < 3041/2500, log tail < 4/1000, full Z < 611/500')
    A=(2,4,1)
    B=rows((1,0,1,0))
    check=0
    for m in range(5):
        for bits in product(range(2),repeat=max(m-1,0)):
            host=(1,*bits) if m else ()
            t=sum(bit<<(m-1-i) for i,bit in enumerate(host))
            hr,_=closed(host)
            for U,W in ((A,B),(B,A)):
                hu,qu=profile(U)
                hw,qw=profile(W)
                expected=(hr+1)*hu*hw+(qu-hu)*((1+t)*qw-hw)
                actual=profile(glue(rows(host),W,U))[0]
                if actual!=expected: raise AssertionError((host,actual,expected))
                check+=1
    print('Literal variable-host glue checks:',check)
if __name__=='__main__': main()
