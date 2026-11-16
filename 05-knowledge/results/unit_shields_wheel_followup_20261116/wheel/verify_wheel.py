"""Exact replay using only stdlib. Run with --6424 for the large cap sieve."""
from functools import cache
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import json, sys

CERT = json.loads(Path(__file__).with_name('prime_certs.json').read_text())
@cache
def certified(p):
    if p < 10000:
        return p >= 2 and all(p % t for t in range(2,isqrt(p)+1))
    a, fs = CERT[str(p)]
    return (prod(q**k for q,k in fs)==p-1 and pow(a,p-1,p)==1 and
            all(certified(q) and gcd(pow(a,(p-1)//q,p)-1,p)==1 for q,k in fs))

def val(q,i):
    return ((1<<q)-1) if i==0 else ((1<<(q<<(i-1)))+1)//((1<<(1<<(i-1)))+1)

layers = {
    11:((23,89),(683,),(397,2113),(353,2931542417),(229153,119782433,43872038849)),
    13:((8191,),(2731,),(53,157,1613),(858001,308761441),(78919881726271091143763623681,),
        (928513,23877647873,18558466369,715668470267111297,21316654212673)),
    17:((131071,),(43691,),(137,953,26317),(17,354689,2879347902817),
        (383521,2368179743873,373200722470799764577),
        (5441,335631827046798245410603730138717057,63406006407727721042109834220642811713))
}

def q_cert():
    result={}
    for q, packets in layers.items():
        counts=[]
        for i, pp in enumerate(packets):
            n=q<<i
            assert len(set(pp))==len(pp) and prod(pp)==val(q,i)
            count=0
            for p in pp:
                assert certified(p)
                if pow(2,n,p)==1 and pow(2,n//q,p)!=1 and (i==0 or pow(2,n//2,p)!=1):
                    count+=1
            counts.append(count)
        result[q]=counts
    assert result=={11:[2,1,2,2,3],13:[1,1,3,2,1,5],17:[1,1,3,2,3,3]}
    # Presence of one primitive prime at order 352 is enough for q=11; 5304641 is such a prime.
    p=5304641
    assert all(p%j for j in range(2,isqrt(p)+1))
    assert val(11,5)%p==0 and pow(2,352,p)==1 and pow(2,176,p)!=1 and pow(2,32,p)!=1
    print('Certified primitive counts:',result)

def actual_covers():
    # Dyadic phases nested, then place each mixed-order clock on a distinct
    # odd CRT slot of the remaining dyadic fiber.
    dyadic=(3,5,17,257,65537)
    for q in (11,13):
        labels=[(p,1<<i,(1<<(i-1))-1) for i,p in enumerate(dyadic,1)]
        mixed=[(p,q<<i) for i,pp in enumerate(layers[q]) for p in pp]
        if q==11:mixed.append((5304641,352))
        for j,(p,d) in enumerate(mixed):
            t=j+q*((31-j)*pow(q,-1,d//q)%(d//q)) if d!=q else j
            labels.append((p,d,t))
        M=prod(p for p,d,t in labels)
        K=sum(pow(2,-t,p)*(M//p)*pow(M//p,-1,p) for p,d,t in labels)%M
        if not K%2:K+=M
        assert (K&1) and K>max(p for p,d,t in labels)
        assert all(any((K*pow(2,e,p)-1)%p==0 for p,d,t in labels) for e in range(32*q))
        print('Absolute CRT cover',32*q,'coefficient',K)

def brier_120():
    orders={3:2,5:4,7:3,11:10,13:12,17:8,31:5,41:20,61:60,151:15,
            241:24,331:30,1321:60,61681:40,4562284561:120}
    # Factorization and primality. Largest trial divisor is below 68 thousand.
    assert 15*prod(orders)==2**120-1
    assert all(all(p%j for j in range(2,isqrt(p)+1)) for p in orders)
    assert all(pow(2,d,p)==1 and all(pow(2,d//q,p)!=1 for q in [2,3,5] if d%q==0)
               for p,d in orders.items())
    F=set(range(8))
    pair=[{0,2},{1,3},{4,6},{5,7}]
    sing=[{i} for i in F]
    sheet=[{0,1,2,3},{4,5,6,7}]
    best=120;cases=0
    for d4,d8,d3,d5,d15,d12,d20,d24,d40 in product(pair,sing,sheet,sheet,sheet,pair,pair,sing,sing):
        alive=F-d4-d8
        need=sum(max(0,(3-(t in d3)-(t in d12)-(t in d24))*
                       (4-(t in d5)-(t in d20)-(t in d40))-1-(t in d15)) for t in alive)
        best=min(best,need);cases+=1
    assert best==11 and cases==262144 and best>2+2+1
    print('Period 120 both-sign obstruction:',cases,'relaxed arrangements; minimum holes',best)

def cap6424():
    cutoff=10**9
    def capacity(primes):
        m=prod(primes);caps=[]
        for i in range(4):
            d=m<<i
            X=2 if not i else -(1<<(1<<(i-1)))
            if len(primes)==1:
                v=(X**m-1)//(X-1)
            else:
                p,q=primes
                v=(X**m-1)*(X-1)//((X**p-1)*(X**q-1))
            g=d*(1+(d%2)); found=[]
            for p in range(g+1,cutoff,g):
                if pow(2,d,p)!=1 or v%p:continue
                if any(pow(2,d//q,p)==1 for q in primes):continue
                if i and pow(2,d//2,p)==1:continue
                if any(p%j==0 for j in range(2,isqrt(p)+1)):continue
                found.append(p)
                while v%p==0:v//=p
            k=0;lo=1;a=cutoff+((1-cutoff)%g)
            while lo*(a+g*k)<=v:
                lo*=a+g*k;k+=1
            caps.append(len(found)+k)
            print('d',d,'known',found,'total capacity',caps[-1],flush=True)
        return caps
    r=7
    c=capacity((73,));s=capacity((11,73))
    assert c==[3,3,6,9] and s==[24,24,48,96]
    assert (11-r)*(73-sum(c))==208>192==sum(s)
    print('Period 6424 has no cover: at least 208 isolated cells required; at most 192 available.')

if __name__=='__main__':
    q_cert(); actual_covers(); brier_120()
    if '--6424' in sys.argv:cap6424()
