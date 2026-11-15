from fractions import Fraction as F
from math import gcd, prod
from pathlib import Path
import csv, sys, hashlib
import check

HERE=Path(__file__).resolve().parent
PINNED_SHA256='5226fb6d6cb71e2c99dc50536744f044c31af55e173e3be9f5f722a39edf3bc2'

def rerun(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==PINNED_SHA256
    strong_profiles=set()
    with path.open() as f:
        for row in csv.DictReader(f):
            if row['strong']=='0': continue
            code=int(row['code'])
            T=tuple((code>>(8*i))&255 for i in range(8))
            h,p=check.profile(T)
            assert h==int(row['H'])
            strong_profiles.add((h,p))
    return sorted(strong_profiles)

def factor(n):
    k=2
    while k*k<=n:
        if n%k==0:
            yield k
            while n%k==0:n//=k
        k+=1
    if n>1:yield n

def main():
    profiles8=[tuple(map(int,l.split())) for l in (HERE/'profiles8.txt').read_text().splitlines()]
    if len(sys.argv)==2:
        assert rerun(Path(sys.argv[1]))==profiles8
        print('pinned order-eight class source replayed')
    ps=[]
    with (HERE/'profiles.txt').open() as f:
        for n,k in map(str.split,iter(f.readline,'')):
            ps.extend(tuple(map(int,f.readline().split())) for _ in range(int(k)))
    ps=sorted(set(ps+profiles8))
    primitive=[(h,p) for h,p in ps if gcd(h,p)==1]
    assert len(profiles8)==1567 and len(ps)==1721 and len(primitive)==982
    assert all(4*h<=p for h,p in ps)
    primes=sorted({k for h,p in ps for k in factor(gcd(h,p))})
    assert 2 not in primes
    M=prod(primes)
    log_upper=sum((F(2,(2*j+1)*7**(2*j+1)) for j in range(3)),F(0))+F(2,7**8)/(1-F(1,49))
    def mass(allowed):
        G=[prod(F(p**i,p**i-h**(i-1)) for h,p in allowed) for i in range(1,4)]
        return log_upper+sum((g-1)/((1<<i)*i) for i,g in enumerate(G,1))+(G[0]-1)/3584
    bound=mass(primitive)
    assert bound<F(486,1000)
    print(f'profiles all/primitive {len(ps)}/{len(primitive)}; certified conditional mass cap {float(bound):.9f}')
    print(f'M={M}, primes {primes}')
    m=3*5*13
    survivors=[(h,p) for h,p in ps if gcd(gcd(h,p),m)==1]
    coarser=mass(survivors)
    assert len(survivors)==1034
    assert coarser<F(4975,10000)
    print(f'm={m}, {len(survivors)} candidate profiles; certified conditional mass cap {float(coarser):.9f}')

if __name__=='__main__':
    main()
