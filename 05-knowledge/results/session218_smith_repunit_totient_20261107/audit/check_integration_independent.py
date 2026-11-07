"""Independent table audit using prime suppliers and trial division; no producer imports."""
import json, sys
from pathlib import Path
from math import isqrt

def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def factors(n):
    out=[];p=2
    while p*p<=n:
        while n%p==0: out.append(p); n//=p
        p+=1
    if n>1:out.append(n)
    return out

def ds(n,b):
    a=[]
    while n:a.append(n%b);n//=b
    return a[::-1]
root=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else root.parent/'totient/integration_receipt.json'
if not source.exists(): source=root.parent/'session218-scaling/integration_receipt.json'
rows=json.loads(source.read_text())['rows']
expected=[(14,1,29),(26,1,53),(34,2,137),(38,5,1217),(50,1,101),(62,7,7937),(68,1,137),(74,1,149),(76,4,1217)]
assert [(r['n'],r['dyadic_repair']['a'],r['dyadic_repair']['preimage']) for r in rows]==expected
for r in rows:
    n=r['n'];fs=factors(n);F=sum(fs);b=n-F+1
    assert r['F']==F and r['base']==b and r['digits']==ds(n,b)
    assert sum(ds(n,b))==sum(sum(ds(p,b)) for p in fs)==r['factor_digit_weight']
    phi=n
    for p in set(fs): phi=phi//p*(p-1)
    assert phi==r['phi_n'] and r['prime_half']==prime(n//2)
    a=r['dyadic_repair']['a'];m=n*2**a;x=r['dyadic_repair']['preimage']
    assert x==m+1 and prime(x) # phi(x)=m; phi(y)<=y−1 proves this is least preimage.
    if a>1:
        odd=n;v=0
        while odd%2==0:odd//=2;v+=1
        assert prime(odd)
        # The proved supplier criterion excludes every earlier multiplier, including a=0.
        for t in range(a):
            assert not any(prime(2**j*odd+1) or odd-1==2**j for j in range(1,v+t+1))
print('PASS all9 rows, exact earliest positive repairs, digit/factor weights, and least preimages.')
for q in (19,31):
    excluded=[];j=1
    while not prime(2**j*q+1):
        excluded.append((j,factors(2**j*q+1)));j+=1
    print('q=',q,'excluded suppliers=',excluded,'first prime supplier j=',j,'p=',2**j*q+1)
