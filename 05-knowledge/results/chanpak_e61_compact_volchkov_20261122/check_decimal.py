from collections import Counter
from math import gcd

def e(p,k=4):
    m=10**k
    return pow(2,p-1,m)*(pow(2,p,m)-1)%m

def parity(x):
    return ''.join('EO'[int(d)%2] for d in f'{x:04d}')

def v(n,p):
    a=0
    while n%p==0:
        n//=p;a+=1
    return a

assert (e(61),parity(e(61)))==(2176,'EOOE')
for p in range(7,507,2):
    assert e(p)==e(p+500)
    if p%4==1: assert e(p)%20==16
    else:
        assert e(p)%100==28
        E=(1<<(p-1))*((1<<p)-1)
        assert v(8*E+1,5)==2*(1+v((p+1)//4,5))
assert len({e(p) for p in range(7,507,2) if p%4==3})==11
assert len({e(p) for p in range(7,507,2) if p%4==1})==125
# There are 125 residues on the first sheet, because modulo 5^4 one
# unit class fixes the first of four base-5 digits.
for p in range(7,107,2): assert parity(e(p))==parity(e(p+100))
for d in (20,50):
    p=next(p for p in range(7,507,2) if parity(e(p))!=parity(e(p+d)))
    print('period witness',p,d,e(p),e(p+d),parity(e(p)),parity(e(p+d)))
a=Counter(parity(e(p)) for p in range(7,107,2))
assert len(a)==6 and sum(a.values())==50
print('odd-exponent parity period histogram',dict(sorted(a.items())))
print('prime-residue filter', dict(sorted(Counter(
    parity(e(p)) for p in range(7,107,2) if gcd(p,100)==1).items())))
