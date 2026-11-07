from verify import inv
from sympy import isprime

for q in range(3,101):
    if not isprime(q):continue
    for a in range(1,7):
        cond=any(isprime(2**j*q+1) or q-1==2**j for j in range(1,a+1))
        assert bool(inv(2**a*q))==cond,(q,a)
assert inv(68)==[] and 137 in inv(136)
assert 17**2 in inv(272) and not isprime(273)
print('Dyadic repair criterion q<=97, a=1..6 with exceptional supplier: PASS')
