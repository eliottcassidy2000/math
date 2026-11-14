from fractions import Fraction
from math import gcd, lcm

def check(ok):
    if not ok:
        raise AssertionError('sieve control failed')

for s in range(1,5):
    L=2**(2**s)-1
    for a in range(1,101):
        W=2**a+1
        if a%(2**s):
            j=(a&-a).bit_length()-1
            F=2**(2**j)+1
            check(L%F==0 and W%F==0 and F>=3)
        else:
            check(gcd(W,2*L)==1)

# Literal availability needs no sampling here: transitive TT_a, cycle C3
# with singleton pads, and order-four H=5/Q=32 with singleton pads.
rays = [(1,2),(3,12),(5,32)]
prefix=7
L=lcm(*(h+(q<<m) for h,q in rays for m in range(prefix)))
T=sum(Fraction(h+2,q<<(prefix-1)) for h,q in rays)
check(T<Fraction(1,2))
for h,q in rays:
    for m in range(25):
        W=h+(q<<m)
        g=gcd(2*L,W)
        if m<prefix:
            check(g==W and 1+h<g)
        else:
            actual = Fraction((g==1)+g*((h+1)%g==0),W)
            check(actual<=Fraction(h+2,W))
print('Fermat phase and literal SCC-ray sieve controls pass.')
