#!/usr/bin/env python3
"""
The a=d rigidity: consecutive case proved + the consecutive-loneliness formula (boxeph-S117)

PROVED: consecutive AP {a,..,a+11} (d=1), a>=2 => M >= a/(2a+11) > 1/13 (witness t=1/(2a+11)).
So d=1 tight => a=1 => {1..12}. Closed form (witness = maximizer, verified):
  M({a,a+1,...,a+n-1}) = a/(2a+n-1);  at a=1 gives 1/(n+1) (LRC tight value for {1..n}).
Necessary conditions (i)-(iv) narrow tight primitive APs to d=1 or d>=17; d>=17 spread APs are
loose (a=2,d=41: M=222/455) but escape (i)-(iv) -- general a=d open.
"""
from math import gcd
from fractions import Fraction as Fr

def Mstar(V, QMAX=800):
    b = Fr(0)
    for q in range(2, QMAX+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = min(min((a*v) % q, q - ((a*v) % q)) for v in V)
            if Fr(m, q) > b: b = Fr(m, q)
    return b

print('closed form M({a,..,a+11}) = a/(2a+11) (n=12):')
for a in range(1, 9):
    C = list(range(a, a+12)); M = Mstar(C); f = Fr(a, 2*a+11)
    print(f'  a={a}: M={M!s:>7}  a/(2a+11)={f!s:>7}  match={M==f}')
print()
print('general n form M({a,..,a+n-1}) = a/(2a+n-1):')
for n in [6, 8, 10]:
    for a in [1, 3]:
        C = list(range(a, a+n)); M = Mstar(C); f = Fr(a, 2*a+n-1)
        print(f'  n={n},a={a}: M={M!s:>7}  a/(2a+n-1)={f!s:>7}  match={M==f}  (a=1 => 1/(n+1)={Fr(1,n+1)})')
print()
print('=> {1,..,n} is the UNIQUE tightest consecutive block (a=1 minimizes a/(2a+n-1)).')
print('d>=17 residual (loose but escapes (i)-(iv)): a=2,d=41 => M =', Mstar([2+41*k for k in range(12)], 900))
