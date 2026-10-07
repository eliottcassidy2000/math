# For D = -4n (n <= 3000): #primitive reduced diagonal forms (a,0,c) vs #genera (Cox Prop 3.11: 2^(mu-1)),
# cross-checked with #ambiguous reduced primitive forms (= #ambiguous classes = #genera, Gauss).
from math import gcd, isqrt
from collections import defaultdict
def odd_primes(n):
    s=set(); m=n
    while m%2==0: m//=2
    q=3
    while q*q<=m:
        while m%q==0: s.add(q); m//=q
        q+=2
    if m>1: s.add(m)
    return s
res = defaultdict(lambda: [0,0])
for n in range(1, 3001):
    D = -4*n
    diag = 0; amb = 0
    for a in range(1, isqrt(4*n//3)+2):
        for B in range(0, a+1, 2):          # B even for D = -4n
            if (B*B - D) % (4*a): continue
            c = (B*B - D)//(4*a)
            if c < a: continue
            if gcd(gcd(a,B),c) != 1: continue
            if B == 0 or B == a or a == c:
                amb += 1
                if B == 0: diag += 1
    r = len(odd_primes(n))
    if n % 4 == 3: mu = r
    elif n % 4 in (1,2): mu = r+1
    elif n % 8 == 4: mu = r+1
    else: mu = r+2
    gen = 2**(mu-1)
    assert amb == gen, (n, amb, gen)
    key = n % 8
    res[key][0] += 1
    res[key][1] += (diag == gen)
for k in sorted(res): print(f"n = {k} mod 8: diag == #genera in {res[k][1]} of {res[k][0]} cases")
