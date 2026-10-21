import gzip, os
def _open_census(n):
    p=f'census_n{n}.txt'
    return open(p) if os.path.exists(p) else gzip.open(p+'.gz','rt')
# labelled multiplicities N_n(h) = sum n!/|Aut| over classes with H=h; check n!/gcd(n!,h) | N_n(h); C-spectrum facts
import math
from collections import Counter
for n in range(3,10):
    R=[tuple(map(int,l.split())) for l in _open_census(n)]
    f=math.factorial(n); N=Counter()
    for H,C,s,A,AA,c3 in R: N[H]+=f//A
    bad=[h for h in N if N[h]%(f//math.gcd(f,h))]
    exact=sum(1 for h in N if N[h]==f//math.gcd(f,h))
    Cs=Counter(r[1] for r in R)
    print(f'n={n}: H values {len(N)}; divisibility n!/gcd(n!,h) | N(h) violations: {len(bad)}; values with N(h) equal to that bound: {exact}; '
          f'C range 0..{max(Cs)}, C values missing: {[c for c in range(max(Cs)+1) if c not in Cs]}')
