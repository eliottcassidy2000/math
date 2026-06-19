import sys
from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a

def setgcd(S): return reduce(gcd, S)

def run(k, V):
    floor = Fraction(1, k+1); best=None
    for S in combinations(range(1,V+1), k):
        if setgcd(S)!=1: continue
        M,t=M_exact_fast(S)
        if M<=floor: continue
        if best is None or M<best[0]: best=(M,S,t)
    return best, floor

import sympy
for k,V in [(9,20),(10,22),(11,22)]:
    best,floor=run(k,V)
    M,S,t=best; g=M-floor
    om=len(sympy.factorint(k-1))
    print(f"k={k} V={V} omega(k-1)={om} sqrt(k)={k**0.5:.2f}: S={S}")
    print(f"   M={M} (~{float(M):.7f}) g*k^2={float(g*k*k):.5f} a={float(level_a(M,k)):.4f} q={M.denominator}/{(level_a(M,k)*(k+1)-1) if False else ''} t*={t}", flush=True)
