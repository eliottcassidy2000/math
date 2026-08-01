"""Does FC(2)'s mechanism transfer to LRC?

FC's engine: a counterexample pins a PERIOD to a RATIONAL value (Vol(Delta)), and
ARITHMETIC RIGIDITY of that value (transcendence) kills it.

LRC(n+1) analogue.  For integer speeds v_1..v_n and threshold 1/(n+1), the SAFE set is
   Safe = [0,1) \ union_i union_k ( k/v_i - 1/((n+1)v_i), k/v_i + 1/((n+1)v_i) ).
Every endpoint is a rational with denominator dividing (n+1) v_i, so mu(Safe) is a
RATIONAL with denominator dividing N := (n+1)*lcm(v_1..v_n).  Hence a GAP PRINCIPLE:

        either mu(Safe) = 0   (covering case)   or   mu(Safe) >= 1/N.

Test: (a) is mu always in (1/N)Z?  (b) how big is the gap in practice -- is it usable?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

def lcm(a, b): return a*b//gcd(a, b)

def safe_measure(vs, n1):
    """exact mu(Safe) for speeds vs, threshold 1/n1, via interval arithmetic on Fractions"""
    bad = []
    for v in vs:
        w = F(1, n1*v)
        for k in range(0, v+1):
            c = F(k, v)
            bad.append((max(F(0), c-w), min(F(1), c+w)))
    bad = [b for b in bad if b[0] < b[1]]
    bad.sort()
    merged = []
    for s, e in bad:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    cov = sum(e-s for s, e in merged)
    return F(1) - cov

print("LRC quantisation test:  mu(Safe) in (1/N)Z with N = (n+1)*lcm(v)?")
tests = [ (1,2,3), (1,2,3,4), (1,3,4,5,7), (1,2,3,4,5,6), (2,3,5,7), (1,4,5,6,7),
          (1,2,3,4,5,6,7), (3,5,7,11,13) ]
allok = True
for vs in tests:
    n = len(vs); n1 = n+1
    L = reduce(lcm, vs); N = n1*L
    mu = safe_measure(vs, n1)
    q = mu*N
    ok = (q.denominator == 1)
    allok &= ok
    print(f"  v={str(vs):22s} n={n} thr=1/{n1}: mu = {str(mu):>16}  N={N:7d}  mu*N = {str(q):>12}"
          f"  integral: {ok}   {'COVERED' if mu==0 else ''}")
print(f"\n  quantisation holds on all samples: {allok}")

print("\nIs the gap usable?  compare mu against the quantum 1/N:")
for vs in tests:
    n = len(vs); n1 = n+1
    L = reduce(lcm, vs); N = n1*L
    mu = safe_measure(vs, n1)
    if mu == 0:
        print(f"  v={str(vs):22s}: mu = 0 (covering) -- gap principle says nothing")
    else:
        print(f"  v={str(vs):22s}: mu = {float(mu):.6f}, quantum 1/N = {1/N:.3e}, "
              f"ratio mu/(1/N) = {float(mu*N):.0f}")
