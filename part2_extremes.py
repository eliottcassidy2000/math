from fractions import Fraction as F
from math import gcd
from lrc import measure_intersection

# The overlap depends on coprime (a,b)=(vi/g, vj/g). Let's nail the EXTREMES exactly.

# EXTREME 1: b = 2a (i.e. vj = 2 vi). Then B_{2a} arcs sit so each B_a arc fully...
# Compute the maximal-overlap pair ratio as function of n.
print("MAXIMAL overlap: pair (a, 2a)  [vj = 2 vi], ratio = measure * n^2/4")
for n in range(3, 16):
    r = measure_intersection([1,2], n) * F(n*n,4)
    o = measure_intersection([1,2], n)
    print(f"  n={n:2d}: overlap={str(o):>8}  ratio={float(r):.4f}  (n+? : (n+1)/2 /1? )  (n-1)/2/(...)= ")
print()
# Conjecture ratio formula for (1,2): r = (n+?)/(...).  Fit:
print("Fit ratio(1,2): values n=3..15:")
vals=[]
for n in range(3,16):
    r = measure_intersection([1,2], n) * F(n*n,4)
    vals.append((n,r))
    print(f"   n={n}: r={r}")
print()

# EXTREME 2: the MINIMUM overlap pairs from the n=7 scan were (1,6),(1,5) etc ->
# b = n-1 (i.e. vj == -vi mod n) ... actually (1,6) with n=7: 6 = -1 mod 7.
# And (1,5): 5 ?  Let's examine pairs (1, n-1) and (1, n-2)... general (1,b).
print("Overlap of pair (1, b) as a function of b, n=7:")
n=7
for b in range(2, 25):
    if gcd(1,b)!=1: continue
    o = measure_intersection([1,b], n)
    r = o*F(n*n,4)
    print(f"   (1,{b:2d}) b mod 7={b%7}: overlap={str(o):>8} ratio={float(r):.4f}")
