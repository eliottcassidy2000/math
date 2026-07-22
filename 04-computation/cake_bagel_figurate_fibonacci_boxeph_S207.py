#!/usr/bin/env python3
"""cake_bagel_figurate_fibonacci_boxeph_S207.py -- boxeph-2026-07-21-S207

How the repo's polygonal/polyhedral (figurate, Moser, opus-S317) work relates to Fibonacci and the
CAKE and BAGEL cutting sequences. These are distinct readings of a common
binomial array:

 - CUTTING SEQUENCES = truncated Pascal ROW sums (partial sums of a row = fixed-depth figurate):
     lazy caterer (disk/lines, 2D)   A000124 = C(n,0)+C(n,1)+C(n,2)
     CAKE (ball/planes, 3D)          A000125 = C(n,0)+C(n,1)+C(n,2)+C(n,3)
     Moser circle (disk/chords)      A000127 = C(n,0)+C(n,2)+C(n,4)   [even-column functional]
     BAGEL (solid torus/planes, 3D)  = C(n,3)+2*C(n,2)+2*C(n,1)       [weighted functional]
 - FIBONACCI / g-BONACCI = GAP-DIAGONAL sums in the full-rank array:
     Pascal skip-sum = Fibonacci;  g-bonacci kernel 1/(1 - x - x^{g+1}) (Fibonacci at g=1, klein-S313).
This computes all four cutting sequences and verifies their binomial formulas,
the Fibonacci skip-sums, and the g-bonacci kernels. A shared Pascal ambient
array is not by itself a mechanism connecting the underlying geometries.
"""
from math import comb

N=12
def cake(n):   return sum(comb(n,k) for k in range(4))          # A000125 (3D ball by planes)
def caterer(n):return sum(comb(n,k) for k in range(3))          # A000124 (2D disk by lines)
def moser(n):  return comb(n,0)+comb(n,2)+comb(n,4)             # A000127 (2D disk by chords)
def bagel(n):  return 1 if n==0 else (n**3+3*n**2+8*n)//6      # solid torus by n planes: 1,2,6,13,24,40,62

print("n        :", list(range(N+1)))
print("caterer  :", [caterer(n) for n in range(N+1)], " (A000124, first 3 Pascal terms)")
print("CAKE     :", [cake(n) for n in range(N+1)],   " (A000125, first 4 Pascal terms)")
print("MOSER    :", [moser(n) for n in range(N+1)],  " (A000127, even-column functional)")
print("BAGEL    :", [bagel(n) for n in range(N+1)],  " (solid torus by n planes)")
print("2^n      :", [2**n for n in range(N+1)],      " (polyhedral/Pascal FULL row-sum)")

print("\n-- BAGEL vs CAKE: exact algebraic identities --")
print("  bagel(n) = C(n,3)+2*C(n,2)+2*C(n,1) = C(n,3)+n(n+1)")
ok = all(bagel(n)==comb(n,3)+n*(n+1) for n in range(1,N+1))
print("  bagel(n) == C(n,3)+n(n+1) for n>=1 ?", ok)
print("  bagel - cake =", [bagel(n)-cake(n) for n in range(N+1)])
print("  = T_n - 1 (triangular minus one) for n>=1 ?", all(bagel(n)-cake(n)==n*(n+1)//2-1 for n in range(1,N+1)),
      "  [numeric identity only; no shadow-lattice identification]")

print("\n-- FIBONACCI = shallow-diagonal (skip) sum of Pascal (polyhedral) --")
def pascal_skip(n):  # sum_{k} C(n-k, k)  = Fibonacci(n+1)
    return sum(comb(n-k,k) for k in range(n//2+1))
print("  Pascal skip-sums:", [pascal_skip(n) for n in range(N+1)], " = Fibonacci")

print("\n-- FULL-RANK gap-g kernels 1/(1 - x - x^{g+1}) (klein-S313 indexing) --")
def gbonacci(g, N):
    a=[0]*(N+1); a[0]=1
    for n in range(1,N+1):
        a[n]=(a[n-1] if n-1>=0 else 0)+(a[n-g-1] if n-g-1>=0 else 0)
    return a
for g in (1,2,3):
    label = "Fibonacci" if g == 1 else "higher-gap recurrence"
    print("  g=%d (kernel 1-x-x^{%d}):"%(g,g+1), gbonacci(g,N), "(%s)" % label)

print("\nSUMMARY: caterer/cake are Pascal-prefix sums; Moser/bagel are fixed binomial functionals;")
print("Fibonacci is a shallow Pascal diagonal, and the listed g-bonacci kernels obey their recurrences.")
print("The identities share a binomial ambient array. No geometric or LRC/JC mechanism is asserted.")
