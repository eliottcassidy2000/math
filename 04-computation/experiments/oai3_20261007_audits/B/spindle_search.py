# For each of the 19 planes with sqrt3 in F: smallest Loeschian n (any parity) with sqrt(-(12n-1)) in L = Q(i, sqrt p: p|m),
# i.e. squarefree part of 12n-1 is a product of primes dividing m.  (THM-4558 generalized spindle, N = 3n.)
from sympy import factorint
def sqf(n):
    s = 1
    for p, e in factorint(n).items():
        if e % 2: s *= p
    return s
def loeschian(n):
    for p, e in factorint(n).items():
        if p % 3 == 2 and e % 2: return False
    return True
ms = [1,5,13,21,33,37,57,85,93,105,133,165,177,253,273,345,357,385,1365]
for m in ms:
    if m % 3: continue
    ps = set(factorint(m))
    hit = None
    for n in range(1, 20001):
        if not loeschian(n): continue
        s = sqf(12*n - 1)
        if set(factorint(s)) <= ps:
            hit = (n, 3*n, 12*n - 1, factorint(12*n - 1)); break
    print(f"m={m:5d}: smallest spindle n, N=3n, 4N-1 = {hit}")
