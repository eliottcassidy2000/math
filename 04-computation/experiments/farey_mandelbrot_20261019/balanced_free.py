"""HYP-9334 test: for every Christoffel (balanced) necklace w of shape (K, X) = (q, p), gcd(2^q - 3^p, B_r) = 1 for every
rotation r of w, where B_r is the carry of F_r(x) = (3^p x + B_r)/2^q.  Consequence: a balanced parity word is an integer
cycle of 3x+d only if (2^q - 3^p) | d (THM-4484's free case): no sporadic balanced cycles."""
from math import gcd, floor
import sys
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 80
def christoffel(p, q): return [floor(k*p/q) - floor((k-1)*p/q) for k in range(1, q+1)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
worst = (1, None); count = 0; nontriv = []
for q in range(2, QMAX+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        D = abs(2**q - 3**p); w = christoffel(p, q)
        for i in range(q):
            r = w[i:] + w[:i]; g = gcd(D, carry(r)); count += 1
            if g > worst[0]: worst = (g, (p, q, i))
            if g > 1: nontriv.append((p, q, i, g))
print(f"checked {count} rotations of Christoffel words with q <= {QMAX}; max gcd(|2^q-3^p|, B_r) = {worst[0]} at {worst[1]}")
print("rotations with gcd > 1 (first 20):", nontriv[:20], "count", len(nontriv))
# also: random unbalanced necklaces of the same shapes: how often is gcd > 1 (control)
import random
random.seed(2); ctrl = 0; tot = 0
for q in range(5, 40):
    for p in range(2, q-1):
        if gcd(p, q) != 1: continue
        D = abs(2**q - 3**p)
        for _ in range(3):
            w = [1]*p + [0]*(q-p); random.shuffle(w); tot += 1
            if gcd(D, carry(w)) > 1: ctrl += 1
print(f"control (random necklaces of the same shapes): gcd > 1 in {ctrl}/{tot}")
