# The two-parameter Legendre triple product behind H18 (three cone vertices over three triples of a 4-point base):
#   Y(q) = sum_{s != t, s,t not in {0,1}} phi(s-1) phi(t-1) L(s) L(t) L(t/s),   L(l) = sum_w phi(w(w-1)(w-l)),
# computed in O(q^2).  Test whether Y lies in base + new directions and whether it carries a new direction.
import sys
import numpy as np
from sympy import primerange
def L_table(q):
    chi = np.full(q, -1, dtype=np.int64); chi[0] = 0; chi[[(x * x) % q for x in range(1, q)]] = 1
    w = np.arange(q); base = chi[(w * (w - 1)) % q]
    L = np.array([int((base * chi[(w - l) % q]).sum()) for l in range(q)], dtype=np.int64)
    return chi, L
def Y(q):
    chi, L = L_table(q)
    inv = np.zeros(q, dtype=np.int64); inv[1:] = [pow(int(x), q - 2, q) for x in range(1, q)]
    s = np.arange(2, q)
    tot = 0
    for t in range(2, q):
        r = (t * inv[s]) % q                    # t/s
        v = chi[(s - 1) % q] * L[s] * L[r]
        v[s == t] = 0
        tot += int(chi[(t - 1) % q] * L[t] * v.sum())
    return tot
if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    for q in primerange(lo, hi):
        if q % 4 == 3: print(f"q={q} Y={Y(q)}", flush=True)
