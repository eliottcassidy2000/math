"""THM-4716 check: lower mechanical words s_k(phi) = floor((k+1)p/q + phi) - floor(kp/q + phi), k = 0..q-1, for the q
intercepts phi = t/q.  Claims: (i) consecutive intercepts differ by exactly one adjacent transposition 01 <-> 10;
(ii) the carries differ by exactly 2^i 3^j with i = the new position of the moved one and j = the number of ones to its
right; (iii) the carry is strictly monotone in the intercept, so sorted order = intercept order; (iv) hence
gcd(2^q - 3^p, carry) = 1 for every rotation (a prime dividing the clock and one element divides all elements, hence a
3-smooth difference)."""
from math import gcd, floor
import sys
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
def mech(p, q, t):   # intercept phi = t/q, exact integer arithmetic
    return [ (( (k+1)*p + t)//q) - ((k*p + t)//q) for k in range(q) ]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
checks = 0; nonmono = 0
for q in range(2, QMAX+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        D = 2**q - 3**p
        words = [mech(p, q, t) for t in range(q)]
        # all q words are distinct rotations of one necklace
        assert len({tuple(w) for w in words}) == q
        rots = {tuple(words[0][i:] + words[0][:i]) for i in range(q)}
        assert all(tuple(w) in rots for w in words)
        carries = [carry(w) for w in words]
        for t in range(q-1):
            a, b = words[t], words[t+1]
            diff = [k for k in range(q) if a[k] != b[k]]
            assert len(diff) == 2 and diff[1] == diff[0] + 1, (p, q, t, diff)
            k0 = diff[0]
            assert (a[k0], a[k0+1]) == (0, 1) and (b[k0], b[k0+1]) == (1, 0), (p, q, t)   # the one moves LEFT as phi grows
            i = k0; j = sum(a[k0+2:])
            assert carries[t] - carries[t+1] == 2**i * 3**j, (p, q, t, carries[t]-carries[t+1], i, j)
            checks += 1
        # monotone: carries strictly decreasing in t
        if any(carries[t] <= carries[t+1] for t in range(q-1)): nonmono += 1
        # coprimality for every rotation
        for c in carries: assert gcd(abs(D), c) == 1
        checks += 1
print(f"q <= {QMAX}: {checks} checks passed; non-monotone cases: {nonmono}; every rotation's carry coprime to the clock")
