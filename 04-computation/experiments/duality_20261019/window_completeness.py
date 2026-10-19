"""THM (window completeness): the ladder L(q,p) = {2^i 3^(p-1-floor((i+1)p/q)) : 0 <= i <= q-2} is the full set of 3-smooth
numbers in [m, M] (m = min, M = max) iff M/m < 3 and M < min(2^(q-1), 3^p).  Test the criterion against brute force for all
coprime (q,p) with q <= QMAX, and record |q log_3 2 - p| for the complete cases."""
from math import gcd, log
import sys
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
theta = log(2)/log(3)
def smooth_in(lo, hi):
    S = []; a = 0
    while 2**a <= hi:
        b = 0
        while 2**a * 3**b <= hi:
            if 2**a * 3**b >= lo: S.append(2**a * 3**b)
            b += 1
        a += 1
    return sorted(S)
complete = []; mism = 0; total = 0
for q in range(3, QMAX+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        L = sorted(2**i * 3**(p-1-((i+1)*p)//q) for i in range(q-1))
        m, M = L[0], L[-1]
        full = (L == smooth_in(m, M))
        crit = (M < 3*m) and (M < min(2**(q-1), 3**p))
        total += 1
        if full != crit: mism += 1; print("MISMATCH", q, p, full, crit)
        if full: complete.append((q, p, round(q*theta - p, 4), round(M/m, 4)))
print(f"{total} coprime shapes q <= {QMAX}: criterion matches brute force in all but {mism}; complete cases: {len(complete)}")
for c in complete: print("  (q,p) =", c[:2], " q*log_3(2) - p =", c[2], " M/m =", c[3])
# sufficient condition: |q theta - p| < ?  list the max |q theta - p| among complete and min among incomplete with M/m<3
inc = [(q, p, round(q*theta - p, 4)) for q in range(3, QMAX+1) for p in range(1, q) if gcd(p,q) == 1 and (q, p) not in [(c[0], c[1]) for c in complete]]
print("max |q theta - p| among complete:", max(abs(c[2]) for c in complete))
print("incomplete shapes with |q theta - p| < 0.1:", [c for c in inc if abs(c[2]) < 0.1])
