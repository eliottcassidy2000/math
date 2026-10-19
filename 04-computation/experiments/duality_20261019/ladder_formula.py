"""THM-4717 check: the consecutive differences of the sorted balanced free cycle of 3x + (2^q - 3^p) are exactly
d_i = 2^i * 3^(p - 1 - floor((i+1) p / q)),  i = 0..q-2  (each power of 2 once; the 3-exponents follow the Beatty sequence).
Also: Fibonacci-shape minima and Fibonacci clocks."""
from math import gcd, floor
import sys, json
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
def mech(p, q, t): return [(((k+1)*p + t)//q) - ((k*p + t)//q) for k in range(q)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
ok = 0
for q in range(2, QMAX+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        elems = sorted(carry(mech(p, q, t)) for t in range(q))
        diffs = sorted(elems[i+1] - elems[i] for i in range(q-1))
        formula = sorted(2**i * 3**(p - 1 - ((i+1)*p)//q) for i in range(q-1))
        assert diffs == formula, (q, p, diffs[:5], formula[:5])
        # the intercept order: carries decrease with t; difference at step t is d_(i) with i = k_1 - 1, k_1 = (q-1-t) p^-1 mod q
        pinv = pow(p, -1, q)
        for t in range(q-1):
            k1 = ((q-1-t)*pinv) % q
            i = k1 - 1
            assert carry(mech(p,q,t)) - carry(mech(p,q,t+1)) == 2**i * 3**(p-1-((i+1)*p)//q), (q,p,t)
        ok += 1
print(f"ladder formula verified for all {ok} coprime (q,p), q <= {QMAX}: differences = {{2^i 3^(p-1-floor((i+1)p/q)) : 0 <= i <= q-2}}, step t moves the one at k_1 = (q-1-t) p^-1 mod q")
# Fibonacci shapes (q,p) = (F_(n+1), F_n): clocks and minima of the balanced (bulb) cycle
F = [1, 1]
for _ in range(14): F.append(F[-1] + F[-2])
mins = []; clocks = []
for n in range(1, 11):
    p, q = F[n], F[n+1]
    if q > 400: break
    D = 2**q - 3**p
    elems = [carry(mech(p, q, t)) for t in range(q)]
    mins.append(min(elems)); clocks.append(D)
    print(f"(q,p)=({q},{p}): clock 2^{q}-3^{p} = {D}; balanced cycle min {min(elems)}, max {max(elems)}, max-min = {max(elems)-min(elems)}")
print("Fibonacci clocks:", clocks)
print("golden Collatz minima:", mins)
json.dump(dict(clocks=[str(c) for c in clocks], minima=[str(m) for m in mins]), open("ladder_formula.json","w"))
