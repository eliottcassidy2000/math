"""Audit: pre-registered P2 control clause - 'non-survivor controls m in {2,3,5,7,11,100,-100,999} have density < 0.10'
(density = fraction of primes 7<=p<=500, p not | m, with F_p(1/m) == 0 mod p)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in
P = primes_in(7, 500)
CTRL = [2, 3, 5, 7, 11, 100, -100, 999]
for fam in ['A', 'B', 'C', 'D', 'E', 'F', 'Apery', 'Domb', 's7', 's10']:
    u = seq(fam, 500)
    out = {}
    for m in CTRL:
        z = t = 0
        for p in P:
            if m % p == 0:
                continue
            t += 1
            xi = pow(m, -1, p)
            F = 0
            for k in range(p - 1, -1, -1):
                F = (F * xi + u[k]) % p
            z += (F == 0)
        out[m] = round(z / t, 3)
    bad = {m: d for m, d in out.items() if d >= 0.10}
    print(f'{fam:6s} control densities {out}  ->  >= 0.10: {bad}')
