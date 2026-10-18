import sys
from fractions import Fraction
from small_n4 import SmallN4
from lie_voa_ss import *
c = Fraction(sys.argv[1]) if len(sys.argv) > 1 else Fraction(-9)
Dmax = Fraction(sys.argv[2]) if len(sys.argv) > 2 else Fraction(3)
P = 2147483647
alg = SmallN4(c, P)
basis = alg.pbw_basis(int(Dmax) + 1)
keys = sorted(k for k in basis if 0 < k[0] <= Dmax)
for k in keys:
    vecs = find_singular(alg, k[0], k[1])
    if vecs:
        print(f"weight {k[0]}, charge {k[1]}: dim singular = {len(vecs)} (block dim {len(basis[k])})")
        for v in vecs[:2]:
            items = sorted(v.items(), key=lambda kv: (alg.li(kv[0]), kv[0]))
            print("   Li-min part:", [(mo, co) for mo, co in items if alg.li(mo) == alg.li(items[0][0])][:8])
