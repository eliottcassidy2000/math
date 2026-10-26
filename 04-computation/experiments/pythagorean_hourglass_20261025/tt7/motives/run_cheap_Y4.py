# Y02, Y05, Y09, Y18 only (the independent cheap triple products) at primes q = 3 (mod 4) in [lo, hi).
import json, sys
from sympy import primerange
from cone_reduce import Y_H, tables
G = json.load(open('genuine37.json'))
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, 'a') as fh:
    for q in primerange(lo, hi):
        if q % 4 != 3: continue
        T = tables(q)
        v = [Y_H(q, G[i]['edges'], T=T)[0] for i in (2, 5, 9, 18)]
        fh.write(f"{q} {v[0]} {-v[0]} {v[1]} {v[1]} {v[2]} {v[3]}\n"); fh.flush()
