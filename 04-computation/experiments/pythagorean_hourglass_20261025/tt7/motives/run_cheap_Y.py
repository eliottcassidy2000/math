# The six cheap cone-reduced triple products Y_H (H = 2, 4, 5, 6, 9, 18; base of four points, O(q^2)) at primes
# q = 3 (mod 4) in [lo, hi).  One line per prime: q Y02 Y04 Y05 Y06 Y09 Y18.  mac-mini-2026-10-25-pyth.
import json, sys, time
from sympy import primerange
from cone_reduce import Y_H, tables
G = json.load(open('genuine37.json'))
cheap = [2, 4, 5, 6, 9, 18]
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, 'a') as fh:
    for q in primerange(lo, hi):
        if q % 4 != 3: continue
        T = tables(q)
        vals = [Y_H(q, G[i]['edges'], T=T)[0] for i in cheap]
        fh.write(f"{q} " + " ".join(map(str, vals)) + "\n"); fh.flush()
