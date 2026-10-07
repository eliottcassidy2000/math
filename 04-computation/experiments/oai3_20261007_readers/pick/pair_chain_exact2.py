"""Exact (dyadic) Haar values for the lag-1 pair chain, n <= NMAX:
P(merged by T-step n) and the exact flip density P(e_n odd | alive at n)."""
import sys
from fractions import Fraction
from pair_chain_verify import step
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
dist = {(0, 1): 1}; merged = 0; dexp = 0
for n in range(NMAX):
    if n == NMAX - 1 or n % 4 == 0 or n >= 24:
        alive_w = sum(dist.values()); odd_w = sum(w for (k, num), w in dist.items() if num & 1)
        flipdens = Fraction(odd_w, alive_w)
    new = {}
    bs = (0,) if n < 4 else (0, 1)
    if n >= 4:
        merged *= 2; dexp += 1
    for (k, num), w in dist.items():
        for b in bs:
            s2 = step(k, num, b)[:2]
            new[s2] = new.get(s2, 0) + w
    merged += new.pop((0, 0), 0)
    dist = new
    if n + 1 >= 24 or (n + 1) % 4 == 0:
        P = Fraction(merged, 2 ** dexp)
        print(f"n={n+1:3d}  P(merged by n) = {P}  = {float(P):.7f}   flip density at step {n}: {float(flipdens):.5f}   live states {len(dist)}", flush=True)
