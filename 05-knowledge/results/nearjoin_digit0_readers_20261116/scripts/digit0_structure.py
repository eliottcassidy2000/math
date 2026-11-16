#!/usr/bin/env python3
"""digit0_structure.py: structure of digit-0 vertices (B_u = H) in STRONG tournaments of orders 3-8:
   out-degree distribution; whether the first-order condition "no Hamiltonian path of S - x starts at u for every
   out-neighbour x of u" is equivalent to digit 0; the k = 0 family: counts h(B') + B_t(B') over bases B' with e_t(B') = 0
   (the rooted append at t is then a strong-digit-0 tournament) -- which odd numbers below 125 it covers from bases <= 7."""
import subprocess
exec(open('nearjoin_check.py').read().split('# (1)')[0])   # helpers: parse, classes, ham_all, starts, B_of, strong, sub
from collections import Counter
outdeg = Counter(); total = 0; firstorder_only = 0; d0 = 0; mism = 0
for n in (3, 4, 5, 6, 7, 8):
    for S in classes(n):
        if not strong(S): continue
        dp, H = ham_all(S); N = 1 << n
        for u in range(n):
            B, h, e = B_of(S, u); isd0 = (B == h)
            outs = [x for x in range(n) if S[u][x]]
            fo = all(starts(sub(S, [v for v in range(n) if v != x]), [v for v in range(n) if v != x].index(u))[0][(1 << (n-1)) - 1] == 0 for x in outs)
            total += 1
            if isd0: d0 += 1; outdeg[len(outs)] += 1
            if fo and not isd0: firstorder_only += 1
            if isd0 and not fo: mism += 1
print(f"strong rooted (S,u), orders 3-8: {total}; digit-0 vertices: {d0}; out-degree distribution of digit-0 vertices: {dict(sorted(outdeg.items()))}")
print(f"first-order condition (no Hamiltonian path of S - x from u, all out-neighbours x) without digit 0: {firstorder_only}; digit 0 without the first-order condition: {mism}")
# k = 0 family
vals = set(); valsk = set()
for n in (2, 3, 4, 5, 6, 7):
    for Bp in classes(n):
        dp, Hs = ham_all(Bp); h = Hs[-1]
        for t in range(n):
            if dp[(1 << n) - 1][t] != 0: continue        # e_t(B') = 0
            B, _, _ = B_of(Bp, t)
            vals.add(h + B)
            for k in range(1, 12): valsk.add(h + (2**k) * B)
odd = [v for v in range(1, 125, 2)]
print("k = 0 family (rooted append at a port with no path ending there), bases <= 7: odd numbers < 125 missed:", [v for v in odd if v not in vals])
print("near-join (k >= 1) from the same bases: odd numbers < 125 missed:", [v for v in odd if v not in valsk])
print("union: missed:", [v for v in odd if v not in vals and v not in valsk])
