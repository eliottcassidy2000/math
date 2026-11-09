#!/usr/bin/env python3
"""spec3_small.py NMAX -- base-b packing counts H_b = I(Omega(T), b) of every tournament of order <= NMAX (gentourng),
by the set-partition brute force of ribbon_packing.py; checks H_2 = H, prints the base-3 value set and its misses
among values 1 mod 3, and the base-3 values of STRONG tournaments; writes spec3_orders_le_{NMAX}.txt."""
import sys, subprocess
from collections import defaultdict
from ribbon_packing import packing_poly_bruteforce, evalpoly
sys.path.insert(0, '../ribbon_spectrum_20261104')
from ribbon import ham_count_bruteforce, is_strong
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
def gentourng(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for bits in out:
        adj = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[k] == '1': adj[i] |= 1 << j
                else: adj[j] |= 1 << i
                k += 1
        res.append(adj)
    return res
vals = defaultdict(set); strong_vals = defaultdict(set); born = {}
polys_by_value3 = defaultdict(set)
for n in range(1, NMAX + 1):
    for adj in gentourng(n):
        p = packing_poly_bruteforce(adj)
        assert evalpoly(p, 2) == ham_count_bruteforce(adj)
        st = is_strong(adj)
        for b in (2, 3, 4, 5):
            v = evalpoly(p, b); vals[b].add(v)
            if st: strong_vals[b].add(v)
        v3 = evalpoly(p, 3)
        polys_by_value3[v3].add(tuple(p))
        born.setdefault(v3, n)
    print(f"order {n} done: base-3 values so far {len(vals[3])}", flush=True)
with open(f"spec3_orders_le_{NMAX}.txt", "w") as f:
    for b in (2, 3, 4, 5):
        V = sorted(vals[b]); S = sorted(strong_vals[b])
        top = max(V)
        misses = [v for v in range(1, top + 1) if v % b == 1 and v not in vals[b]]
        f.write(f"# base {b}: {len(V)} values at orders <= {NMAX}, max {top}; misses (1 mod {b}) below the max: {misses[:60]}{' ...' if len(misses) > 60 else ''}\n")
        f.write(f"# base {b}: strong values {len(S)}; values not attained by a strong tournament of order <= {NMAX}: {[v for v in V if v not in strong_vals[b]][:40]}\n")
    f.write("# base-3 value: birth order, alpha-vectors realising it (orders <= %d)\n" % NMAX)
    for v in sorted(polys_by_value3):
        if v <= 400: f.write(f"{v} {born[v]} {sorted(polys_by_value3[v])}\n")
print(open(f"spec3_orders_le_{NMAX}.txt").read()[:3000])
