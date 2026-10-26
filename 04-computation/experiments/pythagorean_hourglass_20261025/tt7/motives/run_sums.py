# Compute all 38 seven-vertex survivor sums S_H(q) (37 genuine + K7 - C7) at primes q = 3 (mod 4) in [lo, hi);
# one JSON line per prime: {"q": q, "S": [S_H for H in genuine37.json order]}.  Check: sum over the 37 genuine
# m_H S_H must equal g(q) from tt7_genuine_g.out.  Usage: run_sums.py lo hi out.jsonl.  mac-mini-2026-10-25-pyth.
import json, sys, time, re
from sympy import primerange
from sum7 import S7
G = json.load(open('genuine37.json'))
g = {}
for l in open('../tt7_genuine_g.out'):
    m = re.match(r'q=(\d+) g=(-?\d+)', l)
    if m: g[int(m.group(1))] = int(m.group(2))
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, 'a') as fh:
    for p in primerange(lo, hi):
        if p % 4 != 3: continue
        t0 = time.time()
        S = [S7(p, H['edges']) for H in G]
        comb = sum(H['m'] * s for H, s in zip(G, S) if not H['all_even'])
        chk = (comb == g[p]) if p in g else None
        fh.write(json.dumps({'q': p, 'S': S, 'check_g': chk}) + '\n'); fh.flush()
        print(f"q = {p}: done in {time.time() - t0:.0f}s, sum m_H S_H = g(q): {chk}", flush=True)
