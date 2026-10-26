# All genuine sums with a degree-3 vertex (36 of 37) at primes q = 3 (mod 4) in [lo, hi) by the O(q^4) evaluator.
# One JSON line per prime: {"q": q, "S": {index: S_H}}.  Usage: run_fast.py lo hi out.jsonl.  mac-mini-2026-10-25-pyth.
import json, sys, time
from sympy import primerange
from sum7 import chi_table
from sum7_fast import S7_fast, L_table
G = json.load(open('genuine37.json'))
def has3(H):
    E = {tuple(sorted(e)) for e in H['edges']}
    return any(sum(1 for y in range(7) if tuple(sorted((x, y))) in E) == 3 for x in range(7))
todo = [i for i, H in enumerate(G) if not H['all_even'] and has3(H)]
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, 'a') as fh:
    for p in primerange(lo, hi):
        if p % 4 != 3: continue
        t0 = time.time(); chi = chi_table(p); L = L_table(p, chi)
        S = {i: S7_fast(p, G[i]['edges'], chi=chi, L=L) for i in todo}
        fh.write(json.dumps({'q': p, 'S': S}) + '\n'); fh.flush()
        print(f"q = {p}: {len(todo)} sums in {time.time() - t0:.0f}s", flush=True)
