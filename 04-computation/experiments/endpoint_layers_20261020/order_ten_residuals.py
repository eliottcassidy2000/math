#!/usr/bin/env python3
"""Search the three order-ten residuals of order_ten_bank.py with the packet's labelled three-K search (probe_search.find_one)
and replay each certificate with probe_complete.independent_replay. Together with order_ten_bank.json this closes the
(6,3) endpoint layer at n = 10 FINITE-EXACT (every closed standard-blue presentation has three black K moves to K_P(T))."""
import sys, json
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK))
import probe as P, probe_search as S, probe_complete as C
data = json.load(open('order_ten_bank.json'))
out = []
for r in data['residuals']:
    t = tuple(r['rows']); u = tuple(r['incoming']); v = tuple(r['outgoing'])
    P.source_checks(t, u, v)
    other = P.complement(t, tuple(range(10)))
    res = S.find_one(t)
    print('residual', r['kind'], 'I =', u, 'J =', v, 'strong partner', P.strong(other), '->', res['scope'], 'pairs', res['pairs_examined'], flush=True)
    if res['words'] is None:
        out.append({**r, 'certificate': None, 'search': {k: res[k] for k in ('scope', 'first_black_paths', 'last_black_paths', 'pairs_examined')}}); continue
    wit = C.independent_replay(t, res['words'])
    print('   words', res['words'], '\n   independent anti-failures', wit, 'excesses', [s['endpoint_excess'] for s in res['evaluation']['steps']], flush=True)
    out.append({**r, 'certificate': [list(q) for q in res['words']], 'independent_anti_failures': wit,
                'endpoint_excesses': [s['endpoint_excess'] for s in res['evaluation']['steps']],
                'search': {k: res[k] for k in ('scope', 'first_black_paths', 'last_black_paths', 'pairs_examined')}})
json.dump(out, open('order_ten_residuals.json', 'w'), indent=1)
print('done; unresolved:', sum(1 for o in out if o['certificate'] is None))
