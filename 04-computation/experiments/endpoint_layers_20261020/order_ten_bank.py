#!/usr/bin/env python3
"""Order-ten endpoint-three layer: run the packet's own template bank (probe.bank + probe_refined.new_bank) over the
complete marked universe C(7,2)*2^12 = 86,016 closed standard-blue presentations with endpoint scores (6,3) at n = 10
(05-knowledge/results/endpoint_three_repairs_20261020h), with NO mixedness filter and the nonstrong-partner filter
recorded only. Residuals (no template gives three legal black K moves) are saved for the search stage."""
import sys, json, time
from collections import Counter, defaultdict
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK))
import probe as P
import probe_refined as R
n = 10
cands = list(P.bank(n)) + list(R.new_bank(n))
print('templates', len(cands), flush=True)
stats = Counter(); marks = defaultdict(Counter); residuals = []; selected = Counter()
t0 = time.time()
for k, (t, u, v) in enumerate(P.endpoint_three(n)):
    stats['universe'] += 1
    other = P.complement(t, tuple(range(n)))
    strong_partner = P.strong(other)
    if not strong_partner: stats['nonstrong_partner'] += 1
    kind = 'open' if P.arc(t, 1, n - 2) else 'closed'
    key = kind + ' I=' + ','.join(map(str, u))
    marks[key]['total'] += 1
    hit = None
    for name, words in cands:
        if P.evaluate(t, words)['success']:
            hit = name; break
    if hit:
        stats['covered'] += 1; marks[key]['covered'] += 1; selected[hit] += 1
    else:
        stats['residual'] += 1; marks[key]['residual'] += 1
        residuals.append({'rows': list(t), 'incoming': list(u), 'outgoing': list(v), 'kind': kind, 'strong_partner': strong_partner})
    if (k + 1) % 5000 == 0: print(k + 1, dict(stats), f'{time.time()-t0:.0f}s', flush=True)
print(json.dumps({'n': n, **dict(stats)}), flush=True)
print('residual patterns', {k: v['residual'] for k, v in sorted(marks.items()) if v.get('residual')}, flush=True)
print('selected templates', dict(selected.most_common()), flush=True)
Path('order_ten_bank.json').write_text(json.dumps({'stats': dict(stats), 'patterns': {k: dict(v) for k, v in sorted(marks.items())},
    'selected': dict(selected), 'residuals': residuals}, indent=1))
print('saved', len(residuals), 'residuals', flush=True)
