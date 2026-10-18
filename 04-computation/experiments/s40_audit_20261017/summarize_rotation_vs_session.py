"""Compare the audit's best rotation orders (allp500/audit_rotation_witnesses.json) and audit search (audit_search_*.json)
with the session's local-search bests (05-knowledge/results/paley_ranking_local_search*_20261017.json), p <= 500."""
import json, glob, os, math
H = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(H, '..', '..', '..', '05-knowledge', 'results')
rot = {r['p']: r for r in json.load(open(os.path.join(H, 'allp500', 'audit_rotation_witnesses.json')))}
ses = {}
for fn in ('paley_ranking_local_search_20261017.json', 'paley_ranking_local_search_p132_500_20261017.json'):
    for r in json.load(open(os.path.join(R, fn))):
        ses[r['p']] = r
srch = {}
for fn in glob.glob(os.path.join(H, 'audit_search_*.json')):
    for r in json.load(open(fn)):
        srch[r['p']] = max(srch.get(r['p'], -10**9), r['best'])
print(' p    h    L     natural  session-best  rotation-best  audit-search  best-known  natural-optimal?')
unbeaten = []
for p in sorted(rot):
    r = rot[p]; s = ses.get(p, {}).get('best'); a = srch.get(p)
    best = max(x for x in (s, r['best_rotation_excess'], a) if x is not None)
    L = 2 * math.pi * r['h'] / ((6 if p == 3 else 2) * math.sqrt(p))
    flag = 'beaten' if best > r['natural'] else ('exact optimum' if p in (3, 7, 23) else 'unbeaten')
    if best <= r['natural']:
        unbeaten.append((p, round(L, 3)))
    print('%4d %3d %6.3f %8d %12s %13d %13s %10d   %s' % (p, r['h'], L, r['natural'], s, r['best_rotation_excess'], a, best, flag))
print('not beaten by any known order:', unbeaten)
