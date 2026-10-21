#!/usr/bin/env python3
"""The (n-5,4) endpoint layer (|I| = 3) at n = 9, 10: complete marked universe C(n-3,3) * 2^floor((n-3)^2/4) of closed
standard-blue presentations with d(a) = n-5, d(z) = 4; run the inherited endpoint-two/three template bank (probe.bank +
probe_refined.new_bank, words only, all guards re-tested on the actual source), then the labelled three-K search
(probe_search.find_one) on residuals, with independent replay. No mixedness filter; nonstrong partners recorded."""
import sys, json, time
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK))
import probe as P, probe_refined as R, probe_search as S, probe_complete as C
def endpoint_k(n, k):
    for middle in P.OLD.blue_presentations(n - 2):
        for incoming in combinations(range(2, n - 1), k - 1):
            outgoing = tuple(sorted(n - 1 - v for v in incoming))
            rows = [sum(1 << v for v in range(1, n - 1))]
            rows += [(row << 1) | (1 << (n - 1)) for row in middle]
            rows += [1]
            for u in incoming:
                rows[0] ^= 1 << u; rows[u] ^= 1
            for v in outgoing:
                rows[v] ^= 1 << (n - 1); rows[-1] ^= 1 << v
            yield tuple(rows), incoming, outgoing
def checks(t, u, v, k):
    n = len(t); p = tuple(range(n))
    P.need(P.path(t, p) and P.blue(t, p) and P.arc(t, n - 1, 0), 'closed standard blue source')
    P.need(t[0].bit_count() == n - 1 - k and t[-1].bit_count() == k, 'endpoint scores')
    P.need(tuple(w for w in range(1, n - 1) if P.arc(t, w, 0)) == u, 'incoming mark exact')
    P.need(tuple(w for w in range(1, n - 1) if P.arc(t, n - 1, w)) == v, 'outgoing mark exact')
k = 4
def main():
  for n in [int(a) for a in sys.argv[1:]] or [9, 10]:
        cands = list(P.bank(n)) + list(R.new_bank(n))
        stats = Counter(); marks = defaultdict(Counter); residuals = []; selected = Counter(); certs = []
        t0 = time.time()
        for i, (t, u, v) in enumerate(endpoint_k(n, k)):
            checks(t, u, v, k); stats['universe'] += 1
            other = P.complement(t, tuple(range(n)))
            sp = P.strong(other)
            if not sp: stats['nonstrong_partner'] += 1
            kind = 'open' if P.arc(t, 1, n - 2) else 'closed'
            key = kind + ' I=' + ','.join(map(str, u)); marks[key]['total'] += 1
            hit = None
            for name, words in cands:
                if P.evaluate(t, words)['success']: hit = name; break
            if hit: stats['covered'] += 1; marks[key]['covered'] += 1; selected[hit] += 1
            else:
                res = S.find_one(t)
                if res['words'] is None:
                    stats['unresolved'] += 1; marks[key]['unresolved'] += 1
                    residuals.append({'rows': list(t), 'incoming': list(u), 'outgoing': list(v), 'kind': kind, 'strong_partner': sp,
                                      'search': {kk: res[kk] for kk in ('scope', 'first_black_paths', 'last_black_paths', 'pairs_examined')}})
                else:
                    wit = C.independent_replay(t, res['words'])
                    stats['search_certified'] += 1; marks[key]['searched'] += 1
                    certs.append({'rows': list(t), 'incoming': list(u), 'outgoing': list(v), 'kind': kind, 'strong_partner': sp,
                                  'words': [list(q) for q in res['words']], 'independent_anti_failures': wit, 'pairs_examined': res['pairs_examined']})
            if (i + 1) % 5000 == 0: print(n, i + 1, dict(stats), f'{time.time()-t0:.0f}s', flush=True)
        print(json.dumps({'n': n, 'k': k, **dict(stats)}), flush=True)
        print('patterns needing search/unresolved', {kk: (vv.get('searched', 0), vv.get('unresolved', 0)) for kk, vv in sorted(marks.items()) if vv.get('searched') or vv.get('unresolved')}, flush=True)
        print('selected', dict(selected.most_common()), flush=True)
        Path(f'endpoint_four_n{n}.json').write_text(json.dumps({'stats': dict(stats), 'patterns': {kk: dict(vv) for kk, vv in sorted(marks.items())},
            'selected': dict(selected), 'search_certificates': certs, 'unresolved': residuals}, indent=1))

if __name__ == '__main__':
    main()
