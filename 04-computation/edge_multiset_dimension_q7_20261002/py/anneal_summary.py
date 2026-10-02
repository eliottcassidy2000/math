"""anneal_summary.py -- summarise annealing logs; re-verify every reported set with the pure-Python checker.

usage: python3 anneal_summary.py logfile [logfile ...]
For each log: restarts, moves, best number of colliding edge pairs, and the pure-Python defect
(#edges - #distinct histograms) of the best set; every RESOLVING set is re-verified with pycheck.py and
canonicalised with src/canon (distinct orbits are counted).
"""
import sys, os, re, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import pycheck
from collections import Counter

def pairs(d, S):
    c = Counter(pycheck.histograms(d, S))
    return sum(m * (m - 1) // 2 for m in c.values())

for path in sys.argv[1:]:
    txt = open(path).read()
    summ = [l for l in txt.splitlines() if l.startswith('SUMMARY')]
    res = [l for l in txt.splitlines() if l.startswith('RESOLVING')]
    if not summ:
        print('%s: no SUMMARY line (still running?)' % path); continue
    f = dict(t.split('=', 1) for t in summ[-1].split()[1:])
    d, k = int(f['d']), int(f['k'])
    best = [int(x) for x in f['best_set'].split(',')]
    assert len(best) == k and len(set(best)) == k
    dfct, prs = pycheck.defect(d, best), pairs(d, best)
    ok = prs == int(f['best_colliding_pairs'])
    sets = [sorted(int(x) for x in l.split('set=')[1].split(',')) for l in res]
    pyok = all(len(S) == k and pycheck.is_resolving(d, S) for S in sets)
    orbits = set()
    if sets:
        inp = '\n'.join(' '.join(map(str, S)) for S in sets) + '\n'
        out = subprocess.run([os.path.join(ROOT, 'src', 'canon'), str(d)], input=inp, capture_output=True, text=True, check=True).stdout
        orbits = set(l.split()[0] for l in out.splitlines() if l.strip())
    print('ANNEAL %s d=%d k=%d restarts=%s moves=%s best_colliding_pairs=%s (python: pairs=%d defect=%d, agree=%s) '
          'resolving_found=%d python_verified=%s orbits=%d best_set=%s' % (os.path.basename(path), d, k, f['restarts'], f['moves'],
          f['best_colliding_pairs'], prs, dfct, ok, len(sets), pyok, len(orbits), ','.join(map(str, sorted(best)))))
