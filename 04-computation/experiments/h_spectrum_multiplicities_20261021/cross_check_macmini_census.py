"""Independent cross-check of the S50 class census against the concurrent mac-mini-2026-10-21-fanopsp census
(session opus-2026-10-21-S50).  Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md.

The S50 census (h_multiplicities.py) builds the isomorphism classes of n-tournaments, n <= 9, from the S36 generator
with a refined-invariant canonical form.  The mac-mini census (04-computation/experiments/tournament_pseudoprimes_20261021/
census_n*.txt, from gentourng and tour_census.c) lists, for every class, H, the Hamiltonian-cycle count C, strongness,
|Aut|, the number of anti-automorphisms and c3.  Different generator, different canonisation, different path counter.
Check, for every n = 3..9: the multisets of (H, strong, self-converse, |Aut|) over classes coincide exactly.
Prints ALL CHECKS PASSED.  Runtime about 30 seconds."""
import sys, os, gzip, time
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h_multiplicities as HM                                            # noqa: E402

FAIL = []
T0 = time.time()
MM = os.path.join(HERE, '..', 'tournament_pseudoprimes_20261021')


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def read_mm(n):
    path = os.path.join(MM, f'census_n{n}.txt')
    if os.path.exists(path):
        lines = open(path).read().split('\n')
    else:
        lines = gzip.open(path + '.gz', 'rt').read().split('\n')
    out = Counter()
    for ln in lines:
        if not ln.strip():
            continue
        H, C, s, A, AA, c3 = ln.split()
        out[(int(H), int(s) == 1, int(AA) > 0, int(A))] += 1
    return out


data = HM.classes_upto(9)
for n in range(3, 10):
    A = data[n]['A']
    H, st = HM.ham_all(n, A)
    mine = Counter((int(h), bool(s), bool(c), int(a)) for h, s, c, a in zip(H, st, data[n]['sc'], data[n]['aut']))
    theirs = read_mm(n)
    check(mine == theirs, f'n = {n}: {sum(mine.values())} classes; the (H, strong, self-converse, |Aut|) multisets agree with the '
          f'mac-mini census ({sum(theirs.values())} classes)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
