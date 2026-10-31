"""S69 audit: apply session124's paired digit operators to MANY weak B pairs at m = 42..52 produced by the audit's own
enumerator (ss_enum.exe pairdump), plus the eleven threshold42 SMALL_BASES pairs, and check every output twice:
  (1) with an independent checker written here (permutation of [1,n] / [1,n+1], prefix 1,3, ends 8/10 by parity of
      the order, square adjacencies, equal one-based positional parity of every common vertex);
  (2) with session124's own Core.check_pair (imported read-only from core.py, as the session's replay does).

usage: python replay_own_pairs.py PATH_TO_ss_enum.exe [MAX_PAIRS_PER_m]
"""
import importlib.util
import json
import subprocess
import sys
import time
from math import isqrt
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
RES = REPO / '05-knowledge' / 'results'
S124 = RES / 'session124_recursive_square_prefix_20261031'
spec = importlib.util.spec_from_file_location('s124core', S124 / 'core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
C = core.Core.__new__(core.Core)
C.checks = 0

exe = sys.argv[1]
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 60
OPS = {(p, d): json.loads((S124 / 'paired' / f'B_p{p}_d{d}.json').read_text()) for p in (0, 1) for d in range(24, 73)}


def sq(x):
    return isqrt(x) ** 2 == x


def own_check(S, L):
    n = len(S)
    if len(L) != n + 1:
        raise ValueError('lengths')
    pos = []
    for P in (S, L):
        if sorted(P) != list(range(1, len(P) + 1)):
            raise ValueError('not a permutation')
        if P[0] != 1 or P[1] != 3 or P[-1] != (8 if len(P) % 2 == 0 else 10):
            raise ValueError('prefix/end')
        for a, b in zip(P, P[1:]):
            if not sq(a + b):
                raise ValueError('edge')
        q = [0] * (len(P) + 1)
        for i, v in enumerate(P):
            q[v] = i & 1
        pos.append(q)
    for v in range(1, n + 1):
        if pos[0][v] != pos[1][v]:
            raise ValueError('parity')


def pairs_at(m, printed=20000):
    """the first `printed` compatible pairs in the enumerator's order, and the total number of pairs"""
    out = subprocess.run([exe, 'pairdump', str(m), str(printed)], capture_output=True, text=True, check=True).stdout
    lines = [l.strip() for l in out.splitlines()]
    prs = []
    for i in range(0, len(lines) - 1):
        if lines[i].startswith('S: ') and lines[i + 1].startswith('L: '):
            prs.append((list(map(int, lines[i][3:].split(','))), list(map(int, lines[i + 1][3:].split(',')))))
    total_pairs = int(lines[-1].split('=')[-1])
    return prs, total_pairs


small = json.loads((RES / 'session124_square_prefix_threshold42_20261031' / 'SMALL_BASES.json').read_text())
total = 0
for m in range(42, 53):
    prs, kt = pairs_at(m)
    k = len(prs)
    if k > cap:          # deterministic even sample, always keeping the first and the last pair
        idx = sorted({round(i * (k - 1) / (cap - 1)) for i in range(cap)})
        sample = [prs[i] for i in idx]
    else:
        sample = prs
    sample.append(tuple(small[str(m)]))
    for S, L in sample:
        own_check(S, L)
        C.check_pair([S, L])
        for d in range(24, 73):
            a = OPS[(m % 2, d)]
            res = [list(core.vertices([S, L], w)) for w in a['words']]
            n = 49 * m + d
            if len(res[0]) != n or len(res[1]) != n + 1:
                raise ValueError(f'output length at m={m}, d={d}')
            own_check(res[0], res[1])
            C.check_pair(res)
            total += 1
    print(f'm={m}: {kt} compatible pairs ({k} printed); replayed {len(sample)} inputs '
          f'({len(sample) - 1} own + threshold42 base) x 49 digits; all outputs pass both checkers', flush=True)
print(f'ALL {total} operator outputs verified by the independent checker and by session124 check_pair '
      f'(session124 checks = {C.checks}); {time.time() - T0:.0f} s')
