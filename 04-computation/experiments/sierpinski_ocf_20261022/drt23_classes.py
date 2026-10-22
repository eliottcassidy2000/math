# The doubly regular 23-tournaments of drt23_codes_all.txt plus P_23: nauty classes (labelg), exact H (hk64.c,
# Held-Karp, 1.5 GB), 5-vertex census (key = sorted scores + H of the 5-subtournament) and TT5 = number of transitive
# 5-sets.  Same class order as drt_census_hierarchy.py (first appearance; P_23 last).  Output: drt23_classes.json.
# mac-mini-2026-10-22-loop4.
import sys, subprocess, itertools, json
from collections import Counter
import numpy as np
sys.path.insert(0, '.')
from drt_switch import code_to_S, d6
codes = [l.strip() for l in open('drt23_codes_all.txt') if l.strip() and not l.startswith('#')]
q = 23; sq = {(t * t) % q for t in range(1, q)}
pcode = ''.join('1' if (j - i) % q in sq else '0' for i in range(q) for j in range(i + 1, q))
codes.append(pcode)
canon = subprocess.run(['labelg', '-q'], input='\n'.join(d6(code_to_S(c)) for c in codes) + '\n', capture_output=True, text=True).stdout.split()
uniq = {}
for c, cf in zip(codes, canon):
    if cf not in uniq: uniq[cf] = c
reps = list(uniq.items())
perms5 = list(itertools.permutations(range(5)))
def census5(S):
    n = S.shape[0]; cnt = Counter()
    for X in itertools.combinations(range(n), 5):
        sub = S[np.ix_(X, X)]
        sc = tuple(sorted(int((sub[i] == 1).sum()) for i in range(5)))
        h = sum(1 for p in perms5 if all(sub[p[i], p[i + 1]] == 1 for i in range(4)))
        cnt[(sc, h)] += 1
    return cnt
# converse pairs: class of the converse tournament
conv_canon = subprocess.run(['labelg', '-q'], input='\n'.join(d6(-code_to_S(c)) for _, c in reps) + '\n', capture_output=True, text=True).stdout.split()
out = []
for k, (cf, c) in enumerate(reps):
    S = code_to_S(c); n = S.shape[0]
    assert (S @ S.T == n * np.eye(n, dtype=int) - np.ones((n, n), dtype=int)).all()
    H = int(subprocess.run(['./hk64'], input=c + '\n', capture_output=True, text=True).stdout.split()[0])
    cen = census5(S)
    conv = [j for j, (cg, _) in enumerate(reps) if cg == conv_canon[k]]
    rec = {'class': k, 'label': 'P_23' if c == pcode else f'D23_{k}', 'code': c, 'canon': cf, 'H': H,
           'TT5': cen[((0, 1, 2, 3, 4), 1)], 'census5': {f'{sc}|{h}': v for (sc, h), v in sorted(cen.items())},
           'converse_class': conv[0] if conv else None}
    out.append(rec)
    print(k, rec['label'], H, 'TT5', rec['TT5'], 'converse class', rec['converse_class'], flush=True)
json.dump(out, open('drt23_classes.json', 'w'), indent=1)
