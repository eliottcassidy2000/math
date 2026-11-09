"""AUDIT2 for HYP-9610 (A) evidence: re-verify every graph printed by e0search (n = 10, 11) with independent code,
and compare the set of values that actually have a printed witness with the target lists and the 'realised' lines.
"""
import os, sys, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
EXP = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'egcp_20261109')
misses = {}
for line in open(os.path.join(EXP, 'e0_closure_misses.txt')):
    f = line.split()
    if f: misses[int(f[0])] = [int(x) for x in f[1:]]
ok = True
for n in (10, 11):
    found = {}; claimed = set()
    for line in open(os.path.join(EXP, 'e0search_n%d.out' % n)):
        m = re.match(r'FOUND d=(\d+) edges:(.*)', line)
        if m:
            d = int(m.group(1))
            E = [tuple(map(int, t.split('-'))) for t in m.group(2).split()]
            found[d] = from_edges(n, E)
        m = re.match(r'target (\d+) realised', line)
        if m: claimed.add(int(m.group(1)))
    items = [(found[d], -1) for d in sorted(found)]
    res = run_check(items)
    nbad = 0
    for d, r, (g, _) in zip(sorted(found), res, items):
        good = r[0] == d and r[1] == 'EVEN' and r[3] == 0 and g[0] == n
        nbad += not good
        if not good: print("  BAD witness n=%d d=%d -> %s" % (n, d, r))
    print("n=%d: targets %d, 'realised' lines %d, FOUND lines with a graph %d; independent re-check: %d bad"
          % (n, len(misses[n]), len(claimed), len(found), nbad))
    print("      targets claimed 'realised' but with NO printed graph: %s" % sorted(claimed - set(found)))
    print("      index of those targets in the command-line target list: %s"
          % [misses[n].index(x) for x in sorted(claimed - set(found))])
    ok &= nbad == 0
print("witness re-check:", "PASS" if ok else "FAIL")
