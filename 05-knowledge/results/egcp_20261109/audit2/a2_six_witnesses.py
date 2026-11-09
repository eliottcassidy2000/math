"""AUDIT2: explicit even, universal-vertex-free witnesses for the six values that e0search_n10/n11.out mark 'realised'
without printing a graph (found[16] overflow): 773, 779 (n = 10), 1533, 1537, 1541, 1547 (n = 11).  Each is re-verified
by a2_check.exe and by a direct subset-by-subset domination count; the evenness is re-checked with networkx VF2."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
want = {10: (773, 779), 11: (1533, 1537, 1541, 1547)}
ok = True
for n, ks in want.items():
    W = {}
    for line in open(os.path.join(HERE, 'a2_e0full_n%d_witnesses.txt' % n)):
        f = list(map(int, line.split())); W[f[0]] = (n, f[1:])
    for k in ks:
        g = W[k]
        r = run_check([(g, -1)])[0]
        nxodd, cnt = nx_sign_odd(g)
        db = d_brute(g)
        good = r[0] == k == db and r[1] == 'EVEN' and r[3] == 0 and not nxodd
        ok &= good
        print("n=%d d=%d: C %s, brute d=%d, universal=%d, networkx: %d automorphisms, odd=%s -> %s\n   edges: %s"
              % (n, k, r[1], db, r[3], cnt, nxodd, 'OK' if good else 'BAD', ' '.join('%d-%d' % e for e in edges(g))))
print("ALL:", "PASS" if ok else "FAIL")
