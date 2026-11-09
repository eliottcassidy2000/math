"""Check the census base used by THM-5040 (all from base_sets.json / base_witnesses.txt, extracted by base_extract.py
from the exhaustive dom_even census of THM-4996):
  (1) E(n) = O(n) for n <= 9;
  (2) every W witness: even, d as stated, order = least order of an even realiser;
  (3) every L witness (m = 8, 9): even, d as stated, exactly m vertices, has a leaf, no universal vertex;
      and the L witnesses cover T(8), T(9) and 2^(m-1) + {1,3,5,7,11,13}.
"""
import json
from egcp_lib import O, Tj, from_hex, dcount, is_even, has_leaf, no_universal

sets = json.load(open('base_sets.json'))
E = {int(n): set(v) for n, v in sets['E'].items()}
ok = True
for n in range(1, 10):
    r = E[n] == O(n)
    ok &= r
    print("E(%d) == O(%d): %s" % (n, n, r))
least = {}
for n in range(1, 10):
    for v in E[n]:
        least.setdefault(v, n)
W = 0; L = {8: set(), 9: set()}
for line in open('base_witnesses.txt'):
    if line.startswith('#'): continue
    tag, n, d, *e = line.split(); n = int(n); d = int(d); g = from_hex(n, e[0] if e else '')
    good = dcount(g) == d and is_even(g)
    if tag == 'W':
        good &= n == least[d]; W += 1
    else:
        good &= has_leaf(g) and no_universal(g); L[n].add(d)
    if not good:
        ok = False; print("BAD witness", line.strip())
for m in (8, 9):
    need = {2**(m - 1) - j for j in Tj(m)} | {2**(m - 1) + c for c in (1, 3, 5, 7, 11, 13)}
    r = need <= L[m]; ok &= r
    print("L(%d) witnesses cover T(%d) and the offsets: %s" % (m, m, r))
print("W witnesses checked:", W, "| all base facts:", "PASS" if ok else "FAIL")
