"""Check THM-5040 Corollary 4, first bullet: "The top bands T(m) are realisable only on m vertices, and only by
graphs that are neither unions nor joins (Lemma 0)."

(1) explicit counterexample: P_3 u G_7 with d(G_7) = 101 has 10 vertices, d = 5 * 101 = 505 in T(10) = (503, 512),
    and is even (verified with the audit's own code);
(2) counts, for each m, the top-band values (T(m) and X(m-1)) that are products d1 * d2 with d1 in O(a), d2 in O(b),
    a + b = m, d1, d2 > 1 (nontrivial disjoint unions; uses E = O, i.e. the theorem itself, for the census of
    factors), and checks that NO top-band value is a join value (joins of two graphs on a, b >= 1 vertices).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_lib import *

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'egcp_20261109', 'base_witnesses.txt')
Wt = {}
for line in open(BASE):
    if line.startswith('#') or not line.strip():
        continue
    p = line.split()
    if p[0] == 'W':
        Wt[int(p[2])] = from_hex(int(p[1]), p[3] if len(p) > 3 else '')

G7 = pad(Wt[101], 7)
U = disjoint_union(Pk(3), G7)
print("P_3 u G_7: n = %d, d = %d, even = %s, connected components >= 2: union; 505 in T(10) = (%d, %d): %s" % (
    U[0], dom_count(U), is_even(U), M(9), 2 ** 9, M(9) < 505 < 2 ** 9))
# a second, larger example
G = disjoint_union(Pk(3), pad(Wt[203] if 203 in Wt else Wt[201], 8))
print("P_3 u G_8 (d(G_8) = %d): n = %d, d = %d, even = %s" % (dom_count(pad(Wt[203], 8)), G[0], dom_count(G), is_even(G)))
print("  T(11) = (%d, %d); X(10) n/a" % (M(10), 2 ** 10))

for m in range(8, 25):
    top = {2 ** (m - 1) - j for j in Tj(m)}
    prods = set()
    for a in range(3, m - 2):
        b = m - a
        A = [x for x in O(a) if x > 1]
        lo = min(top)
        for d1 in A:
            # d2 range with d1 * d2 in [lo, 2^(m-1))
            for d2 in range(max(3, (lo + d1 - 1) // d1) | 1, 2 ** (m - 1) // d1 + 2, 2):
                if d1 * d2 in top and d2 in O(b):
                    prods.add(d1 * d2)
    print("m=%2d: |top-band values| = %d, realisable as nontrivial disjoint unions: %d %s" % (
        m, len(top), len(prods), sorted(2 ** (m - 1) - v for v in prods)[:8]))
# joins: d(A v B) = 2^m - D(cA) - D(cB) + 1 with D(cA) <= 2^a - 1 ... show join values >= 2^(m-1) + 1 for a, b >= 1
print("join values: d(A v B) = 2^m - (2^a - d(A)) - (2^b - d(B)) + 1 >= 2^m - 2^a - 2^b + 3 > 2^(m-1) for a, b >= 1:",
      all(2 ** m - 2 ** a - 2 ** (m - a) + 3 > 2 ** (m - 1) for m in range(2, 60) for a in range(1, m)))
