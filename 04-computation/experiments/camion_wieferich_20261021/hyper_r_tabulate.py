# Tabulate second-order root counts of r-uniform hypergraph zero classes by class type and compare with random degree-r
# polynomials over F_p (root count = fixed points of a uniform permutation in S_r) and with "forced root + random degree r-1".
# mac-mini-2026-10-21-fanopsp.  Usage: python3 hyper_r_tabulate.py r file...
import sys, re
from math import factorial
from itertools import permutations
from collections import Counter
r = int(sys.argv[1]); files = sys.argv[2:]
def fixdist(k):
    c = Counter(sum(1 for i, x in enumerate(s) if i == x) for s in permutations(range(k))); n = factorial(k)
    return {j: c[j] / n for j in range(k + 1)}
by = Counter(); primes = 0
for f in files:
    for line in open(f):
        if not line.startswith('p='): continue
        primes += 1
        for t in line.split('zero-classes:')[1].split():
            m = re.match(r'(\d+):(\d+)', t); c = int(m.group(1)); k = int(m.group(2))
            by[(c if c in (1, 2, 3) else 'sporadic', k)] += 1
print(f"r = {r}: {primes} primes")
rand = fixdist(r); forced = {j + 1: v for j, v in fixdist(r - 1).items()}
for typ in (1, 2, 3, 'sporadic'):
    tot = sum(v for (t, k), v in by.items() if t == typ)
    if not tot: continue
    obs = {k: by[(typ, k)] for k in range(r + 1)}
    print(f"  class {typ}: n = {tot}  observed {obs}")
    print(f"      random degree-{r}: {{" + ", ".join(f"{k}: {tot * rand.get(k, 0):.1f}" for k in range(r + 1)) + "}")
    print(f"      forced root + random degree-{r-1}: {{" + ", ".join(f"{k}: {tot * forced.get(k, 0):.1f}" for k in range(r + 1)) + "}")
