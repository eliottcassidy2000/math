# Is the A-step parity an XOR (product of +-1 characters) of simple arithmetic features of (m, n)?
# Exhaustive search over all XOR-combinations of a feature list, as bitsets over all primitive triples, m <= MMAX.
import sys, itertools
from math import gcd
from sympy import jacobi_symbol
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 150
def word(m, n):
    w = []
    while (m, n) != (2, 1):
        if m > 3 * n: w.append('D'); m -= 2 * n
        elif m > 2 * n: w.append('A'); m, n = n, m - 2 * n
        else: w.append('U'); m, n = n, 2 * n - m
    return ''.join(reversed(w))
def J(a, b): return jacobi_symbol(a % b, b) if b > 1 else 1
rows = []
for m in range(2, MMAX + 1):
    for n in range(1, m):
        if (m - n) % 2 == 0 or gcd(m, n) != 1: continue
        rows.append((m, n, word(m, n)))
o = lambda m, n: m if m % 2 else n
e = lambda m, n: n if m % 2 else m
feat = {
 'J(e|o)': lambda m, n: J(e(m, n), o(m, n)) < 0,
 'J(m-n|m+n)': lambda m, n: J(m - n, m + n) < 0,
 'J(2|o)': lambda m, n: J(2, o(m, n)) < 0,
 'J(-1|o)': lambda m, n: J(-1, o(m, n)) < 0,
 'J(2|m+n)': lambda m, n: J(2, m + n) < 0,
 'J(-1|m+n)': lambda m, n: J(-1, m + n) < 0,
 'J(2|m-n)': lambda m, n: J(2, m - n) < 0,
 'J(-1|m-n)': lambda m, n: J(-1, m - n) < 0,
 'm odd': lambda m, n: m % 2 == 1,
 'e/2 odd': lambda m, n: (e(m, n) // 2) % 2 == 1,
 'e/4 odd': lambda m, n: (e(m, n) // 4) % 2 == 1,
 'a>b': lambda m, n: m * m - n * n > 2 * m * n,
 'm>3n': lambda m, n: m > 3 * n,
 'm>2n': lambda m, n: m > 2 * n,
 'J(m|n)odd?': lambda m, n: (J(m, n) < 0) if n % 2 else False,
 'J(n|m)odd?': lambda m, n: (J(n, m) < 0) if m % 2 else False,
 'inv n mod m < m/2': lambda m, n: (pow(n, -1, m) < m / 2) if m > 1 else False,
 'inv m mod n < n/2': lambda m, n: (pow(m, -1, n) < n / 2) if n > 1 else False,
}
names = list(feat)
masks = {k: sum(1 << i for i, (m, n, w) in enumerate(rows) if feat[k](m, n)) for k in names}
target = sum(1 << i for i, (m, n, w) in enumerate(rows) if w.count('A') % 2)
full = (1 << len(rows)) - 1
hits = []
for r in range(0, 6):
    for combo in itertools.combinations(names, r):
        x = 0
        for k in combo: x ^= masks[k]
        if x == target or x ^ full == target: hits.append((combo, x != target))
print(len(rows), "triples;", "hits:", hits[:20])
# how far is the best single/double feature?
best = []
for r in (1, 2):
    for combo in itertools.combinations(names, r):
        x = 0
        for k in combo: x ^= masks[k]
        agree = bin(~(x ^ target) & full).count('1')
        best.append((max(agree, len(rows) - agree) / len(rows), combo))
best.sort(reverse=True); print(best[:6])
