# F8 branch parity = parity of the number of middle (A) steps in the Berggren word of a primitive triple.
# Search for a closed form in the Euclid parameters (m, n): descent t = m/n:
#   t > 3 -> D^-1 : (m - 2n, n);   2 < t < 3 -> A^-1 : (n, m - 2n);   1 < t < 2 -> U^-1 : (n, 2n - m).
# Candidates tested on all primitive triples with m <= MMAX.
import sys
from math import gcd
from sympy import jacobi_symbol
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300

def word(m, n):
    w = []
    while (m, n) != (2, 1):
        if m > 3 * n: w.append('D'); m -= 2 * n
        elif m > 2 * n: w.append('A'); m, n = n, m - 2 * n
        else: w.append('U'); m, n = n, 2 * n - m
    return ''.join(reversed(w))

def jac(a, b):           # Jacobi symbol (a/b), b odd positive
    return jacobi_symbol(a % b, b)

rows = []
for m in range(2, MMAX + 1):
    for n in range(1, m):
        if (m - n) % 2 == 0 or gcd(m, n) != 1: continue
        w = word(m, n)
        par = w.count('A') % 2
        a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
        rows.append((m, n, a, b, c, w, par))
print(len(rows), "primitive triples")

def test(name, f):
    vals = [(f(*r[:6]), r[6]) for r in rows]
    # does the candidate determine the parity?  (as a function of its value)
    table = {}
    ok = True
    for v, p in vals:
        if v in table and table[v] != p: ok = False; break
        table[v] = p
    if ok:
        print(f"DETERMINES parity: {name}; value->parity table (first 12): {dict(list(sorted(table.items(), key=str))[:12])}")
    return ok

odd = lambda m, n: m if m % 2 else n
even = lambda m, n: n if m % 2 else m
cands = {
 'jacobi(even/odd)': lambda m, n, a, b, c, w: jac(even(m, n), odd(m, n)),
 'jacobi(n/odd) etc': lambda m, n, a, b, c, w: (jac(n, m) if m % 2 else jac(m, n)),
 'jacobi(m-n / m+n)': lambda m, n, a, b, c, w: jac(m - n, m + n),
 'jacobi(m+n / m-n)': lambda m, n, a, b, c, w: jac(m + n, m - n) if m - n > 1 else 1,
 'jacobi(2/c)': lambda m, n, a, b, c, w: jac(2, c),
 'jacobi(2/a)': lambda m, n, a, b, c, w: jac(2, a),
 'sign(a-b)': lambda m, n, a, b, c, w: (a > b),
 '(a-b) mod 16': lambda m, n, a, b, c, w: (a - b) % 16,
 '(a+b) mod 16': lambda m, n, a, b, c, w: (a + b) % 16,
 'c mod 16': lambda m, n, a, b, c, w: c % 16,
 'a mod 16': lambda m, n, a, b, c, w: a % 16,
 'b mod 32': lambda m, n, a, b, c, w: b % 32,
 '(m mod 8, n mod 8)': lambda m, n, a, b, c, w: (m % 8, n % 8),
 '(m mod 16, n mod 16)': lambda m, n, a, b, c, w: (m % 16, n % 16),
 'len(word) mod 2': lambda m, n, a, b, c, w: len(w) % 2,
 '#U mod 2': lambda m, n, a, b, c, w: w.count('U') % 2,
 '#D mod 2': lambda m, n, a, b, c, w: w.count('D') % 2,
 '(a-b sign, (a-b) mod 8)': lambda m, n, a, b, c, w: (a > b, (a - b) % 8),
 'jacobi(m/ (m^2+n^2)) ': lambda m, n, a, b, c, w: jac(m, c),
 'jacobi(n/c)': lambda m, n, a, b, c, w: jac(n, c),
 'jacobi(m-n / c)': lambda m, n, a, b, c, w: jac(m - n, c),
 'jacobi(m+n / c)': lambda m, n, a, b, c, w: jac(m + n, c),
}
for k, f in cands.items(): test(k, f)
# parity along rays
for ray in 'UAD':
    print(ray, [r[6] for r in rows if set(r[5]) <= {ray}][:10])
