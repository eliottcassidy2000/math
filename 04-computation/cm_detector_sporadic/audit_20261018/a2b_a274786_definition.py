"""Audit: OEIS A274786's own formula a(n) = sum_{j=0}^{2n} (-1)^j C(2n,j) C(j,n)^3 equals C(2n,n) b_n (Bala's product
formula) - checked for n < NMAX - and the first 9 terms match the OEIS %S line."""
from math import comb
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
oeis = [1, 6, 114, 2940, 87570, 2835756, 96982116, 3446781624, 126047377170]
bad = []
for n in range(NMAX):
    a = sum((-1) ** j * comb(2 * n, j) * comb(j, n) ** 3 for j in range(2 * n + 1))
    b = comb(2 * n, n) * sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1))
    if a != b:
        bad.append(n)
    if n < len(oeis) and a != oeis[n]:
        bad.append(('oeis', n))
print(f'A274786 definition == C(2n,n) b_n for 0 <= n < {NMAX}: {not bad}  {bad[:5]}')
