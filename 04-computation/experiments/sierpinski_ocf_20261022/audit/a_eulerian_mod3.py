# Audit: A_n(x) mod 3 = (x-1)^(n-1) (n odd), (x-1)^(n-2)(x+1) (n even); A_n mod 2 = (1+x)^(n-1).
from math import comb
def eulerian_rows(N):
    rows = {1: [1]}
    for m in range(2, N + 1):
        prev = rows[m - 1]
        row = [0] * m
        for k in range(m):
            a = (k + 1) * prev[k] if k < m - 1 else 0
            b = (m - k) * prev[k - 1] if k > 0 else 0
            row[k] = a + b
        rows[m] = row
    return rows
def polymul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return r
def xm1pow(e, p):
    return [(comb(e, k) * (-1) ** (e - k)) % p for k in range(e + 1)]
rows = eulerian_rows(80)
bad = []
for n in range(1, 81):
    A = [a % 3 for a in rows[n]]
    if n % 2 == 1:
        want = xm1pow(n - 1, 3)
    else:
        want = polymul(xm1pow(n - 2, 3), [1, 1], 3)
    if A != want: bad.append(('mod3', n))
    if [a % 2 for a in rows[n]] != [comb(n - 1, k) % 2 for k in range(n)]: bad.append(('mod2', n))
print("Eulerian A_n mod 3 / mod 2 pattern, n = 1..80: failures =", bad)
# collapse at n = 2*3^j+1 and 'flat plus endpoints' at n = 3^j+1
for j in range(0, 4):
    n = 2 * 3 ** j + 1
    nz = [k for k in range(n) if rows[n][k] % 3]
    print(f"n={n}: nonzero A(n,k) mod 3 at k={nz}, values {[rows[n][k] % 3 for k in nz]}")
for j in range(1, 4):
    n = 3 ** j + 1
    print(f"n={n}: A(n,k) mod 3 = {[rows[n][k] % 3 for k in range(n)]}")
