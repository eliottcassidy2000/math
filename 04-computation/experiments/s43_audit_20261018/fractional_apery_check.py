# S43 referee: Apery A and truncated sum_k ((1/3)_k (2/3)_k / k!^2)^2 at n == -1/3 mod p versus eta(3tau)^8 mod p
from math import comb
from fractions import Fraction as Fr
def eta_product(shape, N):
    c = [0]*(N+1); c[0] = 1
    for a, m in shape.items():
        for _ in range(m):
            for n in range(1, N//a+1):
                st = a*n
                for k in range(N, st-1, -1):
                    c[k] -= c[k-st]
    return c
N = 400
c3B = eta_product({3: 8}, N)
a3B = lambda m: c3B[m-1]
ps = [p for p in range(7, N) if all(p % d for d in range(2, int(p**.5)+1))]
A = lambda n: sum((comb(n,k)*comb(n+k,k))**2 for k in range(n+1))
rows = []
for p in ps:
    n = (-pow(3, -1, p)) % p
    b = [Fr(1)]
    for k in range(1, p):
        b.append(b[-1]*(Fr(1,3)+k-1)*(Fr(2,3)+k-1)/(k*k))
    full = sum(x*x for x in b)
    fullp = full.numerator*pow(full.denominator, -1, p) % p
    rows.append((p, p % 3, n, A(n) % p, a3B(p) % p, fullp, (-a3B(p)) % p))
okA = all(r[3] == r[4] for r in rows)
okAm = all(r[3] == r[6] for r in rows)
okT = all(r[5] == r[4] for r in rows)
print('A(n_{1/3}(p)) == a_3B(p) mod p for all primes 7..%d: %s; == -a_3B(p): %s; truncated 4F3(1/3,2/3,1/3,2/3) to p-1 == a_3B(p) mod p: %s' % (N, okA, okAm, okT))
print('first rows (p, p mod 3, n, A mod p, a3B mod p, trunc mod p, -a3B mod p):', rows[:10])
okAT = all(r[3] == r[5] for r in rows)
print('A(n_{1/3}(p)) == truncated sum_{k<p} ((1/3)_k (2/3)_k/k!^2)^2 mod p for all primes 7..%d: %s' % (N, okAT))
