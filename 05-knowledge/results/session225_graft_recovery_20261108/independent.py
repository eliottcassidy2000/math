"""Fresh independent permutation/bin-assignment check; no producer imports."""
from itertools import permutations, product
from fractions import Fraction

B = [8, 1, 3, 6]

def paths(out, vertices):
    return [p for p in permutations(vertices)
            if all(out[u] >> v & 1 for u, v in zip(p, p[1:]))]

def covers(out, k):
    n = len(out)
    for bins in product(range(k), repeat=n):
        lists = [paths(out, [i for i in range(n) if bins[i] == j])
                 for j in range(k)]
        yield from product(*lists)

def profile(out):
    ps = paths(out, range(len(out)))
    return [len(ps), sum(p[-1] == 0 for p in ps),
            sum(bool(c[1]) and c[1][0] == 0 for c in covers(out, 2)),
            len(paths(out, range(1, len(out))))]

ps = paths(B, range(4))
c2 = list(covers(B, 2))
c3 = list(covers(B, 3))
N = [0] * 3
for i in range(1, 4):
    for j in range(1, 4):
        if B[i] >> j & 1:
            N[i - 1] |= 1 << (j - 1)

stats = [len(ps), len(c2), sum(p[0] == 0 for p in ps),
         sum(p[-1] == 0 for p in ps),
         sum(bool(c[0]) and c[0][0] == 0 for c in c2),
         sum(bool(c[0]) and c[0][-1] == 0 for c in c2),
         sum(bool(c[0]) and c[0][0] == 0 for c in c3),
         len(paths(N, range(3))), len(list(covers(N, 2)))]
H, Q, s, e, a, de, T, q, b = stats
K = [[H, -H, Q - H, 2 * H - Q],
     [e, -e, de - e, e - de],
     [a, s - a, 2 * T - a, 3 * a - s - 2 * T],
     [q, -q, b - q, 2 * q - b]]

def det(M):
    n = len(M)
    return sum(
        (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        * __import__('functools').reduce(
            lambda x, i: x * M[i][p[i]], range(n), 1)
        for p in permutations(range(n)))

def mul(A, B):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*B)]
            for row in A]

I = [[int(i == j) for j in range(4)] for i in range(4)]
powers = [I]
for i in range(4):
    powers.append(mul(powers[-1], K))
res = [[sum(c * powers[k][i][j]
            for k, c in enumerate([-8, 31, -37, -27, 1]))
        for j in range(4)] for i in range(4)]
basis = [profile(A) for A in ([0], [2, 0], [2, 4, 1], [6, 4, 0])]
v = [[1], [0], [2], [1]]
seq = [1]
for i in range(4):
    v = mul(K, v)
    seq.append(v[0][0])
print({'stats': stats, 'K': K, 'detK': det(K), 'basis': basis,
       'detBasis': det(basis), 'characteristic_residual': res, 'sequence': seq})
if (stats != [5, 32, 1, 1, 5, 7, 17, 1, 8] or det(K) != -8
        or det(basis) != -4 or any(any(row) for row in res)
        or seq != [1, 37, 1113, 31349, 886465]):
    raise ValueError('mismatch')
