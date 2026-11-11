"""Exact new-energy controls, independent of the inherited compiler engine."""
import json
from pathlib import Path
import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def f(n):
    total = 0
    while n:
        q, b = divmod(n, 2)
        total += b * (q + (q & 1))
        n = q
    return total


def g(n):
    total = 1
    while n:
        q, b = divmod(n, 2)
        total += (1 - b) * (q + (q & 1))
        n = q
    return total


def coeff(m):
    return [1-m, 0, m-3, m-5, -4, -m-1][m % 6]


def weight(d):
    return 2 if d == 1 else 2 ** (d - 1)


matrices = []
for m in range(1, 38):
    W = sp.Matrix(m, m, lambda i, j: 0 if i == j else weight(abs(i-j)))
    D = sp.eye(m)
    C = sp.zeros(m)
    for i in range(1, m):
        D[i, i-1] = -2
        C[i, i] = 4
        C[i, i-1] = C[i-1, i] = -1 if i == 1 else -3
    for i in range(m-2):
        C[i, i+2] = C[i+2, i] = 1
    require(D * W * D.T == -2 * C, f'banded form at {m}')
    det = int(W.det(method='domain-ge'))
    require(det == (-2)**m * coeff(m), f'determinant at {m}')
    if m % 6 == 1:
        v = sp.Matrix([[1, 0, -2, -3, -2, 0][i % 6] for i in range(m)])
        require(W * v == sp.zeros(m, 1), f'kernel at {m}')
        require(sum(2**i * int(v[i]) for i in range(m)) == 1, 'value shift')
        require(len(W.nullspace()) == 1, 'nullity')
        if m > 1:
            require(max(abs(int(z)) for z in v) == 3, 'binary hostile')
    matrices.append({'digits': m, 'determinant': det})

words = 0
for m in range(1, 13):
    for n in range(2**m):
        bits = [(n >> i) & 1 for i in range(m)]
        rev = sum(bits[i] << (m-1-i) for i in range(m))
        pair = sum(weight(j-i)*bits[i]*bits[j]
                   for i in range(m) for j in range(i+1, m))
        require(pair == f(n), 'pair energy')
        require(g(n)+f(n) == n+1-(n & 1), 'g/f complement')
        require(f(rev) == pair, 'reversal')
        require(g((1 << m)-1-n) == g(rev), 'complement-reversal')
        y = [n >> i for i in range(m)]
        local = -4*sum(z*z for z in y[1:])
        local += sum((2 if i == 0 else 6)*y[i]*y[i+1] for i in range(m-1))
        local -= 2*sum(y[i]*y[i+2] for i in range(m-2))
        require(local == pair, 'local prefix energy')
        words += 1

out = {'status': 'FINITE-EXACT controls of proved identities; no coverage claim',
       'matrix_universe': 'all digit lengths 1..37',
       'word_universe': 'all padded binary words of lengths 1..12',
       'words': words, 'matrices': matrices,
       'hostile': {'seven_digit_radical': [1,0,-2,-3,-2,0,1],
                   'encoded_value': 1, 'not_a_binary_difference': True}}
path = Path(__file__).with_name('energy_checks.json')
path.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({'matrices': len(matrices), 'words': words, 'passed': True}))
