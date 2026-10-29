"""Exact characteristic-two exterior-square transport, including F64/F256.

All field units and square roots are exhaustive. Matrix probes are a stated
deterministic finite bank, not an exhaustive Sp4 census beyond F2.
"""
import argparse
import json
import random
from pathlib import Path


def need(ok, why):
    if not ok:
        raise ValueError(why)


class Field:
    def __init__(self, m, polynomial):
        self.m, self.q, self.polynomial = m, 1 << m, polynomial
        self.table = [[self.raw(a, b) for b in range(self.q)] for a in range(self.q)]
        for a in range(1, self.q):
            need(self.power(a, self.q-1) == 1, 'every nonzero field element is a unit')
        self.sqrt = [self.power(a, self.q//2) for a in range(self.q)]
        need(all(self.mul(r, r) == a for a, r in enumerate(self.sqrt)), 'all square roots')

    def raw(self, a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.polynomial
        return result

    def mul(self, a, b):
        return self.table[a][b]

    def power(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            n >>= 1
        return out


def identity():
    return [[int(i == j) for j in range(4)] for i in range(4)]


def xor_sum(values):
    out = 0
    for v in values:
        out ^= v
    return out


def product(A, B, f):
    return [[xor_sum(f.mul(A[i][k], B[k][j]) for k in range(4))
             for j in range(4)] for i in range(4)]


def frobenius(A, f, k=1):
    return [[f.power(x, 1 << k) for x in row] for row in A]


def transvection(v, c, f):
    w = [v[2], v[3], v[0], v[1]]
    return [[int(i == j) ^ f.mul(c, f.mul(v[i], w[j]))
             for j in range(4)] for i in range(4)]


def symplectic(A, f):
    for i in range(4):
        for j in range(4):
            value = xor_sum([f.mul(A[0][i], A[2][j]), f.mul(A[2][i], A[0][j]),
                             f.mul(A[1][i], A[3][j]), f.mul(A[3][i], A[1][j])])
            if value != int((i, j) in [(0, 2), (2, 0), (1, 3), (3, 1)]):
                return False
    return True


def psi(A, f):
    pairs = [(0, 1), (0, 3), (1, 2), (2, 3)]
    def minor(rows, cols):
        i, j = rows
        k, l = cols
        return f.mul(A[i][k], A[j][l]) ^ f.mul(A[i][l], A[j][k])
    need(all(minor(r, (0, 2)) == minor(r, (1, 3)) for r in pairs), 'nucleus quotient')
    N = [[minor(r, c) for c in pairs] for r in pairs]
    perm = [0, 1, 3, 2]  # (a,b,c,d) -> (a,b,d,c), restoring standard omega.
    return [[N[i][j] for j in perm] for i in perm]


def run_field(m, polynomial):
    f = Field(m, polynomial)
    rng = random.Random(29102026 + m)
    bank = [identity()]
    for _ in range(128):
        v = [rng.randrange(f.q) for _ in range(4)]
        T = transvection(v, rng.randrange(f.q), f)
        bank.append(product(bank[-1], T, f))
    for A in bank:
        need(symplectic(A, f), 'input symplectic')
        B = psi(A, f)
        need(symplectic(B, f), 'output symplectic')
        need(psi(B, f) == frobenius(A, f), 'exceptional square equals Frobenius')
        need(psi(frobenius(A, f), f) == frobenius(B, f), 'Frobenius commutation')
    for A, B in zip(bank, bank[1:]):
        need(psi(product(A, B, f), f) == product(psi(A, f), psi(B, f), f), 'homomorphism')
    witness = transvection([1, 0, 0, 0], 1 if m == 1 else 2, f)
    orbit, A = [], witness
    while A not in orbit:
        orbit.append(A)
        A = psi(A, f)
    need(A == witness and len(orbit) == 2*m, 'full phase witness')
    twists = [k for k in range(m) if (2*k+1) % m == 0]
    for k in range(m):
        shifted = lambda B: psi(frobenius(B, f, k), f)
        need((shifted(shifted(witness)) == witness) == (k in twists), 'twist obstruction')
    return dict(degree=m, field_size=f.q, modulus=polynomial,
                lagrangian_count=(f.q+1)*(f.q*f.q+1), matrix_probe_count=len(bank),
                full_phase_period=len(orbit), involutive_frobenius_twists=twists,
                phase_witness=witness, phase_orbit=orbit)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    result = [run_field(m, poly) for m, poly in [(1, 3), (2, 7), (3, 11), (6, 67), (8, 283)]]
    path = Path(__file__).with_suffix('.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == result, 'stored field receipt')
    print(json.dumps([{k: v for k, v in row.items() if k != 'phase_orbit' and k != 'phase_witness'}
                      for row in result]))


if __name__ == '__main__':
    main()
