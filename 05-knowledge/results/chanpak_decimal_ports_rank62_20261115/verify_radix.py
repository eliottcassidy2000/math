#!/usr/bin/env python3
"""Verify exclusion of a literal repeated 13-syllable motif for 1000^k+1."""
from subprocess import run
from sys import argv
from sympy import Matrix, eye
from verify_syllable62 import NAMES, check, mask, mult, read

def main(path, executable):
    v, h, M = read(path)
    v = M[0] * v
    U = {i: mult([mask(a) for a in w.split('.')], M) for i, w in enumerate(NAMES)}
    P = Matrix.hstack(v, *(U[i] * v for i in (1,0,2,3,4,7,9)))
    Q = Matrix.vstack(h, *(h * U[i] for i in (0,1,2,3,5,6,7)))
    Xi = (Q * P).inv()
    R = {i: Q * A * P * Xi for i, A in U.items()}
    check(17699**13 > 1000 * 10000**13, 'growth cutoff')
    p = 101
    flat = []
    for A in R.values():
        for x in A:
            check(int(x.q) % p != 0, 'denominator prime to test modulus')
            flat.append(str(int(x.p) * pow(int(x.q), -1, p) % p))
    out = run([executable], input=' '.join(flat), text=True, capture_output=True, check=True)
    check(out.stderr.split() == ['99420', '945', '21'], 'word count and modular filters')
    candidates = [tuple(map(int, line.split())) for line in out.stdout.splitlines()]
    check(len(candidates) == len(set(candidates)) == 21, 'modular candidates')
    I = eye(8)
    pairs = set()
    for w in candidates:
        T = mult(w, R)
        pair = ((1000 * I - T).det(), (I - T).det())
        check(pair != (0,0), ('exact candidate survived', w))
        pairs.add(pair)
    check(pairs == {(-4324134008000000000000000,-6464),(-3231798000000000000000000,-3030)}, 'determinant controls')
    print('99420 weighted motifs of cost ≤12; modulo 101 filters 945 -> 21; no exact survivor')

if __name__ == '__main__':
    main(argv[1], argv[2])
