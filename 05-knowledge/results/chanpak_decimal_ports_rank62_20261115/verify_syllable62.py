#!/usr/bin/env python3
"""Exact kernel- and padding-bridge replay for the pinned rank-62 alphabet."""
from hashlib import sha256
from pathlib import Path
from sys import argv

from sympy import Matrix, eye, zeros

EXTRAS = '0.12 0.0.3 0.0.13 0.0.123 0.13.12'.split()
NAMES = '0.2 1.2 0.123 2.13 0.1.3 0.3.1 0.13.2 0.1.123 0.3.12 0.12.3 1.2.3 1.12.3 0.1.23'.split()
D = 62
HASH = '0e1e72c81d24f7117f3d0763830e38c20181fafaedd18d5087ba3c7762a9cdea'

def check(ok, message):
    if not ok:
        raise AssertionError(message)

def read(path):
    raw = Path(path).read_bytes()
    check(sha256(raw).hexdigest() == HASH, 'pinned rank-62 data')
    t = iter(map(int, raw.split()))
    check(next(t) == D and next(t) == 8, 'shape')
    v = Matrix(D, 1, [next(t) for _ in range(D)])
    h = Matrix(1, D, [next(t) for _ in range(D)])
    mats = {}
    for _ in range(8):
        k = next(t)
        mats[k] = Matrix(D, D, [next(t) for _ in range(D * D)])
    check(set(mats) == set(range(8)) and next(t, None) is None, 'all letters')
    return v, h, mats

def mult(seq, A):
    dim = A[seq[0]].rows if seq else next(iter(A.values())).rows
    X = eye(dim)
    for i in seq:
        X = A[i] * X
    return X

def mask(a):
    return sum(1 << (int(c)-1) for c in a if c != '0')

def basis(S):
    return S[:, list(S.rref()[1])]

def main(path):
    v, h, M = read(path)
    v = M[0] * v
    U = {i: mult([mask(a) for a in w.split('.')], M) for i, w in enumerate(NAMES)}
    P = Matrix.hstack(v, *(U[i] * v for i in (1,0,2,3,4,7,9)))
    Q = Matrix.vstack(h, *(h * U[i] for i in (0,1,2,3,5,6,7)))
    X = Q * P
    check(X.det() != 0, 'Hankel block')
    R = {i: Q * A * P * X.inv() for i, A in U.items()}

    E = basis(Matrix.hstack(*U.values()))
    check(E.cols == 46 and (Q * E).rank() == 8, 'carrier and readout dimension')
    K = E * Matrix.hstack(*(Q * E).nullspace())
    check(K.cols == 38, 'hidden dimension')
    universally_capped = []
    for i, A in U.items():
        F = Q * A - R[i] * Q
        check(F * E == zeros(8,46), ('intertwining', i))
        check(F.rank() == 1, ('ambient residual rank', i))
        check(Matrix.hstack(R[i],F).rank() == R[i].rank(), ('no left padding', i))
        N = Matrix.hstack(*(Q * A).nullspace())
        B = A * N
        ok = [a for a in range(8) if M[a] * B == zeros(D,N.cols)]
        check(ok == list(range(4)) or ok == list(range(8)), ('literal cap list', i))
        if len(ok) == 8:
            universally_capped.append(i)
    check(universally_capped == [1,4,6,7,8,9,11,12], 'all-mask caps')
    K1 = basis(Matrix.hstack(*(A * K for A in U.values())))
    check(K1.cols == 12 and Q * K1 == zeros(8,12), 'one-unit hidden band')
    for i, A in U.items():
        check(A * K1 == zeros(D,12), ('kills hidden band', i))

    # Full rank never sees a zero matrix: h * M_W * v_all is a true path count.
    w = (5,5)
    check(mult(w, R).rank() == 1, 'w reduced rank')
    check(mult(w, U).rank() == 2, 'zero trailing pads insufficient')
    check((M[0] * mult(w,U)).rank() == 1, 'one literal cap suffices')

    # Explore the full 49 phases from the known short floor resets.
    good = set()
    phases = set()
    for label in NAMES + EXTRAS:
        w = label.split('.')
        for i in range(len(w)):
            c = w[i:] + w[:i]
            label2 = '.'.join(c)
            phases.add(label2)
            A = mult([mask(k) for k in c], M)
            B = Q * A * P * X.inv()
            F = Q * A - B * Q
            N = Matrix.hstack(*(Q * A).nullspace())
            AN = A * N
            if (F * E == zeros(8,46)
                    and Matrix.hstack(E,A).rank() == 46
                    and Matrix.hstack(B,F).rank() == B.rank()
                    and all(M[a] * AN == zeros(D,N.cols) for a in range(4))):
                good.add(label2)
    check(len(phases) == 49, 'phase bank')
    check(good == set(NAMES + EXTRAS + ['2.1','1.0.3','2.3.1','0.3.0','0.13.0']), 'exact 23-phase bridge')
    print('rank-62 to rank-8 bridge: E=46, visible=8, hidden=38; no-left-pad and literal right-cap identities verified')
    print('sharp literal one-cap control: 0.3.1 / 0.3.1 then 0; full ranks 2,1')
    print('23 compatible phases out of 49, including all 18 unrotated representatives')

if __name__ == '__main__':
    main(argv[1])
