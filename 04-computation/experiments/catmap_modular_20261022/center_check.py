"""Is the cat map detected faithfully modulo the centre?  Session opus-2026-10-22-S55.
For each modular data of the census with projective conductor N >= 2: find the least j such that A^j is a central
matrix lambda*I in SL2(Z/N) (lambda^2 = 1 mod N), and compare with the projective order of C = T S T^-1 S^-1.
Faithfulness modulo the centre on <A> means: if C^j is scalar then A^j is central (every scalar power of C is a central
power of A).  The converse is false in general (a central power of A need not act by a scalar: -I acts as the charge
conjugation S^2); it is reported separately.  The lattices and the Deligne product were added after audit_S55.
Usage: python center_check.py"""
import json, math
import numpy as np
import catmap_mtc as cm

CAT = [[2, 1], [1, 1]]


def powers(N, m):
    X = [[1 % N, 0], [0, 1 % N]]; out = [X]
    for _ in range(m):
        X = [[(X[0][0] * 2 + X[0][1]) % N, (X[0][0] + X[0][1]) % N], [(X[1][0] * 2 + X[1][1]) % N, (X[1][0] + X[1][1]) % N]]
        out.append(X)
    return out


def builders():
    yield cm.fibonacci()
    for p in range(2, 7):
        for pp in range(p + 1, 13):
            if math.gcd(p, pp) == 1 and p * pp <= 66:
                yield cm.minimal(p, pp)
    for k in range(1, 61):
        yield cm.su2(k)
    for k in range(1, 11):
        yield cm.su3(k)
    for name, G in cm.LATTICES:
        yield cm.pointed(name, G)
    for N in range(2, 31):
        yield cm.double_ZN(N)
    yield cm.double_S3()
    yield cm.deligne(cm.fibonacci(), cm.deligne(cm.double_ZN(2), cm.double_ZN(3)))


bad, conv_bad, rows = [], [], 0
for md in builders():
    C = md.catmap(); N = md.Nproj
    po = cm.proj_order(C, 10 ** 4)
    if N < 2:
        continue
    m = cm.mat_order(CAT, N)
    P = powers(N, m)
    central = [j for j in range(1, m + 1) if P[j][0][1] == 0 and P[j][1][0] == 0 and P[j][0][0] == P[j][1][1]]
    jc = central[0]
    lam = P[jc][0][0]
    # projective triviality of every central power: C^j scalar for j in central?
    scal = [j for j in range(1, m + 1) if np.allclose(np.linalg.matrix_power(C, j), np.linalg.matrix_power(C, j)[0, 0] * np.eye(len(C)), atol=1e-7)]
    ok = set(scal) <= set(central)
    conv = set(central) <= set(scal)
    rows += 1
    if not ok:
        bad.append(md.name)
    if not conv:
        conv_bad.append(md.name)
    print(f'{md.name:30s} N {N:4d} ord_N A {m:4d}  first central power j = {jc:4d} (A^j = {lam} I)  proj ord C {po:4d}  '
          f'scalar powers are central powers: {ok}; central powers act by scalars: {conv}')
print(f'\n{rows} modular data; faithful modulo the centre on <A>: {rows - len(bad)}; exceptions: {bad}')
print(f'converse (every central power of A acts by a scalar) fails in {len(conv_bad)}: {conv_bad}')
