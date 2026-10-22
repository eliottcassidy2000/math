"""Audit E (HYP-9372): faithfulness of the cat map modulo the centre, in BOTH directions; the Ising spin-structure
monomial structure; Ising-type Deligne factorisation of M(4,p').  Auditor for opus-2026-10-22-S55.
Uses the independent builders of c_modular_data.py.
Usage: python e_center.py"""
import math, cmath, itertools
from fractions import Fraction
import numpy as np
import mpmath as mp
import c_modular_data as cd


def Apowers(N, m):
    X = (1 % N, 0, 0, 1 % N); out = [X]
    for _ in range(m):
        X = ((2 * X[0] + X[1]) % N, (X[0] + X[1]) % N, (2 * X[2] + X[3]) % N, (X[2] + X[3]) % N)
        out.append(X)
    return out


def faithful(md):
    o, po, C, kind = cd.orders(md)
    N = md.N; m = cd.ordA(N)
    P = Apowers(N, m)
    central = [j for j in range(1, m + 1) if P[j][1] == 0 and P[j][2] == 0 and P[j][0] == P[j][3]]
    Cn = cd.to_np(C) if kind == 'mp' else C
    scal, X = [], np.eye(len(Cn), dtype=complex)
    for j in range(1, m + 1):
        X = X @ Cn
        if np.abs(X - X[0, 0] * np.eye(len(Cn))).max() < 1e-8:
            scal.append(j)
    fwd = set(scal) <= set(central)          # C^j scalar  =>  A^j central   (what center_check.py tests)
    bwd = set(central) <= set(scal)          # A^j central =>  C^j scalar    (the converse, asserted by the HYP title)
    selfdual = np.allclose(cd.to_np(md.S) @ cd.to_np(md.S), np.eye(md.r), atol=1e-9)
    print(f'{md.name:32s} N {N:4d} ord_N A {m:4d} central j: {[(j, P[j][0]) for j in central]} scalar j: {scal[:6]}'
          f'{"..." if len(scal) > 6 else ""}  po {po}  [C^j scalar => A^j central: {fwd}] [A^j central => C^j scalar: {bwd}]'
          f'  self-dual: {selfdual}')
    return fwd, bwd


print('=== faithfulness modulo the centre, both directions ===')
cases = [cd.minimal(3, 4), cd.su2(2), cd.minimal(4, 5), cd.minimal(4, 7), cd.minimal(4, 9), cd.minimal(4, 11),
         cd.dzn(3), cd.dzn(5), cd.su3(1), cd.su3(2), cd.double_S3(), cd.fib(), cd.su2(4)]
lat = {'A1 (definite)': [[2]], 'A1 + U (Lorentzian)': [[2, 0, 0], [0, 0, 1], [0, 1, 0]], 'A1(-1)': [[-2]],
       'A2 (definite)': [[2, -1], [-1, 2]], 'A2 + U (Lorentzian)': [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
       'A4 (definite)': [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]],
       'D4 (definite)': [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]],
       'golden [[2,1],[1,-2]]': [[2, 1], [1, -2]]}
cases += [cd.pointed(nm, G) for nm, G in lat.items()] + [cd.deligne(cd.fib(), cd.dzn(2), cd.dzn(3))]
res = {md.name: faithful(md) for md in cases}

# ---------------------------------------------------------------- Ising spin-structure basis
print('\n=== Ising in the spin-structure basis u = (chi_0 + chi_psi, chi_0 - chi_psi, sqrt2 chi_sigma) ~ sqrt(theta_3/eta),'
      ' sqrt(theta_4/eta), sqrt(theta_2/eta) ===')
for md in (cd.minimal(3, 4), cd.su2(2)):
    S = cd.to_np(md.S); h = [float(x) for x in md.h]; c = float(md.c)
    # identify labels: unit (h = 0), psi (h = 1/2 mod 1), sigma (the rest)
    i0 = [i for i, x in enumerate(md.h) if x % 1 == 0][0]
    ipsi = [i for i, x in enumerate(md.h) if x % 1 == Fraction(1, 2)][0]
    isg = [i for i in range(3) if i not in (i0, ipsi)][0]
    sgn = np.sign(S[i0, isg].real)
    B = np.zeros((3, 3)); B[0, i0] = B[0, ipsi] = 1; B[1, i0] = 1; B[1, ipsi] = -1; B[2, isg] = math.sqrt(2) * sgn
    T = np.diag([cmath.exp(2j * math.pi * (x - c / 24)) for x in h])          # linear T on characters
    # characters transform as chi(gamma tau) = rho(gamma) chi(tau), rho(s) = S, rho(t) = T; in the u basis: B rho B^-1
    Su = B @ S @ np.linalg.inv(B); Tu = B @ T @ np.linalg.inv(B)
    Cu = Tu @ Su @ np.linalg.inv(Tu) @ np.linalg.inv(Su)

    def monomial(M):
        nz = np.abs(M) > 1e-12
        return bool((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()), [int(np.argmax(np.abs(M[:, j]))) for j in range(3)]
    print(f'{md.name}: S in u-basis monomial/perm {monomial(Su)}, |entries| {np.round(np.abs(Su[np.abs(Su) > 1e-12]), 12)}; '
          f'T in u-basis monomial/perm {monomial(Tu)}; C monomial/perm {monomial(Cu)}; nonzero entries of C: '
          f'{np.round(Cu[np.abs(Cu) > 1e-12], 6)}; C^3 = I: {np.allclose(np.linalg.matrix_power(Cu, 3), np.eye(3))}')
# mod 2 action of A on the three nonzero vectors of (Z/2)^2 (= the three even spin structures, up to the usual shift)
A2 = np.array([[0, 1], [1, 1]])
vecs = [(1, 0), (0, 1), (1, 1)]
print('A mod 2 on nonzero vectors of (Z/2)^2:', {v: tuple(int(x) for x in (A2 @ np.array(v)) % 2) for v in vecs}, '(a 3-cycle)')

# ---------------------------------------------------------------- Deligne factorisation of M(4,p')
print('\n=== Is M(4,p\') an Ising-type x (rank (p\'-1)/2) Deligne product?  test: permutation + diagonal signs matching S ===')


def ising(nu, dsign=1):
    r2 = math.sqrt(2) * dsign
    S = [[0.5, 0.5, r2 / 2], [0.5, 0.5, -r2 / 2], [r2 / 2, -r2 / 2, 0]]
    return cd.MD(f'Ising_{nu}', [[mp.mpf(x) for x in row] for row in S], [0, Fraction(1, 2), Fraction(nu, 16)], Fraction(nu, 2))


def galois_minimal2(pp, sigma):
    """Galois conjugate of M(2,pp) modular data by zeta_{4pp}^sigma acting on T (S conjugated accordingly) -- found by search"""
    pass


def match(big, A_, B_):
    """find a bijection big-labels -> pairs with T equal and S equal up to a diagonal sign conjugation and global sign."""
    Sb = cd.to_np(big.S); hb = [x % 1 for x in big.h]
    P = cd.deligne(A_, B_); Sp = cd.to_np(P.S); hp = [x % 1 for x in P.h]
    if sorted(hb) != sorted(hp) or len(set(hb)) != len(hb):
        return False
    perm = [hp.index(x) for x in hb]
    Sq = Sp[np.ix_(perm, perm)]
    # S_big = eps * D Sq D for diagonal signs D (unit row fixes D up to global): D_i = sign(Sb[0,i] / Sq[0,i]) * eps ...
    for eps in (1, -1):
        D = np.sign((Sb[0] / (eps * Sq[0])).real)
        if np.allclose(Sb, eps * np.outer(D, D) * Sq, atol=1e-9):
            return True
    return False


for pp in (5, 7, 9, 11):
    big = cd.minimal(4, pp)
    others = [cd.minimal(2, pp)] + ([cd.fib()] if pp == 5 else [])
    found = []
    for nu in range(1, 16, 2):
        for ds in (1, -1):
            I = ising(nu, ds)
            for X in others:
                # also try complex-conjugate (Galois) partner: conj S, conj T
                for conjX in (False, True):
                    Y = X if not conjX else cd.MD(X.name + '-bar', mp.matrix([[mp.conj(X.S[i, j]) for j in range(X.r)] for i in range(X.r)]),
                                                     [(-x) % 1 for x in X.h], -X.c)
                    if match(big, I, Y):
                        found.append((I.name, 'd_sigma sign', ds, Y.name))
    print(f'M(4,{pp}): factorisations found {found[:4]}{" ..." if len(found) > 4 else ""}')


# ---------------------------------------------------------------- better: extract the partner X from M(4,p') itself
print('\n=== M(4,p\') = Ising_nu x X with X = {(r,1): r odd} extracted from M(4,p\') (S_X = 2 S restricted) ===')
for pp in (5, 7, 9, 11):
    big = cd.minimal(4, pp)
    reps = big.reps
    idx = {rs: i for i, rs in enumerate(reps)}

    def find(rs):
        r, s = rs
        return idx[(r, s)] if (r, s) in idx else idx[(pp - r, 4 - s)]
    Xl = [find((r, 1)) for r in range(1, pp, 2)]
    Sb = cd.to_np(big.S)
    SX = 2 * Sb[np.ix_(Xl, Xl)]
    hX = [big.h[i] for i in Xl]
    nu = int(big.h[find((1, 2))] % 1 * 16)
    SXmp = mp.matrix([[2 * big.S[i, j] for j in Xl] for i in Xl])
    X = cd.MD(f'X(4,{pp})', SXmp, hX, big.c - Fraction(nu, 2))
    chkX = cd.check_md(X)
    ok = False
    for ds in (1, -1):
        if match(big, ising(nu, ds), X):
            ok = True
    oI, poI, _, _ = cd.orders(ising(nu))
    oX, poX, _, _ = cd.orders(X)
    oB, poB, _, _ = cd.orders(big)
    print(f'M(4,{pp}) = Ising_{nu} x X(rank {X.r}): {ok}; X relations: sym {chkX["sym"]:.0e} unit {chkX["unit"]:.0e} '
          f'ST3 {chkX["ST3"]:.0e} Verlinde int {chkX.get("verlinde_int"):.0e} nonneg {chkX.get("verlinde_nonneg")}; '
          f'po(Ising) {poI}, po(X) {poX}, lcm {math.lcm(poI, poX)}, po(M) {poB}; ord: {oI}, {oX}, ord(M) {oB}; '
          f'X twists {[str(x % 1) for x in hX]}')
