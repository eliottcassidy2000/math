"""Audit item 3 (independent): the chain x -> (p x + 1) 2^-a on Z/p^J, P(a) = 2^-a.

Exact integer matrix M = (2^L - 1) P, L = ord(2 mod p^J), built by summing the geometric weights a = 1..L directly.
Checks (python-flint, exact):
  * rows of M^n depend only on x mod p^(J-n); M^J has equal rows; M^(J-1) does not
  * rank P^n = p^(J-n) (modular lower bound = exact distinct-row upper bound)
  * ranks of N^n = P^n - Pi exactly -> Jordan type on the whole space vs the note's
  * minimal polynomial of M = (x - (2^L - 1)) x^J (exact, small cases)
  * units-only restriction (the recurrent class): ranks and Jordan type (scope remark)
  * stationary law three ways (chain row of P^J, nullspace of P^T - I, Tao's sum by DP) and the note's displayed
    formula sum_j p^(j-1) 2^-(a_j + ... + a_J), which is tested against them
  * the circulant Circ(w): exact determinant +-(2^L - 1)^(L-1) != 0
"""
import math
import sys
import time
from fractions import Fraction

import flint

Q1 = 2305843009213693951
Q2 = 4611686018427387847


def ordmod(b, m):
    if m == 1:
        return 1
    k, x = 1, b % m
    while x != 1:
        x = x * b % m
        k += 1
    return k


def build_M(p, J):
    Mod = p ** J
    L = ordmod(2, Mod)
    h = pow(2, -1, Mod)
    rows = [[0] * Mod for _ in range(Mod)]
    for x in range(Mod):
        y = (p * x + 1) % Mod
        for a in range(1, L + 1):
            y = y * h % Mod
            rows[x][y] += 1 << (L - a)
    for r in rows:
        assert sum(r) == (1 << L) - 1
    return rows, L


def nrank(M, q):
    return flint.nmod_mat(M.nrows(), M.ncols(), [int(x) % q for x in M.entries()], q).rank()


def distinct(rowsets):
    return len(set(tuple(int(x) for x in r) for r in rowsets))


def analyse(p, J, minpoly=False, units_scope=False):
    t0 = time.time()
    rows, L = build_M(p, J)
    Mod = p ** J
    M = flint.fmpz_mat(rows)
    D = (1 << L) - 1
    I = flint.fmpz_mat(Mod, Mod, [1 if i == j else 0 for i in range(Mod) for j in range(Mod)])
    Mn = [I]
    for s in range(1, J + 1):
        Mn.append(Mn[-1] * M)
    res = {'L': L}
    # rows of M^n depend only on x mod p^(J-n)
    for s in range(1, J + 1):
        Ls = Mn[s].tolist()
        md = p ** (J - s)
        for x in range(Mod):
            assert Ls[x] == Ls[x % md], (p, J, s, x)
    lastrows = Mn[J].tolist()
    assert all(r == lastrows[0] for r in lastrows)
    prev = Mn[J - 1].tolist()
    assert any(r != prev[0] for r in prev)
    # ranks of P^n and of N^n
    piD = [int(v) for v in lastrows[0]]               # D^J * pi
    PiD = flint.fmpz_mat([piD] * Mod)                 # D^J * Pi
    rP, rN = [], []
    for s in range(0, J + 1):
        low = max(nrank(Mn[s], Q1), nrank(Mn[s], Q2))
        up = distinct(Mn[s].tolist())
        assert low == up == p ** (J - s), (p, J, s, low, up)
        rP.append(low)
        if s == 0:
            rN.append(Mod)                            # N^0 = I
            continue
        Ns = Mn[s] * (D ** (J - s)) - PiD               # D^J (P^s - Pi)
        lowN = max(nrank(Ns, Q1), nrank(Ns, Q2))
        assert lowN == low - 1, (p, J, s, lowN, low)     # upper bound rank P^s - 1 (see audit text)
        rN.append(lowN)
    # N = P - Pi nilpotent of index exactly J on the whole space
    assert rN[J] == 0 and rN[J - 1] > 0
    rN.append(0)
    ge = {s: rN[s - 1] - rN[s] for s in range(1, J + 1)}
    blocks = {s: ge[s] - ge.get(s + 1, 0) for s in range(1, J + 1)}
    claim = {s: (p - 1) ** 2 * p ** (J - 1 - s) for s in range(1, J)}
    claim[J] = p - 1
    claim_full = dict(claim)
    claim_full[1] = claim_full.get(1, 0) + 1          # + the constants (N 1 = 0) on the whole space
    if J == 1:
        claim_full = {1: p - 1 + 1}
    assert blocks == claim_full, (p, J, blocks, claim_full)
    res['blocks_full'] = blocks
    if minpoly:
        mp = M.minpoly()
        x = flint.fmpz_poly([0, 1])
        assert mp == (x - D) * x ** J, (p, J, mp)
        res['minpoly'] = True
    if units_scope:
        units = [x for x in range(Mod) if x % p]

        def sub(A):
            T = A.tolist()
            return flint.fmpz_mat([[T[i][j] for j in units] for i in units])
        rU = []
        for s in range(0, J + 1):
            S = sub(Mn[s])
            rU.append(max(nrank(S, Q1), nrank(S, Q2)))
            assert rU[-1] == distinct(S.tolist())
        res['units_ranks'] = rU
        rNu = [r - 1 for r in rU[1:]]
        rNu = [len(units)] + rNu
        rNu[-1] = 0
        rNu.append(0)
        geu = {s: rNu[s - 1] - rNu[s] for s in range(1, J + 1)}
        res['units_blocks'] = {s: geu[s] - geu.get(s + 1, 0) for s in range(1, J + 1)}
    res['pi'] = {y: Fraction(int(piD[y]), D ** J) for y in range(Mod) if piD[y]}
    res['secs'] = time.time() - t0
    return res


def tao_sum_law(p, J, reversed_note=False):
    """Law of sum_{m=1..J} p^(m-1) 2^-(a_1+...+a_m) mod p^J (Tao's Syrac form) by exact DP over (s mod L, sum).
    With reversed_note=True: the note's displayed formula sum_{j=1..J} p^(j-1) 2^-(a_j+...+a_J)."""
    Mod = p ** J
    L = ordmod(2, Mod)
    h = pow(2, -1, Mod)
    hp = [pow(h, s, Mod) for s in range(L)]
    W = [(1 << (L - r)) for r in range(1, L + 1)]     # weights for a = r mod L (over D = 2^L - 1)
    D = (1 << L) - 1
    st = {(0, 0): 1}
    order = list(range(1, J + 1)) if not reversed_note else list(range(J, 0, -1))
    for m in order:
        nst = {}
        for (s, tot), w in st.items():
            for r in range(1, L + 1):
                s2 = (s + r) % L
                t2 = (tot + p ** (m - 1) * hp[s2]) % Mod
                k = (s2, t2)
                nst[k] = nst.get(k, 0) + w * W[r - 1]
        st = nst
    law = {}
    for (s, tot), w in st.items():
        law[tot] = law.get(tot, 0) + w
    return {y: Fraction(w, D ** J) for y, w in law.items()}


def stationary_nullspace(p, J):
    rows, L = build_M(p, J)
    Mod = p ** J
    D = (1 << L) - 1
    A = flint.fmpz_mat(rows).transpose()
    B = A - flint.fmpz_mat(Mod, Mod, [D if i == j else 0 for i in range(Mod) for j in range(Mod)])
    X, nul = B.nullspace()
    assert nul == 1
    v = [int(X[i, 0]) for i in range(Mod)]
    s = sum(v)
    return {y: Fraction(v[y], s) for y in range(Mod) if v[y]}


def circulant_det(L):
    w = [1 << (L - j) for j in range(1, L + 1)]       # first row: w(1..L)
    C = flint.fmpz_mat([[w[(j - i) % L] for j in range(L)] for i in range(L)])
    d = int(C.det())
    return d, abs(d) == ((1 << L) - 1) ** (L - 1)


def main():
    t0 = time.time()
    print('B1: Tao chain x -> (p x + 1) 2^-a on Z/p^J, exact (python-flint)')
    cases = [(3, J) for J in range(1, 6)] + [(5, 1), (5, 2), (5, 3), (7, 1), (7, 2), (7, 3), (9, 1), (9, 2), (11, 2),
                                             (13, 2), (15, 2), (17, 2), (23, 2), (31, 2), (21, 2)]
    for p, J in cases:
        mp = (p ** J <= 125)
        us = (p ** J <= 343)
        r = analyse(p, J, minpoly=mp, units_scope=us)
        Mod = p ** J
        line = (f'p = {p:2d} J = {J}: L = ord(2 mod p^J) = {r["L"]} (index of <2> in units = '
                f'{(Mod - Mod // p) // r["L"]}); rows of P^n depend only on x mod p^(J-n); P^J = Pi, P^(J-1) != Pi; '
                f'rank P^n = p^(J-n) exact; Jordan type of N on Z/p^J (whole space) = {r["blocks_full"]} '
                f'(= note + one 1x1 block from constants)')
        if mp:
            line += '; minpoly(M) = (x - (2^L - 1)) x^J exactly'
        if us:
            line += f'; units-only chain: ranks {r["units_ranks"]}, Jordan type {r["units_blocks"]}'
        line += f'  [{r["secs"]:.1f}s]'
        print(line)
        sys.stdout.flush()
        if p == 3 and J == 2:
            law = r['pi']
            got = [law.get(y, Fraction(0)) for y in range(9)]
            tao = [Fraction(x, 63) for x in (0, 8, 16, 0, 11, 4, 0, 2, 22)]
            print(f'    stationary law mod 9 = {[str(g) for g in got]}; equals Tao (0,8,16,0,11,4,0,2,22)/63: {got == tao}')
            assert got == tao
    # stationary law three ways, and the note's displayed formula
    for J in (1, 2, 3, 4):
        chain = analyse(3, J)['pi']
        nul = stationary_nullspace(3, J)
        tao = tao_sum_law(3, J)
        note = tao_sum_law(3, J, reversed_note=True)
        print(f'p = 3 J = {J}: stationary law: chain row of P^J == nullspace of (P^T - I): {chain == nul}; '
              f'== Tao sum_m 3^(m-1) 2^-(a_1+..+a_m): {chain == tao}; '
              f'== note formula sum_j 3^(j-1) 2^-(a_j+..+a_J): {chain == note}')
        assert chain == nul == tao
        if J == 2:
            m3 = lambda law: [sum((w for y, w in law.items() if y % 3 == c), Fraction(0)) for c in range(3)]
            print(f'    mod-3 marginal: Syrac {[str(x) for x in m3(tao)]} vs note formula {[str(x) for x in m3(note)]}; '
                  f'note formula law mod 9 = {[str(note.get(y, 0)) for y in range(9)]}')
    for p, J in ((5, 2), (7, 2)):
        chain = analyse(p, J)['pi']
        tao = tao_sum_law(p, J)
        note = tao_sum_law(p, J, reversed_note=True)
        print(f'p = {p} J = {J}: chain stationary law == sum_m p^(m-1) 2^-(a_1+..+a_m): {chain == tao}; '
              f'== note formula: {chain == note}')
        assert chain == tao
    for L in (2, 4, 6, 18, 54, 3, 21):
        d, ok = circulant_det(L)
        print(f'Circ(w), L = {L}: |det| = (2^L - 1)^(L-1): {ok} (sign {"+" if d > 0 else "-"})')
        assert ok
    print(f'B1 ALL CHECKS PASSED  [{time.time() - t0:.0f}s]')


if __name__ == '__main__':
    main()
