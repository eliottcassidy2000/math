#!/usr/bin/env python3
"""THM-475 (claudebox-2026-06-11-S1): the DRT flag construction attains the tournament
Barba value 2(n-1)^((n-1)/2) at n == 1 (mod 4). Self-contained, pure python, exact (Bareiss).
Also verifies the doubling corollary det(I + S_{D(T)}) = 2^n det(I+S_T)^2 (THM-447 law)."""

def det_bareiss(mat):
    a = [row[:] for row in mat]; n = len(a); sign = 1; prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for i in range(k + 1, n):
                if a[i][k] != 0: a[k], a[i] = a[i], a[k]; sign = -sign; break
            else: return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]

def paley_S(q):
    QR = {(x * x) % q for x in range(1, q)}
    return [[0 if i == j else (1 if (j - i) % q in QR else -1) for j in range(q)] for i in range(q)]

def gf27_paley_S():
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    def mul(u, v):
        coef = [0] * 5
        for i, ui in enumerate(u):
            for j, vj in enumerate(v): coef[i + j] = (coef[i + j] + ui * vj) % 3
        for k in (4, 3):  # x^3 = x + 1
            coef[k - 2] = (coef[k - 2] + coef[k]) % 3
            coef[k - 3] = (coef[k - 3] + coef[k]) % 3
            coef[k] = 0
        return (coef[0], coef[1], coef[2])
    sq = {mul(e, e) for e in els if e != (0, 0, 0)}
    idx = {e: i for i, e in enumerate(els)}
    sub = lambda u, v: tuple((ui - vi) % 3 for ui, vi in zip(u, v))
    return [[0 if u == v else (1 if sub(v, u) in sq else -1) for v in els] for u in els]

def skewH_from_drt(S0):
    m = len(S0); H = [[0] * (m + 1) for _ in range(m + 1)]; H[0][0] = 1
    for j in range(m): H[0][j + 1] = 1; H[j + 1][0] = -1
    for i in range(m):
        for j in range(m): H[i + 1][j + 1] = (1 if i == j else 0) + S0[i][j]
    return H

def double_skewH(H):
    nn = len(H); H2 = [[0] * (2 * nn) for _ in range(2 * nn)]
    for i in range(nn):
        for j in range(nn):
            H2[i][j] = H[i][j]; H2[i][j + nn] = H[i][j]
            H2[i + nn][j] = -H[j][i]; H2[i + nn][j + nn] = H[j][i]
    return H2

def drt_from_skewH(H):
    n = len(H); D = [H[0][j] for j in range(n)]
    G = [[D[i] * H[i][j] * D[j] for j in range(n)] for i in range(n)]
    return [[G[i + 1][j + 1] - (1 if i == j else 0) for j in range(n - 1)] for i in range(n - 1)]

def check_drt(S0):
    m = len(S0)
    S2 = [[sum(S0[i][k] * S0[k][j] for k in range(m)) for j in range(m)] for i in range(m)]
    return all(S2[i][j] == (1 - (m if i == j else 0)) for i in range(m) for j in range(m))

def flag(S0):
    """Flag(DRT): u,v beat the whole DRT; u beats v. THM-475."""
    m = len(S0); n = m + 2
    S = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m): S[i][j] = S0[i][j]
        S[m][i] = 1; S[i][m] = -1; S[m + 1][i] = 1; S[i][m + 1] = -1
    S[m][m + 1] = 1; S[m + 1][m] = -1
    return S

def M_of_S(S):
    return [[(1 if i == j else 0) + S[i][j] for j in range(len(S))] for i in range(len(S))]

def D_double(S):
    """THM-447 skew-Sylvester double on the S level: M' = [[S, S+I],[S-I, -S]]."""
    m = len(S); n = 2 * m
    S2 = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            S2[i][j] = S[i][j]
            S2[i][j + m] = S[i][j] + (1 if i == j else 0)
            S2[i + m][j] = S[i][j] - (1 if i == j else 0)
            S2[i + m][j + m] = -S[i][j]
    return S2

if __name__ == '__main__':
    S15 = drt_from_skewH(double_skewH(skewH_from_drt(paley_S(7))))
    drts = [(9, 'Paley QR_7', paley_S(7)), (13, 'Paley QR_11', paley_S(11)),
            (17, 'doubling-tower DRT(15)', S15), (25, 'Paley QR_23', paley_S(23)),
            (29, 'GF(27)-Paley', gf27_paley_S())]
    print("== THM-475: flag(DRT(n-2)) attains 2(n-1)^((n-1)/2) ==")
    for n, name, S0 in drts:
        assert check_drt(S0), name
        d = det_bareiss(M_of_S(flag(S0)))
        t = 2 * (n - 1) ** ((n - 1) // 2)
        print(f"n={n:3d} via {name:24s}: det(I+S) = {d} ; target = {t} ; {'ATTAINED' if d == t else 'FAIL'}")
        assert d == t
    print("\n== doubling corollary: det(I+S_D(T)) = 2^n det(I+S_T)^2 ==")
    for n, name, S0 in drts[:2]:
        S = flag(S0); m = len(S)
        d1 = det_bareiss(M_of_S(S)); d2 = det_bareiss(M_of_S(D_double(S)))
        print(f"n={m}: det = {d1}; double(2n={2*m}): det = {d2}; 2^n*det^2 = {2**m * d1 * d1}; {'OK' if d2 == 2**m * d1 * d1 else 'FAIL'}")
        assert d2 == 2 ** m * d1 * d1
    print("\nall checks passed")
