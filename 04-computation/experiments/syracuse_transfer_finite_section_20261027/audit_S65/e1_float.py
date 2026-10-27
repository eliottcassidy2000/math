"""Audit item 5: floating-point claims (a guide only).

(1) numpy |lambda_2| of P_K (K = 8..12) in the natural ordering and after random permutation similarities (the
    spurious values are artefacts and depend on ordering); compare with eps^(1/(K-1)).
(2) the full singular spectrum of P_K (float SVD) against the auditor's closed form: 1 once, sqrt(1 - 2^(m-K)) with
    multiplicity 2^(m-2) for m = 2..K-1, and 0 with multiplicity 2^(K-2).
(3) Tao's chain: numpy |lambda_2| (J = 2..6) and the l^2(pi) norm of P_J - Pi restricted to the units, J = 2..6
    (numpy), and J = 2, 3 at 40 digits (mpmath) as a certificate of the printed digits.
"""
import numpy as np
import mpmath


def v2(n):
    return (n & -n).bit_length() - 1


def P_K(K, p=3):
    mod = 1 << K
    n = mod // 2
    P = np.zeros((n, n))
    for i in range(n):
        a = 2 * i + 1
        t = (p * a + 1) % mod
        if t == 0:
            P[i, :] = 1.0 / n
            continue
        v = v2(t)
        r = (t >> v) % (1 << (K - v))
        for b in range(r, mod, 1 << (K - v)):
            P[i, (b - 1) // 2] = 2.0 ** (-v)
    return P


def ordmod(b, m):
    k, x = 1, b % m
    while x != 1:
        x = x * b % m
        k += 1
    return k


def tao_P(J, p=3, prec=None):
    Mod = p ** J
    L = ordmod(2, Mod)
    h = pow(2, -1, Mod)
    if prec:
        P = mpmath.zeros(Mod, Mod)
        D = mpmath.mpf(2) ** L - 1
    else:
        P = np.zeros((Mod, Mod))
        D = 2.0 ** L - 1
    for x in range(Mod):
        y = (p * x + 1) % Mod
        for a in range(1, L + 1):
            y = y * h % Mod
            P[x, y] += (mpmath.mpf(2) ** (L - a) if prec else 2.0 ** (L - a)) / D
    return P


def main():
    rng = np.random.default_rng(20261027)
    print('E1 (1) numpy eigvals of P_K: |lambda_2| (exact value 0)')
    for K in range(8, 13):
        P = P_K(K)
        ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
        perms = []
        for _ in range(3):
            pi = rng.permutation(len(P))
            Q = P[np.ix_(pi, pi)]
            perms.append(np.sort(np.abs(np.linalg.eigvals(Q)))[::-1][1])
        print(f'  K = {K}: |lambda_2| natural order = {ev[1]:.4f}; three random relabellings: '
              f'{", ".join(f"{x:.4f}" for x in perms)}; (1.1e-16)^(1/(K-1)) = {1.1e-16 ** (1 / (K - 1)):.4f}; '
              f'number of computed |lambda| > 1e-6 besides 1: {int(np.sum(ev[1:] > 1e-6))}')
    print('E1 (2) singular spectrum of P_K vs closed form')
    for K in range(4, 13):
        P = P_K(K)
        s = np.sort(np.linalg.svd(P, compute_uv=False))[::-1]
        pred = [1.0]
        for m in range(2, K):
            pred += [np.sqrt(1 - 2.0 ** (m - K))] * (2 ** (m - 2))
        pred += [0.0] * (2 ** (K - 2))
        pred = np.sort(np.array(pred))[::-1]
        print(f'  K = {K:2d}: max |svd - closed form| = {np.max(np.abs(s - pred)):.2e}; sigma_2 = {s[1]:.10f} vs '
              f'sqrt(1 - 2^-(K-2)) = {np.sqrt(1 - 2.0 ** -(K - 2)):.10f}')
    print('E1 (3) Tao chain: numpy |lambda_2| and ||P - Pi|| in l^2(pi) on the units')
    for J in range(2, 7):
        P = tao_P(J)
        ev = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
        # stationary law: a row of P^J
        PJ = np.linalg.matrix_power(P, J)
        pi = PJ[0]
        units = [x for x in range(3 ** J) if x % 3]
        pu = pi[units]
        Q = P[np.ix_(units, units)]
        Nq = Q - pu[None, :]
        S = np.sqrt(pu)[:, None] * Nq / np.sqrt(pu)[None, :]
        nrm = np.linalg.norm(S, 2)
        print(f'  J = {J}: numpy |lambda_2| = {ev[1]:.4f}; ||P_J - Pi||_(l^2(pi), units) = {nrm:.6f}')
    for J in (2, 3):
        mpmath.mp.dps = 40
        P = tao_P(J, prec=True)
        Mod = 3 ** J
        PJ = P ** J
        pi = [PJ[0, y] for y in range(Mod)]
        units = [x for x in range(Mod) if x % 3]
        n = len(units)
        S = mpmath.zeros(n, n)
        for i, x in enumerate(units):
            for j, y in enumerate(units):
                S[i, j] = mpmath.sqrt(pi[x]) * (P[x, y] - pi[y]) / mpmath.sqrt(pi[y])
        G = S.T * S
        E = mpmath.eigsy(G)[0]
        lam = max(E[i] for i in range(n))
        print(f'  J = {J} (40 digits): ||P_J - Pi||_(l^2(pi), units) = {mpmath.nstr(mpmath.sqrt(lam), 15)}')
    print('E1 DONE')


if __name__ == '__main__':
    main()
