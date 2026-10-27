"""Audit item 3: Tao's chain mod 3^6, all 729 starts, exact integers (flint): M = (2^486 - 1) P.
Checks: rows of M^n depend only on x mod 3^(6-n) (n = 1..6), M^6 has equal rows, M^5 does not; rank P^n = 3^(6-n)
(modular lower bound mod a 61-bit prime = exact distinct-row upper bound)."""
import time

import flint

Q1 = 2305843009213693951


def main():
    t0 = time.time()
    J, p = 6, 3
    Mod = p ** J
    L = 2 * 3 ** (J - 1)
    h = pow(2, -1, Mod)
    rows = [[0] * Mod for _ in range(Mod)]
    for x in range(Mod):
        y = (p * x + 1) % Mod
        for a in range(1, L + 1):
            y = y * h % Mod
            rows[x][y] += 1 << (L - a)
    M = flint.fmpz_mat(rows)
    P = M
    for n in range(1, J + 1):
        if n > 1:
            P = P * M
        T = P.tolist()
        md = p ** (J - n)
        assert all(T[x] == T[x % md] for x in range(Mod)), n
        distinct = len(set(tuple(int(v) for v in r) for r in T))
        r = flint.nmod_mat(Mod, Mod, [int(v) % Q1 for v in P.entries()], Q1).rank()
        assert distinct == r == p ** (J - n), (n, distinct, r)
        print(f'n = {n}: rows of P^n depend only on x mod 3^{J - n}; rank P^n = {r} = 3^(J-n)  [{time.time() - t0:.0f}s]',
              flush=True)
    print('J = 6: P^6 = Pi exactly for all 729 starts; P^5 != Pi; ranks exact')


if __name__ == '__main__':
    main()
