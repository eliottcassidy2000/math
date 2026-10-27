"""Audit item 4 (independent): the T-map contrasts.

(a) T-map random chain on Z/p^J: x -> x/2 or (p x + 1)/2, probability 1/2 each.  The integer matrix 2P has (claimed,
    after the substitution mu = 2 lambda) characteristic polynomial (mu - 2) prod_d (mu^(L_d) - 1)^(c_d), L_d = ord(2 mod
    p^d), c_d = phi(p^d)/L_d; for p = 3 this is the note's (lambda - 1) prod (lambda^(L_d) - 2^-(L_d)).  Exact (flint).
    Also: P^n != Pi for every n (eigenvalues of modulus 1/2), and the identity with THM-4631 on the line u = v = z.
(b) T-map Ulam (de Bruijn) matrix on Z/2^K: row a uniform on the two lifts of T(a) mod 2^(K-1).  Exact: P^K = Pi,
    P^(K-1) != Pi, rank P^n = 2^(K-n), Jordan type (one block K, 2^(K-1-s) blocks of size s on the mean-zero space).
"""
import sys
import time

import flint

Q1 = 2305843009213693951


def ordmod(b, m):
    k, x = 1, b % m
    while x != 1:
        x = x * b % m
        k += 1
    return k


def tmap_chain_2P(p, J):
    Mod = p ** J
    h = pow(2, -1, Mod)
    rows = [[0] * Mod for _ in range(Mod)]
    for x in range(Mod):
        rows[x][x * h % Mod] += 1
        rows[x][(p * x + 1) * h % Mod] += 1
    return flint.fmpz_mat(rows)


def part_a():
    mu = flint.fmpz_poly([0, 1])
    for p, Js in ((3, (1, 2, 3, 4, 5)), (5, (1, 2, 3)), (7, (1, 2, 3)), (17, (1, 2)), (31, (1, 2))):
        for J in Js:
            t0 = time.time()
            A = tmap_chain_2P(p, J)
            cp = A.charpoly()
            target = mu - 2
            for d in range(1, J + 1):
                Ld = ordmod(2, p ** d)
                cd = (p - 1) * p ** (d - 1) // Ld
                target *= (mu ** Ld - 1) ** cd
            ok = (cp == target)
            # THM-4631 at j = 0 on the line u = v = z: det(I - z U) = (1 - z) prod (1 - (z/2)^(L_d))^(c_d)
            # <=> det(lambda - U) = (lambda - 1) prod (lambda^(L_d) - 2^-(L_d))^(c_d)  (same statement as above)
            print(f'  p = {p:2d} J = {J}: charpoly(2P) == (mu - 2) prod_d (mu^(L_d) - 1)^(c_d): {ok}; '
                  f'L_d = {[ordmod(2, p ** d) for d in range(1, J + 1)]}, c_d = '
                  f'{[(p - 1) * p ** (d - 1) // ordmod(2, p ** d) for d in range(1, J + 1)]}  [{time.time() - t0:.1f}s]')
            assert ok
            sys.stdout.flush()
    # 1 is a simple root, all others modulus 1/2 (in lambda); hence P^n != Pi for all n and SLEM = 1/2 exactly.
    print('  -> eigenvalues of P: 1 (simple) and 2^-1 x (L_d-th roots of unity): SLEM exactly 1/2, never exactly mixed')


def debruijn_2P(K):
    n = 1 << K
    rows = [[0] * n for _ in range(n)]
    for a in range(n):
        t = a // 2 if a % 2 == 0 else (3 * a + 1) // 2
        base = t % (n // 2)
        rows[a][base] += 1
        rows[a][base + n // 2] += 1
    return flint.fmpz_mat(rows)


def part_b():
    for K in range(2, 12):
        t0 = time.time()
        A = debruijn_2P(K)
        n = 1 << K
        Jm = flint.fmpz_mat(n, n, [1] * (n * n))
        pw = [flint.fmpz_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)])]
        for s in range(1, K + 1):
            pw.append(pw[-1] * A)
        okK = (pw[K] == Jm)                         # (2P)^K = 2^K Pi = J
        okK1 = (pw[K - 1] * 2 != Jm)                # (2P)^(K-1) = 2^(K-1) Pi  <=> 2 (2P)^(K-1) = J
        ranks = []
        for s in range(K + 1):
            r = flint.nmod_mat(n, n, [int(x) % Q1 for x in pw[s].entries()], Q1).rank()
            # upper bound: rows constant on classes mod 2^(K-s)
            Ls = pw[s].tolist()
            st = 1 << (K - s)
            assert all(row[j] == row[(j + st) % n] for row in Ls for j in range(n))
            assert r == (1 << (K - s)), (K, s, r)
            ranks.append(r)
        rN = [r - 1 for r in ranks]
        rN[0] = n
        ge = {s: rN[s - 1] - rN[s] for s in range(1, K + 1)}
        blocks = {s: ge[s] - ge.get(s + 1, 0) for s in range(1, K + 1)}
        claim = {s: 2 ** (K - 1 - s) for s in range(1, K)}
        claim[K] = 1
        claim[1] = claim.get(1, 0) + 1              # + constants on the whole space
        if K == 1:
            claim = {1: 2}
        print(f'  K = {K:2d}: (2P)^K = J (P^K = Pi): {okK}; P^(K-1) != Pi: {okK1}; rank P^n = 2^(K-n); '
              f'Jordan type of N on the whole space {blocks} == claim: {blocks == claim}  [{time.time() - t0:.1f}s]')
        assert okK and okK1 and blocks == claim
        sys.stdout.flush()


def main():
    print('C1 (a) T-map random chain on Z/p^J, exact characteristic polynomials')
    part_a()
    print('C1 (b) T-map Ulam (de Bruijn) matrix on Z/2^K, exact')
    part_b()
    print('C1 ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
