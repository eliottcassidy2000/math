"""Scope remark for THM-4822: Jordan type of N = P - Pi restricted to the recurrent class R = supp(pi)
= {y : y mod p in <2 mod p>} (size r p^(J-1), r = |<2 mod p>|), on the pi-mean-zero space.
Auditor's formula (from rank P^n|_R = r p^(J-1-n) for n <= J-1): r - 1 blocks of size J, r(p-2) + 1 blocks of size J-1,
r (p-1)^2 p^(J-2-s) blocks of size s <= J-2  (J >= 2).  For p = 3 (r = 2): 1, 3, 8 3^(J-2-s).  Exact (flint)."""
import flint
from b1_tao_chain import build_M, nrank, Q1, Q2


def subgroup2(p):
    s, x = set(), 1
    while True:
        x = x * 2 % p
        s.add(x)
        if x == 1:
            return s


def main():
    for p, J in ((3, 2), (3, 3), (3, 4), (3, 5), (5, 2), (5, 3), (7, 2), (7, 3), (9, 2), (15, 2), (17, 2), (31, 2)):
        rows, L = build_M(p, J)
        Mod = p ** J
        H = subgroup2(p)
        r = len(H)
        R = [x for x in range(Mod) if x % p in H]
        M = flint.fmpz_mat(rows)
        D = (1 << L) - 1
        Mn = M
        pows = [None]
        for n in range(1, J + 1):
            if n > 1:
                Mn = Mn * M
            pows.append(Mn.tolist())
        pi = pows[J][0]
        assert all(pi[y] > 0 for y in R) and all(pi[y] == 0 for y in range(Mod) if y not in set(R))  # supp pi = R
        ranksN = [len(R) - 1]                       # N^0 on the pi-mean-zero space of R
        for n in range(1, J + 1):
            T = pows[n]
            Nn = flint.fmpz_mat([[T[x][y] * D ** (J - n) - pi[y] for y in R] for x in R])
            ranksN.append(max(nrank(Nn, Q1), nrank(Nn, Q2)))
        assert ranksN[J] == 0
        ranksN.append(0)
        ge = {s: ranksN[s - 1] - ranksN[s] for s in range(1, J + 1)}
        blocks = {s: ge[s] - ge.get(s + 1, 0) for s in range(1, J + 1)}
        # ranksN[0] = |R| - 1 is the dimension of the pi-mean-zero space of R and N kills the constants, so rank N^n on
        # the R x R matrix equals rank on the mean-zero space: 'blocks' is already the type on the mean-zero space
        pred = {J: r - 1, J - 1: r * (p - 2) + 1}
        for s in range(1, J - 1):
            pred[s] = r * (p - 1) ** 2 * p ** (J - 2 - s)
        pred = {s: v for s, v in pred.items() if v}
        got = {s: v for s, v in blocks.items() if v}
        print(f'p = {p:2d} J = {J}: r = |<2 mod p>| = {r}, |R| = {len(R)}; Jordan type on R (pi-mean-zero) = {got}; '
              f'formula {pred}: {got == pred}')
        assert got == pred
    print('B4 ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
