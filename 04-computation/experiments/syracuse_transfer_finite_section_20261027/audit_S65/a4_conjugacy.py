"""Audit item 2(6)/7: the universality in p is an isomorphism of Markov chains, inherited from the classical
parity-vector conjugacy (Terras 1976; Lagarias 1985 for p = 3; the same proof for every odd p).

Q_p(a) = (parities of T_p^i(a), i = 0..K-1), T_p(x) = x/2 or (p x + 1)/2, is a bijection Z/2^K -> {0,1}^K (Terras);
psi = Q_3^(-1) o Q_p maps odd residues to odd residues.  Claim: P_K^(p)(a, b) = P_K^(3)(psi(a), psi(b)) for all odd a, b
(the Ulam matrices are permutation-conjugate), and psi is a tree automorphism (preserves every class mod 2^m).
Also: the T-map Ulam (de Bruijn) matrices are permutation-conjugate across p by the same psi on all residues.
"""
from fractions import Fraction


def v2(n):
    return (n & -n).bit_length() - 1


def parity_vector(a, K, p):
    x = a
    bits = []
    for _ in range(K):
        e = x % 2
        bits.append(e)
        x = x // 2 if e == 0 else (p * x + 1) // 2
    return tuple(bits)


def ulam_rows(K, p):
    mod = 1 << K
    rows = {}
    for a in range(1, mod, 2):
        t = (p * a + 1) % mod
        if t == 0:
            rows[a] = frozenset(range(1, mod, 2))
            continue
        v = v2(t)
        r = (t >> v) % (1 << (K - v))
        rows[a] = frozenset(range(r, mod, 1 << (K - v)))
    return rows          # row a is uniform on rows[a]


def tmap_rows(K, p):
    mod = 1 << K
    rows = {}
    for a in range(mod):
        t = a // 2 if a % 2 == 0 else (p * a + 1) // 2
        base = t % (mod // 2)
        rows[a] = frozenset((base, base + mod // 2))
    return rows


def main():
    for p in (5, 7, 1, -1, 9, 15, -3, 1023):
        oks = []
        for K in range(2, 14):
            mod = 1 << K
            Q3 = {parity_vector(a, K, 3): a for a in range(mod)}
            assert len(Q3) == mod                                   # Terras bijection, p = 3
            Qp = {a: parity_vector(a, K, p) for a in range(mod)}
            assert len(set(Qp.values())) == mod                     # Terras bijection, p
            psi = {a: Q3[Qp[a]] for a in range(mod)}
            # tree automorphism: a = b mod 2^m  <=>  psi(a) = psi(b) mod 2^m  (check via images of classes)
            for m in range(1, K + 1):
                img = {}
                for a in range(mod):
                    img.setdefault(a % (1 << m), set()).add(psi[a] % (1 << m))
                assert all(len(s) == 1 for s in img.values())
            R3, Rp = ulam_rows(K, 3), ulam_rows(K, p)
            ok_ulam = all(frozenset(psi[b] for b in Rp[a]) == R3[psi[a]] for a in range(1, mod, 2))
            T3, Tp = tmap_rows(K, 3), tmap_rows(K, p)
            ok_t = all(frozenset(psi[b] for b in Tp[a]) == T3[psi[a]] for a in range(mod))
            assert ok_ulam and ok_t, (p, K)
            oks.append(K)
        print(f'p = {p:5d}: psi = Q_3^-1 o Q_p is a tree automorphism of Z/2^K and conjugates the Syracuse Ulam matrix '
              f'P_K^(p) onto P_K^(3) and the T-map de Bruijn matrix onto that of p = 3, exactly, for K = {oks[0]}..{oks[-1]}')
    print('A4 ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
