"""Audit of THM-4761 (independent code).  Usage: python chaos_check.py > chaos_check.out
(1) E kappa^2 three ways: (a) brute force over all tournaments (ekappa2_brute.out, n = 3..7);
    (b) a permutation-pair sum: E kappa^2 = (1/n!) sum_{tau in S_n, no step a -> a-1} 2^{#steps a -> a+1}  (n <= 10);
    (c) THM-4761 (2) literally: sum over multisets lambda of even parts of 2^p (n)_(L+p) / prod m! / ((n)_L)^2.
(2) exact Var kappa, n Var kappa, the n^-3 coefficient and the score share, exact rationals.
(3) the sign in N_L = sum_F 2^p(F) chi_F: compare with exact N_L (forest_audit.exe) on small tournaments.
(4) N_2 = n(n-1) - ||S 1||^2 on random tournaments."""
import itertools, math, random, sys, os
from fractions import Fraction as Fr
from functools import lru_cache
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn

HERE = os.path.dirname(os.path.abspath(__file__))


def ff(n, k):
    r = 1
    for i in range(k):
        r *= (n - i)
    return r


def even_multisets(L):
    """multisets of positive even parts summing to L, as dicts part -> multiplicity."""
    def parts(rem, maxp):
        if rem == 0:
            yield []
            return
        for p in range(min(rem, maxp), 1, -1):
            if p % 2:
                continue
            for rest in parts(rem - p, p):
                yield [p] + rest
    for lam in parts(L, L):
        m = {}
        for x in lam:
            m[x] = m.get(x, 0) + 1
        yield m


def ENL2_literal(n, L):
    tot = 0
    for m in even_multisets(L):
        p = sum(m.values())
        if L + p > n:
            continue
        den = 1
        for v in m.values():
            den *= math.factorial(v)
        num = 2 ** p * ff(n, L + p)
        assert num % den == 0
        tot += num // den
    return tot


def ENL2_closed(n, L):
    """my closed form: sum_p 2^p (n)_(L+p) C(L/2-1, p-1) / p!   (exponential formula)."""
    if L == 0:
        return 1
    k = L // 2
    tot = Fr(0)
    for p in range(1, k + 1):
        if L + p > n:
            break
        tot += Fr(2 ** p * ff(n, L + p) * math.comb(k - 1, p - 1), math.factorial(p))
    assert tot.denominator == 1
    return int(tot)


def Ekappa2_formula(n, literal=True, Lmax=None):
    tot = Fr(1)
    for L in range(2, n, 2):
        if Lmax is not None and L > Lmax:
            break
        e = ENL2_literal(n, L) if literal else ENL2_closed(n, L)
        tot += Fr(e, ff(n, L) ** 2)
    return tot


def Ekappa2_pairs(n):
    """(1/n!) sum over tau with no fall-by-one of 2^(rises by one); tau compared with the identity path."""
    tot = 0
    for tau in itertools.permutations(range(n)):
        ok = True; s = 0
        for a, b in zip(tau, tau[1:]):
            if b == a - 1:
                ok = False; break
            if b == a + 1:
                s += 1
        if ok:
            tot += 1 << s
    return Fr(tot, math.factorial(n))


def main():
    print('== (1) E kappa^2: brute force over all tournaments vs permutation pairs vs THM-4761 (2) ==')
    brute = {}
    for line in open(os.path.join(HERE, 'ekappa2_brute.out')):
        t = line.split()
        n = int(t[0][2:]); E = n * (n - 1) // 2
        sH = int(t[2].split('=')[1]); sH2 = int(t[3].split('=')[1])
        EH = Fr(sH, 2 ** E); EH2 = Fr(sH2, 2 ** E)
        k = Fr(2 ** (n - 1), math.factorial(n))
        brute[n] = (EH * k, EH2 * k * k)
    for n in range(3, 11):
        lit = Ekappa2_formula(n, literal=True)
        clo = Ekappa2_formula(n, literal=False)
        pairs = Ekappa2_pairs(n) if n <= 10 else None
        b = brute.get(n)
        bs = f'brute E kappa = {b[0]}, E kappa^2 = {b[1]}' if b else 'brute: -'
        print(f'n={n}: formula(literal) {lit} = {float(lit):.6f}; closed form equal {clo == lit}; '
              f'pairs {pairs} equal {pairs == lit}; {bs}' + (f'; equal {b[1] == lit}' if b else ''))
    print()
    print('== (2) exact Var kappa, n Var, n^3 (Var - 2/n), score share, n^2 (share - 1 + 1/n), E(kappa-2+W/n(n-1))^2 ==')
    for n in (8, 12, 20, 30, 50, 100, 200, 400):
        V = Ekappa2_formula(n, literal=False) - 1
        V2 = Fr(2 * (n - 2), n * (n - 1))
        share = V2 / V
        print(f'n={n}: n Var = {float(n * V):.6f}; n^3(Var - 2/n) = {float(n ** 3 * (V - Fr(2, n))):.5f}; '
              f'share = {float(share):.6f}; n^2(share - 1 + 1/n) = {float(n * n * (share - 1 + Fr(1, n))):.5f}; '
              f'n^2 E(kappa - 1 - kappa_2)^2 = {float(n * n * (V - V2)):.5f}')
    # large n with truncation L <= 24 (omitted terms are O(n^-12))
    for n in (10 ** 3, 10 ** 4, 10 ** 5):
        V = Ekappa2_formula(n, literal=False, Lmax=24) - 1
        V2 = Fr(2 * (n - 2), n * (n - 1))
        share = V2 / V
        c3 = n ** 3 * (V - Fr(2, n))
        print(f'n={n} (L <= 24): n^3(Var - 2/n) = {float(c3):.6f} (target -14/3 = {-14/3:.6f}); '
              f'n*(n^3(Var-2/n) + 14/3) = {float(n * (c3 + Fr(14, 3))):.4f}; '
              f'n^2(share - 1 + 1/n) = {float(n * n * (share - 1 + Fr(1, n))):.6f} (target 4/3); '
              f'n^2 E(kappa-1-kappa_2)^2 = {float(n * n * (V - V2)):.6f} (target 2)')
    # symbolic expansion of the four shapes
    import sympy as sp
    x = sp.symbols('x')  # x = 1/n
    nn = 1 / x
    def ffs(k):
        r = 1
        for i in range(k):
            r *= (nn - i)
        return r
    terms = {'(2)': 2 * ffs(3) / ffs(2) ** 2, '(2,2)': 2 * ffs(6) / ffs(4) ** 2, '(4)': 2 * ffs(5) / ffs(4) ** 2,
             '(2,2,2)': sp.Rational(8, 6) * ffs(9) / ffs(6) ** 2, '(4,2)': 4 * ffs(8) / ffs(6) ** 2}
    tot = 0
    for k, v in terms.items():
        ser = sp.series(sp.simplify(v), x, 0, 5).removeO()
        print(f'shape {k}: {sp.expand(ser)}')
        tot += ser
    print('sum of shapes through x^4 (x = 1/n):', sp.expand(sp.series(tot, x, 0, 4).removeO()))
    V2s = 2 * (nn - 2) / (nn * (nn - 1))
    print('share series:', sp.series(sp.simplify(V2s / tot), x, 0, 3))
    print()
    print('== (3) sign check of N_L = sum_F 2^p chi_F (chi_F = prod over undirected edges ij, i<j, of s_ij) ==')
    rng = random.Random(7)
    for n in (3, 4, 5, 6):
        # all undirected linear forests with all paths of even length: enumerate directed paths via vertex sequences
        def forests(n):
            res = []
            verts = list(range(n))
            # enumerate sets of vertex-disjoint paths (as canonical vertex sequences with seq[0] < seq[-1])
            paths = []
            for k in range(3, n + 1, 2):        # k vertices -> k-1 edges even
                for seq in itertools.permutations(verts, k):
                    if seq[0] < seq[-1]:
                        paths.append(seq)
            def rec(start, used, cur):
                res.append(list(cur))
                for i in range(start, len(paths)):
                    s = set(paths[i])
                    if s & used:
                        continue
                    cur.append(paths[i]); rec(i + 1, used | s, cur); cur.pop()
            rec(0, set(), [])
            return res
        F = forests(n)
        for trial in range(4):
            code = ''.join(rng.choice('01') for _ in range(n * (n - 1) // 2))
            if trial == 0:
                code = '1' * (n * (n - 1) // 2)       # transitive
            A = tourn.from_code(code)
            S = A - A.T
            (H, c, _), = tourn.forest_numbers([code], n)
            N = tourn.N_from_c(c, n)
            uns = [0] * n; sgn = [0] * n
            for fr in F:
                L = sum(len(p) - 1 for p in fr)
                chi = 1; eps = 1
                for p in fr:
                    for a, b in zip(p, p[1:]):
                        i, j = min(a, b), max(a, b)
                        chi *= S[i, j]
                        if a > b:
                            eps = -eps
                uns[L] += 2 ** len(fr) * chi
                sgn[L] += 2 ** len(fr) * eps * chi
            print(f'n={n} code {code}: exact N_L {N}; unsigned sum_F 2^p chi_F {uns}; with descent sign {sgn}; '
                  f'unsigned ok {uns == N}; signed ok {sgn == N}')
    print()
    print('== (4) N_2 = n(n-1) - ||S1||^2 on random tournaments ==')
    for n in (5, 8, 11, 13):
        bad = 0
        for t in range(20):
            code = ''.join(rng.choice('01') for _ in range(n * (n - 1) // 2))
            A = tourn.from_code(code); S = A - A.T
            (H, c, _), = tourn.forest_numbers([code], n)
            N = tourn.N_from_c(c, n)
            W = int(np.sum((S.sum(1)) ** 2))
            if N[2] != n * (n - 1) - W:
                bad += 1
        print(f'n={n}: 20 random tournaments, violations {bad}')


if __name__ == '__main__':
    main()
