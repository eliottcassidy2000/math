# Audit2 (A): independent check of THM-4757 (1), the linear-forest identity
#     2^(n-1) H(T) = sum_F (n - |F|)! prod_(e in F) S_e,   i.e.  kappa = H 2^(n-1)/n! = sum_L N_L(S)/(n)_L,
# by DIRECT enumeration of all directed linear forests of K_n (Lah structures: set partitions into ordered lists,
# isolated vertices = singleton lists), exact integer arithmetic.  H by brute force over permutations.
# Also checks: N_L = 0 for odd L; N_L equals the same sum restricted to forests whose paths all have even arc
# length (reversal of one odd path is a sign-reversing involution); N_(n-1) = (-1)^(n-1) F(T,-1) with
# F(T,x) = sum_sigma x^(#forward steps) (THM-4756 convention).
# Tournaments: every labelled tournament for n <= 4, 40 random ones for n = 5..7, 12 random ones for n = 8,
# plus the regular ones R_5, R_7, P_7, and the non-circulant regular 7-tournament.
import itertools, random, sys
from fractions import Fraction
from math import factorial

def lah_structures(n):
    """yield lists of paths (lists of vertices) covering range(n); every ordered-list partition exactly once"""
    def rec(v, paths):
        if v == n:
            yield paths; return
        # new singleton path
        yield from rec(v + 1, paths + [[v]])
        # insert v into an existing path at any position
        for i, p in enumerate(paths):
            for pos in range(len(p) + 1):
                q = p[:pos] + [v] + p[pos:]
                yield from rec(v + 1, paths[:i] + [q] + paths[i + 1:])
    yield from rec(0, [])

def forest_sums(S, n):
    N = [0] * n
    Neven = [0] * n
    count = 0
    for paths in lah_structures(n):
        count += 1
        L = 0; w = 1; alleven = True
        for p in paths:
            l = len(p) - 1
            L += l
            if l % 2: alleven = False
            for a, b in zip(p, p[1:]): w *= S[a][b]
        N[L] += w
        if alleven: Neven[L] += w
    return N, Neven, count

def ham_and_F(A, n):
    H = 0; F = {}
    for perm in itertools.permutations(range(n)):
        f = sum(A[perm[i]][perm[i + 1]] for i in range(n - 1))
        F[f] = F.get(f, 0) + 1
        if f == n - 1: H += 1
    Fm1 = sum(c * (-1) ** f for f, c in F.items())
    return H, Fm1

def tour_from_bits(n, bits):
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k]: A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A

def circulant(n, conn):
    return [[1 if (j - i) % n in conn else 0 for j in range(n)] for i in range(n)]

def check(A, n, label):
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
    N, Neven, cnt = forest_sums(S, n)
    H, Fm1 = ham_and_F(A, n)
    lhs = 2 ** (n - 1) * H
    rhs = sum(N[L] * factorial(n - L) for L in range(n))
    kappa = Fraction(lhs, factorial(n))
    kap2 = sum(Fraction(N[L], factorial(n) // factorial(n - L)) for L in range(n))
    ok = (lhs == rhs) and (kappa == kap2)
    odd0 = all(N[L] == 0 for L in range(1, n, 2))
    evenok = (N == Neven)
    topok = (N[n - 1] == (-1) ** (n - 1) * Fm1)
    return ok, odd0, evenok, topok, cnt, H, N

if __name__ == '__main__':
    random.seed(20261022)
    tests = []
    for n in range(2, 5):
        m = n * (n - 1) // 2
        for bits in itertools.product([0, 1], repeat=m):
            tests.append((n, tour_from_bits(n, bits), f'all n={n}'))
    for n, reps in ((5, 40), (6, 40), (7, 40), (8, 12)):
        m = n * (n - 1) // 2
        for r in range(reps):
            tests.append((n, tour_from_bits(n, [random.randint(0, 1) for _ in range(m)]), f'random n={n}'))
    tests.append((5, circulant(5, {1, 2}), 'R_5'))
    tests.append((7, circulant(7, {1, 2, 3}), 'R_7'))
    tests.append((7, circulant(7, {1, 2, 4}), 'P_7'))
    # the non-circulant regular 7-tournament: reverse a directed 3-cycle of R_7 (keeps regularity); the three regular
    # 7-tournaments are told apart by the number of transitive 4-subtournaments (P_7: 0, R_7: 7)
    import copy
    A = circulant(7, {1, 2, 3})
    # census of TT4 distinguishes the three regular 7-tournaments (P_7: 0? R_7: 14?); print and use all distinct
    seen = {}
    for tri in itertools.combinations(range(7), 3):
        for (x, y, z) in (tri, (tri[0], tri[2], tri[1])):
            if A[x][y] and A[y][z] and A[z][x]:
                B = copy.deepcopy(A)
                B[x][y] = B[y][z] = B[z][x] = 0; B[y][x] = B[z][y] = B[x][z] = 1
                tt4 = sum(1 for Q in itertools.combinations(range(7), 4)
                          if sorted(sum(B[u][v] for v in Q) for u in Q) == [0, 1, 2, 3])
                seen.setdefault(tt4, B)
    tt4 = lambda B: sum(1 for Q in itertools.combinations(range(7), 4) if sorted(sum(B[u][v] for v in Q) for u in Q) == [0, 1, 2, 3])
    print('TT4 counts: R_7 =', tt4(circulant(7, {1, 2, 3})), ' P_7 =', tt4(circulant(7, {1, 2, 4})),
          ' triangle-reversals of R_7 give TT4 in', sorted(seen))
    for t, B in seen.items():
        if t not in (tt4(circulant(7, {1, 2, 3})), tt4(circulant(7, {1, 2, 4}))):
            tests.append((7, B, f'non-circulant regular 7 (TT4={t})'))
            print('non-circulant regular 7-tournament rows:', [''.join(map(str, r)) for r in B])
    bad = 0; summary = {}
    for n, A, label in tests:
        ok, odd0, evenok, topok, cnt, H, N = check(A, n, label)
        key = label if not label.startswith(('all', 'random')) else label
        s = summary.setdefault(key, [0, 0, 0, 0, 0, cnt])
        s[0] += 1; s[1] += ok; s[2] += odd0; s[3] += evenok; s[4] += topok
        if not (ok and odd0 and evenok and topok):
            bad += 1; print('FAIL', label, ok, odd0, evenok, topok, H, N)
        if not label.startswith(('all', 'random')):
            print(f'{label}: H = {H}, N_L = {N}, kappa = {Fraction(2 ** (n - 1) * H, factorial(n))}')
    print('\nlabel: #tournaments, identity ok, N_odd = 0, even-path restriction ok, top term ok, #linear forests of K_n')
    for k, s in summary.items(): print(f'  {k}: {s}')
    print('ALL PASSED' if bad == 0 else f'{bad} FAILURES')
