"""AUDIT: Theorem D (section 5.1; THM-5137 (D)), recomputed independently.

usage: python3 aud_theoremD.py /path/to/aud_cube [big]
 - brute-force cube sums G_35 (C program aud_cube) for P_3, P_7, P_11, P_19, P_23, P_27, CD_3, CD_4 (and P_31, CD_5
   with the argument 'big');
 - rank of A over F_2, kernel of A^T, weights mod 4, solutions of A^T u = 1, and the value given by the kernel
   formula of THM-5132 (3) for m = 3 (mod 4): 2^((m+1)/2) (n_00 + n_10);
 - Theorem D (1): parity of c_5 on regular tournaments (all circulants of orders 5..13, random regular ones),
   and THM-5132 (1) on random tournaments.
"""
import sys, subprocess, random, itertools
from aud_common import paley, is_doubly_regular, check_tournament, cd_tournament

CUBE = sys.argv[1]
BIG = len(sys.argv) > 2 and sys.argv[2] == "big"


def cube(A):
    n = len(A)
    text = str(n) + "\n" + "\n".join("".join(str(x) for x in row) for row in A) + "\n"
    return subprocess.run([CUBE], input=text, capture_output=True, text=True, check=True).stdout.strip().split("\n")


def gf2_kernel(rows, ncols):
    """rows: list of ints (bit j = column j). Returns (rank, basis of the right kernel {x : M x = 0})."""
    rows = rows[:]
    piv = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(rows)):
            if rows[i] >> c & 1:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> c & 1:
                rows[i] ^= rows[r]
        piv.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        x = 1 << f
        for i, c in enumerate(piv):
            if rows[i] >> f & 1:
                x |= 1 << c
        basis.append(x)
    return r, basis


def solve_gf2(rows, rhs, ncols):
    """One solution of M x = rhs over F_2 or None. rows: ints; rhs: list of bits."""
    aug = [rows[i] | (rhs[i] << ncols) for i in range(len(rows))]
    piv = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(aug)):
            if aug[i] >> c & 1:
                p = i
                break
        if p is None:
            continue
        aug[r], aug[p] = aug[p], aug[r]
        for i in range(len(aug)):
            if i != r and aug[i] >> c & 1:
                aug[i] ^= aug[r]
        piv.append(c)
        r += 1
    for i in range(r, len(aug)):
        if aug[i] >> ncols & 1:
            return None
    x = 0
    for i, c in enumerate(piv):
        if aug[i] >> ncols & 1:
            x |= 1 << c
    return x


def kernel_data(A):
    n = len(A)
    # A^T as rows: row j of A^T has bit i iff A[i][j]
    AT = [sum(1 << i for i in range(n) if A[i][j]) for j in range(n)]
    Arows = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    rank, K0 = gf2_kernel(AT, n)
    rankA, _ = gf2_kernel(Arows, n)
    assert rank == rankA
    u1 = solve_gf2(AT, [1] * n, n)
    d = len(K0)
    w0 = [0] * 4
    w1 = [0] * 4
    if d <= 22:
        # enumerate the kernel by Gray code
        x = 0
        w0[0] += 1
        if u1 is not None:
            w1[bin(u1).count("1") % 4] += 1
        for g in range(1, 1 << d):
            b = (g & -g).bit_length() - 1
            x ^= K0[b]
            w0[bin(x).count("1") % 4] += 1
            if u1 is not None:
                w1[bin(x ^ u1).count("1") % 4] += 1
    # A A^T over F_2
    AAT = [[sum(A[i][k] & A[j][k] for k in range(n)) % 2 for j in range(n)] for i in range(n)]
    isI = all(AAT[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n))
    isJ = all(AAT[i][j] == 1 for i in range(n) for j in range(n))
    return rank, d, w0, w1, u1, isI, isJ


def c5_parity_via_walks(A):
    """c_5 from tr A^5 = 5 c_5 (closed 5-walks in a tournament are 5-cycles)."""
    n = len(A)
    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    A2 = mm(A, A); A3 = mm(A2, A); A5 = mm(A3, A2)
    tr = sum(A5[i][i] for i in range(n))
    assert tr % 5 == 0
    return tr // 5


def c5_direct(A):
    n = len(A)
    c = 0
    for S in itertools.combinations(range(n), 5):
        v0 = S[0]
        for p in itertools.permutations(S[1:]):
            if A[v0][p[0]] and A[p[0]][p[1]] and A[p[1]][p[2]] and A[p[2]][p[3]] and A[p[3]][v0]:
                c += 1
    return c


print("== cube sums")
cases = [("P_3", paley(3)), ("P_7", paley(7)), ("CD_3", cd_tournament(3)), ("P_11", paley(11)), ("CD_4", cd_tournament(4)),
         ("P_19", paley(19)), ("P_23", paley(23)), ("P_27", paley(27))]
if BIG:
    cases += [("P_31", paley(31)), ("CD_5", cd_tournament(5))]
for name, A in cases:
    n = len(A)
    assert is_doubly_regular(A), name
    rank, d, w0, w1, u1, isI, isJ = kernel_data(A)
    formula = 2 ** ((n + 1) // 2) * (w0[0] + w1[0])
    pred = 2 ** n if n % 8 == 7 else 2 ** ((n + 1) // 2)
    lines = cube(A)
    print("%s (order %d = %d mod 8): rank_2 A = %d [(n+1)/2 = %d], dim ker A^T = %d, A A^T = I: %s, = J: %s" % (
        name, n, n % 8, rank, (n + 1) // 2, d, isI, isJ))
    print("   weights mod 4 of ker A^T: %s ; of {u : A^T u = 1}: %s ; kernel formula 2^((n+1)/2)(n_00 + n_10) = %d ; Theorem D predicts %d ; equal: %s" % (
        w0, w1, formula, pred, formula == pred))
    for l in lines:
        print("   " + l)

print("== kernel data only (no brute force)")
more = [("P_31", paley(31)), ("P_43", paley(43)), ("P_47", paley(47)), ("CD_5", cd_tournament(5)), ("CD_6", cd_tournament(6)), ("CD_7", cd_tournament(7))]
for name, A in more:
    n = len(A)
    dr = is_doubly_regular(A)
    rank, d, w0, w1, u1, isI, isJ = kernel_data(A)
    if d <= 22:
        formula = 2 ** ((n + 1) // 2) * (w0[0] + w1[0])
        pred = 2 ** n if n % 8 == 7 else 2 ** ((n + 1) // 2)
        txt = "weights mod 4 of ker A^T: %s ; of K_1: %s ; formula = 2^%d ; predicted 2^%d ; equal: %s" % (
            w0, w1, formula.bit_length() - 1, pred.bit_length() - 1, formula == pred)
    else:
        txt = "kernel too large to enumerate (dim %d); checking doubly-evenness on the basis and pairwise even overlaps" % d
        # doubly even code: basis weights = 0 mod 4 and pairwise intersections even
        AT = [sum(1 << i for i in range(n) if A[i][j]) for j in range(n)]
        _, K0 = gf2_kernel(AT, n)
        ok = all(bin(b).count("1") % 4 == 0 for b in K0) and all(bin(a & b).count("1") % 2 == 0 for a in K0 for b in K0)
        txt += ": doubly even: %s" % ok
    print("%s (order %d = %d mod 8): doubly regular: %s; rank_2 A = %d [(n+1)/2 = %d]; A A^T = I: %s, = J: %s; %s" % (
        name, n, n % 8, dr, rank, (n + 1) // 2, isI, isJ, txt))

print("== Theorem D (1): c_5 of regular tournaments")
bad = 0
tested = 0
for n in (5, 7, 9, 11, 13):
    m = (n - 1) // 2
    cnt = 0
    for signs in itertools.product((0, 1), repeat=m):
        S = set((d if s == 0 else n - d) for d, s in zip(range(1, m + 1), signs))
        A = [[1 if ((j - i) % n) in S else 0 for j in range(n)] for i in range(n)]
        c5 = c5_parity_via_walks(A)
        if n <= 9:
            assert c5 == c5_direct(A)
        tested += 1
        cnt += 1
        if c5 % 2:
            bad += 1
    print("order %d: %d circulant tournaments, odd c_5 so far: %d" % (n, cnt, bad))
random.seed(20261111)
for n in (5, 7, 9, 11, 13, 15, 17, 21):
    m = (n - 1) // 2
    # start from the rotational tournament and reverse random directed triangles (keeps all scores)
    A = [[1 if 1 <= (j - i) % n <= m else 0 for j in range(n)] for i in range(n)]
    odd = 0
    distinct = set()
    for step in range(300):
        for _ in range(30):
            a, b, c = random.sample(range(n), 3)
            if A[a][b] and A[b][c] and A[c][a]:
                A[a][b] = A[b][c] = A[c][a] = 0
                A[b][a] = A[c][b] = A[a][c] = 1
        check_tournament(A)
        assert all(sum(r) == m for r in A)
        c5 = c5_parity_via_walks(A)
        distinct.add(c5)
        tested += 1
        if c5 % 2:
            odd += 1
            bad += 1
    print("order %d: 300 random regular tournaments (triangle reversals), %d distinct values of c_5, odd c_5: %d" % (n, len(distinct), odd))
print("regular tournaments tested: %d, with odd c_5: %d" % (tested, bad))

print("== THM-5132 (1) on random tournaments: c_5 = C(q_1,2) + C(q_2,2) (mod 2)")
viol = 0
N = 0
for n in range(3, 16):
    for _ in range(200):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        s = [sum(r) % 2 for r in A]
        t = [sum(A[i][j] * s[j] for j in range(n)) % 2 for i in range(n)]
        q1, q2 = sum(s), sum(t)
        c5 = c5_parity_via_walks(A)
        N += 1
        if (c5 - (q1 * (q1 - 1) // 2 + q2 * (q2 - 1) // 2)) % 2:
            viol += 1
print("random tournaments of orders 3..15: %d tested, violations: %d" % (N, viol))
