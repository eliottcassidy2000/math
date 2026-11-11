"""AUDIT: section 5.3 (HYP-9667): switching classes, the proved reduction, the refuted local forms.

usage: python3 aud_switch.py /path/to/aud_ham [n24]
"""
import sys, random, subprocess, itertools
from aud_common import (paley, ham_python, check_tournament, is_doubly_regular, from_upper_code, cd_tournament,
                        double_tournament)

HAM = sys.argv[1]
N24 = len(sys.argv) > 2 and sys.argv[2] == "n24"
random.seed(9667)


def ham_mod32(A):
    """H mod 2^32 through the audit C program (one run)."""
    n = len(A)
    text = str(n) + "\n" + "\n".join("".join(str(x) for x in row) for row in A) + "\n"
    return int(subprocess.run([HAM, "0"], input=text, capture_output=True, text=True, check=True).stdout)


def H(A):
    return ham_python(A) if len(A) <= 11 else ham_mod32(A)


def with_source(T):
    n = len(T)
    U = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(1, n + 1):
        U[0][j] = 1
    for a in range(n):
        for b in range(n):
            U[a + 1][b + 1] = T[a][b]
    return U


def switch(U, X):
    """Reverse all arcs between X and its complement."""
    n = len(U)
    inX = [False] * n
    for x in X:
        inX[x] = True
    V = [row[:] for row in U]
    for a in range(n):
        for b in range(n):
            if a != b and inX[a] != inX[b]:
                V[a][b] = 1 - U[a][b]
    return V


def skew(U):
    n = len(U)
    return [[(1 if U[a][b] else -1) if a != b else 0 for b in range(n)] for a in range(n)]


def is_skew_hadamard(U):
    n = len(U)
    S = skew(U)
    for a in range(n):
        for b in range(n):
            d = sum(S[a][x] * S[b][x] for x in range(n))
            if d != ((n - 1) if a == b else 0):
                return False
    return True


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def even_path_endpoint_parity(R):
    """g_R(x) = parity of the number of directed paths of R with an even number of vertices having x as first or last vertex;
    also g'(x) = the same for paths on |R| - 1 vertices only."""
    n = len(R)
    out = [sum(1 << j for j in range(n) if R[i][j]) for i in range(n)]
    full = (1 << n) - 1
    g = [0] * n
    g1 = [0] * n
    # dp[(mask, last)] with the first vertex tracked: do one DP per start vertex
    for s in range(n):
        layer = {(1 << s, s): 1}
        size = 1
        while layer:
            nxt = {}
            for (mask, u), x in layer.items():
                t = out[u] & ~mask & full
                while t:
                    w = (t & -t).bit_length() - 1
                    t &= t - 1
                    key = (mask | (1 << w), w)
                    nxt[key] = nxt.get(key, 0) + x
            size += 1
            layer = nxt
            if size % 2 == 0:
                for (mask, u), x in layer.items():
                    g[s] += x
                    g[u] += x
                    if size == n - 1:
                        g1[s] += x
                        g1[u] += x
    return [v % 2 for v in g], [v % 2 for v in g1]


print("== (ii) switching classes of (source) + T: H mod 4 and mod 8 over members")
for name, T, mode in (("P_3", paley(3), "all"), ("P_7", paley(7), "all"), ("P_11", paley(11), 200),
                      ("CD_4", cd_tournament(4), 100), ("P_19", paley(19), 40)):
    U = with_source(T)
    N = len(U)
    assert is_skew_hadamard(U)
    res4, res8 = {}, {}
    members = set()
    if mode == "all":
        subsets = [[v for v in range(1, N) if (m >> (v - 1)) & 1] for m in range(1 << (N - 1))]
    else:
        subsets = [[v for v in range(1, N) if random.random() < 0.5] for _ in range(mode)]
    for X in subsets:
        V = switch(U, X)
        assert is_skew_hadamard(V) or N > 12  # (checked for the small ones; conjugation by a sign matrix)
        members.add(tuple(tuple(r) for r in V))
        h = H(V)
        res4[h % 4] = res4.get(h % 4, 0) + 1
        res8[h % 8] = res8.get(h % 8, 0) + 1
    print("N = %2d, class of (source) + %s: %d switching sets (%s), %d distinct labelled members; H mod 4: %s; H mod 8: %s; 1 + N/2 mod 4 = %d" % (
        N, name, len(subsets), mode if mode == "all" else "random", len(members), dict(sorted(res4.items())),
        dict(sorted(res8.items())), (1 + N // 2) % 4))

lane = sys.argv[3] if len(sys.argv) > 3 else None
if N24:
    import os
    root = "/Users/e/Documents/GitHub/math-wt-atoms/04-computation/experiments/sierpinski_ocf_20261022"
    code19 = [l.strip() for l in open(os.path.join(root, "drt19_codes.txt")) if l.strip() and not l.startswith("#")][0]
    codes23 = [l.strip() for l in open(os.path.join(root, "drt23_codes_all.txt")) if l.strip() and not l.startswith("#")]
    for name, T, cnt in (("D_19", from_upper_code(code19), 40), ("P_23", paley(23), 12), ("double of P_11", double_tournament(paley(11)), 12),
                         ("annealed D_23 code 0", from_upper_code(codes23[0]), 8), ("annealed D_23 code 2", from_upper_code(codes23[2]), 8)):
        U = with_source(T)
        N = len(U)
        res4, res8 = {}, {}
        for _ in range(cnt):
            X = [v for v in range(1, N) if random.random() < 0.5]
            h = ham_mod32(switch(U, X))
            res4[h % 4] = res4.get(h % 4, 0) + 1
            res8[h % 8] = res8.get(h % 8, 0) + 1
        print("N = %2d, class of (source) + %s: %d random members; H mod 4: %s; H mod 8: %s; 1 + N/2 mod 4 = %d" % (
            N, name, cnt, dict(sorted(res4.items())), dict(sorted(res8.items())), (1 + N // 2) % 4))

print("== H mod 4 is not a switching invariant in general")
for n in range(4, 11):
    diff = 0
    trials = 200
    for _ in range(trials):
        A = random_tournament(n)
        X = [v for v in range(n) if random.random() < 0.5]
        if ham_python(A) % 4 != ham_python(switch(A, X)) % 4:
            diff += 1
    print("order %d: H mod 4 differs for %d of %d random (tournament, switching set) pairs" % (n, diff, trials))

print("== the proved reduction: (H(U) - H(U^v))/2 = sum_{x in N^-(v)} g_R(x) (mod 2); for |R| odd g_R = paths on |R|-1 vertices")
bad1 = bad2 = tested = odd_tested = 0
for n in range(3, 10):
    for _ in range(120):
        U = random_tournament(n)
        v = random.randrange(n)
        Uv = switch(U, [v])
        idx = [x for x in range(n) if x != v]
        R = [[U[a][b] for b in idx] for a in idx]
        g, g1 = even_path_endpoint_parity(R)
        lhs = ((ham_python(U) - ham_python(Uv)) // 2) % 2
        rhs = sum(g[i] for i, x in enumerate(idx) if U[x][v]) % 2
        tested += 1
        if lhs != rhs:
            bad1 += 1
        if len(R) % 2 == 1:
            odd_tested += 1
            if g != g1:
                bad2 += 1
print("random tournaments of orders 3..9: %d tested, violations of the reduction: %d; |R| odd: %d tested, g_R != (paths on |R|-1 vertices) in %d" % (
    tested, bad1, odd_tested, bad2))

print("== refuted local form 1: det S = 4H - 3 (mod 16)")
def det_int(M):
    from fractions import Fraction
    n = len(M)
    M = [[Fraction(x) for x in row] for row in M]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return int(d)

for n in (2, 4, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    viol = 0
    tot = 0
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        tot += 1
        if (det_int(skew(A)) - (4 * ham_python(A) - 3)) % 16:
            viol += 1
    print("order %d: all %d labelled tournaments, violations: %d" % (n, tot, viol))

print("== refuted local form 2: one orthogonal row (exhaustive up to isomorphism)")
print("   U has row v of S orthogonal to all other rows iff U is a switching of (source) + R', R' regular; U^v is then (source) + R' switched at the complementary set.")
def regular_classes(n):
    d = (n - 1) // 2
    outp = subprocess.run(["gentourng", "-q", "-d%d" % d, "-D%d" % d, str(n)], capture_output=True, text=True, check=True).stdout.split()
    return [from_upper_code(c) for c in outp]

for n in (3, 5, 7, 9):
    classes = regular_classes(n)
    N = n + 1
    differ = 0
    total = 0
    cls_with_fail = 0
    for T in classes:
        assert all(sum(r) == (n - 1) // 2 for r in T)
        U = with_source(T)
        hs = {}
        fail_here = 0
        for m in range(1 << n):
            X = [v for v in range(1, N) if (m >> (v - 1)) & 1]
            hs[m] = ham_python(switch(U, X)) % 4
        for m in range(1 << n):
            comp = ((1 << n) - 1) ^ m
            if m < comp:
                total += 1
                if hs[m] != hs[comp]:
                    differ += 1
                    fail_here += 1
        if fail_here:
            cls_with_fail += 1
    print("N = %2d: %d regular tournaments of order %d up to isomorphism; pairs {U, U^v}: %d; H mod 4 differs in %d (classes with a failure: %d)" % (
        N, len(classes), n, total, differ, cls_with_fail))

print("== refuted local form 3: score parities: <s_R, g_R> = 0 (mod 2) for every tournament R")
for n in (3, 4, 5, 6, 7):
    outp = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True, check=True).stdout.split()
    viol = 0
    for c in outp:
        R = from_upper_code(c)
        s = [sum(r) % 2 for r in R]
        g, _ = even_path_endpoint_parity(R)
        if sum(a * b for a, b in zip(s, g)) % 2:
            viol += 1
    print("order %d: %d isomorphism classes, <s, g> odd in %d" % (n, len(outp), viol))

print("== the double of an arbitrary base tournament: residues of H mod 4")
for n in (1, 2, 3, 4, 5):
    if n == 1:
        bases = [[[0]]]
    else:
        outp = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True, check=True).stdout.split()
        bases = [from_upper_code(c) for c in outp]
    res = {}
    for T in bases:
        for TT in (T, [[T[j][i] for j in range(n)] for i in range(n)]):
            h = ham_python(double_tournament(TT))
            res[h % 4] = res.get(h % 4, 0) + 1
    print("base order %d (%d classes and their converses): H(double) mod 4: %s" % (n, len(bases), dict(sorted(res.items()))))
