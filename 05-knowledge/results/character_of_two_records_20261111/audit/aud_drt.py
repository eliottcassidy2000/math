"""AUDIT: doubly regular tournaments on record (orders 15, 19, 23): Theorem E, exact path counts, residues.

usage: python3 aud_drt.py /path/to/aud_ham REPO_ROOT [all23]
"""
import sys, os, re
from math import comb
from aud_common import (paley, from_upper_code, is_doubly_regular, ham_exact, ham_python, anti_automorphisms,
                        fold, paths_ending_in, check_tournament, cd_sign, cd_tournament, double_tournament)

HAM = sys.argv[1]
ROOT = sys.argv[2]
ALL23 = len(sys.argv) > 3 and sys.argv[3] == "all23"
LANE = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")


def double_matrix_tournament(T):
    """The same double from the matrix [[M, M], [-M^T, M^T]], M = I + S the normalised skew Hadamard matrix."""
    n = len(T)
    N = n + 1
    M = [[0] * N for _ in range(N)]
    for i in range(N):
        M[i][i] = 1
    for j in range(1, N):
        M[0][j] = 1
        M[j][0] = -1
    for a in range(n):
        for b in range(n):
            if a != b:
                M[a + 1][b + 1] = 1 if T[a][b] else -1
    D = [[0] * (2 * N) for _ in range(2 * N)]
    for i in range(N):
        for j in range(N):
            D[i][j] = M[i][j]
            D[i][N + j] = M[i][j]
            D[N + i][j] = -M[j][i]
            D[N + i][N + j] = M[j][i]
    # Hadamard and skew
    for i in range(2 * N):
        for j in range(2 * N):
            dot = sum(D[i][x] * D[j][x] for x in range(2 * N))
            assert dot == (2 * N if i == j else 0)
            if i != j:
                assert D[i][j] == -D[j][i]
    # row 0 is all +1: vertex 0 is a source of the tournament of D - I; delete it
    assert all(D[0][j] == 1 for j in range(2 * N))
    A = [[1 if (i != j and D[i][j] == 1) else 0 for j in range(1, 2 * N)] for i in range(1, 2 * N)]
    check_tournament(A)
    return A


def report_sc(name, A):
    """Anti-involutions, fold, Theorem E."""
    n = len(A)
    t = (n - 3) // 4
    antis = anti_automorphisms(A, limit=2000)
    invs = [p for p in antis if all(p[p[i]] == i for i in range(n))]
    print("%s: order %d, doubly regular: %s, anti-automorphisms found: %d, anti-involutions: %d" % (
        name, n, is_doubly_regular(A), len(antis), len(invs)))
    out = None
    for psi in invs[:3]:
        Z, D1, plus, minus, v0 = fold(A, psi)
        m = len(Z)
        zp = paths_ending_in(Z, D1)
        cyc = sum(1 for i in range(m) if A[plus[i]][minus[i]])  # triple {v0, a+, a-} cyclic iff a+ -> a-
        # Theorem E count: #{x : a -> x -> psi(a)} for a = a^+
        Ecounts = sorted(set(sum(1 for x in range(n) if A[plus[i]][x] and A[x][minus[i]]) for i in range(m)))
        pred = [comb(m, k) % 2 if t % 2 == 0 else 0 for k in range(1, m + 1)]
        print("   fixed vertex %d: invariant triples cyclic: %d of %d (Theorem E predicts %s); |E| values %s; "
              "paths of Z ending in D_1 mod 2 by k: %s; predicted C(m,k)[t even]: %s; fold formula H mod 4 = %d" % (
                  v0, cyc, m, "all" if t % 2 == 0 else "none", Ecounts,
                  "".join(str(z % 2) for z in zp[1:]), "".join(map(str, pred)), (1 + 2 * sum(zp)) % 4))
        out = (1 + 2 * sum(zp)) % 4
    return len(antis), out


print("== order 15: Cayley-Dickson sign tournament CD_4 (own construction)")
for k in (2, 3, 4):
    A = cd_tournament(k)
    n = len(A)
    H = ham_exact(A, HAM)
    print("CD_%d: order %d, doubly regular: %s, H = %d, H mod 8 = %d, n(n-1)/2 mod 8 = %d" % (
        k, n, is_doubly_regular(A), H, H % 8, (n * (n - 1) // 2) % 8))
    na, _ = report_sc("CD_%d" % k, A)
A7 = cd_tournament(3)
print("CD_3 isomorphic to P_7 (both doubly regular of order 7, unique class): H(CD_3) = %d, H(P_7) = %d" % (
    ham_python(A7), ham_python(paley(7))))
# double of P_7 against CD_4
D15 = double_tournament(paley(7))
D15m = double_matrix_tournament(paley(7))
print("double of P_7: tournament form equals matrix form: %s; doubly regular: %s; H = %d (CD_4: 198335025)" % (
    D15 == D15m, is_doubly_regular(D15), ham_exact(D15, HAM)))
for q in (3, 11):
    print("double of P_%d: tournament form equals matrix form: %s" % (q, double_tournament(paley(q)) == double_matrix_tournament(paley(q))))

print("== order 19")
code19 = None
for line in open(os.path.join(LANE, "drt19_codes.txt")):
    line = line.strip()
    if line and not line.startswith("#"):
        code19 = line
D19 = from_upper_code(code19)
H19 = ham_exact(D19, HAM)
print("D_19: doubly regular: %s; H = %d (record 1172732961771: %s); H mod 4 = %d, mod 8 = %d, mod 16 = %d; n(n-1)/2 mod 8 = %d" % (
    is_doubly_regular(D19), H19, H19 == 1172732961771, H19 % 4, H19 % 8, H19 % 16, (19 * 9) % 8))
report_sc("D_19", D19)
P19 = paley(19)
report_sc("P_19", P19)
# D_19 is not Paley: compare a simple invariant (number of transitive 5-subsets is heavy; use H)
print("D_19 not isomorphic to P_19: H differs: %s" % (H19 != 1172695746915))

print("== order 23")
recorded = {}
for line in open(os.path.join(LANE, "drt23_classes.out")):
    parts = line.split()
    if len(parts) >= 3 and parts[2].isdigit():
        recorded[parts[1]] = int(parts[2])
print("recorded H (drt23_classes.out): %d classes; residues mod 8: %s; mod 16: %s; mod 32: %s" % (
    len(recorded), sorted(set(h % 8 for h in recorded.values())), sorted(set(h % 16 for h in recorded.values())),
    sorted(set(h % 32 for h in recorded.values()))))
rec_vals = sorted(recorded.values())
codes = [l.strip() for l in open(os.path.join(LANE, "drt23_codes_all.txt")) if l.strip() and not l.startswith("#")]
print("codes in drt23_codes_all.txt: %d" % len(codes))
todo = codes if ALL23 else codes[:3]
seen = {}
for i, c in enumerate(todo):
    A = from_upper_code(c)
    H = ham_exact(A, HAM)
    na = len(anti_automorphisms(A, limit=1))
    seen.setdefault(H, 0)
    seen[H] += 1
    print("  code %2d: doubly regular: %s; H = %d; in recorded list: %s; H mod 8 = %d; mod 16 = %d; self-converse: %s" % (
        i, is_doubly_regular(A), H, H in rec_vals, H % 8, H % 16, na > 0))
print("distinct H among recomputed codes: %d" % len(seen))
# the doubled Paley class
D23 = double_tournament(paley(11))
H = ham_exact(D23, HAM)
print("double of P_11: doubly regular: %s; H = %d (record 15759752857737405: %s); H mod 8 = %d; mod 16 = %d; self-converse: %s" % (
    is_doubly_regular(D23), H, H == 15759752857737405, H % 8, H % 16, len(anti_automorphisms(D23, limit=1)) > 0))
P23 = paley(23)
print("P_23: self-converse: %s" % (len(anti_automorphisms(P23, limit=1)) > 0))
