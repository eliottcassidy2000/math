"""ker_verify.py - checks for the two-port kernel statements of S23 (independent of kernels.c: direct recounts).

V1. H(K[x <- A, y <- B]) = sum_{i,j} N_ij P_i(A) P_j(B) for random kernels from ker_m.txt and random modules.
V2. the strong 4-tournament (code 34, ports 1, 3):   H = Pi_2(A) Pi_2(B) + H(A) H(B).
V3. order 5, code 209, ports (0, 3):                 H = 3 Pi_2(A) Pi_2(B) - H(A) H(B).
V4. order 8, code 17082497, ports (0, 7), transitive modules: H = -1 - 16 2^a + 62 2^(a+b)  (a at port 0).
V5. factorisation N22 H = L_A L'_B + det * H(A) H(B)."""
import os, random, sys
from math import comb
here = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")


def from_code(m, code):
    A = [[0] * m for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def hs_table(A):
    n = len(A)
    full = (1 << n) - 1
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    cnt = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        cnt[1 << v][v] = 1
    hs = [0] * (full + 1)
    hs[0] = 1
    for s in range(1, full + 1):
        row = cnt[s]
        t = 0
        for v in range(n):
            c = row[v]
            if c:
                t += c
                nx = out[v] & ~s
                while nx:
                    u = (nx & -nx).bit_length() - 1
                    nx &= nx - 1
                    cnt[s | 1 << u][u] += c
        hs[s] = t
    return hs


def H(A):
    return hs_table(A)[(1 << len(A)) - 1] if A else 1


def P_values(A, K):
    """P_1..P_K: ordered i-tuples of disjoint NONEMPTY paths covering A."""
    n = len(A)
    full = (1 << n) - 1
    hs = hs_table(A)
    f = [[0] * (full + 1) for _ in range(K + 1)]
    f[0][0] = 1
    for k in range(1, K + 1):
        for S in range(full + 1):
            t = 0
            sub = S
            while True:
                t += hs[sub] * f[k - 1][S ^ sub]
                if sub == 0:
                    break
                sub = (sub - 1) & S
            f[k][S] = t
    F = [f[k][full] for k in range(K + 1)]
    return [None] + [sum((-1) ** (i - j) * comb(i, j) * F[j] for j in range(i + 1)) for i in range(1, K + 1)]


def substitute(K, x, y, A, B):
    m = len(K)
    verts = []
    for v in range(m):
        if v == x:
            verts += [("A", i) for i in range(len(A))]
        elif v == y:
            verts += [("B", i) for i in range(len(B))]
        else:
            verts.append(("K", v))
    N = len(verts)
    T = [[0] * N for _ in range(N)]

    def kv(t):
        return x if t[0] == "A" else (y if t[0] == "B" else t[1])
    for p in range(N):
        for q in range(N):
            if p == q:
                continue
            a, b = verts[p], verts[q]
            if a[0] == b[0] == "A":
                T[p][q] = A[a[1]][b[1]]
            elif a[0] == b[0] == "B":
                T[p][q] = B[a[1]][b[1]]
            else:
                T[p][q] = K[kv(a)][kv(b)]
    return T


def rand_t(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def TT(n):
    return [[1 if i < j else 0 for j in range(n)] for i in range(n)]


def read_ker(m):
    rows = []
    for l in open(os.path.join(here, "ker_%d.txt" % m)):
        parts = l.split("|")
        head = parts[0].split()
        mm, code, x, y, maxi, maxj, n11, n12, n21, n22, det, two = map(int, head)
        if len(parts) > 1:
            mat = [[int(t) for t in r.split()] for r in parts[1].split(";") if r.strip()]
        else:
            mat = [[n11, n12], [n21, n22]]
        rows.append((code, x, y, mat))
    return rows


if __name__ == "__main__":
    rng = random.Random(2311)
    bad = tot = 0
    for m in (3, 4, 5, 6):
        rows = read_ker(m)
        for _ in range(60 if m < 6 else 40):
            code, x, y, mat = rng.choice(rows)
            K = from_code(m, code)
            A = rand_t(rng.randint(1, 4), rng)
            B = rand_t(rng.randint(1, 4), rng)
            kk = max(len(mat), max(len(r) for r in mat))
            PA = P_values(A, kk)
            PB = P_values(B, kk)
            pred = sum(mat[i][j] * PA[i + 1] * PB[j + 1] for i in range(len(mat)) for j in range(len(mat[i])))
            tot += 1
            bad += pred != H(substitute(K, x, y, A, B))
    print("V1. visit-matrix formula: cases %d, failures %d" % (tot, bad))
    bad = tot = 0
    K4 = from_code(4, 34)
    K5 = from_code(5, 209)
    for _ in range(60):
        A = rand_t(rng.randint(1, 5), rng)
        B = rand_t(rng.randint(1, 5), rng)
        hA, hB = H(A), H(B)
        qA, qB = P_values(A, 2)[2] + 2 * hA, P_values(B, 2)[2] + 2 * hB       # Pi_2 = P_2 + 2 H
        tot += 2
        bad += H(substitute(K4, 1, 3, A, B)) != qA * qB + hA * hB
        bad += H(substitute(K5, 0, 3, A, B)) != 3 * qA * qB - hA * hB
    print("V2/V3. Gaussian kernel (H = Pi Pi' + H H') and first hyperbolic kernel (H = 3 Pi Pi' - H H'): cases %d, failures %d" % (tot, bad))
    print("    kernel 34 arcs:", [(i, j) for i in range(4) for j in range(4) if K4[i][j]], " kernel 209 arcs:", [(i, j) for i in range(5) for j in range(5) if K5[i][j]])
    K8 = from_code(8, 17082497)
    bad = tot = 0
    vals = []
    for a in range(1, 5):
        for b in range(1, 5):
            h = H(substitute(K8, 0, 7, TT(a), TT(b)))
            tot += 1
            ok = (h == -1 - 16 * 2 ** a + 62 * 2 ** (a + b))
            bad += not ok
            vals.append(h)
    print("V4. anti-carrier (order 8, code 17082497, ports 0, 7): H = -1 - 2^(a+4) + 31 * 2^(a+b+1): cases %d, failures %d; values %s" % (tot, bad, vals[:8]))
    for (code, form, name) in ((172278145, (1, 16, 0, 24), "d = 3"), (71647617, (1, 4, 0, 42), "d = 21")):
        K = from_code(8, code)
        res = []
        for (pa, pb) in ((0, 7), (7, 0)):
            okc = all(H(substitute(K, pa, pb, TT(a), TT(b))) == form[0] + form[1] * 2 ** a + form[2] * 2 ** b + form[3] * 2 ** (a + b) for a in range(1, 4) for b in range(1, 4))
            okd = all(H(substitute(K, pa, pb, TT(a), TT(b))) == form[0] + form[1] * 2 ** b + form[2] * 2 ** a + form[3] * 2 ** (a + b) for a in range(1, 4) for b in range(1, 4))
            res.append((pa, pb, okc, okd))
        print("    known carrier %s (code %d) form %s: %s" % (name, code, form, res))
    # V5
    bad = tot = 0
    for m in (3, 4, 5, 6):
        rows = [r for r in read_ker(m) if len(r[3]) == 2 and len(r[3][0]) == 2 and r[3][1][1] > 0]
        for _ in range(40):
            code, x, y, mat = rng.choice(rows)
            K = from_code(m, code)
            A = rand_t(rng.randint(1, 4), rng)
            B = rand_t(rng.randint(1, 4), rng)
            hA, hB = H(A), H(B)
            pA, pB = P_values(A, 2)[2], P_values(B, 2)[2]
            (n11, n12), (n21, n22) = mat
            LA = n22 * pA + n12 * hA
            LB = n22 * pB + n21 * hB
            tot += 1
            bad += n22 * H(substitute(K, x, y, A, B)) != LA * LB + (n11 * n22 - n12 * n21) * hA * hB
    print("V5. N22 H = L_A L'_B + det H(A) H(B): cases %d, failures %d" % (tot, bad))
