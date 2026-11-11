#!/usr/bin/env python3
"""Run G: doubly regular tournaments that are not Paley tournaments, and their path counts mod 4 and mod 8.

  drt_zoo.py controls EXACT LAYERED     H(P_7), H(P_11), H(P_19), H(CD_4) exact; H(P_23) mod 256
  drt_zoo.py g1 EXACT [stride]          order 19: Goethals-Seidel skew Hadamard matrices of order 20 (circulants of
                                        order 5, the first of skew type), all descendants, exact H
  drt_zoo.py g2 LAYERED                 order 23: the double [[M, M], [-M^T, M^T]] of the Paley matrix of order 12,
                                        its 24 descendants, H mod 256
  drt_zoo.py g3 [max_matrices [stride]] order 27: Goethals-Seidel with circulants of order 7; classes of descendants;
                                        for the self-converse ones H mod 4 by the folding formula (Theorem Q)

EXACT = path of the compiled hampaths_exact.c, LAYERED = path of the compiled hampaths_layered.c.
A descendant of a skew Hadamard matrix M (M + M^T = 2I) at index c: conjugate by the diagonal matrix that makes
row c all +1 (then column c is -1 off the diagonal), delete row and column c; a -> b iff the entry is +1."""
import sys, itertools, subprocess
from paley_tools import paley_out
from fold import fold, paths_ending_in
sys.setrecursionlimit(100000)

def circ(row):
    m = len(row); return [[row[(j - i) % m] for j in range(m)] for i in range(m)]
def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def transpose(A): return [list(r) for r in zip(*A)]
def neg(A): return [[-x for x in r] for r in A]
def block(rows):
    out = []
    for brow in rows:
        m = len(brow[0])
        for i in range(m):
            r = []
            for B in brow: r += B[i]
            out.append(r)
    return out
def is_skew_hadamard(M):
    N = len(M)
    if any(M[a][a] != 1 for a in range(N)): return False
    if any(M[a][b] != -M[b][a] for a in range(N) for b in range(a + 1, N)): return False
    return all(sum(M[a][c] * M[b][c] for c in range(N)) == 0 for a in range(N) for b in range(a + 1, N))
def goethals_seidel(a, b, c, d):
    m = len(a); A, B, C, D = circ(a), circ(b), circ(c), circ(d)
    R = [[1 if (i + j) % m == 0 else 0 for j in range(m)] for i in range(m)]
    BR, CR, DR = matmul(B, R), matmul(C, R), matmul(D, R)
    BtR, CtR, DtR = matmul(transpose(B), R), matmul(transpose(C), R), matmul(transpose(D), R)
    return block([[A, BR, CR, DR], [neg(BR), A, DtR, neg(CtR)], [neg(CR), neg(DtR), A, BtR], [neg(DR), CtR, neg(BtR), A]])
def paf(row):
    m = len(row); return tuple(sum(row[i] * row[(i + s) % m] for i in range(m)) for s in range(1, (m - 1) // 2 + 1))
def gs_solutions(m):
    """all (a, b, c, d): a of skew type (a_0 = 1, a_{m-j} = -a_j), b, c, d arbitrary +-1 rows, PAF sums zero."""
    h = (m - 1) // 2
    skew = []
    for bits in itertools.product((1, -1), repeat=h):
        a = [1] + list(bits) + [0] * h
        for j in range(1, h + 1): a[m - j] = -a[j]
        skew.append(tuple(a))
    allrows = [tuple(r) for r in itertools.product((1, -1), repeat=m)]
    pafs = {r: paf(r) for r in allrows}
    pair = {}
    for c in allrows:
        for d in allrows:
            key = tuple(x + y for x, y in zip(pafs[c], pafs[d]))
            pair.setdefault(key, []).append((c, d))
    for a in skew:
        pa = paf(a)
        for b in allrows:
            key = tuple(-(x + y) for x, y in zip(pa, pafs[b]))
            for (c, d) in pair.get(key, ()):
                yield a, b, c, d
def descendant(M, c):
    N = len(M); d = [M[c][i] if i != c else 1 for i in range(N)]
    idx = [i for i in range(N) if i != c]; n = N - 1; out = [0] * n
    for x, a in enumerate(idx):
        for y, b in enumerate(idx):
            if a != b and d[a] * M[a][b] * d[b] == 1: out[x] |= 1 << y
    return out
def is_drt(out):
    n = len(out)
    return (set(bin(o).count('1') for o in out) == {(n - 1) // 2}
            and set(bin(out[a] & out[b]).count('1') for a in range(n) for b in range(a + 1, n)) == {(n - 3) // 4})
def paley_matrix(q):
    out = paley_out(q); N = q + 1
    M = [[1] * N] + [[-1] + [1 if a == b else (1 if (out[a] >> b) & 1 else -1) for b in range(q)] for a in range(q)]
    assert is_skew_hadamard(M); return M
def double(M):
    Mt = transpose(M); K = block([[M, M], [neg(Mt), Mt]]); assert is_skew_hadamard(K); return K
def invariant(out):
    """cheap isomorphism invariant: sorted numbers of 3-cycles inside the out-neighbourhoods, and of 3-cycles inside
    the common out-neighbourhoods of the arcs."""
    n = len(out)
    def c3(mask):
        vs = [v for v in range(n) if (mask >> v) & 1]; k = len(vs)
        sc = [bin(out[v] & mask).count('1') for v in vs]
        return k * (k - 1) * (k - 2) // 6 - sum(s * (s - 1) // 2 for s in sc)
    a = tuple(sorted(c3(out[v]) for v in range(n)))
    b = tuple(sorted(c3(out[u] & out[v]) for u in range(n) for v in range(n) if (out[u] >> v) & 1))
    return a, b
def find_map(A, B, anti=False, count=False, limit=1):
    """an isomorphism A -> B (or, with anti, a map f with a -> b in A iff f(b) -> f(a) in B); count=True counts."""
    n = len(A)
    arcA = [[(A[a] >> b) & 1 for b in range(n)] for a in range(n)]
    arcB = [[(B[a] >> b) & 1 for b in range(n)] for a in range(n)]
    perm = [-1] * n; used = [False] * n; res = []; cnt = [0]
    def rec(i):
        if i == n:
            cnt[0] += 1
            if not count: res.append(perm[:])
            return (not count) or (limit is not None and cnt[0] >= limit)
        for c in range(n):
            if used[c]: continue
            ok = True
            for j in range(i):
                if anti:
                    if arcA[i][j] != arcB[perm[j]][c]: ok = False; break
                else:
                    if arcA[i][j] != arcB[c][perm[j]]: ok = False; break
            if not ok: continue
            perm[i] = c; used[c] = True
            if rec(i + 1): return True
            perm[i] = -1; used[c] = False
        return False
    rec(0)
    return cnt[0] if count else (res[0] if res else None)
def anti_involution(out):
    """an anti-involution of the tournament, or None: an odd power of any anti-automorphism."""
    f = find_map(out, out, anti=True)
    if f is None: return None
    n = len(out); order = 1; g = f[:]
    while g != list(range(n)):
        g = [f[x] for x in g]; order += 1
    b = order
    while b % 2 == 0: b //= 2
    psi = list(range(n))
    for _ in range(b): psi = [f[x] for x in psi]
    assert psi != list(range(n)) and all(psi[psi[x]] == x for x in range(n))
    return psi
def exact_H(exe, tours):
    inp = "".join("%d\n%s\n" % (len(o), "\n".join(str(x) for x in o)) for o in tours)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    return [int(l.split()[1]) for l in res if l.strip()]
def layered_H(exe, out):
    inp = "%d\n" % len(out) + "".join("%d 1\n" % o for o in out)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout
    return int(res.strip().split("\n")[-1].split(":")[1])
def classes(tours, reps=None):
    """partition into isomorphism classes; returns list of (representative, [indices])."""
    reps = [] if reps is None else reps
    for i, o in enumerate(tours):
        inv = invariant(o)
        for r in reps:
            if r[0] == inv and find_map(o, r[1]) is not None: r[2].append(i); break
        else: reps.append((inv, o, [i]))
    return reps

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "controls":
        exact, layered = sys.argv[2], sys.argv[3]
        from cayley_dickson import tournament
        ts = [paley_out(7), paley_out(11), paley_out(19), tournament(4)[0]]
        hs = exact_H(exact, ts)
        want = [189, 95095, 1172695746915, 198335025]
        for name, h, w in zip(["P_7", "P_11", "P_19", "CD_4"], hs, want):
            print("control %s: H = %d (%s); doubly regular: %s" % (name, h, "ok" if h == w else "MISMATCH", is_drt(ts[["P_7", "P_11", "P_19", "CD_4"].index(name)])))
        h23 = layered_H(layered, paley_out(23))
        print("control P_23: H mod 256 = %d, mod 8 = %d (%s)" % (h23, h23 % 8, "ok" if h23 % 8 == 5 and h23 == (23 * 11 * 62293308207033) % 256 else "MISMATCH"))
    elif mode == "g1":
        exact = sys.argv[2]; stride = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        m = 5; seen = set(); nsol = 0; nmat = 0; reps = []
        for a, b, c, d in gs_solutions(m):
            nsol += 1
            if (nsol - 1) % stride: continue
            M = goethals_seidel(a, b, c, d); assert is_skew_hadamard(M); nmat += 1
            new = []
            for idx in range(4 * m):
                o = descendant(M, idx); assert is_drt(o)
                if tuple(o) not in seen: seen.add(tuple(o)); new.append(o)
            classes(new, reps)
        print("G1: %d solutions (a, b, c, d) of the Goethals-Seidel equations with circulants of order 5; %d skew Hadamard matrices of order 20 examined (every %d-th solution); %d descendants, %d distinct labelled tournaments, all doubly regular; %d isomorphism classes" % (
            nsol, nmat, stride, 20 * nmat, len(seen), len(reps)))
        P = paley_out(19); hs = exact_H(exact, [r[1] for r in reps])
        for (inv, rep, members), h in zip(reps, hs):
            aut = find_map(rep, rep, count=True, limit=None); psi = anti_involution(rep)
            line = "  class with %d labelled members: H = %d = %d (mod 4) = %d (mod 8); |Aut| = %d; isomorphic to P_19: %s; self-converse: %s" % (
                len(members), h, h % 4, h % 8, aut, "yes" if find_map(rep, P) is not None else "no", "yes" if psi else "no")
            if psi:
                Z, D1, plus = fold(rep, 19, psi); t = paths_ending_in(Z, D1)
                line += "; fold: |D_1| = %d of %d, paths ending in D_1 = %d, 1 + 2*paths = %d (mod 4)" % (bin(D1).count('1'), len(plus), t, (1 + 2 * t) % 4)
            print(line); sys.stdout.flush()
        print("R9 (G1: H = 3 mod 4): %s; R10 (G1: H = 3 mod 8): %s" % (
            "HOLDS" if all(h % 4 == 3 for h in hs) else "FAILS", "HOLDS" if all(h % 8 == 3 for h in hs) else "FAILS"))
    elif mode == "g2":
        layered = sys.argv[2]
        K = double(paley_matrix(11)); tours = [descendant(K, c) for c in range(24)]
        assert all(is_drt(o) for o in tours)
        cl = classes(tours); P = paley_out(23)
        print("G2: the double of the Paley matrix of order 12; 24 descendants, all doubly regular; %d isomorphism classes" % len(cl))
        res = []
        for inv, rep, members in cl:
            h = layered_H(layered, rep); aut = find_map(rep, rep, count=True, limit=None); psi = anti_involution(rep)
            line = "  class of descendants %s: H = %d (mod 256), %d (mod 4), %d (mod 8); |Aut| = %d; isomorphic to P_23: %s; self-converse: %s" % (
                members, h, h % 4, h % 8, aut, "yes" if find_map(rep, P) is not None else "no", "yes" if psi else "no")
            if psi:
                Z, D1, plus = fold(rep, 23, psi); t = paths_ending_in(Z, D1)
                line += "; fold: |D_1| = %d of %d, 1 + 2*paths = %d (mod 4)" % (bin(D1).count('1'), len(plus), (1 + 2 * t) % 4)
            print(line); sys.stdout.flush(); res.append(h)
        # converse relations
        for i in range(len(cl)):
            for j in range(i, len(cl)):
                if find_map(cl[i][1], cl[j][1], anti=True) is not None:
                    print("  classes %d and %d are converse to each other" % (i, j))
        print("R9 (G2: H = 1 mod 4): %s; R10 (G2: H = 5 mod 8): %s" % (
            "HOLDS" if all(h % 4 == 1 for h in res) else "FAILS", "HOLDS" if all(h % 8 == 5 for h in res) else "FAILS"))
    elif mode == "g3":
        maxmat = int(sys.argv[2]) if len(sys.argv) > 2 else 200
        m = 7; reps = []; nmat = 0; ntour = 0; seenM = set()
        stride = int(sys.argv[3]) if len(sys.argv) > 3 else 1; nsol = 0
        for a, b, c, d in gs_solutions(m):
            nsol += 1
            if (nsol - 1) % stride or nmat >= maxmat: continue
            M = goethals_seidel(a, b, c, d); assert is_skew_hadamard(M); nmat += 1
            tours = [descendant(M, idx) for idx in range(28)]
            assert all(is_drt(o) for o in tours); ntour += len(tours)
            classes(tours, reps)
        P = paley_out(27)
        print("G3: %d solutions (a, b, c, d) of the Goethals-Seidel equations with circulants of order 7; %d skew Hadamard matrices of order 28 examined (every %d-th solution, at most %d); %d descendants, all doubly regular; %d isomorphism classes" % (nsol, nmat, stride, maxmat, ntour, len(reps)))
        ok = True; computed = 0
        for inv, rep, members in reps:
            psi = anti_involution(rep); isP = find_map(rep, P) is not None
            line = "  class with %d members: isomorphic to P_27: %s; self-converse: %s" % (len(members), "yes" if isP else "no", "yes" if psi else "no")
            if psi:
                Z, D1, plus = fold(rep, 27, psi); t = paths_ending_in(Z, D1); r = (1 + 2 * t) % 4; computed += 1
                line += "; fold: |D_1| = %d of %d, H = %d (mod 4)" % (bin(D1).count('1'), len(plus), r)
                if r != 3: ok = False
            print(line); sys.stdout.flush()
        print("R9 (G3: H = 3 mod 4 for the %d classes with a converse symmetry): %s" % (computed, "HOLDS" if ok else "FAILS"))
