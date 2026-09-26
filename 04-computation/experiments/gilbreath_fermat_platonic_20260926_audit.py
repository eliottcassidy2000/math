#!/usr/bin/env python3
"""gilbreath_fermat_platonic_20260926_audit.py -- independent audit of the S10 session
(note gilbreath_fermat_platonic_20260926.md, THM-4511, HYP-9162, scripts gilbreath_fermat_tower_20260926*.py,
gilbreath_fermat_platonic_20260926_groups.py, gilbreath_extinction_20260926.py).  Auditor subagent, 2026-09-26.

Everything below is re-implemented from scratch (nothing is imported from the audited scripts).
 A. sea kernel theorem: Lucas; the binary value identity; the light-cone/window condition tested on random rows
    that contain defects; the prime triangle (primes < 200000): first all-0/2 row, frontier values, the 69 kernel
    checks from row 65, the 1653 in-cone cells of row 10.
 B. single-seed diagram: EVERY interior zero of rows < 2^K lies in exactly one inverted triangle topped by a
    maximal run of ones (vertical left edge of ones / diagonal right edge of ones), all sides Mersenne; counts vs
    the sum formula and the closed form 3^(K-1-m); the identity itself checked to K = 16.
 C. the skew-Hadamard tower: skewness, orthogonality, normalization, double regularity, the recursion, the
    diagonal extension of automorphisms; |Aut(T_k)| by (i) invariant-pruned arc-consistency backtracking,
    (ii) networkx VF2 (k <= 5 and P_7, P_11, P_31), (iii) the orbit-stabilizer induction
    |Aut(T_(k+1))| = |Orb(0')| * |Aut(T_k)| (Stab(0') = diagonal Aut(T_k), which needs every arc on a 3-cycle;
    Orb(0') read off a vertex invariant) up to k = 9; element orders and orbits of Aut(T_5); T_5 vs P_31 by a
    vertex invariant; the 4-vertex census.
 D. the groups PSL(2,q), PGL(2,3) from ALL matrices (not generators), identified with sympy / by hand.
 E. von Staudt-Clausen denominators of B_(2^r) (sympy) vs 2 prod_(2^j <= r) F_j.
 F. THM-4511: exhaustive small contexts, random right contexts of length 30 and 60, column-sequence shapes,
    and explicit failures outside the hypothesis (non-sea right context; two defects; size 6).
 G. extinction table p_d(F) re-computed; dependence on the right context length R for d = 6, 8; the Monte Carlo
    replicated with the same seed.
 H. prime frontier: fresh fronts, sizes, distances, risk sums, largest terms, first-two-rows survival.
 J. miscellaneous numbers quoted in the note (drifts, constructible rows, non-result column).
Usage: python3 gilbreath_fermat_platonic_20260926_audit.py
"""
import math, sys, time, itertools
import numpy as np

sys.setrecursionlimit(20000)
T_START = time.time()


def hdr(s):
    print("\n== %s ==  (t = %.0fs)" % (s, time.time() - T_START))


def FERM(i):
    return 2 ** (2 ** i) + 1


def primes_below(P):
    s = np.ones(P, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


# ----------------------------------------------------------------------------------------------------------- A
def kernel_xor(bits, t):
    """XOR_(j subset t) bits[i + j]; the subset loop is a plain range filter (different from the audited code)"""
    n = len(bits) - t
    out = np.zeros(n, dtype=np.uint8)
    for j in range(t + 1):
        if (j & ~t) == 0:
            out ^= bits[j:j + n]
    return out


def part_A():
    hdr("A. sea kernel theorem")
    for n in range(256):
        for j in range(n + 1):
            assert (math.comb(n, j) & 1) == (1 if (j & ~n) == 0 else 0)
    print(" A1 Lucas: C(n, j) odd iff j subset n, all n < 256: OK")
    for t in range(1024):
        s = sum(1 << j for j in range(t + 1) if (j & ~t) == 0)
        p = 1
        for i in range(t.bit_length()):
            if (t >> i) & 1:
                p *= FERM(i)
        assert s == p
    print(" A2 sum_(j subset t) 2^j = prod_(i in bits t) F_i for all t < 1024: OK (expand the product: one term per subset)")
    # A3: light cone / window condition on a random row WITH defects
    rng = np.random.default_rng(1)
    L, T = 300, 100
    row = rng.choice(np.array([0, 2], dtype=np.int64), size=L)
    defects = rng.choice(L, size=25, replace=False)
    row[defects] = rng.choice(np.array([4, 6, 8, 10]), size=25)
    rows = [row.copy()]
    for r in range(T):
        row = np.abs(np.diff(row)); rows.append(row.copy())
    b0 = rows[0]
    ok_in = bad_in = bad_out = n_out = 0
    for r in range(1, T + 1):
        subs = [j for j in range(r + 1) if (j & ~r) == 0]
        sea = (b0 == 0) | (b0 == 2)
        for i in range(L - r):
            pred = 2 * (sum(int(b0[i + j] == 2) for j in subs) & 1)
            if sea[i:i + r + 1].all():
                if rows[r][i] == pred:
                    ok_in += 1
                else:
                    bad_in += 1
            else:
                n_out += 1
                if rows[r][i] != pred:
                    bad_out += 1
    print(" A3 random row (L = 300, 25 defects in {4,6,8,10}, 100 rows): cells whose row-0 window is all 0/2: %d correct, %d wrong;"
          " cells whose window contains a defect: %d, XOR formula wrong for %d of them (so the window condition is needed)" % (ok_in, bad_in, n_out, bad_out))
    assert bad_in == 0
    # A4: the prime triangle
    P = 200000
    row = primes_below(P)
    keep, frontier, first_all = {}, {}, None
    want = set(range(1, 66)) | {65 + t for t in range(1, 65)} | {65 + t for t in (127, 255, 511, 1023, 2047)} | set(range(10, 70))
    for r in range(1, 65 + 2047 + 1):
        row = np.abs(np.diff(row))
        assert row[0] == 1
        big = np.nonzero(row[1:] >= 4)[0]
        frontier[r] = int(big[0]) + 1 if len(big) else None
        if first_all is None and len(big) == 0:
            first_all = r
        if r in want:
            keep[r] = row.copy()
    print(" A4 primes < %d: first row that is 1 followed by 0/2 only: %d (note says 65); rows 1..64 each contain an entry >= 4: %s"
          % (P, first_all, all(frontier[r] is not None for r in range(1, 65))))
    fr = [frontier[r] for r in (1, 2, 5, 10, 20, 30, 40, 50, 60, 64)]
    print("    frontier F(r) at r = 1,2,5,10,20,30,40,50,60,64:", fr, "(note: [3, 8, 25, 59, 291, 870, 2770, 2763, 5942, 5940])")
    assert fr == [3, 8, 25, 59, 291, 870, 2770, 2763, 5942, 5940]
    w = (keep[65][1:] // 2).astype(np.uint8)
    assert set(np.unique(keep[65][1:]).tolist()) <= {0, 2}
    nchk = 0
    for t in list(range(1, 65)) + [127, 255, 511, 1023, 2047]:
        act = (keep[65 + t][1:] // 2).astype(np.uint8)
        k = kernel_xor(w, t)
        assert len(act) == len(k) and np.array_equal(act, k) and keep[65 + t][0] == 1, t
        nchk += 1
    print("    rows 65+t = Lucas-kernel image of row 65 for t = 1..64, 127, 255, 511, 1023, 2047: %d checks OK, leading 1 throughout (note: 69)" % nchk)
    F10 = frontier[10]; d10 = int(keep[10][F10])
    w = (keep[10][1:F10] // 2).astype(np.uint8)
    tot = 0
    for t in range(1, F10 - 1):
        k = kernel_xor(w, t)
        assert np.array_equal(keep[10 + t][1:1 + len(k)], 2 * k)
        tot += len(k)
    print("    row 10: F = %d, defect %d, in-cone cells over t = 1..%d all equal the kernel image: %d cells (note: 1653; 57*58/2 = %d)" % (F10, d10, F10 - 2, tot, 57 * 58 // 2))
    assert tot == 1653
    return frontier


# ----------------------------------------------------------------------------------------------------------- B
def part_B():
    hdr("B. single-seed diagram: zero triangles")
    for K in (5, 7, 9):
        p = 2 ** K
        a = np.zeros(2 ** (K + 1) + 2, dtype=np.uint8); a[p] = 1
        rows = [a.copy()]
        for r in range(2 ** (K + 1)):
            a = a[:-1] ^ a[1:]; rows.append(a.copy())
        cover = {}; hist = {}; tid = 0
        for r in range(2 ** K):
            row = rows[r]; n = len(row); i = 0
            while i < n:
                if row[i]:
                    j = i
                    while j < n and row[j]:
                        j += 1
                    Lr = j - i
                    assert Lr & (Lr - 1) == 0, ("run length not a power of two", r, i, Lr)
                    assert i > 0 and j < n and row[i - 1] == 0 and row[j] == 0
                    t = Lr - 1
                    if t >= 1:
                        hist[t] = hist.get(t, 0) + 1
                        for s in range(1, t + 1):
                            seg = rows[r + s][i:i + t - s + 1]
                            assert len(seg) == t - s + 1 and not seg.any()
                            assert rows[r + s][i - 1] == 1 and rows[r + s][i + t - s + 1] == 1   # ones on both edges
                            for c in range(i, i + t - s + 1):
                                assert (r + s, c) not in cover, ("zero covered twice", r + s, c)
                                cover[(r + s, c)] = tid
                        tid += 1
                    i = j
                else:
                    i += 1
        interior = missing = 0
        for r in range(1, 2 ** K):
            row = rows[r]
            for c in range(p - r, p + 1):          # the diagram occupies columns p-r..p in row r
                if row[c] == 0:
                    interior += 1
                    if (r, c) not in cover:
                        missing += 1
        sides = sorted(hist)
        pred = {}
        for n in range(2 ** K):
            m = 0
            while (n >> m) & 1:
                m += 1
            if m >= 1:
                pred[2 ** m - 1] = pred.get(2 ** m - 1, 0) + 2 ** (bin(n).count('1') - m)
        closed = {2 ** m - 1: (3 ** (K - 1 - m) if m < K else 1) for m in range(1, K + 1)}
        print(" K=%d: sides %s; counts %s; sum formula %s; closed form 3^(K-1-m) (+1 for the all-ones row) %s; interior zeros in rows 1..%d: %d, not covered by a triangle: %d"
              % (K, sides, [hist[t] for t in sides], [pred[t] for t in sides], [closed[t] for t in sides], 2 ** K - 1, interior, missing))
        assert hist == pred == closed and missing == 0
    for K in range(1, 17):
        for m in range(1, K):
            s = sum(2 ** (bin(n).count('1') - m) for n in range(2 ** K) if (n & (2 ** m - 1)) == 2 ** m - 1 and not (n >> m) & 1)
            assert s == 3 ** (K - 1 - m), (K, m, s)
    print(" identity  sum_(n < 2^K, exactly m trailing ones) 2^(popcount(n) - m) = 3^(K-1-m)  for K <= 16, 1 <= m <= K-1: OK"
          " (proof: bits 0..m-1 are 1, bit m is 0, the K-1-m bits above are free and contribute (1+2)^(K-1-m))")


# ----------------------------------------------------------------------------------------------------------- C
def tower_H(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1):
        H = np.block([[H, H], [-H.T, H.T]])
    return H


def tour_from_H(H):
    A = (H[1:, 1:] == 1).copy()
    np.fill_diagonal(A, False)
    return A


def paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    A = np.zeros((q, q), dtype=bool)
    for i in range(q):
        for j in range(q):
            if i != j and (j - i) % q in sq:
                A[i, j] = True
    return A


def is_drt(A):
    n = A.shape[0]; Ai = A.astype(np.int64)
    if not np.array_equal(Ai + Ai.T + np.eye(n, dtype=np.int64), np.ones((n, n), dtype=np.int64)):
        return False
    if not np.all(Ai.sum(1) == (n - 1) // 2):
        return False
    C = Ai @ Ai.T
    return bool(np.all(C[~np.eye(n, dtype=bool)] == (n - 3) // 4))


def check_recursion(Ak, Ak1):
    n = Ak.shape[0]; assert Ak1.shape[0] == 2 * n + 1
    for i in range(n):
        for j in range(n):
            if i != j:
                assert Ak1[i, j] == Ak[i, j] and Ak1[n + 1 + i, n + 1 + j] == Ak[j, i]
                assert Ak1[i, n + 1 + j] == Ak[i, j] and Ak1[n + 1 + i, j] == Ak[i, j]
        assert Ak1[i, n + 1 + i] and not Ak1[n + 1 + i, i]
        assert Ak1[n, i] and not Ak1[i, n]
        assert Ak1[n + 1 + i, n] and not Ak1[n, n + 1 + i]
    return True


def outmasks(A):
    n = A.shape[0]
    return [int(sum(1 << j for j in np.nonzero(A[i])[0].tolist())) for i in range(n)]


def vertex_invariant(A):
    """for each v: histograms of |N(u) cap N(w) cap N| over pairs u < w in N, for N = N^+(v) and N = N^-(v)
    (N(u) = out-neighbourhood). An isomorphism invariant of the vertex."""
    n = A.shape[0]; Af = A.astype(np.float32); inv = []
    for v in range(n):
        res = []
        for N in (A[v], ~A[v] & ~(np.arange(n) == v)):
            idx = np.nonzero(N)[0]
            M = Af[np.ix_(idx, idx)]
            C = np.rint(M @ M.T).astype(np.int64)
            iu = np.triu_indices(len(idx), 1)
            vals, cnts = np.unique(C[iu], return_counts=True)
            res.append(tuple(zip(vals.tolist(), cnts.tolist())))
        inv.append(tuple(res))
    return inv


class Cap(Exception):
    pass


def aut_backtrack(out, n, classes, fixed=None, collect=False, node_cap=3_000_000, tcap=900, first_only=False):
    """Count (or collect) bijections phi with u -> v iff phi(u) -> phi(v), phi respecting the vertex classes.
    Plain arc-consistency backtracking on bitmasks (no refinement). fixed = {v: w} pre-assigns images."""
    csize = {}
    for v in range(n):
        csize[classes[v]] = csize.get(classes[v], 0) + 1
    cmask = {}
    for v in range(n):
        cmask[classes[v]] = cmask.get(classes[v], 0) | (1 << v)
    order = sorted(range(n), key=lambda v: (csize[classes[v]], v))
    if fixed:
        order = list(fixed) + [v for v in order if v not in fixed]
    phi = [-1] * n; res = {'count': 0, 'nodes': 0}; auts = []; t0 = time.time()

    def rec(d, used):
        res['nodes'] += 1
        if res['nodes'] > node_cap or time.time() - t0 > tcap:
            raise Cap()
        if d == n:
            res['count'] += 1
            if collect:
                auts.append(tuple(phi))
            if first_only:
                raise Cap()
            return
        v = order[d]
        cand = cmask[classes[v]] & ~used
        if fixed and v in fixed:
            cand &= (1 << fixed[v])
        for u in order[:d]:
            pu = phi[u]
            if (out[u] >> v) & 1:
                cand &= out[pu]
            else:
                cand &= ~out[pu]
            if not cand:
                return
        while cand:
            w = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            phi[v] = w
            rec(d + 1, used | (1 << w))
            phi[v] = -1

    status = 'done'
    try:
        rec(0, 0)
    except Cap:
        status = 'found' if (first_only and res['count']) else 'cap'
    return status, res['count'], res['nodes'], auts


def aut_count_vf2(A):
    import networkx as nx
    from networkx.algorithms.isomorphism import DiGraphMatcher
    n = A.shape[0]
    G = nx.DiGraph(); G.add_nodes_from(range(n))
    G.add_edges_from((i, j) for i in range(n) for j in range(n) if A[i, j])
    return sum(1 for _ in DiGraphMatcher(G, G).isomorphisms_iter())


def census4(out, n):
    types = {}
    for S in itertools.combinations(range(n), 4):
        m = sum(1 << v for v in S)
        sc = tuple(sorted((out[v] & m).bit_count() for v in S))
        types[sc] = types.get(sc, 0) + 1
    return types


def part_C():
    hdr("C. the skew-Hadamard doubling tower and its tournaments")
    Ts = {}; prev = None
    for k in range(1, 10):
        H = tower_H(k); N = H.shape[0]
        assert np.array_equal(H + H.T, 2 * np.eye(N, dtype=np.int64)), k          # skew: H - I antisymmetric
        assert np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64)), k          # Hadamard
        assert np.all(H[0, :] == 1) and np.all(H[1:, 0] == -1)                   # normalized
        A = tour_from_H(H)
        if k >= 2:
            assert is_drt(A), k
        if prev is not None and k >= 2:
            check_recursion(prev, A)
        Ts[k] = A; prev = A
    print(" H_(2^k) skew Hadamard + normalized for k = 1..9; T_k doubly regular for k = 2..9 (n = 3..511); recursion T_(k+1) = T_k + {0'} + T_k' verified for k = 1..8")
    # every arc on a 3-cycle (needed for Stab(0') = diagonal Aut(T_k)); count 2-paths y -> z -> x for each arc x -> y
    for k in range(2, 10):
        A = Ts[k]; Ai = A.astype(np.int64); n = A.shape[0]
        two = Ai @ Ai                    # two[y, x] = # z with y -> z -> x
        assert np.all(two.T[A] >= 1)     # for each arc x -> y: two[y, x] >= 1
        assert np.all(two.T[A] == (n + 1) // 4)
    print(" every arc x -> y of T_k (k = 2..9) lies on exactly (n+1)/4 three-cycles, so the closed out-neighbourhoods N^+[i] separate points")
    # diagonal extension of automorphisms: collect Aut(T_k) for k = 2..6 by backtracking (own code)
    auts = {}
    print(" |Aut(T_k)| by invariant-pruned arc-consistency backtracking (own code):")
    for k in range(2, 10):
        A = Ts[k]; n = A.shape[0]; out = outmasks(A)
        t0 = time.time()
        inv = vertex_invariant(A)
        cls = {}
        classes = [cls.setdefault(inv[v], len(cls)) for v in range(n)]
        st, cnt, nodes, al = aut_backtrack(out, n, classes, collect=(k <= 6), node_cap=20_000_000, tcap=1200)
        sizes = sorted([classes.count(c) for c in set(classes)], reverse=True)
        print("   T_%d (n = %3d): |Aut| = %s [%s, %d nodes, %.1fs]; invariant classes (sizes) %s" % (k, n, cnt if st == 'done' else '>= %d' % cnt, st, nodes, time.time() - t0, sizes))
        if st == 'done' and k <= 6:
            auts[k] = al
    # the diagonal extension is an automorphism; and all automorphisms of T_(k+1) fixing 0' are diagonal
    for k in range(2, 6):
        n = Ts[k].shape[0]; out1 = outmasks(Ts[k + 1]); N = 2 * n + 1
        for s in auts[k]:
            ext = list(s) + [n] + [n + 1 + s[i] for i in range(n)]
            assert all(((out1[i] >> j) & 1) == ((out1[ext[i]] >> ext[j]) & 1) for i in range(N) for j in range(N) if i != j)
        fix = [g for g in auts[k + 1] if g[n] == n]
        diag = [g for g in fix if all(g[n + 1 + i] == n + 1 + g[i] for i in range(n))]
        print("   k=%d: all %d automorphisms of T_%d extend diagonally to T_%d (checked); automorphisms of T_%d fixing 0': %d, all diagonal: %s; |Aut(T_%d)| = %d"
              % (k, len(auts[k]), k, k + 1, k + 1, len(fix), len(fix) == len(diag), k + 1, len(auts[k + 1])))
        assert len(fix) == len(diag) == len(auts[k])
    # element orders and orbits of Aut(T_5); Frobenius identification
    al = auts[5]; n = 31
    orders = {}
    for g in al:
        o = 1; x = g
        while x != tuple(range(n)):
            x = tuple(g[x[i]] for i in range(n)); o += 1
        orders[o] = orders.get(o, 0) + 1
    seen = set(); orbits = []
    for v in range(n):
        if v not in seen:
            orb = sorted({g[v] for g in al}); orbits.append(orb); seen |= set(orb)
    comm = all(tuple(g[h[i]] for i in range(n)) == tuple(h[g[i]] for i in range(n)) for g in al for h in al)
    print(" Aut(T_5): order %d, element orders %s (note: {3: 14, 7: 6}), abelian: %s (the nonabelian group of order 21 is Z_7 x| Z_3), orbit sizes %s (note: 7,1,7,1,7,1,7)"
          % (len(al), dict(sorted(orders.items())), comm, [len(o) for o in orbits]))
    assert len(al) == 21 and orders == {1: 1, 3: 14, 7: 6} and not comm and [len(o) for o in orbits] == [7, 1, 7, 1, 7, 1, 7]
    # VF2 cross-check (networkx), plus Paley groups by both methods
    print(" networkx VF2 automorphism counts (independent implementation):")
    for name, A in (("T_2", Ts[2]), ("T_3", Ts[3]), ("T_4", Ts[4]), ("T_5", Ts[5]), ("T_6", Ts[6]), ("P_7", paley(7)), ("P_11", paley(11)), ("P_31", paley(31))):
        t0 = time.time(); c = aut_count_vf2(A)
        print("   %-4s (n = %2d): |Aut| = %d  (%.1fs)" % (name, A.shape[0], c, time.time() - t0))
        if name.startswith("P_"):
            q = A.shape[0]; assert c == q * (q - 1) // 2, (name, c)
        if name in ("T_4", "T_5", "T_6"):
            assert c == 21
    for q in (7, 11, 31):
        A = paley(q); assert is_drt(A); out = outmasks(A); inv = vertex_invariant(A)
        cls = {}; classes = [cls.setdefault(inv[v], len(cls)) for v in range(q)]
        st, cnt, nodes, _ = aut_backtrack(out, q, classes)
        assert st == 'done' and cnt == q * (q - 1) // 2, (q, st, cnt)
    print("   backtracking on P_7, P_11, P_31: 21, 55, 465 = q(q-1)/2: OK")
    import networkx as nx
    G3 = nx.DiGraph([(i, j) for i in range(7) for j in range(7) if Ts[3][i, j]]); P7 = nx.DiGraph([(i, j) for i in range(7) for j in range(7) if paley(7)[i, j]])
    print(" T_3 isomorphic to the Paley heptagon P_7 (networkx): %s" % nx.is_isomorphic(G3, P7))
    assert nx.is_isomorphic(G3, P7)
    # T_5 vs P_31: vertex invariant multisets
    invT = vertex_invariant(Ts[5]); invP = vertex_invariant(paley(31))
    nT = len(set(invT)); nP = len(set(invP))
    print(" T_5 vs P_31: the vertex invariant takes %d distinct values on T_5 (class sizes %s) and %d on P_31; multisets equal: %s -> NOT isomorphic (also |Aut| 21 vs 465)"
          % (nT, sorted([invT.count(x) for x in set(invT)], reverse=True), nP, sorted(invT) == sorted(invP)))
    assert sorted(invT) != sorted(invP)
    # 4-vertex census
    names = {(0, 1, 2, 3): 'transitive', (1, 1, 1, 3): 'source+3cycle', (0, 2, 2, 2): 'sink+3cycle', (1, 1, 2, 2): 'strong'}
    for lab, A in (("T_3", Ts[3]), ("T_4", Ts[4]), ("T_5", Ts[5]), ("P_7", paley(7)), ("P_31", paley(31))):
        c = census4(outmasks(A), A.shape[0])
        print("   4-vertex census %-4s: %s" % (lab, {names[k]: c[k] for k in sorted(c)}))
    cT = census4(outmasks(Ts[5]), 31); cP = census4(outmasks(paley(31)), 31)
    assert cT == cP and cT[(1, 1, 1, 3)] == 4340 and cT[(1, 1, 2, 2)] == 13020 and cT[(0, 2, 2, 2)] == 4340 and cT[(0, 1, 2, 3)] == 9765
    print("   T_5 and P_31 have the same 4-vertex census (4340 source, 13020 strong, 4340 sink, 9765 transitive; sum 31465 = C(31,4)): OK")
    # orbit-stabilizer induction for larger k: |Aut(T_(k+1))| = |Orb(0')| * |Aut(T_k)|
    print(" induction |Aut(T_(k+1))| = |Orb(0')| * |Aut(T_k)| with Orb(0') from the vertex invariant (0' = vertex index 2^k - 1 of T_(k+1)):")
    order = {3: 21}
    for k in range(3, 9):
        A = Ts[k + 1]; N = A.shape[0]; n0 = (N - 1) // 2
        t0 = time.time(); inv = vertex_invariant(A)
        same = [v for v in range(N) if inv[v] == inv[n0]]
        orb = 1; capped = 0
        if len(same) > 1:
            out = outmasks(A); cls = {}; classes = [cls.setdefault(inv[v], len(cls)) for v in range(N)]
            for v in same:
                if v != n0:
                    st, cnt, nodes, _ = aut_backtrack(out, N, classes, fixed={n0: v}, first_only=True, node_cap=5_000_000, tcap=600)
                    if st == 'found':
                        orb += 1
                    elif st != 'done':
                        capped += 1
        if capped or order.get(k) is None:
            order[k + 1] = None
            print("   T_%d (n = %3d): vertices with the invariant of 0': %d, %d candidate images undecided (search capped): |Aut| NOT determined here" % (k + 1, N, len(same), capped))
        else:
            order[k + 1] = order[k] * orb
            print("   T_%d (n = %3d): vertices with the invariant of 0': %d, |Orb(0')| = %d, |Aut(T_%d)| = %d  (%.1fs)" % (k + 1, N, len(same), orb, k + 1, order[k + 1], time.time() - t0))
    return Ts


# ----------------------------------------------------------------------------------------------------------- D
class GF:
    def __init__(self, q):
        self.q = q
        if q == 4:
            # elements 0,1,2=a,3=a+1 with a^2 = a + 1
            self.mul_t = [[0] * 4 for _ in range(4)]
            for x in range(4):
                for y in range(4):
                    r = 0
                    for i in range(2):
                        if (y >> i) & 1:
                            r ^= x << i
                    if r & 4:
                        r ^= 0b111
                    self.mul_t[x][y] = r
            self.add = lambda x, y: x ^ y
            self.mul = lambda x, y: self.mul_t[x][y]
            self.inv = {x: next(y for y in range(1, 4) if self.mul_t[x][y] == 1) for x in range(1, 4)}
        else:
            self.add = lambda x, y: (x + y) % q
            self.mul = lambda x, y: (x * y) % q
            self.inv = {x: pow(x, q - 2, q) for x in range(1, q)}

    def det(self, a, b, c, d):
        return self.add(self.mul(a, d), (-self.mul(b, c)) % self.q if self.q != 4 else self.mul(b, c))


def mobius(K, a, b, c, d):
    q = K.q; perm = []
    for x in range(q + 1):
        if x == q:
            perm.append(q if c == 0 else K.mul(a, K.inv[c]))
        else:
            den = K.add(K.mul(c, x), d); num = K.add(K.mul(a, x), b)
            perm.append(q if den == 0 else K.mul(num, K.inv[den]))
    return tuple(perm)


def proj_group(q, pgl=False):
    K = GF(q); G = set()
    for a in range(q):
        for b in range(q):
            for c in range(q):
                for d in range(q):
                    dt = K.det(a, b, c, d)
                    if dt == 0:
                        continue
                    if pgl or dt == 1:
                        G.add(mobius(K, a, b, c, d))
    return sorted(G)


def part_D():
    hdr("D. rotation groups of the solids as projective groups (all matrices enumerated; sympy identification)")
    from sympy.combinatorics import Permutation, PermutationGroup
    sph = [(p, q) for p in range(3, 40) for q in range(3, 40) if (p - 2) * (q - 2) < 4]
    euc = [(p, q) for p in range(3, 40) for q in range(3, 40) if (p - 2) * (q - 2) == 4]
    print(" Schlafli: (p-2)(q-2) < 4: %s (5); = 4: %s (3)" % (sph, euc))
    assert len(sph) == 5 and len(euc) == 3
    res = {}
    for q, pgl, name in ((2, False, 'PSL(2,2)'), (3, False, 'PSL(2,3)'), (3, True, 'PGL(2,3)'), (4, False, 'PSL(2,4)=SL(2,4)'), (5, False, 'PSL(2,5)'), (7, False, 'PSL(2,7)'), (17, False, 'PSL(2,17)')):
        G = proj_group(q, pgl)
        expected = q * (q * q - 1) * (1 if (q % 2 == 0 or pgl) else 1) // (1 if (q % 2 == 0 or pgl) else 2)
        PG = PermutationGroup([Permutation(list(g)) for g in G]) if len(G) <= 200 else None
        line = " %-16s order %4d on %2d points (formula %4d)" % (name, len(G), q + 1, expected)
        if PG is not None:
            line += "; sympy: is_alternating=%s is_symmetric=%s" % (PG.is_alternating, PG.is_symmetric)
            if name == 'PSL(2,2)':
                assert PG.is_symmetric and len(G) == 6
            if name == 'PSL(2,3)':
                assert PG.is_alternating and len(G) == 12
            if name == 'PGL(2,3)':
                assert PG.is_symmetric and len(G) == 24
            if name.startswith('PSL(2,4)'):
                assert PG.is_alternating and len(G) == 60
            if name == 'PSL(2,5)':
                ident = tuple(range(6))
                comp = lambda g, h: tuple(g[h[i]] for i in range(6))
                invs = [g for g in G if g != ident and comp(g, g) == ident]
                kleins = set()
                for a, b in itertools.combinations(invs, 2):
                    if comp(a, b) == comp(b, a):
                        kleins.add(frozenset([ident, a, b, comp(a, b)]))
                kleins = sorted(kleins, key=lambda s: sorted(s)); idx = {s: i for i, s in enumerate(kleins)}
                images = set(); kernel = 0
                for g in G:
                    gi = [0] * 6
                    for i, x in enumerate(g):
                        gi[x] = i
                    gi = tuple(gi)
                    img = tuple(idx[frozenset(comp(comp(g, h), gi) for h in s)] for s in kleins)
                    images.add(img); kernel += (img == tuple(range(len(kleins))))
                Im = PermutationGroup([Permutation(list(p)) for p in images])
                stab = sum(1 for g in G if g[0] == 0)
                line += "; involutions %d, Klein four-subgroups %d, conjugation action on them: kernel %d, image order %d, image is A_5: %s; point stabilizer order %d (D_10), transitive on 6 points: %s" % (
                    len(invs), len(kleins), kernel, len(images), Im.is_alternating, stab, len({g[0] for g in G}) == 6)
                assert len(invs) == 15 and len(kleins) == 5 and kernel == 1 and len(images) == 60 and Im.is_alternating and stab == 10
        assert len(G) == expected
        print(line)
    print(" orders 6, 12, 24, 60, 60, 168, 2448 confirmed; finite rotation groups of R^3 have orders n, 2n, 12, 24, 60 (Klein), so PSL(2,7) and PSL(2,17) are not among them")


# ----------------------------------------------------------------------------------------------------------- E
def part_E():
    hdr("E. von Staudt-Clausen: denominator of B_(2^r) = 2 prod_(F_j prime, 2^j <= r) F_j")
    from sympy import bernoulli, Rational
    known = [3, 5, 17, 257, 65537]
    for r in range(1, 10):
        n = 2 ** r
        B = Rational(bernoulli(n)); den = int(B.q)
        pred = 2
        for j, Fj in enumerate(known):
            if 2 ** j <= r:
                pred *= Fj
        vsc = 1
        for p in range(2, n + 2):
            if all(p % d for d in range(2, int(p ** 0.5) + 1)) and n % (p - 1) == 0:
                vsc *= p
        print("  r=%d, n=%3d: denom(B_n) = %d ; 2 prod F_j = %d ; prod_(p-1 | n) p = %d ; equal: %s" % (r, n, den, pred, vsc, den == pred == vsc))
        assert den == pred == vsc


# ----------------------------------------------------------------------------------------------------------- F
def evolve_iter(row0):
    """yields the successive rows (N, W - r), r = 1, 2, ..., down to width 1, without storing them"""
    r = row0
    while r.shape[1] > 1:
        r = np.abs(r[:, :-1].astype(np.int16) - r[:, 1:].astype(np.int16)).astype(np.uint8)
        yield r


def evolve(row0):
    """row0: (N, W) uint8; returns the list of all rows (N, W - r) down to width 1 (small inputs only)"""
    return [row0] + list(evolve_iter(row0))


def loss_fraction(row0):
    lost = np.zeros(row0.shape[0], dtype=bool)
    for r in evolve_iter(row0):
        lost |= (r[:, 0] != 1)
    return float(lost.mean())


def build(left, right, d, F):
    N = left.shape[0]; R = right.shape[1]
    row = np.empty((N, F + R + 1), dtype=np.uint8)
    row[:, 0] = 1; row[:, 1:F] = left; row[:, F] = d; row[:, F + 1:] = right
    return row


def all_words(nb):
    idx = np.arange(1 << nb, dtype=np.int64)
    return (((idx[:, None] >> np.arange(nb, dtype=np.int64)[None, :]) & 1) * 2).astype(np.uint8)


def wall_checks(left, right, F, label):
    N = left.shape[0]
    lost = np.zeros(N, dtype=bool); mincol4 = np.full(N, F, dtype=np.int64); capF1 = capF = None
    for s, r in enumerate(evolve_iter(build(left, right, 4, F)), start=1):
        lost |= (r[:, 0] != 1)
        W = r.shape[1]; cols = np.arange(W)
        has4 = (r >= 4); has4[:, 0] = False
        mincol4 = np.minimum(mincol4, np.where(has4, cols[None, :], 10 ** 6).min(1))
        if s == F - 1:
            capF1 = r[:, 1].copy()
        if s == F:
            capF = r[:, 0].copy()
    all0 = (left == 0).all(1)
    # z = trailing zeros of the left sea (cells F-1, F-2, ...), c = F - 1 - z (0 means: no c)
    z = np.zeros(N, dtype=np.int64); alive = np.ones(N, dtype=bool)
    for col in range(F - 1, 0, -1):
        alive &= (left[:, col - 1] == 0); z += alive
    c = F - 1 - z
    ok1 = np.array_equal(lost, all0)
    ok2 = (bool(np.all(capF1[all0] == 4)) and bool(np.all(capF[all0] == 3))) if all0.any() else True
    # no entry >= 4 in columns 1..c whenever c >= 1; the 4 reaches exactly column F - z
    ok3 = bool(np.all(mincol4[c >= 1] > c[c >= 1])) and bool(np.all(mincol4 == F - z))
    print("  %s: N = %d configurations: (leading 1 destroyed) == (left cells all 0): %s; if destroyed, a_(F-1)(1) = 4 and a_F(0) = 3: %s;"
          " min column of any entry >= 4 is exactly F - z (> c): %s" % (label, N, ok1, ok2, ok3))
    assert ok1 and ok2 and ok3


def column_shapes(left, right, F):
    """per configuration, verify the column-sequence shapes of the proof: columns c+1..F are {0,4}* (2 {0,2}*)?,
    column c is 2 until column c+1 first shows 2, then 0, then {0,2}; columns 1..c-1 in {0,2}"""
    N = left.shape[0]; rows = evolve(build(left, right, 4, F)); bad = 0
    for cfg in range(N):
        z = 0
        while z < F - 1 and left[cfg, F - 2 - z] == 0:
            z += 1
        c = F - 1 - z
        cols = {}
        for i in range(0, F + 1):
            cols[i] = [int(r[cfg, i]) for r in rows if r.shape[1] > i]
        first2 = {}
        for i in range(F, c, -1):
            seq = cols[i]
            k = next((t for t, x in enumerate(seq) if x not in (0, 4)), None)
            first2[i] = k
            if k is not None and (seq[k] != 2 or any(x not in (0, 2) for x in seq[k:])):
                bad += 1
            if i < F and first2[i + 1] is not None and k != first2[i + 1] + 1:
                bad += 1                      # column i shows its 2 exactly one row after column i+1
        if c >= 1:
            seq = cols[c]; k = first2[c + 1]
            if k is None:
                if any(x != 2 for x in seq):
                    bad += 1
            else:
                if any(x != 2 for x in seq[:k + 1]) or (len(seq) > k + 1 and seq[k + 1] != 0) or any(x not in (0, 2) for x in seq[k + 1:]):
                    bad += 1
            for i in range(1, c):
                if any(x not in (0, 2) for x in cols[i]):
                    bad += 1
        if any(x != 1 for x in cols[0]) != (c == 0):
            bad += 1
    return bad


def part_F():
    hdr("F. THM-4511 (wall theorem) brute force")
    for F in range(2, 11):
        R = 10
        words = all_words(F - 1 + R)
        wall_checks(words[:, :F - 1], words[:, F - 1:], F, "exhaustive F=%2d R=%2d" % (F, R))
    rng = np.random.default_rng(20260926)
    for F in range(2, 11):
        for R, reps in ((30, 64), (60, 8)):
            left = np.repeat(all_words(F - 1), reps, axis=0)
            right = (rng.integers(0, 2, size=(left.shape[0], R)) * 2).astype(np.uint8)
            wall_checks(left, right, F, "all left words x %d random right contexts, F=%2d R=%2d" % (reps, F, R))
    for F in range(2, 8):
        R = 7
        words = all_words(F - 1 + R)
        bad = column_shapes(words[:, :F - 1], words[:, F - 1:], F)
        print("  column-sequence shapes of the proof (columns c+1..F are {0,4}* then 2 then {0,2}*, each one row after its right neighbour;"
              " column c is 2 until then, then 0; columns < c in {0,2}): F=%d R=%d, violations: %d" % (F, R, bad))
        assert bad == 0
    # front-only counts: C(F-1, w) sea words give a front path of weight w
    for F in range(3, 13):
        cnt = {}
        for bits in range(1 << (F - 1)):
            b = [(bits >> i) & 1 for i in range(F - 1)]   # b[i] = cell at column i+1
            wgt = 0
            for s in range(F - 1):
                col = F - 1 - s
                v = 0
                for j in range(s + 1):
                    if (j & ~s) == 0:
                        v ^= b[col - 1 + j]
                wgt += v
            cnt[wgt] = cnt.get(wgt, 0) + 1
        assert all(cnt.get(w, 0) == math.comb(F - 1, w) for w in range(F)), (F, cnt)
    print("  front path weights: exactly C(F-1, w) sea words give weight w, F = 3..12 (own subset loop): OK")
    # outside the hypothesis
    print("  outside the hypothesis (these are NOT covered by the theorem, and the theorem's conclusion fails):")
    row = np.array([[1, 2, 4, 10]], dtype=np.uint8)
    rows = evolve(row)
    print("    non-sea right context 1 2 4 10 (c = 1 exists): rows =", [r[0].tolist() for r in rows], "-> leading entry", int(rows[-1][0, 0]))
    rows = evolve(np.array([[1, 2, 0, 6, 0, 0, 0, 0]], dtype=np.uint8))
    print("    a lone 6 (c = 1): rows =", [r[0].tolist() for r in rows], "-> the wall column becomes |6-2| = 4 and the leading entry ends as", int(rows[-1][0, 0]))
    # SEVERAL size-4 defects: the note/THM say a second 4 can 'reopen the wall'. Test the multi-defect statement:
    # cells in {0,2,4} (any number of 4s), F1 = first 4, c = largest 2-column below F1 (if any). Claim: no entry >= 4
    # ever reaches a column <= c, and the leading 1 is destroyed iff a 4 exists and no 2 precedes the first 4.
    print("  SEVERAL size-4 defects (all cells in {0,2,4}, any number of 4s; c = largest 2 before the FIRST 4):")
    for L in range(2, 14):
        N = 3 ** L
        idx = np.arange(N, dtype=np.int64); cells = np.empty((N, L), dtype=np.uint8); t = idx.copy()
        for j in range(L):
            cells[:, j] = (t % 3) * 2; t //= 3
        row0 = np.concatenate([np.ones((N, 1), dtype=np.uint8), cells], axis=1)
        is4 = cells == 4; has4 = is4.any(1)
        F1 = np.where(has4, is4.argmax(1) + 1, L + 1)                     # column of the first 4 (L+1 if none)
        is2 = (cells == 2) & (np.arange(1, L + 1)[None, :] < F1[:, None])
        c = np.where(is2.any(1), L - 1 - np.argmax(is2[:, ::-1], axis=1) + 1, 0)   # largest 2-column below F1
        lost = np.zeros(N, dtype=bool); mincol4 = np.full(N, L + 1, dtype=np.int64)
        for r in evolve_iter(row0):
            lost |= (r[:, 0] != 1)
            W = r.shape[1]; big = (r >= 4); big[:, 0] = False
            mincol4 = np.minimum(mincol4, np.where(big, np.arange(W)[None, :], 10 ** 6).min(1))
        pred_lost = has4 & (c == 0)
        ok_a = np.array_equal(lost, pred_lost)
        ok_b = bool(np.all(mincol4[c >= 1] > c[c >= 1]))
        n_multi = int((is4.sum(1) >= 2).sum())
        print("    L=%2d: all %7d rows (%7d with >= 2 fours): destroyed iff (a 4 exists and no 2 before the first 4): %s; no entry >= 4 ever at a column <= c: %s" % (L, N, n_multi, ok_a, ok_b))
        assert ok_a and ok_b
    rng2 = np.random.default_rng(7)
    L = 40; N = 200000
    cells = rng2.choice(np.array([0, 2, 4], dtype=np.uint8), size=(N, L), p=[0.45, 0.35, 0.2])
    row0 = np.concatenate([np.ones((N, 1), dtype=np.uint8), cells], axis=1)
    is4 = cells == 4; has4 = is4.any(1); F1 = np.where(has4, is4.argmax(1) + 1, L + 1)
    is2 = (cells == 2) & (np.arange(1, L + 1)[None, :] < F1[:, None])
    c = np.where(is2.any(1), L - 1 - np.argmax(is2[:, ::-1], axis=1) + 1, 0)
    lost = np.zeros(N, dtype=bool); mincol4 = np.full(N, L + 1, dtype=np.int64)
    for r in evolve_iter(row0):
        lost |= (r[:, 0] != 1)
        W = r.shape[1]; big = (r >= 4); big[:, 0] = False
        mincol4 = np.minimum(mincol4, np.where(big, np.arange(W)[None, :], 10 ** 6).min(1))
    ok_a = np.array_equal(lost, has4 & (c == 0)); ok_b = bool(np.all(mincol4[c >= 1] > c[c >= 1]))
    print("    L=40, %d random rows with P(4) = 0.2 (mean %.1f fours per row): destroyed iff (4 exists, no 2 before it): %s; no entry >= 4 at a column <= c: %s" % (N, is4.sum(1).mean(), ok_a, ok_b))
    assert ok_a and ok_b
    print("    => the wall holds for ANY number of size-4 defects (proof: the 'good column' induction only needs each column right of c to START in {0,2,4};"
          " start the induction at the rightmost 4). The note's 'a second 4 ... can reopen the wall' does not happen.")


# ----------------------------------------------------------------------------------------------------------- G
def exact_p(d, F, R):
    words = all_words(F - 1 + R)
    return loss_fraction(build(words[:, :F - 1], words[:, F - 1:], d, F))


def front_only(d, F):
    j = d // 2
    return sum(math.comb(F - 1, t) for t in range(0, j - 1)) / 2 ** (F - 1)


def part_G():
    hdr("G. extinction table p_d(F) (own code), R-dependence, Monte Carlo replication")
    quoted = {(4, 8): 0.007812500, (4, 12): 0.000488281, (6, 8): 0.067993164, (6, 12): 0.006338120, (8, 8): 0.244003296, (8, 12): 0.034801483}
    print("  d   F   R   p_d(F)         front-only     ratio")
    ratios = {4: [], 6: [], 8: []}
    for d in (4, 6, 8):
        for F in range(3, 13):
            R = 10 if F <= 11 else 9
            p = exact_p(d, F, R); fo = front_only(d, F)
            ratios[d].append(p / fo)
            flag = ''
            if (d, F) in quoted:
                flag = '  (quoted %.9f: %s)' % (quoted[(d, F)], 'match' if abs(p - quoted[(d, F)]) < 1e-9 else 'MISMATCH')
            if d == 4:
                assert abs(p - 2.0 ** (1 - F)) < 1e-15
            print("  %d  %2d  %2d   %.9f    %.9f    %.4f%s" % (d, F, R, p, fo, p / fo, flag))
    print("  ratio ranges over F = 3..12: d=6: %.4f..%.4f ; d=8 (F >= 4): %.4f..%.4f ; d=4: all exactly 1 (p_4(F) = 2^(1-F))"
          % (min(ratios[6]), max(ratios[6]), min(ratios[8][1:]), max(ratios[8][1:])))
    for d in (6, 8):
        vals = [(R, exact_p(d, 8, R)) for R in (6, 8, 10, 12, 14)]
        print("  R-dependence of the 'exact' value, d=%d, F=8: %s  (the quoted 6-digit values are the R=10 truncation)" % (d, ['R=%d: %.7f' % v for v in vals]))
    # Monte Carlo replication: same seed and the same draw shapes as the audited script (its exact table draws nothing)
    rng = np.random.default_rng(20260926)
    out = []
    for F in (14, 16):
        nb = F - 1 + 16; samples = 4000000
        try:
            raw = rng.integers(0, 2, size=(samples, nb))            # same call as the audited script (int64 stream)
            bits = raw.astype(np.uint8); del raw; bits *= 2
            p = loss_fraction(build(bits[:, :F - 1], bits[:, F - 1:], 4, F))
            out.append("F=%d: %.3e (quoted %s) vs 2^(1-F) = %.3e" % (F, p, '1.267e-04' if F == 14 else '2.875e-05', 2.0 ** (1 - F)))
            del bits
        except MemoryError:
            out.append("F=%d: MemoryError at 4e6 samples (skipped)" % F)
    print("  Monte Carlo size 4, 4e6 samples, R = 16, seed 20260926:", out)


# ----------------------------------------------------------------------------------------------------------- H
def part_H():
    hdr("H. prime frontier (primes < 200000)")
    row = primes_below(200000); fronts = []; prevF = None; nrows = 0; risk_rows = 0.0; per_row = []
    for r in range(1, 20000):
        row = np.abs(np.diff(row))
        assert row[0] == 1
        big = np.nonzero(row[1:] >= 4)[0]
        if len(big) == 0:
            break
        nrows += 1
        F = int(big[0]) + 1; d = int(row[F]); per_row.append((r, F, d))
        risk_rows += front_only(d, F)
        if prevF is None or F != prevF - 1:
            fronts.append((r, F, d))
        prevF = F
    sizes = {}
    for _, F, d in fronts:
        sizes[d] = sizes.get(d, 0) + 1
    print("  rows with a defect: %d (note: 64); fresh fronts: %d (note: 29); sizes: %s (note: 23 of size 4)" % (nrows, len(fronts), dict(sorted(sizes.items()))))
    print("  size-4 fresh-front distances:", [F for _, F, d in fronts if d == 4])
    print("  fresh fronts (row, F, d):", fronts)
    risk_fresh = sum(front_only(d, F) for _, F, d in fronts)
    top = sorted(((front_only(d, F), r, F, d) for r, F, d in fronts), reverse=True)[:6]
    print("  front-only risk sum over fresh fronts: %.6f (note: 0.258); over all rows: %.6f (note: 0.258); largest terms: %s" % (risk_fresh, risk_rows, top))
    print("  (1 - 1/4)(1 - 1/128) = %.4f (note: 'about 0.74')" % ((1 - 0.25) * (1 - 1 / 128)))
    assert nrows == 64 and len(fronts) == 29 and sizes[4] == 23
    # rows 1 and 2 of the primes are NOT lone-defect rows: count entries >= 4 in them
    row = primes_below(200000)
    for r in (1, 2):
        row = np.abs(np.diff(row))
        print("  row %d: first 40 entries %s; entries >= 4 in the row: %d (the wall theorem's hypothesis 'every other cell in {0,2}' fails)" % (r, row[:40].tolist(), int((row[1:] >= 4).sum())))


# ----------------------------------------------------------------------------------------------------------- J
def part_J():
    hdr("J. miscellaneous numbers quoted in the note")
    print("  drifts log2(F_k) - 2:", ["%d: %+.4f" % (p, math.log2(p) - 2) for p in (3, 5, 17, 257, 65537)])

    def row_value(n):
        return sum(1 << j for j in range(n + 1) if (j & ~n) == 0)
    known = [3, 5, 17, 257, 65537]
    prods = sorted({math.prod(c) for r in range(6) for c in itertools.combinations(known, r)})
    rows = sorted(row_value(n) for n in range(32))
    print("  rows 0..31 of the single-seed diagram as binary numbers == the 32 products of distinct known Fermat primes: %s; row 32 = %d = 641 * 6700417: %s"
          % (prods == rows, row_value(32), row_value(32) == 641 * 6700417 == FERM(5)))
    assert prods == rows and row_value(32) == 641 * 6700417
    print("  rows 2^k - 1 all ones, value 2^(2^k) - 1 = prod_(i<k) F_i:", [(2 ** k - 1, row_value(2 ** k - 1), math.prod(FERM(i) for i in range(k))) for k in range(1, 6)])
    s = [row_value(n) for n in range(64)]; lead = []; cur = s
    for r in range(11):
        cur = [abs(cur[i + 1] - cur[i]) for i in range(len(cur) - 1)]; lead.append(cur[0])
    print("  leading column of the absolute difference triangle of s_n = prod_(i in bits n) F_i, rows 1..11:", lead, "(note: 2, 0, 8, 8, 16, 0, 80, 64, 32, 32, 896)")
    assert lead == [2, 0, 8, 8, 16, 0, 80, 64, 32, 32, 896]
    print("  F_k in binary is 1 0^(2^k - 1) 1: the zeros of row 2^k ARE the top row of the side-(2^k - 1) zero triangle (the note says the triangle 'opens under row 2^k'; it opens under row 2^k - 1, IN row 2^k)")


def main():
    print("gilbreath_fermat_platonic_20260926_audit.py -- independent auditor run, numpy %s" % np.__version__)
    part_A(); part_B(); part_C(); part_D(); part_E(); part_F(); part_G(); part_H(); part_J()
    print("\nALL AUDIT ASSERTIONS PASSED (t = %.0fs)" % (time.time() - T_START))


if __name__ == '__main__':
    main()
