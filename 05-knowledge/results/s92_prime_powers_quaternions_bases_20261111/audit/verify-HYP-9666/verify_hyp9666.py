#!/usr/bin/env python3
# verify_hyp9666.py - verifier of HYP-9666 after the audit corrections (opus-2026-11-11-S92).
# Own code, written from the definitions; it imports none of the auditors' scripts.
#
#   [Y] front matter (YAML), path tokens, ASCII, stale tokens of the document
#   [N] the fractions quoted in the document
#   [L] class laws for uniform states: the fusion corner in both coordinate systems, rank 3, six squares
#   [F] fusion states of all labelled pointed tournaments of order <= 7 (root = vertex 0): definitions against an
#       explicit path enumeration, the law (3) of Session264, parities, path-mass coordinates, lattices,
#       affine relations modulo 3, 5, 7, residue vectors, class-0 factors of H(A * B)
#   [K] two-port kernels: visit matrices of orders 3, 4, 5, the census of rank-3 kernels, the substitution
#       formula of THM-5133 (A), the module congruences of THM-4756, the collapse modulo 3 with its explicit
#       constant, the Gaussian kernel
#
# Run:  python -B verify_hyp9666.py > verify_hyp9666.out      (about one minute; numba, no cache files)
import sys
import re
import glob
import os
import itertools
from fractions import Fraction

import numpy as np
import yaml
from numba import njit

ROOT = "C:/tmp/wt-s92"
DOC = (ROOT + "/05-knowledge/hypotheses/HYP-9666-seams-of-rank-four-and-the-octonion-law-equidistributed-"
       "families-rank-three-kernels-degenerate-primes.md")
AUD = ROOT + "/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit"

FAILS = []


def P(s=""):
    print(s, flush=True)


def check(name, ok):
    if not ok:
        FAILS.append(name)
    return "OK" if ok else "FAIL"


# ----------------------------------------------------------------------------------------------- numba kernels
@njit
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit
def bit_index(b):
    w = 0
    while b > 1:
        b >>= 1
        w += 1
    return w


@njit
def ham_all(adj, n, Hs):
    # adj[v] = bitmask of the out-neighbours of v.  Hs[S] = number of Hamiltonian paths of T[S]; Hs[0] = 1.
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full):
        for v in range(n):
            c = dp[S, v]
            if c != 0:
                o = adj[v] & ~S
                while o != 0:
                    b = o & -o
                    dp[S | b, bit_index(b)] += c
                    o ^= b
    Hs[0] = 1
    for S in range(1, full):
        t = 0
        for v in range(n):
            t += dp[S, v]
        Hs[S] = t


@njit
def ham_count(adj, n):
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full):
        for v in range(n):
            c = dp[S, v]
            if c != 0:
                o = adj[v] & ~S
                while o != 0:
                    b = o & -o
                    dp[S | b, bit_index(b)] += c
                    o ^= b
    t = 0
    for v in range(n):
        t += dp[full - 1, v]
    return t


@njit
def masses(adj, n, Hs, res):
    # Paths R starting at vertex 0 along adj.  res[0] = sum of H(T - V(R)) over the paths with an odd number of
    # vertices (the one-vertex path included), res[1] = the same sum over the paths with an even number.
    full = 1 << n
    f = np.zeros((full, n), dtype=np.int64)
    f[1, 0] = 1
    odd = 0
    even = 0
    for S in range(1, full, 2):
        cnt = 0
        for v in range(n):
            c = f[S, v]
            if c != 0:
                cnt += c
                o = adj[v] & ~S
                while o != 0:
                    b = o & -o
                    f[S | b, bit_index(b)] += c
                    o ^= b
        if cnt != 0:
            m = cnt * Hs[(full - 1) ^ S]
            if popcount(S) % 2 == 1:
                odd += m
            else:
                even += m
    res[0] = odd
    res[1] = even


@njit
def all_states(n, st):
    # st[code] = (H, h, U, V, d, e, a_o, a_e, b_o, b_e); code bit k (pairs i < j in lexicographic order) = 1: i -> j
    npairs = n * (n - 1) // 2
    full = 1 << n
    Hs = np.zeros(full, dtype=np.int64)
    adj = np.zeros(n, dtype=np.int64)
    radj = np.zeros(n, dtype=np.int64)
    r1 = np.zeros(2, dtype=np.int64)
    r2 = np.zeros(2, dtype=np.int64)
    for code in range(1 << npairs):
        for v in range(n):
            adj[v] = 0
            radj[v] = 0
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> k) & 1:
                    adj[i] |= 1 << j
                    radj[j] |= 1 << i
                else:
                    adj[j] |= 1 << i
                    radj[i] |= 1 << j
                k += 1
        ham_all(adj, n, Hs)
        masses(adj, n, Hs, r1)
        masses(radj, n, Hs, r2)
        H = Hs[full - 1]
        h = Hs[(full - 1) ^ 1]
        st[code, 0] = H
        st[code, 1] = h
        st[code, 2] = r1[0] + r1[1]
        st[code, 3] = r2[0] + r2[1]
        st[code, 4] = r1[0] - r1[1]
        st[code, 5] = r2[0] - r2[1]
        st[code, 6] = r1[0] - h
        st[code, 7] = r1[1]
        st[code, 8] = r2[0] - h
        st[code, 9] = r2[1]


@njit
def fused_hist(n, hist):
    # histogram of H(A * B) modulo 105 over ALL ordered pairs of labelled pointed tournaments of order n,
    # by direct path counting on the fused tournament (order 2n - 1); the fusion law is not used.
    npairs = n * (n - 1) // 2
    N = 2 * n - 1
    adj = np.zeros(N, dtype=np.int64)
    for ca in range(1 << npairs):
        for cb in range(1 << npairs):
            for v in range(N):
                adj[v] = 0
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (ca >> k) & 1:
                        adj[i] |= 1 << j
                    else:
                        adj[j] |= 1 << i
                    bi = 0
                    if i != 0:
                        bi = n + i - 1
                    bj = n + j - 1
                    if (cb >> k) & 1:
                        adj[bi] |= 1 << bj
                    else:
                        adj[bj] |= 1 << bi
                    k += 1
            for i in range(1, n):
                for j in range(1, n):
                    adj[i] |= 1 << (n + j - 1)
            hist[ham_count(adj, N) % 105] += 1


@njit
def ham_sub(K, m, x, y, A, na, B, nb):
    # H(K[A, B]): the vertex x of K replaced by A, the vertex y by B (0/1 adjacency matrices)
    N = m - 2 + na + nb
    owner = np.zeros(N, dtype=np.int64)
    idx = np.zeros(N, dtype=np.int64)
    t = 0
    for v in range(m):
        if v == x:
            for i in range(na):
                owner[t] = v
                idx[t] = i
                t += 1
        elif v == y:
            for i in range(nb):
                owner[t] = v
                idx[t] = i
                t += 1
        else:
            owner[t] = v
            idx[t] = 0
            t += 1
    adj = np.zeros(N, dtype=np.int64)
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            va = owner[a]
            vb = owner[b]
            if va == vb:
                if va == x:
                    arc = A[idx[a], idx[b]]
                else:
                    arc = B[idx[a], idx[b]]
            else:
                arc = K[va, vb]
            if arc != 0:
                adj[a] |= 1 << b
    return ham_count(adj, N)


@njit
def module_congruences(n, do_p3, res):
    # over all labelled tournaments of order n:
    # res[0] = their number; res[1] = violations of Pi_2 = (-1)^n H (mod 3); res[2] = violations of
    # P_2 = ((-1)^n + 1) H (mod 3); res[3] = P_2 odd; res[4] = P_3 not divisible by 6; res[5] = H even
    npairs = n * (n - 1) // 2
    full = 1 << n
    Hs = np.zeros(full, dtype=np.int64)
    adj = np.zeros(n, dtype=np.int64)
    sign = -1
    if n % 2 == 0:
        sign = 1
    for code in range(1 << npairs):
        for v in range(n):
            adj[v] = 0
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> k) & 1:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                k += 1
        ham_all(adj, n, Hs)
        H = Hs[full - 1]
        pi2 = 0
        for S in range(full):
            pi2 += Hs[S] * Hs[(full - 1) ^ S]
        P2 = pi2 - 2 * H
        res[0] += 1
        if (((pi2 - sign * H) % 3) + 3) % 3 != 0:
            res[1] += 1
        if (((P2 - (sign + 1) * H) % 3) + 3) % 3 != 0:
            res[2] += 1
        if P2 % 2 != 0:
            res[3] += 1
        if H % 2 == 0:
            res[5] += 1
        if do_p3:
            P3 = 0
            for S1 in range(1, full - 1):
                R = (full - 1) ^ S1
                p2r = 0
                sub = (R - 1) & R
                while sub != 0:
                    p2r += Hs[sub] * Hs[R ^ sub]
                    sub = (sub - 1) & R
                P3 += Hs[S1] * p2r
            if P3 % 6 != 0:
                res[4] += 1


# ----------------------------------------------------------------------------------------------- python helpers
def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        M = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (code >> k) & 1:
                M[i][j] = 1
            else:
                M[j][i] = 1
        yield code, M


def masks_of(M):
    n = len(M)
    return np.array([sum(1 << j for j in range(n) if M[i][j]) for i in range(n)], dtype=np.int64)


def H_perm(M, verts):
    if len(verts) == 0:
        return 1
    c = 0
    for perm in itertools.permutations(verts):
        ok = True
        for a, b in zip(perm, perm[1:]):
            if not M[a][b]:
                ok = False
                break
        if ok:
            c += 1
    return c


def brute_state(M):
    # the ten numbers from the definitions, by explicit enumeration of the rooted paths (root = vertex 0)
    n = len(M)
    others = list(range(1, n))
    H = H_perm(M, list(range(n)))
    h = H_perm(M, others)
    U = V = d = e = ao = ae = bo = be = 0
    for k in range(0, n):
        for tail in itertools.permutations(others, k):
            rest = [v for v in others if v not in tail]
            # path starting at the root: 0, tail
            seq = (0,) + tail
            if all(M[a][b] for a, b in zip(seq, seq[1:])):
                w = H_perm(M, rest)
                U += w
                d += w if (k + 1) % 2 == 1 else -w
                if k >= 1:
                    if (k + 1) % 2 == 1:
                        ao += w
                    else:
                        ae += w
            # path ending at the root: tail, 0
            seq = tail + (0,)
            if all(M[a][b] for a, b in zip(seq, seq[1:])):
                w = H_perm(M, rest)
                V += w
                e += w if (k + 1) % 2 == 1 else -w
                if k >= 1:
                    if (k + 1) % 2 == 1:
                        bo += w
                    else:
                        be += w
    return (H, h, U, V, d, e, ao, ae, bo, be)


def fuse(A, B):
    na, nb = len(A), len(B)
    N = na + nb - 1
    M = [[0] * N for _ in range(N)]
    for i in range(na):
        for j in range(na):
            M[i][j] = A[i][j]

    def ib(j):
        return 0 if j == 0 else na + j - 1
    for i in range(nb):
        for j in range(nb):
            if i != j and B[i][j]:
                M[ib(i)][ib(j)] = 1
    for i in range(1, na):
        for j in range(1, nb):
            M[i][na + j - 1] = 1
    return M


def hnf_basis(vectors, dim):
    # integer row reduction: a basis (in echelon form) of the lattice spanned by the vectors
    rows = [list(v) for v in vectors]
    basis = []
    for col in range(dim):
        nz = [r for r in rows if r[col] != 0]
        rest = [r for r in rows if r[col] == 0]
        if nz:
            piv = nz[0]
            for r in nz[1:]:
                while r[col] != 0:
                    q = piv[col] // r[col]
                    piv, r = r, [a - q * b for a, b in zip(piv, r)]
                rest.append(r)
            basis.append(piv)
        rows = rest
    return basis


def det_int(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    d = Fraction(1)
    for i in range(n):
        piv = None
        for r in range(i, n):
            if A[r][i] != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return int(d)


def lattice_index(vectors, dim):
    b = hnf_basis(vectors, dim)
    if len(b) < dim:
        return None
    return abs(det_int(b))


def affine_rank(vecs, p):
    vecs = list(vecs)
    v0 = vecs[0]
    basis = []
    for v in vecs[1:]:
        r = [(a - b) % p for a, b in zip(v, v0)]
        for piv, row in basis:
            if r[piv]:
                f = r[piv]
                r = [(x - f * y) % p for x, y in zip(r, row)]
        nz = [i for i, x in enumerate(r) if x]
        if nz:
            piv = nz[0]
            inv = pow(r[piv], -1, p)
            basis.append((piv, [(x * inv) % p for x in r]))
            if len(basis) == len(v0):
                break
    return len(basis)


def rank_law_counts(p, k, r):
    # numbers of (u, w) in (Z/p^k)^r x (Z/p^k)^r with u . w = c, from statement (1) of THM-5135
    q = p ** k
    tot = Fraction(q) ** (2 * r)
    out = []
    for c in range(q):
        if c == 0:
            fac = (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(k)) \
                + Fraction(1, p ** ((r - 1) * k))
        else:
            v = 0
            cc = c
            while cc % p == 0:
                cc //= p
                v += 1
            fac = (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(v + 1))
        cnt = fac * tot / q
        assert cnt.denominator == 1
        out.append(int(cnt))
    return out


def brute_law(G, m):
    r = len(G)
    Gm = np.array(G, dtype=np.int64) % m
    grid = np.array(list(itertools.product(range(m), repeat=r)), dtype=np.int64)
    Gw = (grid @ Gm.T) % m
    cnt = np.zeros(m, dtype=np.int64)
    for u in grid:
        cnt += np.bincount((Gw @ u) % m, minlength=m)
    return [int(x) for x in cnt]


def conv(a, b, q):
    out = [0] * q
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[(i + j) % q] += ai * bj
    return out


def visit_matrix(K, m, x, y):
    # N[i][j] = number of walks of K visiting x exactly i times, y exactly j times and every other vertex once.
    # No cap is imposed: a gap between two visits of a port contains another vertex, so the walks are finite.
    size = m + 2
    N = [[0] * size for _ in range(size)]
    others = ((1 << m) - 1) & ~(1 << x) & ~(1 << y)
    out = [[w for w in range(m) if K[v][w]] for v in range(m)]

    def rec(v, mask, cx, cy):
        if mask == others:
            N[cx][cy] += 1
        for w in out[v]:
            if w == x:
                rec(w, mask, cx + 1, cy)
            elif w == y:
                rec(w, mask, cx, cy + 1)
            elif not (mask >> w) & 1:
                rec(w, mask | (1 << w), cx, cy)

    for v in range(m):
        if v == x:
            rec(v, 0, 1, 0)
        elif v == y:
            rec(v, 0, 0, 1)
        else:
            rec(v, 1 << v, 0, 0)
    return N


def cover_numbers(M, imax):
    # P_i = number of ordered i-tuples of disjoint non-empty paths covering the tournament
    n = len(M)
    full = (1 << n) - 1
    Hsub = {}
    for S in range(1, full + 1):
        verts = [v for v in range(n) if (S >> v) & 1]
        Hsub[S] = H_perm(M, verts)
    cur = {0: 1}
    Pn = [0] * (imax + 1)
    for i in range(1, imax + 1):
        nxt = {}
        for mask, c in cur.items():
            rest = full & ~mask
            sub = rest
            while sub:
                nxt[mask | sub] = nxt.get(mask | sub, 0) + c * Hsub[sub]
                sub = (sub - 1) & rest
        cur = nxt
        Pn[i] = cur.get(full, 0)
    return Pn


def class_reps(n):
    seen = {}
    for code, T in tournaments(n):
        best = min(tuple(T[p[i]][p[j]] for i in range(n) for j in range(n))
                   for p in itertools.permutations(range(n)))
        if best not in seen:
            seen[best] = T
    return list(seen.values())


def np_adj(M):
    return np.array(M, dtype=np.int8)


# ----------------------------------------------------------------------------------------------- [Y]
def section_Y():
    P("[Y] front matter and text of the document")
    text = open(DOC, encoding="utf-8").read()
    ok_start = text.startswith("---\n")
    end = text.index("\n---\n", 4)
    fm = yaml.safe_load(text[4:end])
    body = text[end + 5:]
    keys = list(fm.keys())
    P("  YAML front matter parses: %s; keys: %s" % (check("yaml", ok_start and isinstance(fm, dict)), keys))
    need = ["id", "title", "status", "session", "date", "source", "scripts", "related"]
    P("  required keys present: %s; id = %s; title and status are strings: %s; scripts %d items, related %d items"
      % (check("keys", all(k in fm for k in need)), fm.get("id"),
         isinstance(fm.get("title"), str) and isinstance(fm.get("status"), str),
         len(fm.get("scripts", [])), len(fm.get("related", []))))
    bad = []
    npath = 0
    for entry in [fm["source"]] + list(fm["scripts"]) + list(fm["related"]):
        if not isinstance(entry, str):
            bad.append(repr(entry))
            continue
        tok = entry.split(" ")[0].rstrip(",")
        npath += 1
        if "*" in tok:
            hit = glob.glob(ROOT + "/" + tok)
            if len(hit) != 1:
                bad.append(tok + " (%d matches)" % len(hit))
        else:
            if not os.path.exists(ROOT + "/" + tok):
                bad.append(tok)
            if "(+ .out" in entry:
                npath += 1
                if not os.path.exists(ROOT + "/" + tok[:-3] + ".out"):
                    bad.append(tok[:-3] + ".out")
    P("  path tokens of source / scripts / related resolving to existing files: %d of %d %s %s"
      % (npath - len(bad), npath, check("paths", not bad), bad if bad else ""))
    toks = sorted(set(re.findall(r"audit/[A-Za-z0-9_\-]+/[A-Za-z0-9_\-\.]+\.(?:out|py|txt|md)", body)))
    badb = [t for t in toks if not os.path.exists(AUD + "/" + t[len("audit/"):])]
    P("  receipts named audit/... in the body: %d distinct, missing: %s %s" % (len(toks), badb, check("body paths", not badb)))
    non_ascii = sorted(set(ch for ch in text if ord(ch) > 127))
    P("  non-ASCII characters: %s %s" % (non_ascii, check("ascii", not non_ascii)))
    stale = {}
    for pat in ["PRE-AUDIT", "sigma(", "THM-5134 (H)", "falsifier", "Falsifier", "(4)(d)", "X1", "X2"]:
        stale[pat] = text.count(pat)
    P("  occurrences of stale tokens: %s %s" % (stale, check("stale", not any(stale.values()))))
    P("  mentions of Chan-Pak: %d, of Collatz: %d" % (text.count("Chan"), text.count("Collatz")))
    labs = sorted(set(re.findall(r"THM-5135 \([0-9]\)(?:\([a-z]\))?", text)))
    P("  labels of THM-5135 cited: %s" % labs)
    labs2 = sorted(set(re.findall(r"THM-5134 \([A-Z][0-9]?\)", text)))
    P("  labels of THM-5134 cited: %s" % labs2)
    long_lines = [i + 1 for i, ln in enumerate(body.split("\n")) if len(ln) > 120]
    P("  body lines longer than 120 characters: %s" % long_lines)
    fnd = re.search(r"findings ([F0-9, to]+) mention this file", text)
    P("  findings line: %s" % (fnd.group(1) if fnd else "absent"))
    P()


# ----------------------------------------------------------------------------------------------- [N]
def section_N():
    P("[N] fractions quoted in the document")
    for p, want, bias in ((3, Fraction(83, 81), Fraction(2, 81)), (5, Fraction(629, 625), Fraction(4, 625))):
        f0 = 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)
        P("  p = %d: 1 + p^-3 - p^-4 = %s (expected %s), bias %s (expected %s) %s"
          % (p, f0, want, f0 - 1, bias, check("f0 %d" % p, f0 == want and f0 - 1 == bias)))
    vals = [float(1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)) for p in (3, 5, 7)]
    P("  uniform values at 3, 5, 7: %.4f %.4f %.4f (document: 1.0247, 1.0064, 1.0025) %s"
      % (vals[0], vals[1], vals[2],
         check("uniform decimals", ["%.4f" % v for v in vals] == ["1.0247", "1.0064", "1.0025"])))
    for p, want in ((3, Fraction(34, 27)), (5, Fraction(151, 125))):
        t = 1 + Fraction((3 * p - 2) * (p - 1), 2 * p ** 3)
        P("  THM-4916, two root changes, p = %d: 1 + (3p - 2)(p - 1)/(2 p^3) = %s (expected %s) %s; p (L(0) - 1) = %s"
          % (p, t, want, check("4916 %d" % p, t == want), (t - 1) * p))
    P("  p (L(0) - 1) at p = 10007: %.5f (tends to 3/2)"
      % float(Fraction((3 * 10007 - 2) * (10007 - 1), 2 * 10007 ** 2)))
    c = rank_law_counts(2, 1, 4)
    P("  rank 4 modulo 2: share of the class 0 = %s (expected 17/32) %s"
      % (Fraction(c[0], 2 ** 8), check("17/32", Fraction(c[0], 2 ** 8) == Fraction(17, 32))))
    P("  0.0004 against the rank-4 biases p^-3 - p^-4 at 5, 7, 11, 13: %s"
      % ", ".join("%.5f" % float(Fraction(1, p ** 3) - Fraction(1, p ** 4)) for p in (5, 7, 11, 13)))
    P()


# ----------------------------------------------------------------------------------------------- [L]
def section_L():
    P("[L] class laws for uniform states")
    G = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    Gm = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]
    P("  fusion corner in (h, H-h, U-h, d-h): det = %d; in the path-mass coordinates (h, H-h, a_o, a_e): det = %d %s"
      % (det_int(G), det_int(Gm), check("dets", det_int(G) == -1 and det_int(Gm) == -4)))
    res = []
    for m, (p, k) in ((2, (2, 1)), (3, (3, 1)), (4, (2, 2)), (5, (5, 1)), (9, (3, 2))):
        same = brute_law(G, m) == rank_law_counts(p, k, 4)
        res.append("%d:%s" % (m, same))
        check("corner law %d" % m, same)
    P("  corner, uniform states, enumeration of all pairs, equal to the rank-4 law (1): %s" % "  ".join(res))
    res = []
    for m, (p, k) in ((2, (2, 1)), (3, (3, 1)), (4, (2, 2)), (5, (5, 1)), (9, (3, 2))):
        same = brute_law(Gm, m) == rank_law_counts(p, k, 4)
        res.append("%d:%s" % (m, same))
        check("mass law %d" % m, same == (p != 2))
    P("  the same form in path-mass coordinates (determinant -4): %s (unimodular at the odd primes only)"
      % "  ".join(res))
    # K-form of THM-4916 (1d) in the coordinates (U, E, d, a): U U' + E E' - d d' + a d' + a' d
    GK = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 1], [0, 0, 1, 0]]
    res = []
    for m, (p, k) in ((2, (2, 1)), (3, (3, 1)), (4, (2, 2)), (5, (5, 1)), (9, (3, 2))):
        same = brute_law(GK, m) == rank_law_counts(p, k, 4)
        res.append("%d:%s" % (m, same))
        check("K-form law %d" % m, same)
    P("  K-form in (U, E, d, a): det = %d %s; uniform states, equal to the rank-4 law (1): %s"
      % (det_int(GK), check("K-form det", det_int(GK) == -1), "  ".join(res)))
    # rank 3 and six squares
    for p, k in ((3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2), (13, 1)):
        q = p ** k
        tab = [0] * q
        sq = [0] * q
        for a in range(q):
            sq[(a * a) % q] += 1
            for b in range(q):
                tab[(a * b) % q] += 1
        law3 = conv(conv(tab, tab, q), tab, q)
        okf = law3 == rank_law_counts(p, k, 3)
        s2 = conv(sq, sq, q)
        s6 = conv(conv(s2, s2, q), s2, q)
        same6 = s6 == law3
        P("  rank 3 modulo %2d: factor (1 - p^-3) sigma_2(p^v)/p^(2v), with the class 0: %s; six squares have the same law: %s (expected %s) %s"
          % (q, okf, same6, p % 4 == 1, check("rank3 %d" % q, okf and same6 == (p % 4 == 1))))
    P()


# ----------------------------------------------------------------------------------------------- [F]
def section_F():
    P("[F] fusion states of the labelled pointed tournaments (root = vertex 0), orders 1 to 7")
    st = {}
    for n in range(1, 8):
        a = np.zeros((1 << (n * (n - 1) // 2), 10), dtype=np.int64)
        all_states(n, a)
        st[n] = a
    P("  numbers of tournaments: %s" % [len(st[n]) for n in range(1, 8)])
    # definitions against explicit path enumeration, n <= 5
    bad = 0
    tot = 0
    for n in range(1, 6):
        for code, M in tournaments(n):
            tot += 1
            if tuple(int(x) for x in st[n][code]) != brute_state(M):
                bad += 1
    P("  (H, h, U, V, d, e, a_o, a_e, b_o, b_e) against explicit enumeration of rooted paths, orders <= 5: %d tournaments, %d differences %s"
      % (tot, bad, check("brute states", bad == 0)))
    tri = sorted(set(tuple(int(x) for x in st[3][c][:6]) for c in range(8)))
    P("  order 3, (H, h, U, V, d, e): %s (Session264: cyclic triangle (3, 1, 3, 3, 1, 1))" % tri)
    # identities of the path-mass coordinates
    bad = 0
    for n in range(1, 8):
        H, h, U, V, d, e, ao, ae, bo, be = (st[n][:, i] for i in range(10))
        bad += int(np.sum(2 * ao != U + d - 2 * h)) + int(np.sum(2 * ae != U - d))
        bad += int(np.sum(2 * bo != V + e - 2 * h)) + int(np.sum(2 * be != V - e))
        bad += int(np.sum(U - h != ao + ae)) + int(np.sum(d - h != ao - ae))
        bad += int(np.sum(ao < 0)) + int(np.sum(ae < 0)) + int(np.sum(bo < 0)) + int(np.sum(be < 0))
    P("  a_o = (U + d)/2 - h, a_e = (U - d)/2, b_o = (V + e)/2 - h, b_e = (V - e)/2 with the path masses computed directly, orders <= 7: %d violations %s"
      % (bad, check("mass identities", bad == 0)))
    # the fusion law (3) by direct counting
    tot = 0
    bad = 0
    bad2 = 0
    for na in range(1, 5):
        for ca, A in tournaments(na):
            sa = [int(x) for x in st[na][ca]]
            for nb in range(1, 5):
                for cb, B in tournaments(nb):
                    sb = [int(x) for x in st[nb][cb]]
                    M = fuse(A, B)
                    true = int(ham_count(masks_of(M), len(M)))
                    H, h, U, V, d, e, ao, ae = sa[:8]
                    H2, h2, U2, V2, d2, e2, _, _, bo2, be2 = sb
                    pred = h * h2 + (H - h) * h2 + h * (H2 - h2) + (U - h) * (V2 - h2) + (d - h) * (e2 - h2)
                    pred2 = H * h2 + h * H2 - h * h2 + 2 * (ao * bo2 + ae * be2)
                    tot += 1
                    bad += (true != pred)
                    bad2 += (true != pred2)
    P("  fusion law (3) by direct path counts of A * B, all pairs of orders <= 4: %d pairs, %d violations; with the last two terms written 2 (a_o b_o + a_e b_e): %d violations %s"
      % (tot, bad, bad2, check("fusion law", tot == 5625 and bad == 0 and bad2 == 0)))
    # parities
    for n in range(1, 8):
        H, h, U, V, d, e, ao, ae, bo, be = (st[n][:, i] for i in range(10))
        okp = bool(np.all(H % 2 == 1) and np.all(h % 2 == 1) and np.all((U - d) % 2 == 0) and np.all((V - e) % 2 == 0))
        cl = sorted(set(zip((ao % 2).tolist(), (ae % 2).tolist())))
        cr = sorted(set(zip((bo % 2).tolist(), (be % 2).tolist())))
        pr4 = sorted(set(zip(((U - h) % 4).tolist(), ((d - h) % 4).tolist())))
        P("  n = %d: H, h odd and U = d, V = e (mod 2): %s; (a_o, a_e) mod 2 takes %d classes, (b_o, b_e) mod 2 %d; ((U-h) mod 4, (d-h) mod 4): %d pairs"
          % (n, check("parity %d" % n, okp), len(cl), len(cr), len(pr4)))
    # lattices
    left = set()
    right = set()
    for n in range(1, 8):
        H, h, U, V, d, e, ao, ae, bo, be = (st[n][:, i] for i in range(10))
        left |= set(zip(h.tolist(), (H - h).tolist(), (U - h).tolist(), (d - h).tolist()))
        right |= set(zip(h.tolist(), (H - h).tolist(), (V - h).tolist(), (e - h).tolist()))
    P("  distinct left states (orders <= 7): %d; distinct right states: %d; equal sets: %s" % (len(left), len(right), left == right))
    idx2 = lattice_index(set((v[2], v[3]) for v in left), 2)
    idx2m = lattice_index(set(((v[2] + v[3]) // 2, (v[2] - v[3]) // 2) for v in left), 2)
    idx4 = lattice_index(left, 4)
    v0 = (1, 0, 0, 0)
    idx4d = lattice_index(set(tuple(a - b for a, b in zip(v, v0)) for v in left if v != v0), 4)
    P("  Z-span of the pairs (U - h, d - h): index %s in Z^2 (expected 2); Z-span of the pairs (a_o, a_e): index %s (expected 1) %s"
      % (idx2, idx2m, check("index 2", idx2 == 2 and idx2m == 1)))
    P("  Z-span of the left states: index %s in Z^4 (the conditions H - h even and U = d mod 2); Z-span of their differences: index %s (h odd fixes a coset)"
      % (idx4, idx4d))
    # the two witnesses quoted in the introduction of the document
    arc = brute_state([[0, 1], [0, 0]])
    tri = brute_state([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    w1 = (arc[2] - arc[1], arc[4] - arc[1])
    w2 = (tri[2] - tri[1], tri[4] - tri[1])
    P("  arc 0 -> 1 rooted at its tail: (U - h, d - h) = %s, (a_o, a_e) = %s; cyclic triangle: (U - h, d - h) = %s, (a_o, a_e) = %s; index of the lattice they span: %s %s"
      % (w1, (arc[6], arc[7]), w2, (tri[6], tri[7]), lattice_index([w1, w2], 2),
         check("witnesses", w1 == (1, -1) and w2 == (2, 0) and lattice_index([w1, w2], 2) == 2)))
    # affine relations and residue vectors
    G = np.array([[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.int64)
    expect7 = {3: 81, 5: 591, 7: 1384}
    hist = {}
    for p in (3, 5, 7):
        for n in range(4, 8):
            H, h, U, V, d, e, ao, ae, bo, be = (st[n][:, i] for i in range(10))
            L = [x % p for x in (h, H - h, U - h, d - h)]
            R = [x % p for x in (h, H - h, V - h, e - h)]
            kL = ((L[0] * p + L[1]) * p + L[2]) * p + L[3]
            kR = ((R[0] * p + R[1]) * p + R[2]) * p + R[3]
            cL = np.bincount(kL, minlength=p ** 4)
            cR = np.bincount(kR, minlength=p ** 4)
            hist[(p, n)] = (cL, cR)
            out = []
            for name, cnt in (("left", cL), ("right", cR)):
                occ = np.flatnonzero(cnt)
                vecs = [((i // p ** 3) % p, (i // p ** 2) % p, (i // p) % p, i % p) for i in occ.tolist()]
                rk = affine_rank(vecs, p)
                out.append("%s: %d of %d residue vectors, affine rank %d" % (name, len(vecs), p ** 4, rk))
                check("affine %d %d %s" % (p, n, name), rk == 4)
                if n == 7:
                    check("occ %d %s" % (p, name), len(vecs) == expect7[p])
            P("  p = %d, n = %d: %s (rank 4 = no affine relation)" % (p, n, "; ".join(out)))
    # class-0 factors of H(A * B)
    doc = {3: ["1.31", "1.30", "1.27", "1.16"], 5: ["1.26", "1.18", "1.13", "1.05"], 7: ["0.25", "0.71", "0.92", "0.993"]}
    facs = {}
    for p in (3, 5, 7):
        idx = np.arange(p ** 4)
        comp = np.stack([(idx // p ** 3) % p, (idx // p ** 2) % p, (idx // p) % p, idx % p], axis=1)
        row = []
        for n in range(4, 8):
            cL, cR = hist[(p, n)]
            law = np.zeros(p, dtype=np.float64)
            for i in np.flatnonzero(cL).tolist():
                w = comp[i] @ G
                val = (comp @ w) % p
                law += float(cL[i]) * np.bincount(val, weights=cR.astype(np.float64), minlength=p)
            tot = float(cL.sum()) * float(cR.sum())
            fac = p * law / tot
            facs[(p, n)] = fac
            row.append(fac[0])
        shown = ["%.2f" % x for x in row[:3]] + ["%.3f" % row[3] if p == 7 else "%.2f" % row[3]]
        P("  class-0 factor of H(A * B) modulo %d, A and B uniform of order n = 4, 5, 6, 7 (from the state histograms and the law): %s; rounded %s; document %s %s"
          % (p, " ".join("%.4f" % x for x in row), shown, doc[p], check("factors %d" % p, shown == doc[p])))
    for n in (4, 5):
        hh = np.zeros(105, dtype=np.int64)
        fused_hist(n, hh)
        tot = int(hh.sum())
        out = []
        for p in (3, 5, 7):
            c0 = sum(int(hh[r]) for r in range(105) if r % p == 0)
            f = p * c0 / tot
            out.append("%d: %.4f" % (p, f))
            check("direct %d %d" % (n, p), abs(f - facs[(p, n)][0]) < 1e-9)
        P("  the same factors by direct path counts of all %d fused pairs of order %d (no use of the law): %s"
          % (tot, n, ", ".join(out)))
    P()


# ----------------------------------------------------------------------------------------------- [K]
def section_K():
    P("[K] two-port kernels")
    # module congruences
    for n in range(1, 8):
        res = np.zeros(6, dtype=np.int64)
        module_congruences(n, n <= 6, res)
        P("  n = %d: %8d labelled tournaments; violations of Pi_2 = (-1)^n H (mod 3): %d; of P_2 = ((-1)^n + 1) H (mod 3): %d; P_2 odd: %d; H even: %d; P_3 not divisible by 6: %s %s"
          % (n, res[0], res[1], res[2], res[3], res[5], res[4] if n <= 6 else "not computed",
             check("module %d" % n, res[1] == 0 and res[2] == 0 and res[3] == 0 and res[4] == 0 and res[5] == 0)))
    # visit matrices
    VM = {}
    maxvis = {}
    for m in (3, 4, 5):
        mv = 0
        for code, K in tournaments(m):
            for x in range(m):
                for y in range(m):
                    if x != y:
                        N = visit_matrix(K, m, x, y)
                        VM[(m, code, x, y)] = N
                        for i in range(m + 2):
                            for j in range(m + 2):
                                if N[i][j]:
                                    mv = max(mv, i, j)
        maxvis[m] = mv
    P("  kernels (labelled tournament, ordered port pair): %s; largest number of visits of a port: %s (bound m - 1)"
      % ([sum(1 for k in VM if k[0] == m) for m in (3, 4, 5)], [maxvis[m] for m in (3, 4, 5)]))
    # census of the rank-3 kernels
    for m in (4, 5):
        n_all = n_le2 = n_some3 = n_both3 = n_four = n_det = 0
        dets = {}
        for (mm, code, x, y), N in VM.items():
            if mm != m:
                continue
            n_all += 1
            mx = max(i for i in range(m + 2) for j in range(m + 2) if N[i][j])
            my = max(j for i in range(m + 2) for j in range(m + 2) if N[i][j])
            if mx <= 2 and my <= 2:
                n_le2 += 1
            elif mx <= 3 and my <= 3:
                n_some3 += 1
                G3 = [[N[i][j] for j in (1, 2, 3)] for i in (1, 2, 3)]
                dd = det_int(G3)
                if mx == 3 and my == 3:
                    n_both3 += 1
                    if dd != 0:
                        n_det += 1
                    dets[abs(dd)] = dets.get(abs(dd), 0) + 1
                else:
                    if dd != 0:
                        n_det += 10 ** 6
            else:
                n_four += 1
        P("  order %d: %d kernels; at most two visits %d; some port three times and none four %d, of which both ports three times %d, with det N != 0: %d; a port four times: %d"
          % (m, n_all, n_le2, n_some3, n_both3, n_det, n_four))
        P("     |det N| over the kernels with both ports visited three times (count): %s" % sorted(dets.items()))
        if m == 4:
            check("census 4", (n_all, n_le2, n_some3, n_both3, n_det, n_four) == (768, 720, 48, 48, 48, 0)
                  and sorted(dets) == [32])
        else:
            check("census 5", (n_all, n_some3, n_both3, n_det) == (20480, 5760, 3840, 3840)
                  and sorted(dets) == [6, 14, 32, 128, 192, 228, 242, 248])
    # the example of the document: arcs 0->3, 1->0, 2->0, 2->1, 3->1, 3->2, ports (0, 3)
    arcs = [(0, 3), (1, 0), (2, 0), (2, 1), (3, 1), (3, 2)]
    Kex = [[0] * 4 for _ in range(4)]
    for a, b in arcs:
        Kex[a][b] = 1
    okT = all(Kex[i][j] + Kex[j][i] == 1 for i in range(4) for j in range(4) if i != j)
    N = visit_matrix(Kex, 4, 0, 3)
    G3 = [[N[i][j] for j in (1, 2, 3)] for i in (1, 2, 3)]
    P("  example: arcs %s form a tournament: %s; ports (0, 3): N = %s, det = %d %s"
      % (arcs, okT, G3, det_int(G3), check("example", okT and G3 == [[5, 3, 0], [3, 7, 2], [0, 2, 2]] and det_int(G3) == 32)))
    Ng = visit_matrix(Kex, 4, 1, 2)
    G2 = [[Ng[i][j] for j in (1, 2)] for i in (1, 2)]
    mxg = max(max(i, j) for i in range(6) for j in range(6) if Ng[i][j])
    P("  the same tournament with the ports (1, 2): at most %d visits, N = %s, det = %d (the Gaussian kernel [[5, 2], [2, 1]]) %s"
      % (mxg, G2, det_int(G2), check("gauss N", mxg == 2 and G2 == [[5, 2], [2, 1]])))
    # substitution formula of THM-5133 (A)
    small = []
    for n in (1, 2, 3):
        for code, T in tournaments(n):
            small.append((n, np_adj(T), cover_numbers(T, 4)))
    tot = 0
    bad = 0
    T4 = {code: np_adj(T) for code, T in tournaments(4)}
    T5 = {code: np_adj(T) for code, T in tournaments(5)}
    for (mm, code, x, y), N in VM.items():
        if mm != 4:
            continue
        K = T4[code]
        for (na, A, PA) in small:
            for (nb, B, PB) in small:
                pred = sum(N[i][j] * PA[i] * PB[j] for i in range(1, 5) for j in range(1, 5))
                tot += 1
                if int(ham_sub(K, 4, x, y, A, na, B, nb)) != pred:
                    bad += 1
    P("  H(K[A, B]) = sum N_ij P_i(A) P_j(B): all 768 kernels of order 4, A and B over all labelled tournaments of orders 1 to 3: %d cases, %d violations %s"
      % (tot, bad, check("substitution 4", bad == 0 and tot == 768 * 121)))
    reps = []
    for n in (1, 2, 3, 4):
        for T in class_reps(n):
            reps.append((n, np_adj(T), cover_numbers(T, 4), int(ham_count(masks_of(T), n))))
    P("  isomorphism classes of orders 1 to 4: %d, with (order, H) = %s" % (len(reps), [(n, h) for n, _, _, h in reps]))
    tot = 0
    bad = 0
    for (mm, code, x, y), N in VM.items():
        if mm != 5:
            continue
        mx = max(i for i in range(7) for j in range(7) if N[i][j])
        my = max(j for i in range(7) for j in range(7) if N[i][j])
        if not (mx == 3 and my == 3):
            continue
        K = T5[code]
        for (na, A, PA, hA) in reps:
            if na > 3:
                continue
            for (nb, B, PB, hB) in reps:
                if nb > 3:
                    continue
                pred = sum(N[i][j] * PA[i] * PB[j] for i in range(1, 4) for j in range(1, 4))
                tot += 1
                if int(ham_sub(K, 5, x, y, A, na, B, nb)) != pred:
                    bad += 1
    P("  the same with the 3 x 3 matrix only, for the 3840 kernels of order 5 with both ports visited three times, A and B over the classes of orders 1 to 3: %d cases, %d violations %s"
      % (tot, bad, check("substitution 5", bad == 0 and tot == 3840 * 16)))
    # collapse modulo 3 with the explicit constant
    TT = {3: {code: np_adj(T) for code, T in tournaments(3)}, 4: T4, 5: T5}
    for m in (3, 4, 5):
        tot = 0
        bad = 0
        nk = 0
        informative = 0
        for (mm, code, x, y), N in VM.items():
            if mm != m:
                continue
            nk += 1
            K = TT[m][code]
            for (na, A, PA, hA) in reps:
                sA = 2 if na % 2 == 0 else 0
                for (nb, B, PB, hB) in reps:
                    if m - 2 + na + nb > 10:
                        continue
                    sB = 2 if nb % 2 == 0 else 0
                    c = (N[1][1] + N[1][2] * sB + N[2][1] * sA + N[2][2] * sA * sB) % 3
                    tot += 1
                    if (hA * hB) % 3 != 0:
                        informative += 1
                    if int(ham_sub(K, m, x, y, A, na, B, nb)) % 3 != (c * hA * hB) % 3:
                        bad += 1
        P("  order %d: %d kernels; H(K[A, B]) = c H(A) H(B) (mod 3) with c = N_11 + N_12 s_B + N_21 s_A + N_22 s_A s_B, s_X = (-1)^|X| + 1, A and B over the 8 classes of orders 1 to 4, |K[A, B]| <= 10: %d cases (%d with 3 not dividing H(A) H(B)), %d violations %s"
          % (m, nk, tot, informative, bad, check("collapse %d" % m, bad == 0)))
    # the Gaussian kernel on all pairs of labelled tournaments of orders <= 4
    lab = []
    for n in (1, 2, 3, 4):
        for code, T in tournaments(n):
            Pn = cover_numbers(T, 2)
            lab.append((n, np_adj(T), Pn[1], Pn[2] + 2 * Pn[1]))
    Kg = np_adj(Kex)
    tot = bad = badm = odd = odd_div = 0
    for (na, A, hA, piA) in lab:
        for (nb, B, hB, piB) in lab:
            v = int(ham_sub(Kg, 4, 1, 2, A, na, B, nb))
            tot += 1
            bad += (v != piA * piB + hA * hB)
            badm += (v % 3 != (((-1) ** (na + nb) + 1) * hA * hB) % 3)
            if (na + nb) % 2 == 1:
                odd += 1
                odd_div += (v % 3 == 0)
    P("  Gaussian kernel (ports (1, 2)): H(K[A, B]) = Pi_2(A) Pi_2(B) + H(A) H(B) on ALL %d pairs of labelled tournaments of orders <= 4: %d violations; reduction modulo 3: %d violations; |A| + |B| odd: %d pairs, 3 | H(K[A, B]) in %d %s"
      % (tot, bad, badm, odd, odd_div, check("gauss", tot == 5625 and bad == 0 and badm == 0 and odd == 1188 and odd_div == 1188)))
    P("  the audit's sample: 5625 - 4096 + 600 = %d pairs, of which %d with |A| + |B| odd (none of them in the sampled block)"
      % (5625 - 4096 + 600, odd))
    P()


def main():
    P("# verify_hyp9666.py (verifier of HYP-9666; own code)")
    P()
    section_Y()
    section_N()
    section_L()
    section_F()
    section_K()
    P("FAILED CHECKS: %s" % (FAILS if FAILS else "none"))


if __name__ == "__main__":
    main()
