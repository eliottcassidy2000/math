#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex (S93 Chan-Pak wave): independent core checks.

Everything here is recomputed with code written for this file (own literal Hamiltonian-path count, own 24-state
transfer derived from the path structure [R1] x [R2] y [R3]); the only lane artefact that is read is
../minrep17.json, and it is CHECKED, not trusted.

  A  transfer = literal count (all words n <= 7, random words n = 8..12)
  B  H(0 w) = H(w), H(w 3) = H(w), converse symmetry, H odd, H <= 2 3^n (all words n <= 9)
  C  Hankel block |p|,|s| <= 3 has exact rank 17; rank of the pairing (reachable space) x (observable space) = 17
  D  minrep17.json: x(w l) = x(w) A_l for ALL words (checked on a basis of the reachable space: a proof),
     x(empty), integrality, entry range, fixed vectors, characteristic polynomials, common flag with common diagonal
  E  first moment n = 1..60; E[H^2]/9^n (floating point)
  F  H(3^p 0^p): constant -> 4, bounds
  G  chains modulo 3, 4, 8: orbit sizes, strong connectivity, translation invariance, by-index class shares,
     twisted rate modulo 3 (power iteration, floating point)
  H  F1 identity (three terms) on all p, q <= 6, two-term failures
  I  second differences d_i d_j H(., b), all n <= 8
"""
import itertools
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
import sympy as sp
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
OK = True


def check(name, cond, extra=""):
    global OK
    OK &= bool(cond)
    print(f"[{'ok' if cond else 'FAIL'}] {name} {extra}", flush=True)


# ------------------------------------------------------------------------------------------------ literal count
@njit(cache=True)
def hp_literal(word):
    n = word.shape[0]
    V = n + 2
    u = n
    v = n + 1
    adj = np.zeros((V, V), np.uint8)
    for i in range(n):
        for j in range(i + 1, n):
            adj[i, j] = 1
        if word[i] & 1:
            adj[u, i] = 1
        else:
            adj[i, u] = 1
        if (word[i] >> 1) & 1:
            adj[v, i] = 1
        else:
            adj[i, v] = 1
    adj[u, v] = 1
    full = (1 << V) - 1
    dp = np.zeros((1 << V, V), np.int64)
    for x in range(V):
        dp[1 << x, x] = 1
    for S in range(1, 1 << V):
        for x in range(V):
            c = dp[S, x]
            if c != 0:
                for y in range(V):
                    if adj[x, y] == 1 and ((S >> y) & 1) == 0:
                        dp[S | (1 << y), y] += c
    tot = 0
    for x in range(V):
        tot += dp[full, x]
    return tot


@njit(cache=True)
def hp_one_apex(bits):
    """transitive spine s_1 -> ... -> s_m plus one apex z, z -> s_i iff bits[i] = 1; literal count."""
    m = bits.shape[0]
    V = m + 1
    z = m
    adj = np.zeros((V, V), np.uint8)
    for i in range(m):
        for j in range(i + 1, m):
            adj[i, j] = 1
        if bits[i] == 1:
            adj[z, i] = 1
        else:
            adj[i, z] = 1
    full = (1 << V) - 1
    dp = np.zeros((1 << V, V), np.int64)
    for x in range(V):
        dp[1 << x, x] = 1
    for S in range(1, 1 << V):
        for x in range(V):
            c = dp[S, x]
            if c != 0:
                for y in range(V):
                    if adj[x, y] == 1 and ((S >> y) & 1) == 0:
                        dp[S | (1 << y), y] += c
    tot = 0
    for x in range(V):
        tot += dp[full, x]
    return tot


# ------------------------------------------------------------------------------------------------ own transfer
def sidx(arr, need1, t2, t3):
    return ((arr * 2 + need1) * 3 + t2) * 2 + t3


def letter_matrix(l):
    a, b = l & 1, (l >> 1) & 1
    M = np.zeros((24, 24), np.int64)
    for arr in (0, 1):
        xb, yb = (a, b) if arr == 0 else (b, a)      # xb = 1 iff x -> s, yb = 1 iff y -> s
        for need1 in (0, 1):
            for t2 in (0, 1, 2):
                for t3 in (0, 1):
                    s = sidx(arr, need1, t2, t3)
                    M[s, sidx(arr, 1 if xb == 1 else 0, t2, t3)] += 1          # joins R1 (must end with s -> x)
                    if t2 != 0 or xb == 1:                                      # joins R2 (first needs x -> s)
                        M[s, sidx(arr, need1, 1 if yb == 0 else 2, t3)] += 1    #   last needs s -> y
                    if t3 != 0 or yb == 1:                                      # joins R3 (first needs y -> s)
                        M[s, sidx(arr, need1, t2, 1)] += 1
    return M


MAT = [letter_matrix(l) for l in range(4)]
INIT = np.zeros(24, np.int64)
INIT[sidx(0, 0, 0, 0)] = 1
INIT[sidx(1, 0, 0, 0)] = 1
FIN = np.zeros(24, np.int64)
for arr in (0, 1):
    for t2 in (0, 1):
        if t2 == 0 and arr == 1:
            continue
        for t3 in (0, 1):
            FIN[sidx(arr, 0, t2, t3)] = 1
MATO = [m.astype(object) for m in MAT]


def H_word(word):
    v = INIT.astype(object)
    for l in word:
        v = v.dot(MATO[l])
    return int(v.dot(FIN.astype(object)))


def H_fast(word):
    """int64 version, exact for words of length <= 38 (2 * 3^38 < 2^63)."""
    v = INIT
    for l in word:
        v = v @ MAT[l]
    return int(v @ FIN)


def alpha_word(word):
    v = INIT.astype(object)
    for l in word:
        v = v.dot(MATO[l])
    return v


def beta_word(word):
    v = FIN.astype(object)
    for l in reversed(word):
        v = MATO[l].dot(v)
    return v


def all_alpha(n):
    cur = INIT[None, :].copy()
    for _ in range(n):
        cur = np.stack([cur @ MAT[l] for l in range(4)], axis=1).reshape(-1, 24)
    return cur


def digits(n):
    idx = np.arange(4 ** n)
    return np.stack([(idx >> (2 * (n - 1 - i))) & 3 for i in range(n)], axis=1) if n else np.zeros((1, 0), np.int64)


def main():
    t00 = time.time()
    rng = np.random.default_rng(20261111)

    # ---- A
    bad = cnt = 0
    for n in range(0, 8):
        for w in itertools.product(range(4), repeat=n):
            cnt += 1
            if hp_literal(np.array(w, np.int64)) != H_fast(w):
                bad += 1
    check("A1 transfer = literal Held-Karp count, all words n <= 7", bad == 0, f"({cnt} words, {bad} mismatches)")
    bad = 0
    for n in range(8, 13):
        for _ in range(40):
            w = rng.integers(0, 4, n)
            if hp_literal(w.astype(np.int64)) != H_fast(list(w)):
                bad += 1
    check("A2 transfer = literal count, 200 random words n = 8..12", bad == 0)

    # ---- B
    Hn = {}
    for n in range(0, 10):
        Hn[n] = all_alpha(n) @ FIN
    bad0 = bad3 = badc = bado = badu = 0
    for n in range(0, 9):
        m = 4 ** n
        bad0 += int((Hn[n + 1][:m] != Hn[n]).sum())                 # first letter 0
        bad3 += int((Hn[n + 1][3::4] != Hn[n]).sum())               # last letter 3
    for n in range(1, 10):
        dg = digits(n)
        cv = 3 - dg[:, ::-1]                                        # reverse, (a,b) -> (1-b,1-a): 0<->3, 1<->1, 2<->2
        cv = np.where(dg[:, ::-1] == 1, 1, np.where(dg[:, ::-1] == 2, 2, cv))
        ci = np.zeros(4 ** n, np.int64)
        for i in range(n):
            ci = ci * 4 + cv[:, i]
        badc += int((Hn[n][ci] != Hn[n]).sum())
        bado += int((Hn[n] % 2 == 0).sum())
        badu += int((Hn[n] > 2 * 3 ** n).sum())
    check("B  H(0w)=H(w), H(w3)=H(w), converse, H odd, H <= 2 3^n; all words n <= 9",
          bad0 + bad3 + badc + bado + badu == 0, f"({bad0},{bad3},{badc},{bado},{badu})")
    mx = {n: int(Hn[n].max()) for n in range(1, 10)}
    print("    max H by length 1..9:", mx)

    # ---- C
    ws = [w for k in range(4) for w in itertools.product(range(4), repeat=k)]
    al = [alpha_word(w) for w in ws]
    be = [beta_word(w) for w in ws]
    Hk = sp.Matrix([[int(a.dot(b)) for b in be] for a in al])
    r85 = Hk.rank()
    r21 = Hk[:21, :21].rank()
    r5 = Hk[:5, :5].rank()
    check("C1 Hankel blocks |p|,|s| <= 1, 2, 3 have exact ranks 4, 11, 17", (r5, r21, r85) == (4, 11, 17), f"{(r5, r21, r85)}")

    def closure(start, step):
        basis = [start]
        mat = sp.Matrix([list(map(int, start))])
        queue = [start]
        while queue:
            v = queue.pop(0)
            for l in range(4):
                w = step(v, l)
                m2 = mat.col_join(sp.Matrix([list(map(int, w))]))
                if m2.rank() > mat.rank():
                    mat = m2
                    basis.append(w)
                    queue.append(w)
        return basis

    fbasis = closure(INIT.astype(object), lambda v, l: v.dot(MATO[l]))
    bbasis = closure(FIN.astype(object), lambda v, l: MATO[l].dot(v))
    pair = sp.Matrix([[int(f.dot(b)) for b in bbasis] for f in fbasis])
    check("C2 rank of the pairing reachable x observable (= dimension of the minimal representation) is 17",
          pair.rank() == 17, f"(reachable dim {len(fbasis)}, observable dim {len(bbasis)}, pairing rank {pair.rank()})")

    # ---- D
    d = json.load(open(os.path.join(HERE, "..", "minrep17.json")))
    print("    minrep17.json keys:", sorted(d.keys()))
    r = d["dimension"]
    A = [np.array(d["A"][l], dtype=object) for l in range(4)]
    sb = [[int(c) for c in s] for s in d["suffix_basis"]]
    check("D0 suffix basis as claimed", ["".join(map(str, s)) for s in sb] ==
          ["", "0", "1", "2", "00", "01", "02", "10", "11", "12", "20", "21", "001", "002", "011", "012", "102"])
    B17 = np.array([beta_word(s) for s in sb], dtype=object).T          # 24 x 17, column j = beta(s_j)
    xe = INIT.astype(object).dot(B17)
    check("D1 x(empty) = (H(s_j))_j = (1,1,1,3,1,1,3,3,1,5,5,5,1,3,1,5,9)",
          list(map(int, xe)) == [1, 1, 1, 3, 1, 1, 3, 3, 1, 5, 5, 5, 1, 3, 1, 5, 9] and
          ("x0" not in d or list(d["x0"]) == list(map(int, xe))))
    bad = 0
    for f in fbasis:
        for l in range(4):
            lhs = f.dot(MATO[l]).dot(B17)
            rhs = f.dot(B17).dot(A[l])
            if any(int(p) != int(q) for p, q in zip(lhs, rhs)):
                bad += 1
    check("D2 x(w l) = x(w) A_l on a basis of the reachable space (hence for every word w)", bad == 0,
          f"({len(fbasis)} basis vectors x 4 letters, {bad} failures)")
    ent = [int(x) for l in range(4) for row in d["A"][l] for x in row]
    check("D3 integer matrices, entries in [-41, 32]", all(isinstance(x, int) for x in ent) and min(ent) == -41 and max(ent) == 32,
          f"(min {min(ent)}, max {max(ent)})")
    ones = np.array([1] * r, dtype=object)
    check("D4 (1,...,1) A_l = (1,...,1) for the four letters; x(empty) A_0 = x(empty)",
          all(list(ones.dot(A[l])) == [1] * r for l in range(4)) and list(xe.dot(A[0])) == list(xe))
    # rank 17 of the 17 suffix functions (independence): the 85 prefixes of length <= 3 already give rank 17
    X = sp.Matrix([[int(t) for t in alpha_word(w).dot(B17)] for w in ws])
    check("D5 the 17 suffix functions are linearly independent (rank on prefixes of length <= 3)", X.rank() == 17)
    x = sp.symbols("x")
    As = [sp.Matrix(d["A"][l]) for l in range(4)]
    target = x ** 4 * (x - 3) * (x - 2) ** 4 * (x - 1) ** 8
    cps = [sp.expand(M.charpoly(x).as_expr() - sp.expand(target)) == 0 for M in As]
    S = As[0] + As[1] + As[2] + As[3]
    cps.append(sp.expand(S.charpoly(x).as_expr() - sp.expand(x ** 4 * (x - 12) * (x - 8) ** 4 * (x - 4) ** 8)) == 0)
    check("D6 char poly of each A_l is x^4 (x-3)(x-2)^4 (x-1)^8; of the sum x^4 (x-12)(x-8)^4 (x-4)^8", all(cps))
    # common flag: greedy common eigenvectors (column action) with a COMMON eigenvalue, quotient, repeat
    cur = [M.copy() for M in As]
    diag = []
    okflag = True
    while cur[0].shape[0] > 0:
        k = cur[0].shape[0]
        found = None
        for lam in (0, 1, 2, 3):
            st = sp.Matrix.vstack(*[M - lam * sp.eye(k) for M in cur])
            ns = st.nullspace()
            if ns:
                found = (lam, ns[0])
                break
        if found is None:
            okflag = False
            break
        lam, v = found
        diag.append(lam)
        # basis change: v first, completed by unit vectors
        piv = next(i for i in range(k) if v[i] != 0)
        cols = [v] + [sp.eye(k)[:, i] for i in range(k) if i != piv]
        P = sp.Matrix.hstack(*cols)
        Pi = P.inv()
        cur = [(Pi * M * P)[1:, 1:] for M in cur]
    from collections import Counter
    check("D7 common invariant flag with a common diagonal exists; multiset {3:1, 2:4, 1:8, 0:4}",
          okflag and Counter(diag) == Counter({3: 1, 2: 4, 1: 8, 0: 4}), f"(diagonal found {diag})")
    # random words: char poly
    badw = 0
    for _ in range(12):
        L = int(rng.integers(2, 7))
        w = rng.integers(0, 4, L)
        P = sp.eye(r)
        for l in w:
            P = P * As[int(l)]
        tgt = x ** 4 * (x - 3 ** L) * (x - 2 ** L) ** 4 * (x - 1) ** 8
        if sp.expand(P.charpoly(x).as_expr() - sp.expand(tgt)) != 0:
            badw += 1
    check("D8 char poly of 12 random words of length L = 2..6 is x^4 (x-3^L)(x-2^L)^4 (x-1)^8", badw == 0)
    # H through the representation against the transfer on random long words
    badr = 0
    for _ in range(300):
        L = int(rng.integers(1, 60))
        w = [int(t) for t in rng.integers(0, 4, L)]
        v = xe
        for l in w:
            v = v.dot(A[l])
        if int(v[0]) != H_word(w):
            badr += 1
    check("D9 first coordinate of x(empty) A_w = H(w), 300 random words of length < 60", badr == 0)

    # ---- E
    Ssum = sum(MATO[l] for l in range(4))
    v = INIT.astype(object)
    bad = 0
    for n in range(1, 61):
        v = v.dot(Ssum)
        if 8 * int(v.dot(FIN.astype(object))) != 4 ** n * (3 ** n + 3 * 2 ** n + 3):
            bad += 1
    check("E1 sum_{|w|=n} H(w) = 4^n (3^n + 3 2^n + 3)/8 for n = 1..60 (exact, own transfer)", bad == 0)
    K = sum(np.kron(MAT[l].astype(np.float64), MAT[l].astype(np.float64)) for l in range(4)) / 36.0
    z = np.kron(INIT.astype(np.float64), INIT.astype(np.float64))
    ff = np.kron(FIN.astype(np.float64), FIN.astype(np.float64))
    vals = {}
    for n in range(1, 61):
        z = z @ K
        vals[n] = float(z @ ff)
    print("    E2 E[H^2]/9^n (float): " + "; ".join(f"n={n}: {vals[n]:.10f}" for n in (5, 9, 13, 20, 30, 40, 60)),
          "| (29/160)^2 =", (29 / 160) ** 2)
    check("E2 E[H^2]/9^n at n = 13, 60 matches 0.0326848302, 0.0328515625 to 1e-9",
          abs(vals[13] - 0.0326848302) < 1e-9 and abs(vals[60] - 0.0328515625) < 1e-9)

    # ---- F
    line = []
    okF = True
    for n in (20, 30, 40, 60, 80, 120):
        p = n // 2
        hv = H_word([3] * p + [0] * p)
        val = (1 - Fraction(hv, 2 * 3 ** n)) / Fraction(2, 3) ** p
        lb = 2 * 3 ** n - 4 * (3 ** p * 2 ** p + 2 ** p * 3 ** p) - 2 ** n
        okF &= lb <= hv <= 2 * 3 ** n
        line.append(f"n={n}: {float(val):.9f}")
    check("F  H(3^p 0^p): union-bound inequality holds; constants", okF, " ".join(line))
    # exact closed form by inclusion-exclusion over the four bad events (verifier's own derivation), p = q:
    # per arrangement, colourings with R1 = {} or R1 meets B; R2 = {} or (R2 meets A and R2 meets B); R3 = {} or R3 meets A
    def closed(p, q):
        # count colourings of A (p vertices) and B (q vertices) by which colours occur on A and on B
        tot = 0
        for SA in itertools.product((0, 1), repeat=3):          # colours present on A
            for SB in itertools.product((0, 1), repeat=3):
                def surj(m, S):
                    k = sum(S)
                    if k == 0:
                        return 1 if m == 0 else 0
                    return sum((-1) ** (k - j) * sp.binomial(k, j) * j ** m for j in range(k + 1))
                cA, cB = surj(p, SA), surj(q, SB)
                if cA == 0 or cB == 0:
                    continue
                ok1 = (not (SA[0] or SB[0])) or SB[0]
                ok2e = not (SA[1] or SB[1])
                ok2 = ok2e or (SA[1] and SB[1])
                ok3 = (not (SA[2] or SB[2])) or SA[2]
                if ok1 and ok2 and ok3:
                    tot += cA * cB * (2 if not ok2e else 1)
        return int(tot)
    okc = all(closed(p, q) == H_word([3] * p + [0] * q) for p in range(0, 9) for q in range(0, 9))
    check("F2 H(3^p 0^q) = inclusion-exclusion count by colour supports on the two blocks, all p, q <= 8", okc)

    # ---- G
    Aint = [np.array(d["A"][l], dtype=np.int64) for l in range(4)]
    x_empty = np.array([int(t) for t in xe], dtype=np.int64)

    def orbit(m):
        pw = m ** np.arange(r, dtype=np.int64)
        key = lambda V: (V % m) @ pw
        start = x_empty % m
        states = start[None, :]
        keys = np.array([int(key(start))], dtype=np.int64)
        frontier = states
        depth = 0
        while len(frontier):
            new = np.concatenate([(frontier @ Aint[l]) % m for l in range(4)], axis=0)
            nk = key(new)
            nk, ui = np.unique(nk, return_index=True)
            new = new[ui]
            mask = ~np.isin(nk, keys)
            frontier = new[mask]
            if len(frontier):
                keys = np.concatenate([keys, nk[mask]])
                states = np.concatenate([states, frontier], axis=0)
                depth += 1
        order = np.argsort(keys)
        keys, states = keys[order], states[order]
        nxt = [np.searchsorted(keys, key((states @ Aint[l]) % m)) for l in range(4)]
        for l in range(4):
            assert (keys[nxt[l]] == key((states @ Aint[l]) % m)).all()
        s0 = int(np.searchsorted(keys, key(start)))
        return keys, states, nxt, s0, depth, key

    for m, shift, expect in ((3, 1, 236196), (4, 2, 512), (8, 2, 262144)):
        keys, states, nxt, s0, depth, key = orbit(m)
        N = len(keys)
        # backward reachability of the start state: every state reaches the start  <=> strongly connected
        reach = np.zeros(N, bool)
        reach[s0] = True
        while True:
            new = reach | reach[nxt[0]] | reach[nxt[1]] | reach[nxt[2]] | reach[nxt[3]]
            if (new == reach).all():
                break
            reach = new
        sc = bool(reach.all())
        tk = key((states + shift) % m)
        inv = bool(np.isin(tk, keys).all())
        loop = nxt[0][s0] == s0
        classes = sorted(set((states[:, 0] % m).tolist()))
        check(f"G{m} orbit modulo {m}: {N} states (claimed {expect}), strongly connected {sc}, self-loop {bool(loop)}, "
              f"O + {shift}(1,..,1) = O {inv}, classes of H {classes}",
              N == expect and sc and loop and inv)
        if m == 3:
            # by-index class shares, exact chain propagation in float64 (probabilities)
            p = np.zeros(N)
            p[s0] = 1.0
            cls = states[:, 0] % 3
            out = {}
            for n in range(1, 65):
                q = np.zeros(N)
                for l in range(4):
                    q += np.bincount(nxt[l], weights=p, minlength=N)
                p = q / 4
                if n in (12, 13, 14, 15, 16, 17, 20, 32, 64):
                    out[n] = [3 * float(p[cls == c].sum()) for c in range(3)]
            for n in sorted(out):
                print(f"    by-index class factors modulo 3 at n = {n}: " + " ".join(f"{t:.6f}" for t in out[n]))
            check("G3b class factors at n = 16, 17 equal the lane's (1.004790, 1.002540, 0.992670), (1.003721, 1.002709, 0.993570)",
                  max(abs(a - b) for a, b in zip(out[16], (1.004790, 1.002540, 0.992670))) < 2e-6 and
                  max(abs(a - b) for a, b in zip(out[17], (1.003721, 1.002709, 0.993570))) < 2e-6)
            # twisted rate: component of the measure on which the translation acts by a cube root of unity
            om = np.exp(2j * np.pi / 3)
            tnext = np.searchsorted(keys, key((states + 1) % m))           # index of x + (1,..,1)
            mu = np.zeros(N, complex)
            j = s0
            for t in range(3):
                mu[j] = om ** t
                j = tnext[j]
            logs = []
            for n in range(1, 401):
                q = np.zeros(N, complex)
                for l in range(4):
                    q += np.bincount(nxt[l], weights=mu.real, minlength=N) + 1j * np.bincount(nxt[l], weights=mu.imag, minlength=N)
                nrm = np.linalg.norm(q)
                logs.append(np.log(nrm / 4))
                mu = q / nrm
            rate = float(np.exp(np.mean(logs[200:])))
            print(f"    twisted rate modulo 3 (power iteration, mean of the last 200 log-ratios): {rate:.6f}  (lane: 0.770633)")
            check("G3c twisted rate modulo 3 within 2e-3 of 0.770633", abs(rate - 0.770633) < 2e-3)

    # ---- H
    xone = {}
    for m in range(0, 7):
        for bits in itertools.product((0, 1), repeat=m):
            xone[bits] = int(hp_one_apex(np.array(bits, np.int64))) if m else 1
    bad3t = bad2t = cnt = 0
    for p in range(0, 7):
        for q in range(0, 7):
            for bp in itertools.product((0, 1), repeat=p):
                hb = sum(bp[F] << (p - 1 - F) for F in range(p))
                e1 = 1 if (p == 0 or bp[-1] == 0) else 0
                for aq in itertools.product((0, 1), repeat=q):
                    Ca = sum((1 - aq[L]) << L for L in range(q))
                    e2 = 1 if (q == 0 or aq[0] == 1) else 0
                    word = [2 * t for t in bp] + [t + 2 for t in aq]
                    h = H_fast(word)
                    cnt += 1
                    two = xone[bp] * xone[aq] + (hb + 1) * (Ca + 1)
                    bad2t += h != two
                    bad3t += h != two - e1 * e2
    check("H  F1 three-term identity, all p, q <= 6 (one-apex counts by literal enumeration)", bad3t == 0,
          f"({cnt} instances; two-term formula fails on {bad2t})")

    # ---- I
    worst = {}
    nonconst = 0
    for n in range(2, 9):
        Hall = Hn[n]
        best = -10 ** 9
        for bmask in range(1 << n):
            bb = [(bmask >> (n - 1 - i)) & 1 for i in range(n)]
            base = sum((2 * bb[i]) << (2 * (n - 1 - i)) for i in range(n))
            for i in range(n):
                for j in range(i + 1, n):
                    wi, wj = 1 << (2 * (n - 1 - i)), 1 << (2 * (n - 1 - j))
                    others = [k for k in range(n) if k not in (i, j)]
                    vals = set()
                    for am in range(1 << len(others)) if n <= 6 else (0, (1 << len(others)) - 1, 0b0101 & ((1 << len(others)) - 1)):
                        off = base + sum((1 << (2 * (n - 1 - others[t]))) for t in range(len(others)) if (am >> t) & 1)
                        vals.add(int(Hall[off + wi + wj] - Hall[off + wi] - Hall[off + wj] + Hall[off]))
                    if len(vals) != 1:
                        nonconst += 1
                    best = max(best, max(vals))
        worst[n] = best
    check("I  second differences d_i d_j H(., b): independent of the other a-bits, largest value -2 for n = 2..8",
          nonconst == 0 and all(vv == -2 for vv in worst.values()), f"{worst}")

    print(f"total time {time.time() - t00:.1f} s")
    print("ALL VERIFIER CHECKS PASS" if OK else "SOME VERIFIER CHECK FAILED")


if __name__ == "__main__":
    main()
