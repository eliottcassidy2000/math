#!/usr/bin/env python3
"""S92 (note, section 10): an independent recomputation of two statements of THM-5133 (G) (kind-pasteur S23, found by
its audit): the Cayley-Dickson sign tournament CD_n (vertices = the 2^n - 1 imaginary units, a -> b iff
e_a e_b = + e_(a xor b); doubling (a, b)(c, d) = (ac - conj(d) b, da + b conj(c))) IS the core T_(2^n - 1) of the
skew-Hadamard doubling tower H -> [[H, H], [-H^T, H^T]] of THM-448 for n = 2, 3, 4 (explicit isomorphisms), and
is NOT for n = 5 (31 vertices: neither T_31 nor its converse; exhaustive search, with two controls).
Nothing here is new; it is a receipt for the sentence in the note.  Exit status 1 on failure.
"""
import random
from itertools import combinations
import sys
from functools import lru_cache

sys.setrecursionlimit(10000)


def conj(v):
    return [v[0]] + [-c for c in v[1:]]


def cd_mul(x, y):
    """Cayley-Dickson product of coefficient lists of length 2^n: (a, b)(c, d) = (ac - d* b, da + b c*)."""
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    h = n // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    ac = cd_mul(a, c)
    dsb = cd_mul(conj(d), b)
    da = cd_mul(d, a)
    bcs = cd_mul(b, conj(c))
    return [p - q for p, q in zip(ac, dsb)] + [p + q for p, q in zip(da, bcs)]


def cd_tournament(n):
    N = 1 << n
    adj = [[0] * (N - 1) for _ in range(N - 1)]
    for a in range(1, N):
        ea = [0] * N
        ea[a] = 1
        for b in range(1, N):
            if a == b:
                continue
            eb = [0] * N
            eb[b] = 1
            pr = cd_mul(ea, eb)
            nz = [(i, v) for i, v in enumerate(pr) if v]
            assert len(nz) == 1 and nz[0][0] == a ^ b and abs(nz[0][1]) == 1
            if nz[0][1] == 1:
                adj[a - 1][b - 1] = 1
    return adj


def tower_tournament(n):
    H = [[1]]
    for _ in range(n):
        m = len(H)
        Ht = [[H[j][i] for j in range(m)] for i in range(m)]
        top = [H[i] + H[i] for i in range(m)]
        bot = [[-v for v in Ht[i]] + Ht[i] for i in range(m)]
        H = top + bot
    m = len(H)
    assert all(H[0][j] == 1 for j in range(m)) and all(H[i][0] == -1 for i in range(1, m))
    assert all(sum(H[i][k] * H[j][k] for k in range(m)) == (m if i == j else 0) for i in range(m) for j in range(m))
    return [[1 if (i != j and H[i][j] == 1) else 0 for j in range(1, m)] for i in range(1, m)]


def is_tournament(adj):
    n = len(adj)
    return all(adj[i][j] + adj[j][i] == 1 for i in range(n) for j in range(n) if i != j)


def find_iso(A, B):
    n = len(A)
    img = [-1] * n
    used = [False] * n

    def rec(k):
        if k == n:
            return True
        for c in range(n):
            if used[c]:
                continue
            if all(A[k][j] == B[c][img[j]] for j in range(k)):
                img[k] = c
                used[c] = True
                if rec(k + 1):
                    return True
                used[c] = False
        img[k] = -1
        return False
    return list(img) if rec(0) else None


def count_aut(A):
    n = len(A)
    img = [-1] * n
    used = [False] * n
    cnt = 0

    def rec(k):
        nonlocal cnt
        if k == n:
            cnt += 1
            return
        for c in range(n):
            if used[c]:
                continue
            if all(A[k][j] == A[c][img[j]] and A[j][k] == A[img[j]][c] for j in range(k)):
                img[k] = c
                used[c] = True
                rec(k + 1)
                used[c] = False
        img[k] = -1
    rec(0)
    return cnt


def ham(adj):
    n = len(adj)
    out = [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v, c in dp[mask].items():
            nxt = out[v] & ~mask
            while nxt:
                w = (nxt & -nxt).bit_length() - 1
                d = dp[mask | (1 << w)]
                d[w] = d.get(w, 0) + c
                nxt &= nxt - 1
        if mask != full:
            dp[mask] = None
    return sum(dp[full].values())


def main():
    ok = True
    for n in (2, 3, 4):
        A = cd_tournament(n)
        B = tower_tournament(n)
        Bc = [[B[j][i] for j in range(len(B))] for i in range(len(B))]
        assert is_tournament(A) and is_tournament(B)
        iso = find_iso(A, B)
        isoc = find_iso(A, Bc)
        scoresA = sorted(sum(r) for r in A)
        line = (f"n = {n}: {len(A)} vertices; regular: {len(set(scoresA)) == 1}; CD_n isomorphic to the tower core: {iso is not None}; "
                f"to its converse: {isoc is not None}; automorphisms of CD_n: {count_aut(A)}")
        if n <= 4:
            hA, hB = ham(A), ham(B)
            line += f"; H(CD_n) = {hA}, H(tower core) = {hB}"
            ok &= hA == hB
        print(line)
        if iso is not None:
            print(f"    isomorphism (unit a of CD_n -> row of the core, 1-based): {[v + 1 for v in iso]}")
        ok &= (iso is not None) or (isoc is not None)
    # n = 5: not isomorphic, with controls (a relabelling of each is found isomorphic to itself)
    A, B = cd_tournament(5), tower_tournament(5)
    N = len(A)
    Bc = [[B[j][i] for j in range(N)] for i in range(N)]
    perm = list(range(N))
    random.seed(1)
    random.shuffle(perm)
    A2 = [[A[perm[i]][perm[j]] for j in range(N)] for i in range(N)]
    B2 = [[B[perm[i]][perm[j]] for j in range(N)] for i in range(N)]
    c1, c2 = find_iso(A, A2) is not None, find_iso(B, B2) is not None
    i1, i2 = find_iso(A, B) is not None, find_iso(A, Bc) is not None
    autA, autB = count_aut(A), count_aut(B)

    def census(T):
        c4 = {}
        for S in combinations(range(N), 4):
            t = sum(1 for a, b, c in combinations(S, 3)
                    if (T[a][b] and T[b][c] and T[c][a]) or (T[b][a] and T[c][b] and T[a][c]))
            c4[t] = c4.get(t, 0) + 1
        return sorted(c4.items())
    print(f"n = 5: 31 vertices; controls (each against a relabelling of itself): {c1}, {c2}; CD_5 isomorphic to T_31: {i1}; "
          f"to its converse: {i2}; automorphisms: CD_5 {autA}, T_31 {autB}; 4-subsets by number of 3-cycles: "
          f"CD_5 {census(A)}, T_31 {census(B)}")
    ok &= c1 and c2 and (not i1) and (not i2)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
