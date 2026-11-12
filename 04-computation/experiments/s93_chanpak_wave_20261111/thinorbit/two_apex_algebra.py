#!/usr/bin/env python3
"""Lane thinorbit (S93).  The spine families are triangular (solvable): one-apex R(h) and two-apex G(2,0).

(A) Two-apex family (seed/two_apex_probe.py): transitive spine s_1 -> ... -> s_n, apexes u -> v, letter
    (a_i, b_i): a_i = 1 iff u -> s_i, b_i = 1 iff v -> s_i.  H(word) = f . L_{l_n} ... L_{l_1} v0 with four 36 x 36
    nonnegative integer matrices (re-implemented here from the description of the transfer; re-checked against a
    direct Hamiltonian path count for all words of length <= 4 and random words of length 5, 6).
    Structure: for each arrangement (x, y) = (u, v) or (v, u) the letter acts on the state (r1, r2, r3) in
    {0,1,2} x {0,1,2} x {0,1} as a Kronecker sum X1 (+) X2 (+) X3 (the spine vertex joins run R1, R2 or R3), and
    each factor is a set of rank-one "reset" matrices with an explicit common flag.  The script builds the change
    of basis and CHECKS that all four letters become upper triangular with the same diagonal, entries in
    {0, 1, 2, 3}.  Hence every word of length n has eigenvalues in {3^n, 2^n, 1, 0}: the semigroup lies in a
    Borel subgroup (solvable Zariski closure).
(B) One-apex family: state (q, h, g, 1) with h -> 2h + b, g -> h, q -> q + b (h - g)  (g = floor(h/2),
    q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), x(h) = 2h + 1 - 2 q(h)).  Checks the recurrence against a direct
    Hamiltonian path count of R(h) for h < 2^9 and that both digit matrices are upper triangular in the basis
    (e_q, e_g, e_h, e_1), diagonal (1, 0, 2, 1).
"""
import itertools

import numpy as np
import sympy as sp


def step_matrix(a, b):
    L = np.zeros((36, 36), dtype=np.int64)
    for arr in range(2):
        bx = a if arr == 0 else b
        by = b if arr == 0 else a
        for r1 in range(3):
            for r2 in range(3):
                for r3 in range(2):
                    src = arr * 18 + r1 * 6 + r2 * 2 + r3
                    n1 = 1 if bx == 0 else 2
                    L[arr * 18 + n1 * 6 + r2 * 2 + r3, src] += 1
                    if r2 != 0 or bx == 1:
                        n2 = 1 if by == 0 else 2
                        L[arr * 18 + r1 * 6 + n2 * 2 + r3, src] += 1
                    if r3 != 0 or by == 1:
                        L[arr * 18 + r1 * 6 + r2 * 2 + 1, src] += 1
    return L


def final_vec():
    f = np.zeros(36, dtype=np.int64)
    for arr in range(2):
        for r1 in range(2):
            for r2 in range(2):
                if r2 == 0 and arr == 1:
                    continue
                for r3 in range(2):
                    f[arr * 18 + r1 * 6 + r2 * 2 + r3] = 1
    return f


def ham_paths(adj):
    V = len(adj)
    dp = [[0] * V for _ in range(1 << V)]
    for x in range(V):
        dp[1 << x][x] = 1
    for S in range(1, 1 << V):
        for x in range(V):
            c = dp[S][x]
            if c:
                for y in range(V):
                    if adj[x][y] and not (S >> y) & 1:
                        dp[S | (1 << y)][y] += c
    return sum(dp[(1 << V) - 1])


def brute_two_apex(n, a, b):
    V = n + 2
    u, v = n, n + 1
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if a[i]:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
        if b[i]:
            adj[v][i] = 1
        else:
            adj[i][v] = 1
    adj[u][v] = 1
    return ham_paths(adj)


def main():
    letters = [(0, 0), (1, 0), (0, 1), (1, 1)]
    Ls = [step_matrix(a, b) for (a, b) in letters]
    v0 = np.zeros(36, dtype=np.int64)
    v0[0] = 1
    v0[18] = 1
    f = final_vec()
    ok = True
    rng = np.random.default_rng(5)
    nchk = 0
    for n in range(1, 7):
        words = list(itertools.product(range(4), repeat=n)) if n <= 4 else [tuple(int(x) for x in rng.integers(0, 4, n)) for _ in range(60)]
        for wd in words:
            vec = v0.copy()
            for l in wd:
                vec = Ls[l].dot(vec)
            a = [letters[l][0] for l in wd]
            b = [letters[l][1] for l in wd]
            nchk += 1
            if int(f.dot(vec)) != brute_two_apex(n, a, b):
                ok = False
    print(f"(A) two-apex transfer = direct count on {nchk} words (all of length <= 4, 60 random of lengths 5 and 6): {ok}")
    # factor bases (columns)
    P1 = sp.Matrix([[0, 1, 0], [1, -1, 1], [-1, 0, 0]])     # columns e1 - e2, e0 - e1, e1   (rows indexed by r1 = 0,1,2)
    P2 = sp.Matrix([[0, 0, 1], [1, 1, 0], [-1, 0, 0]])      # columns e1 - e2, e1, e0
    P3 = sp.Matrix([[0, 1], [1, 0]])                        # columns e1, e0
    Pblock = sp.kronecker_product(P1, P2, P3)               # index r1 * 6 + r2 * 2 + r3
    P = sp.diag(Pblock, Pblock)
    assert P.det() != 0
    Pinv = P.inv()
    tri = True
    diags = []
    for L in Ls:
        Lm = sp.Matrix(36, 36, lambda i, j: int(L[i, j]))
        U = Pinv * Lm * P
        if any(U[i, j] != 0 for i in range(36) for j in range(i)):
            tri = False
        diags.append(tuple(U[i, i] for i in range(36)))
    print("    all four letters upper triangular in the explicit tensor basis:", tri)
    print("    the four diagonals coincide:", len(set(diags)) == 1)
    from collections import Counter
    print("    diagonal multiset (per letter):", sorted(Counter(int(x) for x in diags[0]).items()))
    S = Ls[0] + Ls[1] + Ls[2] + Ls[3]
    ev = np.linalg.eigvals(S.astype(float))
    rho = max(abs(ev))
    print(f"    spectral radius of the sum of the four letters: {rho:.6f} (= 4 x 3: every letter has top diagonal entry 3)")
    # max and mean of H over words of length n
    for n in range(1, 11):
        vecs = {(): v0}
        mx = 0
        tot = 0
        stack = [((), v0)]
        cur = [v0]
        for _ in range(n):
            cur = [L.dot(vv) for vv in cur for L in Ls]
        vals = [int(f.dot(vv)) for vv in cur]
        print(f"    n = {n:2d}: max H = {max(vals)}, mean H = {sum(vals) / len(vals):.2f}, 2 * 3^n = {2 * 3 ** n}, max^(1/n) = {max(vals) ** (1 / n):.4f}")
    # (B) one-apex
    def x_rec(h):
        q = 0
        bits = bin(h)[2:] if h else ""
        hh = 0
        for c in bits:
            b = int(c)
            if b:
                q += (hh + 1) // 2
            hh = 2 * hh + b
        return 2 * h + 1 - 2 * q

    def one_apex(h, m):
        # spine t_1 -> ... -> t_m transitive (t_i -> t_j for i < j), apex z; arcs of z = binary digits of h
        V = m + 1
        z = m
        adj = [[0] * V for _ in range(V)]
        for i in range(m):
            for j in range(i + 1, m):
                adj[i][j] = 1
        return adj, z

    # find the digit convention by brute force: try the four natural conventions and keep those matching x(h)
    conv_ok = {}
    for conv in range(4):
        good = True
        for m in range(1, 9):
            for h in range(1 << m):
                # the arcs of the apex are the binary digits of 2h on a spine of m + 1 vertices
                adj, z = one_apex(2 * h, m + 1)
                for i in range(m + 1):
                    bit = ((2 * h) >> i) & 1 if conv & 1 else ((2 * h) >> (m - i)) & 1
                    if (bit == 1) == bool(conv & 2):
                        adj[z][i] = 1
                    else:
                        adj[i][z] = 1
                if ham_paths(adj) != x_rec(h):
                    good = False
                    break
            if not good:
                break
        conv_ok[conv] = good
    print("(B) one-apex: x(h) = 2h + 1 - 2 q(h) equals the direct count of R(h) (apex arcs = binary digits of 2h on a spine "
          "of m + 1 vertices) for all h < 2^m, m = 1..8, under digit conventions (lsb at the spine source?, bit 1 = apex beats?):",
          {(bool(c & 1), bool(c & 2)): v for c, v in conv_ok.items()})
    M0 = sp.Matrix([[1, 0, 0, 0], [0, 2, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])      # coordinates (q, h, g, 1), digit 0
    M1 = sp.Matrix([[1, 1, -1, 0], [0, 2, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])     # digit 1
    okm = True
    for h in range(1, 1 << 10):
        st = sp.Matrix([0, 0, 0, 1])
        for c in bin(h)[2:]:
            st = (M1 if c == "1" else M0) * st
        if st[1] != h or st[2] != h // 2 or 2 * h + 1 - 2 * st[0] != x_rec(h):
            okm = False
    print("    digit matrices reproduce (q, h, floor(h/2)) for all h < 2^10:", okm)
    Q = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])       # basis order (e_q, e_g, e_h, e_1)
    for name, M in (("M0", M0), ("M1", M1)):
        U = Q.inv() * M * Q
        print(f"    {name} in the basis (e_q, e_g, e_h, e_1): {U.tolist()}, upper triangular: "
              f"{all(U[i, j] == 0 for i in range(4) for j in range(i))}, diagonal {[U[i, i] for i in range(4)]}")


if __name__ == "__main__":
    main()
