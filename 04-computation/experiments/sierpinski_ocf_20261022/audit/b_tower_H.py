# Audit of THM-4755 statement 4: exact H(S_k), S_k = C3^[k], k <= 4, by code written from scratch.
#  (a) subset DP for path covers of small tournaments (ground truth for <= 12 vertices);
#  (b) "shape-sum" substitution law (derived by the auditor, cf. THM-3121): for T = C3[B0,B1,B2],
#      c_M(T) = sum over ordered M-tuples of C3-walks (start s, length L>=1) of prod_i c_{m_i}(B_i),
#      m_i = total visits to block i; validated against (a) on random C3[B0,B1,B2] with distinct blocks;
#  (c) H(S_3) also by the different decomposition S_3 = S_2[C3] (word DP on S_2 with multiplicities <= 3);
#  (d) palindromy of F(S_3,x) = sum_m c_m (x-1)^(27-m) as a check on all c_m(S_3).
import random, sys
from math import factorial, comb
from fractions import Fraction

def path_covers(n, adj):
    """ordered path-cover numbers c_1..c_n of a digraph on n vertices (adj[i][j] = i->j)."""
    N = 1 << n
    end = [[0] * n for _ in range(N)]
    for v in range(n): end[1 << v][v] = 1
    for S in range(1, N):
        row = end[S]
        for v in range(n):
            c = row[v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and adj[v][w]:
                        end[S | (1 << w)][w] += c
    P = [sum(end[S]) for S in range(N)]
    P[0] = 1
    U = [[0] * N for _ in range(n + 1)]
    U[0][0] = 1
    for mm in range(1, n + 1):
        for S in range(1, N):
            low = S & -S
            rest = S ^ low
            s = 0
            sub = rest
            while True:
                B = sub | low
                s += P[B] * U[mm - 1][S ^ B]
                if sub == 0: break
                sub = (sub - 1) & rest
            U[mm][S] = s
    return [None] + [factorial(mm) * U[mm][N - 1] for mm in range(1, n + 1)]

def lexprod(AQ, AB):
    q, b = len(AQ), len(AB)
    return [[(AQ[r // b][c // b] if r // b != c // b else AB[r % b][c % b]) for c in range(q * b)] for r in range(q * b)]

def subst(AQ, blocks):
    sizes = [len(B) for B in blocks]
    off = [sum(sizes[:i]) for i in range(len(sizes))]
    n = sum(sizes)
    A = [[0] * n for _ in range(n)]
    for i, Bi in enumerate(blocks):
        for j, Bj in enumerate(blocks):
            for x in range(sizes[i]):
                for y in range(sizes[j]):
                    if i == j:
                        A[off[i] + x][off[j] + y] = Bi[x][y]
                    else:
                        A[off[i] + x][off[j] + y] = AQ[i][j]
    return A

C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]

def walk_counts(s, L):
    cnt = [0, 0, 0]
    for t in range(L): cnt[(s + t) % 3] += 1
    return tuple(cnt)

def shape_sum(cs, M_max):
    """cs = [c^(0), c^(1), c^(2)] ordered path-cover lists of the blocks (index 1..b_i).
       returns c_M(C3[B0,B1,B2]) for M = 1..M_max via polynomial powers of the walk generating function."""
    b = [len(c) - 1 for c in cs]
    # G: dict (m0,m1,m2) -> number of walks with those counts (truncated)
    G = {}
    for s in range(3):
        for L in range(1, sum(b) + 1):
            k = walk_counts(s, L)
            if all(k[i] <= b[i] for i in range(3)):
                G[k] = G.get(k, 0) + 1
    res = [None]
    cur = {(0, 0, 0): 1}
    for M in range(1, M_max + 1):
        nxt = {}
        for k1, v1 in cur.items():
            for k2, v2 in G.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                if k[0] <= b[0] and k[1] <= b[1] and k[2] <= b[2]:
                    nxt[k] = nxt.get(k, 0) + v1 * v2
        cur = nxt
        tot = 0
        for k, v in cur.items():
            if min(k) >= 1:
                tot += v * cs[0][k[0]] * cs[1][k[1]] * cs[2][k[2]]
        res.append(tot)
    return res

def rand_tour(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def v3(x):
    v = 0
    while x % 3 == 0: x //= 3; v += 1
    return v

if __name__ == "__main__":
    rng = random.Random(20261022)
    # (b) validation on random substitutions with distinct blocks
    bad = 0
    for trial in range(40):
        blocks = [rand_tour(rng.randint(1, 4), rng) for _ in range(3)]
        T = subst(C3, blocks)
        n = len(T)
        direct = path_covers(n, T)
        cs = [path_covers(len(B), B) for B in blocks]
        ss = shape_sum(cs, n)
        if direct[1:] != ss[1:]: bad += 1
    print("shape-sum law vs direct DP on 40 random C3[B0,B1,B2] (<= 12 vertices): mismatches =", bad)
    # tower
    c1 = path_covers(3, C3)
    S2 = lexprod(C3, C3)
    c2_direct = path_covers(9, S2)
    c2_ss = shape_sum([c1, c1, c1], 9)
    print("c_m(C3) =", c1[1:])
    print("c_m(S_2) direct DP =", c2_direct[1:])
    print("c_m(S_2) shape-sum == direct:", c2_ss[1:] == c2_direct[1:])
    c3 = shape_sum([c2_direct, c2_direct, c2_direct], 27)
    H = {1: c1[1], 2: c2_direct[1], 3: c3[1]}
    # closed form H(C3[B,B,B]) = 3 sum_a c_a (c_a^2 + c_a c_{a+1} + c_{a+1}^2)
    def Hclosed(c):
        b = len(c) - 1
        cc = c + [0]
        return 3 * sum(cc[a] * (cc[a] ** 2 + cc[a] * cc[a + 1] + cc[a + 1] ** 2) for a in range(1, b + 1))
    print("H(S_2) closed form from c(C3):", Hclosed(c1), " H(S_3) closed form from c(S_2):", Hclosed(c2_direct), " shape-sum M=1:", c3[1])
    # (d) palindromy of F(S_3, x)
    n = 27
    F = [0] * n
    for mm in range(1, n + 1):
        j = n - mm
        for k in range(j + 1):
            F[k] += c3[mm] * comb(j, k) * (-1) ** (j - k)
    print("F(S_3,x): palindromic:", F == F[::-1], " all >= 0:", all(f >= 0 for f in F), " sum = 27!:", sum(F) == factorial(27), " F_26 = H:", F[26] == c3[1])
    # H(S_4) from c(S_3)
    H[4] = Hclosed(c3)
    print("H(S_4) via closed form from c(S_3):", H[4])
    for k in range(1, 5):
        nn = 3 ** k
        h = H[k]
        kap = Fraction(h * 2 ** (nn - 1), factorial(nn))
        print(f"k={k} n={nn}: H has {len(str(h))} digits, v3(H)={v3(h)}, (3^k-1)/2={(3**k-1)//2}, kappa=H 2^(n-1)/n! = {float(kap):.7f}")
    print("H(S_3) =", H[3])
