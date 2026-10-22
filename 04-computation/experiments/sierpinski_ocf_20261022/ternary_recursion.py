# Exact path-cover recursion for lexicographic products over C_3, and the ternary (Sierpinski) tournaments
# S_1 = C_3, S_k = C_3[S_(k-1), S_(k-1), S_(k-1)]  (= Cay(Z_(3^k), {x : unit part of x = 1 mod 3})).
#
# c_j(X) = number of ORDERED j-tuples of nonempty vertex-disjoint directed paths covering X (c_1 = H).
# In T = C_3[A, B, C] (A => B => C => A) a path is a run sequence through the blocks in cyclic order, so a path is
# a word (start block s, length l) and an ordered j-cover of T is an ordered j-tuple of words together with ordered
# covers of each block by its total number m_X of runs:
#     c_j(T) = sum over j-tuples of words  prod_X c_(m_X)(X),     m_X = occurrences of block X in the words.
# For j = 1:  H(T) = sum_(s, l) prod_X c_(m_X(s, l))(X).
# mac-mini-2026-10-22-sierp.
import sys
from functools import lru_cache
sys.setrecursionlimit(10000)

def occ(s, l):
    """occurrence vector (m_0, m_1, m_2) of the cyclic word of length l starting at block s"""
    m = [0, 0, 0]
    for i in range(l): m[(s + i) % 3] += 1
    return tuple(m)

def covers_of_product(blocks, J):
    """blocks = [cA, cB, cC] with cX[m] = ordered m-covers of block X (cX[0] = 1 for the empty cover is NOT used:
    every block must be covered, so m_X >= 1).  Returns c[0..J] for C_3[A, B, C]."""
    M = [len(b) - 1 for b in blocks]               # max useful runs per block
    # words: all (s, l) with occurrence <= M componentwise
    words = []
    for s in range(3):
        l = 1
        while True:
            o = occ(s, l)
            if any(o[i] > M[i] for i in range(3)): break
            words.append(o); l += 1
    # N[j] = dict occurrence-vector -> number of ordered j-tuples of words
    res = [0] * (J + 1)
    cur = {(0, 0, 0): 1}
    for j in range(1, J + 1):
        nxt = {}
        for v, cnt in cur.items():
            for w in words:
                u = (v[0] + w[0], v[1] + w[1], v[2] + w[2])
                if u[0] > M[0] or u[1] > M[1] or u[2] > M[2]: continue
                nxt[u] = nxt.get(u, 0) + cnt
        cur = nxt
        tot = 0
        for (a, b, c), cnt in cur.items():
            if a and b and c: tot += cnt * blocks[0][a] * blocks[1][b] * blocks[2][c]
        res[j] = tot
        if not cur: break
    return res

def point():
    return [1, 1]          # one vertex: c_0 (unused) = 1, c_1 = 1

if __name__ == '__main__':
    import math
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    c = point()
    for k in range(1, K + 1):
        n = 3 ** k
        c = covers_of_product([c, c, c], n)
        H = c[1]
        v3 = 0; t = H
        while t % 3 == 0: t //= 3; v3 += 1
        autexp = (3 ** k - 1) // 2
        norm = H * 2 ** (n - 1) / math.factorial(n)
        print(f"k={k} n={n}: H(S_k) = {H}")
        print(f"    v_3(H) = {v3}, v_3|Aut| = {autexp}, H/|Aut| mod 3 = {(H // 3 ** autexp) % 3}, H*2^(n-1)/n! = {norm:.6f}, sum_j c_j = check c_n = n! : {c[n] == math.factorial(n)}")
        sys.stdout.flush()
