# Cayley-Dickson sign tournaments.  Basis units e_0..e_{2^k-1}; e_a e_b = gamma(a,b) e_{a xor b}.
# Tournament CD_k on the 2^k - 1 imaginary units: a -> b iff gamma(a,b) = +1.
import sys, itertools
from functools import lru_cache

def cd_mult(x, y):
    """Cayley-Dickson product of coefficient lists of length 2^k: (a,b)(c,d) = (ac - d* b, d a + b c*)."""
    n = len(x)
    if n == 1: return [x[0] * y[0]]
    h = n // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    def conj(u):
        if len(u) == 1: return u[:]
        hh = len(u) // 2
        return conj(u[:hh]) + [-t for t in u[hh:]]
    def add(u, v): return [s + t for s, t in zip(u, v)]
    def neg(u): return [-t for t in u]
    first = add(cd_mult(a, c), neg(cd_mult(conj(d), b)))
    second = add(cd_mult(d, a), cd_mult(b, conj(c)))
    return first + second

def gamma_table(k):
    n = 1 << k
    g = [[0] * n for _ in range(n)]
    for a in range(n):
        ea = [0] * n; ea[a] = 1
        for b in range(n):
            eb = [0] * n; eb[b] = 1
            p = cd_mult(ea, eb)
            nz = [(i, v) for i, v in enumerate(p) if v]
            assert len(nz) == 1 and nz[0][0] == (a ^ b) and abs(nz[0][1]) == 1, (a, b, nz)
            g[a][b] = nz[0][1]
    return g

def tournament(k):
    g = gamma_table(k); n = (1 << k) - 1
    out = [0] * n
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            if a != b:
                assert g[a][b] == -g[b][a]
                if g[a][b] == 1: out[a - 1] |= 1 << (b - 1)
    return out, g

def ham(out, n):
    f = [dict() for _ in range(1 << n)] if n > 20 else None
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            m = out[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[full]), dp

if __name__ == "__main__":
    for k in range(2, int(sys.argv[1]) + 1):
        out, g = tournament(k); n = (1 << k) - 1
        scores = sorted(bin(o).count('1') for o in out)
        common = sorted(set(bin(out[a] & out[b]).count('1') for a in range(n) for b in range(a + 1, n)))
        # Hadamard check of the full sign table
        N = 1 << k
        had = all(sum(g[a][c] * g[b][c] for c in range(N)) == 0 for a in range(N) for b in range(a + 1, N))
        print("k =", k, "n =", n, "scores:", sorted(set(scores)), "common out-neighbour counts:", common, "sign table Hadamard:", had)
        if n <= 15:
            H, dp = ham(out, n)
            print("   H =", H)
