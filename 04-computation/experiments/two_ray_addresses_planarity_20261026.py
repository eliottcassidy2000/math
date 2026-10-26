"""two_ray_addresses_planarity_20261026.py -- opus-2026-10-26-S61

(1) Two-ray addresses versus the combined table.  In the two-ray carrier (rows D(m) = {2^j m : j >= 0}, columns
    P(m) = {m^k : k >= 1}, m >= 1, with P(1) = {1}) the number of addresses of x (a pair (ray, position) with value x,
    the value itself counted once) is  A(x) = (v2(x) + 1) + tau(g(x)) - 1  for x >= 2, where g(x) is the gcd of the
    exponents of x (tau(g) power addresses (m, k) with k | g; the doubling addresses are (2^b U, v2(x) - b));
    the combined table of aliquot124 (cells 2^j r^k, r >= 2, k >= 1, j >= 0) has  C(x) = sum_(k | g_odd) (floor(A/k) + 1)
    for odd part U > 1 and sum_(k=1)^A floor(A/k) for x = 2^A.  Brute force both to 10^4 and report the difference.
(2) Planarity thresholds of arithmetic unions: on {1..N}, the aliquot functional graph F_N = {{n, s(n)} : s(n) in [1, N]},
    the doubling graph Dbl_N = {{n, 2n}}, the squaring graph Sq_N = {{n, n^2}}, and the square-sum graph Q_N
    ({x, y} : x + y a square).  Smallest N with a non-planar union (networkx check_planarity) and the Kuratowski
    witness type (K_5 or K_(3,3) subdivision, by vertex count of the witness: 5 or 6 branch vertices).
"""
import networkx as nx
from math import isqrt, gcd

def s_aliquot_table(N):
    sig = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d): sig[m] += d
    return [sig[n] - n if n else 0 for n in range(N + 1)]

def exps_gcd(x):
    g = 0; y = x; p = 2
    while p * p <= y:
        if y % p == 0:
            e = 0
            while y % p == 0: y //= p; e += 1
            g = gcd(g, e)
        p += 1
    if y > 1: g = gcd(g, 1)
    return g

def tau(n): return sum(1 for d in range(1, n + 1) if n % d == 0)

def two_ray_brute(x, M):
    cnt = 0
    for m in range(1, M + 1):
        y = m
        while y <= x:                      # row of m
            if y == x: cnt += 1
            y *= 2
        if m >= 2:
            y = m
            while y <= x:                  # column of m
                if y == x: cnt += 1
                y *= m
        else:
            if x == 1: cnt += 1
    # the value x itself sits on its own row (position 0) and own column (position 1): counted twice above; count once
    return cnt - 1 if x >= 2 else cnt - 1

def combined_brute(x, M):
    cnt = 0
    for r in range(2, M + 1):
        rk = r; k = 1
        while rk <= x:
            y = rk
            while y <= x:
                if y == x: cnt += 1
                y *= 2
            rk *= r; k += 1
    return cnt

def main():
    out = []
    # (1)
    N = 10000; bad = 0; diffs = {}
    for x in range(2, N + 1):
        A = 0; y = x
        while y % 2 == 0: y //= 2; A += 1
        g = exps_gcd(x)
        formula = (A + 1) + tau(g) - 1
        br = two_ray_brute(x, x)
        if br != formula: bad += 1
        U = y
        if U > 1:
            godd = 0; z = U; p = 3
            while p * p <= z:
                if z % p == 0:
                    e = 0
                    while z % p == 0: z //= p; e += 1
                    godd = gcd(godd, e)
                p += 2
            if z > 1: godd = gcd(godd, 1)
            comb = sum(A // k + 1 for k in range(1, godd + 1) if godd % k == 0)
        else:
            comb = sum(A // k for k in range(1, A + 1))
        if x <= 100 or x in (36, 216, 1296, 100, 400, 1024, 4096):
            if comb != formula: diffs[x] = (formula, comb)
    out.append(f"(1) two-ray address formula A(x) = v2(x) + tau(g) checked against brute force for 2 <= x <= {N}: mismatches {bad}")
    out.append(f"    two-ray count versus aliquot124's combined-table count, where they differ (x <= 100 and selected): {diffs}")
    out.append(f"    example 36: two-ray {two_ray_brute(36, 36)} (36, 2*18, 4*9 on rows; 6^2 on a column), combined {combined_brute(36, 36)} (adds the mixed cell 4*3^2)")
    # (2)
    M = 3000
    s = s_aliquot_table(M)
    def first_nonplanar(edge_fn, name, limit=M):
        G = nx.Graph()
        for n in range(1, limit + 1):
            G.add_node(n)
            for e in edge_fn(n):
                if 1 <= e[0] <= limit and 1 <= e[1] <= limit and e[0] != e[1]: G.add_edge(*e)
            if n >= 5 and n % 1 == 0:
                ok, cert = nx.check_planarity(G, counterexample=True)
                if not ok:
                    K = cert; branch = [v for v in K.nodes if K.degree(v) >= 3]
                    typ = "K_5" if len(branch) == 5 else ("K_(3,3)" if len(branch) == 6 else f"branch {len(branch)}")
                    return n, typ, sorted(branch)
        return None, None, None
    ali = lambda n: [(n, s[n])] if 1 <= s[n] <= M else []
    dbl = lambda n: [(n, 2 * n)]
    sq = lambda n: [(n, n * n)]
    ss = lambda n: [(n, k * k - n) for k in range(isqrt(n) + 1, isqrt(2 * n) + 1) if 1 <= k * k - n <= M and k * k - n != n]
    for name, fn in [("aliquot F_N", ali), ("aliquot + doubling", lambda n: ali(n) + dbl(n)), ("aliquot + squaring", lambda n: ali(n) + sq(n)),
                     ("doubling + squaring (two-ray graph)", lambda n: dbl(n) + sq(n)), ("aliquot + doubling + squaring", lambda n: ali(n) + dbl(n) + sq(n)),
                     ("square-sum Q_N", ss), ("square-sum + doubling", lambda n: ss(n) + dbl(n))]:
        n, typ, br = first_nonplanar(fn, name, limit=min(M, 400 if "square-sum" in name else M))
        out.append(f"(2) {name}: first non-planar N = {n}, Kuratowski witness {typ}, branch vertices {br}")
    print("\n".join(out))

if __name__ == "__main__":
    main()
