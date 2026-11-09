"""S19 audit: validate tbrute (C subset-DP counter) against a naive pure-Python counter on small tournaments.
Naive: Hamiltonian paths by permutations; odd cycles by permutations of each odd vertex set; packings by explicit
recursion over disjoint cycle sets; port sums by DFS over all paths from u.
Usage: python naive_check.py TBRUTE_EXE NUM SEED
"""
import itertools, random, subprocess, sys
from functools import lru_cache

def naive_stats(n, out, u, xs=(2, 3, 4)):
    A = [[(out[i] >> j) & 1 for j in range(n)] for i in range(n)]
    H = 0; hend_u = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            H += 1
            if p[-1] == u: hend_u += 1
    if n == 0: H = 1
    # odd directed cycles: list of (frozenset, count)
    cycles = []
    for k in range(3, n + 1, 2):
        for U in itertools.combinations(range(n), k):
            first, rest = U[0], U[1:]
            c = 0
            for p in itertools.permutations(rest):
                seq = (first,) + p
                if all(A[seq[i]][seq[(i + 1) % k]] for i in range(k)):
                    c += 1
            if c: cycles.append((frozenset(U), c))
    def packing(V, x):
        # sum over sets of pairwise disjoint cycles inside V of prod(count) x^#
        cyc = [(s, c) for (s, c) in cycles if s <= V]
        def rec(i, used):
            if i == len(cyc): return 1
            s, c = cyc[i]
            tot = rec(i + 1, used)
            if not (s & used): tot += x * c * rec(i + 1, used | s)
            return tot
        return rec(0, frozenset())
    allv = frozenset(range(n))
    res = [H, hend_u]
    for x in xs:
        memo = {}
        def P(V):
            if V not in memo: memo[V] = packing(V, x)
            return memo[V]
        E = O = 0
        stack = [(u, frozenset([u]))]
        while stack:
            v, S = stack.pop()
            w = P(allv - S)
            if len(S) % 2 == 0: E += w
            else: O += w
            for y in range(n):
                if y not in S and A[v][y]: stack.append((y, S | {y}))
        res += [P(allv), E, O, P(allv - {u})]
    return res

def rand_tour(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

def main():
    exe, num, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed)
    cases = []
    for _ in range(num):
        n = rng.randint(1, 8)
        out = rand_tour(n, rng)
        u = rng.randrange(n)
        cases.append((n, out, u))
    inp = "".join(f"{n} {u} " + " ".join(map(str, out)) + "\n" for (n, out, u) in cases)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True).stdout.split("\n")
    bad = 0
    for (n, out, u), line in zip(cases, res):
        got = list(map(int, line.split()))
        exp = naive_stats(n, out, u)
        if got != exp:
            bad += 1
            if bad <= 5: print("MISMATCH", n, out, u, got, exp)
    print(f"naive vs tbrute: {len(cases)} random tournaments of order 1..8 with a random port, {bad} mismatches")
    # identity checks on the same cases: OCF (I_2 = H) and Delta_2 = -hend_u
    ocf = sum(1 for line in res[:len(cases)] if int(line.split()[0]) != int(line.split()[2]))
    dl = 0
    for line in res[:len(cases)]:
        f = list(map(int, line.split()))
        if f[3] - f[4] != -f[1]: dl += 1
    print(f"OCF failures (H != I(Omega,2)): {ocf}; failures of E_2 - O_2 = -hend_u: {dl}")

if __name__ == "__main__":
    main()
