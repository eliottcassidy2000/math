# AUDIT: test of the own pipeline (aud_fwdpoly -> Lagrange -> CRT -> forest_N) and of the THM-4757 section-5
# conversion against a direct enumeration of directed linear forests (set partitions x orderings of each block)
# for random tournaments with n = 5..8.
import itertools, random, sys
from math import factorial
sys.path.insert(0, '.')
from aud_task1 import F_mod, lagrange_coeffs, crt2, forest_N, P1, P2

def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part

random.seed(7)
for n in range(5, 9):
    code = ''.join(random.choice('01') for _ in range(n * (n - 1) // 2))
    S = [[0] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = 1 if code[k] == '1' else -1; S[i][j] = s; S[j][i] = -s; k += 1
    Nbf = [0] * n
    for part in set_partitions(list(range(n))):
        # each block is ordered in all ways (a directed path); product over blocks
        per_block = []
        for b in part:
            opts = {}
            for o in itertools.permutations(b):
                L = len(o) - 1; w = 1
                for a in range(L): w *= S[o[a]][o[a + 1]]
                opts[L] = opts.get(L, 0) + w
            per_block.append(opts)
        tot = {0: 1}
        for opts in per_block:
            new = {}
            for L1, w1 in tot.items():
                for L2, w2 in opts.items(): new[L1 + L2] = new.get(L1 + L2, 0) + w1 * w2
            tot = new
        for L, w in tot.items(): Nbf[L] += w
    v1 = F_mod(code, P1, range(n)); v2 = F_mod(code, P2, range(n))
    f = [crt2(a, b) for a, b in zip(lagrange_coeffs([v1[x] for x in range(n)], P1), lagrange_coeffs([v2[x] for x in range(n)], P2))]
    N = forest_N(f, n)
    print(f"n = {n}: brute-force linear forests N_L = {Nbf}; via own F + conversion {N}; equal: {N == Nbf}; odd L zero: {all(Nbf[L] == 0 for L in range(1, n, 2))}")
