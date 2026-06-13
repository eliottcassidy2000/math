# c3_formula_kpc1_verify.py — ADVERSARIAL VERIFICATION (session kind-pasteur-2026-06-10-S1)
# Independent check of claim A1: c3(T) = C(n,3) - sum_v C(s_v,2).
# Method (fresh, different from worker): direct cyclic-triple test on every triple
# (a triple {a,b,c} is cyclic iff a->b->c->a or a->c->b->a), exhaustive over ALL
# labeled tournaments for n=4,5 AND n=6 (worker only sampled n=6), plus 800 random
# tournaments each at n=7,8,9 with a different seed (999).
import itertools, random
from math import comb

def direct_c3(adj, n):
    c = 0
    for a, b, cc in itertools.combinations(range(n), 3):
        if (adj[a][b] and adj[b][cc] and adj[cc][a]) or \
           (adj[a][cc] and adj[cc][b] and adj[b][a]):
            c += 1
    return c

def formula_c3(adj, n):
    return comb(n, 3) - sum(comb(sum(adj[v]), 2) for v in range(n))

def check_exhaustive(n):
    pairs = list(itertools.combinations(range(n), 2))
    m = len(pairs)
    mismatches = 0
    for mask in range(1 << m):
        adj = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if (mask >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        if direct_c3(adj, n) != formula_c3(adj, n):
            mismatches += 1
    print(f"n={n}: exhaustive over all {1 << m} labeled tournaments, mismatches={mismatches}")
    return mismatches

def check_random(n, trials, rng):
    pairs = list(itertools.combinations(range(n), 2))
    mismatches = 0
    for _ in range(trials):
        adj = [[0] * n for _ in range(n)]
        for (i, j) in pairs:
            if rng.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        if direct_c3(adj, n) != formula_c3(adj, n):
            mismatches += 1
    print(f"n={n}: {trials} random tournaments (seed 999), mismatches={mismatches}")
    return mismatches

total_mism = 0
for n in (4, 5, 6):
    total_mism += check_exhaustive(n)
rng = random.Random(999)
for n in (7, 8, 9):
    total_mism += check_random(n, 800, rng)
print(f"TOTAL MISMATCHES: {total_mism}")
print("VERDICT(A1 formula):", "PASS" if total_mism == 0 else "FAIL")
