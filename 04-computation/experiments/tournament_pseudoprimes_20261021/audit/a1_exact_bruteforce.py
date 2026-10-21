# AUDIT a1 (independent auditor, 2026-10-21): exact strong-tournament counts and brute-force checks.
#  * exact s(m), m <= 200, computed as minus the Hurwitz (binomial-convolution) inverse of T(m) = 2^C(m,2)
#    (b_0 = 1, b_m = -sum_{k=1..m} C(m,k) T(k) b_{m-k}; s(m) = -b_m for m >= 1);
#  * brute force: enumerate every labelled tournament on n <= 6 vertices, test strong connectivity directly,
#    count circulants on Z_n among them, and test Camion's theorem (strong <=> Hamiltonian cycle);
#  * D(r) = s(r) - c(r) for r <= 16 factored with sympy; D(r) > 0 for 4 <= r <= 200.
# Writes a1_exact_bruteforce.out (via stdout) and a1_exact_D.txt (exact D(r), r <= 200, one per line).
import itertools, sys
from math import comb
from sympy import factorint

N = 200
T = [2 ** (m * (m - 1) // 2) for m in range(N + 1)]
b = [0] * (N + 1)
b[0] = 1
for m in range(1, N + 1):
    b[m] = -sum(comb(m, k) * T[k] * b[m - k] for k in range(1, m + 1))
s = [0] + [-b[m] for m in range(1, N + 1)]

def c(m):
    if m == 1:
        return 1
    return 2 ** ((m - 1) // 2) if m % 2 else 0

D = [None] + [s[m] - c(m) for m in range(1, N + 1)]
print("s(1..10) =", s[1:11])
oeis = [1, 0, 2, 24, 544, 22320, 1677488, 236522496, 64026088576, 33832910196480]
print("matches A054946 prefix (from memory/OEIS data):", s[1:11] == oeis)
# Wright's recurrence as an independent consistency check
ok = all(T[m] == sum(comb(m, k) * s[k] * T[m - k] for k in range(1, m + 1)) for m in range(1, 60))
print("Wright recurrence holds for m<60:", ok)

# brute force n <= 6
def strong(n, adj):
    # adj[i] = bitmask of out-neighbours; strong iff vertex 0 reaches all and all reach 0
    full = (1 << n) - 1
    seen = 1; frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if frontier >> v & 1:
                nxt |= adj[v]
        frontier = nxt & ~seen; seen |= nxt
    if seen != full:
        return False
    radj = [0] * n
    for v in range(n):
        for w in range(n):
            if adj[v] >> w & 1:
                radj[w] |= 1 << v
    seen = 1; frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if frontier >> v & 1:
                nxt |= radj[v]
        frontier = nxt & ~seen; seen |= nxt
    return seen == full

def has_ham_cycle(n, adj):
    if n == 1:
        return True
    if n == 2:
        return False
    # DP from vertex 0
    dp = [0] * (1 << n)
    dp[1] = 1  # bitmask of endpoints reachable with path from 0 covering mask
    for mask in range(1, 1 << n):
        if not (mask & 1) or dp[mask] == 0:
            continue
        ends = dp[mask]
        for v in range(n):
            if ends >> v & 1:
                out = adj[v] & ~mask
                w = 0
                while out:
                    if out & 1:
                        dp[mask | (1 << w)] |= 1 << w
                    out >>= 1; w += 1
    full = (1 << n) - 1
    ends = dp[full]
    return any((ends >> v & 1) and (adj[v] & 1) for v in range(n))

for n in range(1, 7):
    pairs = list(itertools.combinations(range(n), 2))
    cnt_strong = 0; cnt_circ = 0; camion_bad = 0; circ_total = 0
    for bits in range(1 << len(pairs)):
        adj = [0] * n
        for idx, (i, j) in enumerate(pairs):
            if bits >> idx & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        st = strong(n, adj)
        if st != has_ham_cycle(n, adj):
            camion_bad += 1
        # circulant on Z_n: i->j iff (j-i) mod n in S, rotation invariant
        circ = all(((adj[i] >> j & 1) == (adj[(i + 1) % n] >> ((j + 1) % n) & 1)) for i in range(n) for j in range(n) if i != j)
        circ_total += circ
        if st:
            cnt_strong += 1
            cnt_circ += circ
    print(f"n={n}: strong labelled = {cnt_strong} (recurrence {s[n]}), circulant tournaments = {circ_total} "
          f"(formula c={c(n)}), strong circulants = {cnt_circ}, Camion violations = {camion_bad}")

print("D(1..7) =", D[1:8])
print("D(r) > 0 for all 4 <= r <= 200:", all(D[r] > 0 for r in range(4, N + 1)))
for r in range(1, 17):
    print(r, s[r], D[r], factorint(D[r]) if D[r] else 0)
with open("a1_exact_D.txt", "w") as f:
    for r in range(1, N + 1):
        f.write(f"{r} {D[r]}\n")
with open("a1_exact_s.txt", "w") as f:
    for r in range(1, N + 1):
        f.write(f"{r} {s[r]}\n")
