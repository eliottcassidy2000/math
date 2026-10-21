# AUDIT a19 (independent auditor, 2026-10-21): THM-4746 Corollary 2, second (rotation) proof, replayed numerically.
# F_k := # strong tournaments on Z_{kp} invariant under the rotation subgroup of order p (k cosets of size p), computed
# from the coset recurrence t_k = sum_j C(k,j) F_j t_{k-j}, t_k = 2^(k(kp-1)/2) (invariant tournaments on k cosets).
# Checks: (i) t_k = eps^k 2^C(k,2) mod p; (ii) F_k = eps^k s(k) mod p; (iii) s(kp) = F_k mod p (direct s(kp) mod p);
# (iv) brute force for p = 3, k = 1, 2: count strong tournaments on Z_6 invariant under x -> x+2 by enumeration.
import subprocess, itertools
from math import comb
def s_direct(n, M):
    out = subprocess.run(["./a2_smod", "single", str(n), str(M)], capture_output=True, text=True).stdout
    return int(out.split("s(n) mod M = ")[1].split()[0])
allok = True
for p in [3, 5, 7, 11, 13]:
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    K = 8
    t = [pow(2, k * (k * p - 1) // 2, p) for k in range(K + 1)]
    F = [0] * (K + 1)
    for k in range(1, K + 1):
        F[k] = (t[k] - sum(comb(k, j) * F[j] * t[k - j] for j in range(1, k))) % p
    s = [0] + [s_direct(k, p) for k in range(1, K + 1)]
    ok1 = all((t[k] - eps ** k * pow(2, k * (k - 1) // 2, p)) % p == 0 for k in range(K + 1))
    ok2 = all((F[k] - eps ** k * s[k]) % p == 0 for k in range(1, K + 1))
    ok3 = all((s_direct(k * p, p) - F[k]) % p == 0 for k in range(1, K + 1))
    allok &= ok1 and ok2 and ok3
    print(f"p={p} eps={eps:+d}: t_k = eps^k T(k) {ok1}; F_k = eps^k s(k) {ok2}; s(kp) = F_k mod p {ok3}  (k <= {K})")
# (iv) brute force: tournaments on Z_6 invariant under x -> x+2 (order 3), strong ones; compare with F_2 (exact, p=3,k=2)
n = 6; pairs = list(itertools.combinations(range(n), 2)); idx = {pq: i for i, pq in enumerate(pairs)}
def strong(adj):
    full = (1 << n) - 1
    def reach(a):
        seen = 1; fr = 1
        while fr:
            nx = 0
            for v in range(n):
                if fr >> v & 1: nx |= a[v]
            fr = nx & ~seen; seen |= nx
        return seen == full
    radj = [0] * n
    for v in range(n):
        for w in range(n):
            if adj[v] >> w & 1: radj[w] |= 1 << v
    return reach(adj) and reach(radj)
inv = 0; inv_strong = 0
for bits in range(1 << len(pairs)):
    adj = [0] * n
    for k, (i, j) in enumerate(pairs):
        if bits >> k & 1: adj[i] |= 1 << j
        else: adj[j] |= 1 << i
    if all(((adj[i] >> j) & 1) == ((adj[(i + 2) % n] >> ((j + 2) % n)) & 1) for i in range(n) for j in range(n) if i != j):
        inv += 1; inv_strong += strong(adj)
t_exact = [2 ** (k * (k * 3 - 1) // 2) for k in range(3)]
F1 = t_exact[1]; F2 = t_exact[2] - 2 * F1 * t_exact[1]
print(f"brute force Z_6, x->x+2: invariant tournaments {inv} (t_2 = {t_exact[2]}), strong invariant {inv_strong} (F_2 exact = {F2})")
print("ALL ROTATION-PROOF CHECKS PASSED" if allok and inv == t_exact[2] and inv_strong == F2 else "CHECK FAILED")
