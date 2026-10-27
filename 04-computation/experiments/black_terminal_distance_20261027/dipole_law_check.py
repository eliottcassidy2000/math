"""dipole_law_check.py -- mac-mini-2026-10-27-potentials

Checks of the native-move dipole law on random tournaments (exact integer arithmetic except the two spectral
Lipschitz bounds, which are floating-point inequalities with margin).

For a tournament T with skew matrix K (K_uv = +1 if u->v, -1 if v->u) and a Hamiltonian path Q = (a, ..., b),
the native complement K_Q(T) keeps the arcs of Q and reverses every other pair.  Claims:
  (L1) K' = -K + 2 A_Q, A_Q the skew matrix of the directed path Q;
  (L2) score deviations d = (K 1)/2 obey d' = -d + e_a - e_b;
  (L3) cut charges c(X) = sum_{x in X} d_x = e(X -> X^c) - |X||X^c|/2 obey c'(X) = -c(X) + [a in X] - [b in X];
       equivalently f_{T'}(X^c) = f_T(X) - [a in X] + [b in X] for f(X) = e(X -> X^c);
  (L4) c3(T') - c3(T) = s(a) - s(b) - 1;
  (L5) | ||K'||_op - ||K||_op | <= 2||A_Q||_op = 4 cos(pi/(n+1)) and | ||K'||_* - ||K||_* | <= 2||A_{P_n}||_*,
       with ||A_{P_n}||_* = sum_k |2 cos(k pi/(n+1))|;
  (L6) the folded L1 score potential Phi(T) = (1/2) sum_v |d_v| is 1-Lipschitz.
"""
import itertools
import random
import numpy as np


def random_tournament(n, rng):
    K = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < 0.5:
                K[u, v], K[v, u] = 1, -1
            else:
                K[u, v], K[v, u] = -1, 1
    return K


def random_ham_path(K, rng, tries=2000):
    n = len(K)
    # randomized insertion (Redei): always succeeds
    order = list(range(n))
    rng.shuffle(order)
    path = [order[0]]
    for v in order[1:]:
        slots = []
        if K[v, path[0]] == 1:
            slots.append(0)
        if K[path[-1], v] == 1:
            slots.append(len(path))
        for i in range(len(path) - 1):
            if K[path[i], v] == 1 and K[v, path[i + 1]] == 1:
                slots.append(i + 1)
        path.insert(rng.choice(slots), v)
    return path


def c3(K):
    n = len(K)
    s = (K.sum(axis=1) + (n - 1)) // 2
    return n * (n - 1) * (n - 2) // 6 - sum(x * (x - 1) // 2 for x in s)


def main():
    rng = random.Random(20261027)
    fails = {k: 0 for k in ("L1", "L2", "L3", "L4", "L5op", "L5tr", "L6")}
    trials = 0
    for n in list(range(4, 13)) + [16, 20, 25, 31]:
        AP_op = 2 * np.cos(np.pi / (n + 1))
        AP_tr = sum(abs(2 * np.cos(k * np.pi / (n + 1))) for k in range(1, n + 1))
        for _ in range(300 if n <= 12 else 60):
            K = random_tournament(n, rng)
            Q = random_ham_path(K, rng)
            a, b = Q[0], Q[-1]
            A = np.zeros((n, n), dtype=int)
            for u, v in zip(Q, Q[1:]):
                A[u, v], A[v, u] = 1, -1
            # native complement
            K2 = -K.copy()
            for u, v in zip(Q, Q[1:]):
                K2[u, v], K2[v, u] = 1, -1
            trials += 1
            if not np.array_equal(K2, -K + 2 * A):
                fails["L1"] += 1
            d, d2 = K.sum(axis=1), K2.sum(axis=1)       # = 2 * deviations
            e = np.zeros(n, dtype=int); e[a] += 2; e[b] -= 2
            if not np.array_equal(d2, -d + e):
                fails["L2"] += 1
            if n <= 12:
                for r in range(1, n):
                    for X in itertools.combinations(range(n), r) if n <= 9 else [tuple(rng.sample(range(n), r)) for _ in range(5)]:
                        Xs = set(X)
                        c = sum(d[x] for x in X)
                        c2 = sum(d2[x] for x in X)
                        if c2 != -c + 2 * ((a in Xs) - (b in Xs)):
                            fails["L3"] += 1
            s = (K.sum(axis=1) + (n - 1)) // 2
            if c3(K2) - c3(K) != s[a] - s[b] - 1:
                fails["L4"] += 1
            sv1 = np.linalg.svd(K.astype(float), compute_uv=False)
            sv2 = np.linalg.svd(K2.astype(float), compute_uv=False)
            if abs(sv2[0] - sv1[0]) > 2 * AP_op + 1e-9:
                fails["L5op"] += 1
            if abs(sv2.sum() - sv1.sum()) > 2 * AP_tr + 1e-9:
                fails["L5tr"] += 1
            if abs(np.abs(d2).sum() - np.abs(d).sum()) > 4:     # d here is 2*deviation, so |change in Phi| <= 1
                fails["L6"] += 1
    print(f"dipole law on {trials} random (tournament, Hamiltonian path) pairs, n in 4..12, 16, 20, 25, 31: failures {fails}")


if __name__ == "__main__":
    main()
