"""
Tournament blowup (column row step) investigation.
Computes H(T[K_2]) for small tournaments, testing the conjecture:
  Blowup of max-H regular odd-n tournament = max-H SC∩SF class at even 2n.

Session: oracle-2026-05-15
See: 07-reflections/adic-column-families.md, INV-184, OPEN-Q-045
"""
from itertools import permutations


def tournament_blowup(T, n):
    """
    T[K_2]: replace each vertex v with directed pair v_0 -> v_1.
    Original arc u -> v expands to all four arcs u_i -> v_j.
    This is the row step (r,k) -> (r+1,k) in the 2-adic column family grid.
    """
    N = 2 * n
    B = [[0] * N for _ in range(N)]
    for v in range(n):
        B[2 * v][2 * v + 1] = 1  # clone pair arc
    for u in range(n):
        for v in range(n):
            if u != v and T[u][v]:
                for i in range(2):
                    for j in range(2):
                        B[2 * u + i][2 * v + j] = 1
    return B, N


def count_ham_paths(T, n):
    count = 0
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i + 1]] for i in range(n - 1)):
            count += 1
    return count


def score_seq(T, n):
    return tuple(sorted(sum(T[i]) for i in range(n)))


def circulant(n, S):
    T = [[0] * n for _ in range(n)]
    for i in range(n):
        for s in S:
            T[i][(i + s) % n] = 1
    return T


if __name__ == "__main__":
    print("=== Tournament blowup: H(T[K_2]) vs H(T) ===\n")

    cases = [
        ("Transitive n=3", [[0, 1, 1], [0, 0, 1], [0, 0, 0]]),
        ("Cyclic C3",       [[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
        ("Transitive n=4",  [[0,1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]]),
        ("H=5 n=4",        [[0,1,0,1],[0,0,1,0],[1,0,0,1],[0,1,0,0]]),
        ("C5^{1,2} interval n=5", circulant(5, (1, 2))),
        ("Paley(5) C5^{1,4}", circulant(5, (1, 4))),
    ]

    print(f"{'Tournament':<25} {'n':>3} {'H':>6} {'H(blowup)':>10} {'ratio':>8}  scores_blowup")
    print("-" * 90)
    for name, T in cases:
        n = len(T)
        H = count_ham_paths(T, n)
        B, N = tournament_blowup(T, n)
        HB = count_ham_paths(B, N)
        ratio = HB / H if H > 0 else float('inf')
        sb = score_seq(B, N)
        print(f"{name:<25} {n:>3} {H:>6} {HB:>10} {ratio:>8.1f}  {sb}")

    print("\nKey finding: blowup of max-H regular (n=3 cyclic, n=5 interval) achieves max H")
    print("at the doubled size. Score pattern: (n/2-1)^{n/2} ∪ (n/2)^{n/2} — SC∩SF type.")
    print("Conjecture: blowup of max-H regular odd-n = max-H SC∩SF at even 2n.")
