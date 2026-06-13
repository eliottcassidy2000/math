"""
Test: does the signed position identity hold for GENERAL digraphs (not just tournaments)?

Identity: sum_S (-1)^|S| [E_i(S)*B_j(R) - E_j(S)*B_i(R)] = 0
where E_v(S) = h_end(S+{v}, v), B_v(R) = h_start({v}+R, v), R = M\S.

If it holds only for tournaments, proof must use T[a][b]+T[b][a]=1.
If it holds for all digraphs, proof might be simpler (pure combinatorics).

Instance: opus-2026-03-05-S4
"""

from itertools import permutations, combinations
import random


def h_end_digraph(G, n, verts, target):
    """# Ham paths in G[verts] ending at target."""
    count = 0
    for perm in permutations(verts):
        if perm[-1] != target:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if G[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def h_start_digraph(G, n, verts, source):
    """# Ham paths in G[verts] starting at source."""
    count = 0
    for perm in permutations(verts):
        if perm[0] != source:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if G[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def test_digraph(n, num_tests=200):
    """Test the identity on random digraphs (not necessarily tournaments)."""
    print(f"=== Digraph test at n={n} ({num_tests} random digraphs) ===")
    passes = 0
    fails = 0

    for trial in range(num_tests):
        # Random digraph: each arc independently present with prob 0.5
        G = [[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(n):
                if a != b:
                    G[a][b] = random.randint(0, 1)

        i, j = 0, 1
        others = [x for x in range(n) if x != i and x != j]
        m = len(others)

        alt_sum = 0
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                S_list = list(S)
                sign = (-1) ** size

                ei = h_end_digraph(G, n, S_list + [i], i)
                bj = h_start_digraph(G, n, [j] + R, j)
                ej = h_end_digraph(G, n, S_list + [j], j)
                bi = h_start_digraph(G, n, [i] + R, i)

                alt_sum += sign * (ei * bj - ej * bi)

        if alt_sum == 0:
            passes += 1
        else:
            fails += 1
            if fails <= 3:
                # Show the digraph
                print(f"  FAIL #{fails}: alt_sum = {alt_sum}")
                for a in range(n):
                    print(f"    G[{a}] = {G[a]}")

    print(f"  Result: {passes}/{num_tests} pass, {fails} fail")
    return fails == 0


def test_tournament_specific(n, num_tests=200):
    """Same test but ONLY on tournaments (to confirm they all pass)."""
    print(f"=== Tournament-only test at n={n} ({num_tests} random) ===")
    passes = 0
    for trial in range(num_tests):
        G = [[0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                if random.random() < 0.5:
                    G[a][b] = 1
                else:
                    G[b][a] = 1

        i, j = 0, 1
        others = [x for x in range(n) if x != i and x != j]
        m = len(others)

        alt_sum = 0
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                S_list = list(S)
                sign = (-1) ** size

                ei = h_end_digraph(G, n, S_list + [i], i)
                bj = h_start_digraph(G, n, [j] + R, j)
                ej = h_end_digraph(G, n, S_list + [j], j)
                bi = h_start_digraph(G, n, [i] + R, i)

                alt_sum += sign * (ei * bj - ej * bi)

        if alt_sum == 0:
            passes += 1

    print(f"  Result: {passes}/{num_tests} pass")


if __name__ == "__main__":
    random.seed(42)
    for n in [4, 5, 6]:
        test_digraph(n, num_tests=500 if n <= 5 else 200)
    print()
    for n in [4, 5, 6]:
        test_tournament_specific(n, num_tests=500 if n <= 5 else 200)
