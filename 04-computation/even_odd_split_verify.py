"""
Even-Odd Split Lemma verification.

From the snippet (opus-S4): The adj decomposition
  delta = sum_S Delta(S, V\{i,j}\S)
satisfies:
  sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)
i.e., sum (-1)^|S| Delta(S,R) = 0.

Where:
  Delta(S, R) = h_end(S+{i}, i) * h_start({j}+R, j) - h_end(S+{j}, j) * h_start({i}+R, i)
and the sum is over all subsets S of V\{i,j}, with R = V\{i,j}\S.

This means delta = 2 * (odd-|S| sum).

Instance: opus-2026-03-05-S4
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import all_tournaments, hamiltonian_path_count
from itertools import permutations, combinations


def h_end(T, verts, target):
    """# Ham paths of T[verts] ending at target."""
    count = 0
    for perm in permutations(verts):
        if perm[-1] != target:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def h_start(T, verts, source):
    """# Ham paths of T[verts] starting at source."""
    count = 0
    for perm in permutations(verts):
        if perm[0] != source:
            continue
        valid = True
        for k in range(len(perm) - 1):
            if T[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def verify_even_odd_split(n, max_tournaments=None):
    print(f"=== Even-Odd Split Lemma verification at n={n} ===")
    total = 0
    passes = 0
    fails = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        if max_tournaments is not None and t_idx >= max_tournaments:
            break
        for i in range(n):
            for j in range(n):
                if i == j or T[i][j] != 1:
                    continue
                others = [x for x in range(n) if x != i and x != j]
                m = len(others)  # n-2

                even_sum = 0
                odd_sum = 0
                alt_sum = 0

                for size in range(m + 1):
                    for S in combinations(others, size):
                        S_set = set(S)
                        R = [x for x in others if x not in S_set]
                        S_list = list(S)

                        he_i = h_end(T, S_list + [i], i)
                        hs_j = h_start(T, [j] + R, j)
                        he_j = h_end(T, S_list + [j], j)
                        hs_i = h_start(T, [i] + R, i)

                        delta_SR = he_i * hs_j - he_j * hs_i

                        if size % 2 == 0:
                            even_sum += delta_SR
                        else:
                            odd_sum += delta_SR
                        alt_sum += ((-1) ** size) * delta_SR

                total += 1
                if alt_sum == 0:
                    passes += 1
                else:
                    fails += 1
                    if fails <= 5:
                        print(f"  FAIL: T#{t_idx} flip {i}->{j}: alt_sum={alt_sum}, even={even_sum}, odd={odd_sum}")

        if t_idx % 50 == 49:
            print(f"  Checked {t_idx+1} tournaments, {total} flips, {fails} failures so far")

    print(f"Result: {passes}/{total} pass, {fails} fail")
    if fails == 0:
        # Also check that delta = 2 * odd_sum for a few cases
        print("Verifying delta = 2 * odd_sum (spot check)...")
        for t_idx, T in enumerate(all_tournaments(n)):
            if t_idx >= 10:
                break
            for i in range(n):
                for j in range(n):
                    if i == j or T[i][j] != 1:
                        continue
                    others = [x for x in range(n) if x != i and x != j]
                    m = len(others)

                    odd_sum = 0
                    total_delta = 0
                    for size in range(m + 1):
                        for S in combinations(others, size):
                            S_set = set(S)
                            R = [x for x in others if x not in S_set]
                            S_list = list(S)

                            he_i = h_end(T, S_list + [i], i)
                            hs_j = h_start(T, [j] + R, j)
                            he_j = h_end(T, S_list + [j], j)
                            hs_i = h_start(T, [i] + R, i)

                            d = he_i * hs_j - he_j * hs_i
                            total_delta += d
                            if size % 2 == 1:
                                odd_sum += d

                    # total_delta should be adj(i,j) - adj'(j,i)
                    # which should equal 2 * odd_sum
                    assert total_delta == 2 * odd_sum, f"delta={total_delta}, 2*odd={2*odd_sum}"
        print("  All spot checks pass: delta = 2 * (odd-S sum)")

    return fails == 0


if __name__ == "__main__":
    for n in [4, 5, 6]:
        ok = verify_even_odd_split(n, max_tournaments=None if n <= 5 else 50)
        print()
        if not ok:
            break
