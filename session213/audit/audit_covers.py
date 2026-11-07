"""Permutation-and-cut oracle for the attachment kernel and comb coefficients."""
import itertools
import json
from independent_counts import check, count_paths
from audit_comb import comb, quadratic


def kernel(a):
    n = len(a)
    start, end = [0] * n, [0] * n
    pair = [[0] * n for _ in range(n)]
    for p in itertools.permutations(range(n)):
        invalid = [i for i in range(n - 1) if not a[p[i]][p[i + 1]]]
        if not invalid:
            start[p[0]] += 1
            end[p[-1]] += 1
        for cut in range(1, n):
            if not invalid or invalid == [cut - 1]:
                pair[p[cut - 1]][p[cut]] += 1
    return start, end, pair


def main():
    counts = {"general_cores": 0, "general_attachments": 0, "coefficient_cores": 0, "ordered_cover_entries": 0}
    for n in range(1, 5):
        pairs = list(itertools.combinations(range(n), 2))
        for encoding in range(1 << len(pairs)):
            a = [[False] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                a[i][j] = bool(encoding >> k & 1)
                a[j][i] = not a[i][j]
            s, e, K = kernel(a)
            counts["general_cores"] += 1
            for attachment in range(1 << n):
                b = [row + [bool(attachment >> i & 1)] for i, row in enumerate(a)]
                b.append([not bool(attachment >> i & 1) for i in range(n)] + [False])
                x = [(attachment >> i) & 1 for i in range(n)]
                expected = sum(s[i] * (1 - x[i]) + e[i] * x[i] for i in range(n)) + sum(K[i][j] * x[i] * (1 - x[j]) for i in range(n) for j in range(n))
                check(count_paths(b) == expected, ("general attachment", n, encoding, attachment))
                counts["general_attachments"] += 1
    for m in range(1, 6):
        b = comb(m, 0, 0, False)
        keep = [i for i in range(len(b)) if i != m + 2]  # Delete coding apex b; last index is a.
        u = [[b[i][j] for j in keep] for i in keep]
        a = m + 2
        s, e, K = kernel(u)
        A = [e[i] - s[i] + sum(K[i]) for i in range(len(u))]
        check(sum(s) == 2 ** m + 1, ("base core", m))
        check(A[a] == 4 * 2 ** m - 2, ("special A_a", m))
        check(A[1] == 2 ** (m + 1) + 2, ("special A_x1", m))
        check(K[1][a] == 0 and K[a][1] == 2 ** m, ("special K_a_x1", m))
        for i in range(m):
            p = i + 2
            first, last = i == 0, i == m - 1
            expected = 2 ** (m + 1 - i) * 3 ** i + 2 ** (i + 2) - 2 * last
            check(A[p] == expected, ("linear A_p", m, i, A[p], expected))
            check(K[p][a] == 2 ** i + last, ("K_p_a", m, i))
            check(K[a][p] == 2 ** (m - 1 - i) - first, ("K_a_p", m, i))
            check(K[p][1] == 2 ** i - last, ("K_p_x1", m, i))
            check(K[1][p] == 2 ** (m - 1 - i) + first, ("K_x1_p", m, i))
            counts["ordered_cover_entries"] += 4
            for j in range(i + 1, m):
                r, d = j + 2, j - i
                check(K[p][r] == 2 ** (m - d) + (d == 1), ("K_p_r", m, i, j))
                check(K[r][p] == 2 ** (m - d) * 3 ** (d - 1) + 2 ** (d - 1), ("K_r_p", m, i, j))
                counts["ordered_cover_entries"] += 2
        counts["coefficient_cores"] += 1
    support = set()
    for m in range(1, 6):
        for q in range(1 << m):
            for tau, forward in itertools.product(range(2), repeat=2):
                h = quadratic(m, q, tau, forward)
                if h <= 35:
                    support.add(h)
    expected = {3, 5, 9, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33}
    check(support == expected, ("support through35", support))
    print(json.dumps({"status": "PASS", "checks": counts, "all_length_comb_support_through_35": sorted(support)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
