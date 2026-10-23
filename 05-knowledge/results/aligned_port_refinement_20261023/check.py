"""Exact controls for joint source/target refinement of authenticated ports.

Run from any directory: python -B <this file>. Writes adjacent verification.json.
No third-party dependencies. All matrices use common, retained vertex labels.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
from math import comb


def tournament(n, mask, free=None):
    pairs = list(combinations(range(n), 2))
    free = pairs if free is None else free
    bits = {edge: (mask >> k) & 1 for k, edge in enumerate(free)}
    a = [[0] * n for _ in range(n)]
    for i, j in pairs:
        direction = bits.get((i, j), 1)
        a[i][j], a[j][i] = direction, 1 - direction
    return a


def partner(a):
    n = len(a)
    assert all(a[i][i + 1] for i in range(n - 1))
    return [[0 if i == j else a[i][j] if abs(i-j) == 1 else a[j][i]
             for j in range(n)] for i in range(n)]


def refine(a, b):
    n = len(a)
    colors = [0] * n
    history = [list(colors)]
    while len(set(colors)) < n:
        cells = [[i for i in range(n) if colors[i] == c]
                 for c in sorted(set(colors))]
        keys = [(colors[i], tuple(sum(a[i][j] for j in cell) for cell in cells),
                 tuple(sum(b[i][j] for j in cell) for cell in cells))
                for i in range(n)]
        palette = {k: c for c, k in enumerate(sorted(set(keys)))}
        new = [palette[k] for k in keys]
        history.append(new)
        if len(set(new)) == len(set(colors)):
            break
        colors = new
    return history


def agreement(a, b):
    return [[int(x and y) for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def trace_path(d):
    n = len(d)
    outgoing = [sum(row) for row in d]
    incoming = [sum(row[j] for row in d) for j in range(n)]
    starts = [i for i in range(n) if (incoming[i], outgoing[i]) == (0, 1)]
    ends = [i for i in range(n) if (incoming[i], outgoing[i]) == (1, 0)]
    if len(starts) != 1 or len(ends) != 1:
        return None, False
    s, t = starts[0], ends[0]
    if any((incoming[i], outgoing[i]) != (1, 1)
           for i in range(n) if i not in (s, t)):
        return None, False
    path, v = [], s
    while v not in path:
        path.append(v)
        if v == t:
            break
        v = d[v].index(1)
    return path, len(path) == n and path[-1] == t


def determinant_int(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            row = next((i for i in range(k + 1, n) if a[i][k]), None)
            if row is None:
                return 0
            a[k], a[row] = a[row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                assert value % previous == 0
                a[i][j] = value // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def all_controls():
    rows = []
    for n in range(2, 8):
        free = [(i,j) for i,j in combinations(range(n), 2) if j != i + 1]
        rounds = Counter()
        maximal_witness = None
        maximum = 0
        for mask in range(1 << len(free)):
            a = tournament(n, mask, free)
            b = partner(a)
            history = refine(a, b)
            steps = len(history) - 1
            assert len(set(history[-1])) == n and steps <= n // 2
            # Independent check of the signed cancellation, including diagonals.
            for i in range(n):
                for j in range(n):
                    assert a[i][j] + b[i][j] - int(i != j) == (
                        int(j == i + 1) - int(i == j + 1))
            d = agreement(a, b)
            path, accepted = trace_path(d)
            assert accepted and path == list(range(n))
            # Converse swaps the actual endpoints, not an arbitrary color gauge.
            ac = list(map(list, zip(*a)))
            bc = list(map(list, zip(*b)))
            dc = agreement(ac, bc)
            reversed_path, ok = trace_path(dc)
            assert ok and reversed_path == list(reversed(range(n)))
            converse_history = refine(ac, bc)
            assert len(converse_history) == len(history)
            for original, converse in zip(history, converse_history):
                assert all((original[i] == original[j]) == (converse[i] == converse[j])
                           for i in range(n) for j in range(n))
            rounds[steps] += 1
            if steps > maximum:
                maximum, maximal_witness = steps, mask
        rows.append(dict(n=n, fixed_path_tournaments=1 << len(free),
                         refinement_rounds=dict(sorted(rounds.items())),
                         maximum_rounds=maximum, maximal_mask=maximal_witness))

    # Hostile: same authentic endpoints and every local path degree, but a
    # detached directed triangle remains. This is minimal with s != t.
    n = 5
    retained = {(0,1), (2,3), (3,4), (4,2)}
    a = tournament(n, (1 << 10) - 1)
    for i,j in retained:
        a[i][j], a[j][i] = 1, 0
    b = [[0] * n for _ in range(n)]
    for i,j in combinations(range(n),2):
        same = (i,j) in retained or (j,i) in retained
        b[i][j] = a[i][j] if same else a[j][i]
        b[j][i] = 1-b[i][j]
    d = agreement(a,b)
    p, ok = trace_path(d)
    assert p == [0,1] and not ok
    matrix = [[int(i == j)-d[i][j] for j in range(n)] for i in range(n)]
    assert determinant_int(matrix) == 0
    hostile = dict(n=n, source=a, target=b, retained_arcs=sorted(retained),
                   traced_path=p, determinant=0,
                   scope="This specified alignment is not native")

    # All aligned pairs at n <= 4: independent permutation-path recognition
    # versus degree/trace recognition; verify disjoint path/cycle determinant.
    from itertools import permutations
    pairs_summary = []
    for n in range(2,5):
        ts = [tournament(n,m) for m in range(1 << (n*(n-1)//2))]
        paths = list(permutations(range(n)))
        native = degree_passes = false_discrete = 0
        for a in ts:
            for b in ts:
                d = agreement(a,b)
                brute = any(all(d[p[i]][p[i+1]] for i in range(n-1)) and
                            sum(map(sum,d)) == n-1 for p in paths)
                path, ok = trace_path(d)
                assert brute == ok
                native += ok
                if path is not None:
                    degree_passes += 1
                    det = determinant_int([[int(i==j)-d[i][j] for j in range(n)]
                                           for i in range(n)])
                    assert det == int(ok)
                if not ok and len(set(refine(a,b)[-1])) == n:
                    false_discrete += 1
        pairs_summary.append(dict(n=n, all_aligned_pairs=len(ts)**2, native=native,
                                  endpoint_degree_passes=degree_passes,
                                  nonnative_but_jointly_discrete=false_discrete))

    alternating = []
    for n in range(2,65):
        a = [[0]*n for _ in range(n)]
        for i,j in combinations(range(n),2):
            a[i][j] = (j-i) % 2
            a[j][i] = 1-a[i][j]
        hist = refine(a,partner(a))
        alternating.append(dict(n=n, rounds=len(hist)-1, cell_sizes=[
            sorted(Counter(h).values()) for h in hist]))

    # An independent chronological update (stay or advance one path arc) vs
    # Lucas's digit condition; the integer time is never reduced modulo two.
    pascal_entries = 0
    clock_periods = []
    for n in range(2,65):
        state = 1
        first_return = None
        for t in range(128):
            expected = sum(1 << k for k in range(n) if k & t == k)
            assert state == expected
            assert all(((state >> k) & 1) == comb(t,k) % 2 if k <= t
                       else not ((state >> k) & 1) for k in range(n))
            pascal_entries += n
            if t > 0 and state == 1 and first_return is None:
                first_return = t
            state = (state ^ (state << 1)) & ((1 << n)-1)
        expected_period = 1 << ((n-1).bit_length())
        assert first_return == expected_period
        clock_periods.append([n,first_return])
    return dict(status="FINITE-EXACT", universe="All fixed-path tournaments n=2..7; all aligned pairs n=2..4",
                fixed_path=rows, aligned_pairs=pairs_summary,
                endpoint_cycle_hostile=hostile, alternating_family=alternating,
                pascal=dict(orders="2..64", clocks="0..127", entries=pascal_entries,
                            first_snapshot_return=clock_periods))


if __name__ == "__main__":
    result = all_controls()
    out = Path(__file__).with_name("verification.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "alternating_family"}, indent=2))
