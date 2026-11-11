"""Exact new-shadow controls; run from any directory with Python 3."""
from itertools import combinations
from random import Random
import json


def tournament(n, mask):
    out = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if mask >> k & 1 else (j, i)
        out[a] |= 1 << b
    return out


def paths(out):
    n = len(out)
    def walk(path, left):
        if not left:
            yield path
            return
        opts = left if not path else left & out[path[-1]]
        while opts:
            bit = opts & -opts
            opts -= bit
            yield from walk(path + (bit.bit_length() - 1,), left ^ bit)
    return list(walk((), (1 << n) - 1))


def shadow(out):
    ps = paths(out)
    n = len(out)
    w = [[0] * n for _ in out]
    for p in ps:
        for i, j in zip(p, p[1:]):
            w[i][j] += 1
            w[j][i] += 1
    return len(ps), w


def determinant(a):
    """Bareiss exact determinant; all tested minors are positive definite."""
    a = [row[:] for row in a]
    last = 1
    for k in range(len(a) - 1):
        pivot = a[k][k]
        assert pivot
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                assert value % last == 0
                a[i][j] = value // last
        last = pivot
    return a[-1][-1] if a else 1


def tree_count(w):
    n = len(w)
    return determinant([[sum(w[i]) if i == j else -w[i][j]
                         for j in range(n - 1)] for i in range(n - 1)])


def run():
    rng = Random(265)
    tested = equality = 0
    examples = {}
    for n in range(2, 8):
        masks = range(1 << (n * (n - 1) // 2)) if n <= 5 else (
            rng.randrange(1 << (n * (n - 1) // 2)) for _ in range(128))
        for mask in masks:
            h, w = shadow(tournament(n, mask))
            tau = tree_count(w)
            assert tau >= h ** (n - 1)
            assert (tau == h ** (n - 1)) == (h == 1)
            assert sum(map(sum, w)) == 2 * (n - 1) * h
            assert all(h <= sum(row) <= 2 * h for row in w)
            tested += 1
            equality += h == 1
            if h in (3, 5, 63) and str(h) not in examples:
                examples[str(h)] = dict(n=n, mask=mask, H=h, trees=tau,
                                        lower_bound=h ** (n - 1), W=w)
    # Exact constructive hostile: all linear path-polytope tests admit 3P_4,
    # but equality rigidity says no tournament can have this shadow.
    fake = [[0,3,0,0],[3,0,3,0],[0,3,0,3],[0,0,3,0]]
    assert tree_count(fake) == 27
    return dict(status="PASS", cases=tested, equality_cases=equality,
                examples=examples, scaled_path_hostile=dict(H=3, n=4, trees=27))


if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True, indent=2))
