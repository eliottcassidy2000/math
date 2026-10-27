"""Small-case ground truth for the THM-4820 recurrence (audit).

1. Explicit enumeration of EVERY binary merge tree on labelled leaves (all (2L-3)!! trees) for every three-unit multiset
   with L <= 7 leaves; the minimum is compared with a count-vector DP written here from scratch.
2. A labelled-subset DP (state = subset of labelled leaves, not a count vector) for multisets with up to 12 leaves,
   compared with the count-vector DP -- this tests 'a subtree's cost depends only on its leaf multiset'.
3. Six-unit leaves: for every six-unit multiset with <= 6 leaves, the explicit minimum over all trees is >= f of the
   three-unit multiset obtained by replacing -1 -> w + w^2, -w -> 1 + w^2, -w^2 -> 1 + w (the (A4) reduction), and the
   root values agree.
"""
import itertools
import sys
from functools import lru_cache

FAIL = []


def check(c, m):
    if not c:
        FAIL.append(m)
        print("FAIL:", m)


def Q(a, b):
    return a * a - a * b + b * b


UNIT = {"1": (1, 0), "w": (0, 1), "w2": (-1, -1), "-1": (-1, 0), "-w": (0, -1), "-w2": (1, 1)}


def merge_cost(y, z):
    x = (y[0] + z[0], y[1] + z[1])
    return x, max(0, Q(*x) - Q(*y) - Q(*z))


def all_trees_min(vals):
    """minimum cost over all binary trees with the given labelled leaf values (explicit enumeration)."""
    idx = tuple(range(len(vals)))

    def trees(S):
        if len(S) == 1:
            yield vals[S[0]], 0
            return
        first, rest = S[0], S[1:]
        for r in range(0, len(rest)):
            for comb in itertools.combinations(rest, r):
                S1 = (first,) + comb
                S2 = tuple(i for i in rest if i not in comb)
                if not S2:
                    continue
                for v1, c1 in trees(S1):
                    for v2, c2 in trees(S2):
                        x, c = merge_cost(v1, v2)
                        yield x, c1 + c2 + c

    best = None
    count = 0
    for v, c in trees(idx):
        count += 1
        best = c if best is None else min(best, c)
    return best, count


@lru_cache(maxsize=None)
def f_counts(a, b, c):
    """count-vector DP, written independently: f(e_i) = 0, f(u) = min over all ordered splits."""
    if a + b + c == 1:
        return 0
    best = None
    for x in range(a + 1):
        for y in range(b + 1):
            for z in range(c + 1):
                s1 = x + y + z
                if s1 == 0 or s1 == a + b + c:
                    continue
                v1 = (x - z, y - z)
                v2 = (a - x - (c - z), b - y - (c - z))
                _, m = merge_cost(v1, v2)
                val = m + f_counts(x, y, z) + f_counts(a - x, b - y, c - z)
                best = val if best is None else min(best, val)
    return best


def subset_dp(vals):
    L = len(vals)
    full = (1 << L) - 1
    val = [None] * (1 << L)
    f = [0] * (1 << L)
    for m in range(1, 1 << L):
        a = b = 0
        for i in range(L):
            if m >> i & 1:
                a += vals[i][0]
                b += vals[i][1]
        val[m] = (a, b)
    for m in range(1, 1 << L):
        if m & (m - 1) == 0:
            f[m] = 0
            continue
        low = m & -m
        best = None
        s = (m - 1) & m
        while s:
            if s & low:  # unordered: the part containing the lowest leaf
                t = m ^ s
                if t:
                    _, c = merge_cost(val[s], val[t])
                    v = f[s] + f[t] + c
                    best = v if best is None else min(best, v)
            s = (s - 1) & m
        f[m] = best
    return f[full]


ntrees = 0
nms = 0
for L in range(1, 8):
    for a in range(L + 1):
        for b in range(L - a + 1):
            c = L - a - b
            vals = [UNIT["1"]] * a + [UNIT["w"]] * b + [UNIT["w2"]] * c
            best, cnt = all_trees_min(vals)
            ntrees += cnt
            nms += 1
            check(best == f_counts(a, b, c), f"explicit vs count DP {(a, b, c)}: {best} {f_counts(a, b, c)}")
print(f"1. explicit enumeration: {nms} three-unit multisets with <= 7 leaves, {ntrees} trees in total; "
      f"minimum == count-vector DP in every case")

nsub = 0
for L in range(8, 13):
    for a in range(L + 1):
        for b in range(L - a + 1):
            c = L - a - b
            if (a + b + c) != L or (L > 10 and not (a >= b >= c)):
                continue
            vals = [UNIT["1"]] * a + [UNIT["w"]] * b + [UNIT["w2"]] * c
            check(subset_dp(vals) == f_counts(a, b, c), f"subset DP vs count DP {(a, b, c)}")
            nsub += 1
print(f"2. labelled-subset DP == count-vector DP on {nsub} multisets with 8..12 leaves")

# symmetry of f under permuting the three unit types (not used by the theorem, a sanity property)
for u in itertools.product(range(5), repeat=3):
    if sum(u) == 0:
        continue
    vs = {f_counts(*p) for p in itertools.permutations(u)}
    check(len(vs) == 1, f"permutation symmetry {u}")
print("   f is invariant under permuting the unit types (all u <= (4,4,4))")

# 3. six-unit reduction
names = list(UNIT)
REPL = {"1": ["1"], "w": ["w"], "w2": ["w2"], "-1": ["w", "w2"], "-w": ["1", "w2"], "-w2": ["1", "w"]}
n6 = 0
for L in range(1, 7):
    for combo in itertools.combinations_with_replacement(names, L):
        vals = [UNIT[nm] for nm in combo]
        best, _ = all_trees_min(vals)
        three = [t for nm in combo for t in REPL[nm]]
        a, b, c = three.count("1"), three.count("w"), three.count("w2")
        root6 = (sum(v[0] for v in vals), sum(v[1] for v in vals))
        root3 = (a - c, b - c)
        check(root6 == root3, f"root value preserved {combo}")
        check(f_counts(a, b, c) <= best, f"six-unit reduction {combo}: f3={f_counts(a, b, c)} > {best}")
        n6 += 1
print(f"3. six-unit leaves: {n6} multisets with <= 6 leaves; replacing negative units by free pairs never increases "
      f"the optimum and preserves the root value")
print("TOTAL FAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
