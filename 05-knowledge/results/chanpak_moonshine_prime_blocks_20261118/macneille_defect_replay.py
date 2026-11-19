from itertools import permutations

PATCH = [(0, 1), (0, 2), (0, 5), (1, 7), (3, 4), (3, 7)]
EXCEPTIONS = [((1, 3),), ((2, 4),), ((0, 2), (4, 6)), ((0, 5), (5, 7))]


def check(cond, detail):
    if not cond:
        raise AssertionError(detail)


def canonical(j):
    rows = [2, 0]
    for b in [0, 0] + [1] * (j + 1):
        n = len(rows)
        rows = [r | ((1 << n) if i or not b else 0)
                for i, r in enumerate(rows)] + [b]
    return rows


def flip(rows, edges):
    rows = rows[:]
    for i, j in edges:
        rows[i] ^= 1 << j
        rows[j] ^= 1 << i
    return rows


def inflate(rows, m):
    mark = len(rows) - 1
    labels = list(range(mark)) + [mark] * m
    return [sum(1 << j for j, y in enumerate(labels)
                if (i < j if x == y else rows[x] >> y & 1))
            for i, x in enumerate(labels)]


def h(rows):
    n = len(rows)
    check(all(not (a >> i & 1) for i, a in enumerate(rows)), rows)
    check(all((rows[i] >> j & 1) + (rows[j] >> i & 1) == 1
              for i in range(n) for j in range(i)), rows)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s, entries in enumerate(dp):
        for i, ways in enumerate(entries):
            if not ways:
                continue
            rest = rows[i] & ~s
            while rest:
                b = rest & -rest
                rest -= b
                dp[s | b][b.bit_length() - 1] += ways
    return sum(dp[-1])


def direct_profile(rows):
    mark = len(rows) - 1
    check(rows[mark] == 1, "mark must beat only 0")
    ans = [0, 0]

    def legal(word):
        return all(rows[i] >> j & 1 for i, j in zip(word, word[1:]))

    for p in permutations(range(mark)):
        for i in range(mark + 1):
            ans[0] += legal(p[:i] + (mark,) + p[i:])
            for j in range(i + 1, mark + 1):
                if p[i] == 0:
                    ans[1] += legal(p[:i] + (mark,) + p[i:j] + (mark,) + p[j:])
    return tuple(ans)


for j, patch in enumerate(EXCEPTIONS):
    a = canonical(j)
    check((h(a), h(flip(a, patch))) == (14 * 2**j - 5, 14 * 2**j + 1), j)

a = canonical(4)
b = flip(a, PATCH)
check(b == [40, 381, 505, 352, 489, 448, 385, 267, 1], b)
check(direct_profile(a) == (219, 112), a)
check(direct_profile(b) == (225, 112), b)
baseline, continuation = set(), set()
for mask in range(64):
    c = flip(a, [edge for i, edge in enumerate(PATCH) if mask >> i & 1])
    hm = h(c)
    if hm == 225:
        baseline.add(mask)
    if h(inflate(c, 2)) - hm == 224:
        continuation.add(mask)
check(baseline == {63} and continuation == {0, 63}, (baseline, continuation))

for m in range(1, 7):
    check(inflate(a, m) == canonical(m + 3), m)
    check(h(inflate(a, m)) == 112 * 2**m - 5, m)
    check(h(inflate(b, m)) == 112 * 2**m + 1, m)
print("PASS: four finite patches, exact marked words, 64 cut faces, six widths")
