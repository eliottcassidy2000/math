"""Exact finite checks for mixed cuts, odd cycles and one-vertex extensions."""
from collections import Counter


def hamilton(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used, last_counts in enumerate(dp):
        for last, count in enumerate(last_counts):
            if count:
                choices = rows[last] & ~used
                while choices:
                    bit = choices & -choices
                    choices -= bit
                    dp[used | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def mixed_cycles(n, m, mask):
    # Canonicalize an oriented cycle by cyclically moving its least vertex to the front.
    cycles = set()
    for i in range(n):
        for j in range(m):
            for s in range(1, n+1):
                a = tuple((i+r) % n for r in range(s))
                for t in range(1, m+1):
                    if (s+t) % 2 == 0:
                        continue
                    b = tuple((j+r) % m for r in range(t))
                    if mask >> (m*a[-1]+b[0]) & 1 and not (mask >> (m*a[0]+b[-1]) & 1):
                        cycle = a + tuple(n+v for v in b)
                        pivot = cycle.index(min(cycle))
                        cycles.add(cycle[pivot:] + cycle[:pivot])
    return cycles


def boundary(n, m, mask):
    shifts = [(-x, y) for x in range(n) for y in range(m) if (x+y) % 2]
    return sum(not (mask >> (((i+x) % n)*m + (j+y) % m) & 1)
               for i in range(n) for j in range(m) if mask >> (m*i+j) & 1
               for x, y in shifts)


def triangles(mask):
    g = [2, 4, 1, 16, 32, 8]
    for i in range(3):
        for j in range(3):
            if mask >> (3*i+j) & 1:
                g[i] |= 1 << (3+j)
            else:
                g[3+j] |= 1 << i
    return g


def mixed_check(n, m):
    values = Counter()
    gains = Counter()
    for mask in range(1, (1 << (n*m))-1):
        count = len(mixed_cycles(n, m, mask))
        assert count == boundary(n, m, mask)
        values[count] += 1
        if n == m == 3:
            gain = hamilton(triangles(mask)) - 9
            assert gain >= 2 * count
            assert (gain == 8) == (mask.bit_count() in (1, 8))
            gains[gain] += 1
    assert min(values) == n*m//2
    print((n,m), "masks",sum(values.values()),"minimum cycles",min(values),
          "extremizers",values[min(values)])
    if gains:
        assert hamilton(triangles(0)) == hamilton(triangles((1 << 9)-1)) == 9
        print("two triangles, least height gain",min(gains),"extremizers",gains[min(gains)])


def extensions():
    # Adjoin a last distinguished vertex for every possible orientation at each size.
    old = [(0,)]
    for n in range(1,6):
        pairs = Counter()
        total = 0
        following = []
        for g in old:
            h = hamilton(g)
            for new_beats in range(1<<n):
                G = tuple(row | (int(not (new_beats >> i & 1)) << n) for i,row in enumerate(g)) + (new_beats,)
                high = hamilton(G)
                assert high >= h
                if high == h+2:
                    pairs[h] += 1
                total += 1
                if n < 5:
                    following.append(G)
        assert set(pairs) <= {1,3}
        print("parent order",n,"labelled extensions",total,"+2 gains by lower count",dict(sorted(pairs.items())))
        old = following


if __name__ == "__main__":
    for n, m in [(1,3), (1,4), (1,5), (3,3)]:
        mixed_check(n, m)
    extensions()
