"""Independent subset path DP and bounded modular square certificates."""
from itertools import combinations
from math import gcd, isqrt


def join(*blocks):
    n = sum(map(len, blocks))
    rows, offset = [], 0
    for block in blocks:
        suffix = (1 << n) - (1 << (offset + len(block)))
        rows.extend((row << offset) | suffix for row in block)
        offset += len(block)
    return rows


def flip(g, u, v):
    h = g.copy()
    h[u] ^= 1 << v
    h[v] ^= 1 << u
    return h


def converse(g):
    return [((1 << len(g)) - 1) ^ (1 << i) ^ r for i, r in enumerate(g)]


def cap(e):
    return flip(join(*([[0]] * (e + 3))), 0, e + 2)


def bits(mask):
    while mask:
        bit = mask & -mask
        mask ^= bit
        yield bit.bit_length() - 1


def reads(g, root=None):
    n = len(g)
    full = (1 << n) - 1
    table = [[0] * n for _ in range(full + 1)]
    all_paths = [1] + [0] * full
    for v in range(n):
        table[1 << v][v] = 1
    for s in range(1, full + 1):
        for v in bits(s):
            val = table[s][v]
            if val:
                all_paths[s] += val
                for w in bits(g[v] & (full ^ s)):
                    table[s | (1 << w)][w] += val
    q = sum(all_paths[s] * all_paths[full ^ s] for s in range(full + 1))
    # 'root' here marks an ending path. Apply to the converse for a start.
    u = None if root is None else sum(
        table[s][root] * all_paths[full ^ s] for s in range(full + 1))
    return all_paths[full], q, u


def precedence(g):
    full = (1 << len(g)) - 1
    precedes = [0] * len(g)

    def walk(path, used):
        if used == full:
            future = full
            for v in path:
                future ^= 1 << v
                precedes[v] ^= future
            return
        for w in bits(g[path[-1]] & (full ^ used)):
            walk(path + (w,), used | (1 << w))

    for v in range(len(g)):
        walk((v,), 1 << v)
    return precedes


def strong(g):
    for v in range(len(g)):
        mask, frontier = 1 << v, 1 << v
        while frontier:
            w = next(bits(frontier))
            frontier ^= 1 << w
            frontier |= g[w] & ~mask
            mask |= g[w]
        if mask != (1 << len(g)) - 1:
            return False
    return True


for e in range(1, 10):
    h, _, u = reads(converse(cap(e)), e + 1)
    assert h == 1 + 2**(e + 1)
    assert u == 2 + 2**(e + 1)

for e, f, middle in [(1, 1, [6, 4, 0]), (1, 2, [2, 4, 1]),
                     (2, 2, [6, 8, 10, 1]), (1, 3, [6, 4, 0])]:
    a, b = cap(e), converse(cap(f))
    old = join(a, middle, b)
    new = flip(old, e + 1, len(a) + len(middle) + f + 1)
    hm, qm, _ = reads(middle)
    expect = (1 + 2**(e + 1)) * hm * (1 + 2**(f + 1))
    expect += (2 + 2**(e + 1)) * qm * (2 + 2**(f + 1))
    assert reads(new)[0] == expect
    assert precedence(new) == precedence(old)
    assert strong(new)
    if (e, f) == (1, 1):
        assert expect == 313

for p, f in [(3, 1), (11, 5)]:
    for b in range(1, 10):
        e = f * p**(b - 1)
        mod = p**b
        assert (2 + pow(2, e + 1, mod)) % mod == 0
        assert gcd(1 + pow(2, e + 1, 33), 33) == 1


def residues(m):
    return {x*x % m for x in range(m//2 + 1)}


def bad(m, end):
    sq = residues(m)
    return [x for x in range(3, end+1) if isqrt(x)**2 != x and x % m in sq]


rest = set(range(2, 144))
for s, count in zip((40, 15, 24, 33, 28, 22, 20, 21, 32),
                    (58, 33, 28, 11, 6, 3, 1, 1, 1)):
    assert 3 <= s <= 47 and isqrt(s)**2 != s
    remove = {m for m in rest if s % m in residues(m)}
    assert len(remove) == count
    rest -= remove
assert not rest
assert bad(144, 52) == [52]
assert bad(1584, 97) == [97]
for m in range(1, 1584):
    k = m
    for p in (2, 3, 11):
        while k % p == 0:
            k //= p
    if k == 1:
        assert bad(m, 95)
assert isqrt(97)**2 != 97 and [97 % m for m in (8, 3, 11)] == [1, 1, 9]
print("strong cap path/precedence identities and sharp square decoder passed")
