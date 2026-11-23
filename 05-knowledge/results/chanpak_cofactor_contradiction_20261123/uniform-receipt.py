from math import gcd, factorial

def all_path_counts(rows):
    n = len(rows)
    d = [[0] * n for _ in range(1 << n)]
    counts = [1] + [0] * ((1 << n) - 1)
    for i, row in enumerate(rows):
        assert not row >> i & 1
        assert all((row >> j & 1) + (rows[j] >> i & 1) == 1
                   for j in range(i))
        d[1 << i][i] = 1
    for s in range(1, 1 << n):
        counts[s] = sum(d[s])
        for v, value in enumerate(d[s]):
            if not value:
                continue
            rem = rows[v] & ~s
            while rem:
                b = rem & -rem
                rem -= b
                d[s | b][b.bit_length() - 1] += value
    return counts

def slots(rows):
    h = all_path_counts(rows)
    prev = [1] + [0] * (len(h) - 1)
    top = []
    for _ in rows:
        cur = [0] * len(h)
        for s in range(1, len(h)):
            a = s
            while a:
                cur[s] += h[a] * prev[s ^ a]
                a = (a - 1) & s
        top.append(cur[-1])
        prev = cur
    return top

def ladder(rows, r):
    n = len(rows)
    full = (1 << (n + 2*r)) - 1
    outside = tuple(range(n, n + 2*r))
    bmask = sum(1 << j for j in outside[::2])
    return [x | bmask for x in rows] + [
        (full ^ ((1 << (j + 1)) - 1)) | (((1 << n) - 1) if (j - n) % 2 else 0)
        for j in outside]

def strong(rows):
    for i in range(len(rows)):
        s = 1 << i
        while True:
            t = s
            for j, row in enumerate(rows):
                if s >> j & 1:
                    t |= row
            if s == t:
                break
            s = t
        assert s == (1 << len(rows)) - 1

fixtures = [
    ((48,49,3,7,12,28),45,432,9),
    ((106,52,121,114,97,64,2),49,812,7),
    ((10,116,121,114,97,65,1),49,872,1),
    ((10,20,25,18,1),13,104,1),
    ((38,124,120,113,97,64,1),37,740,1),
    ((50,116,121,115,96,64,1),29,580,1),
    ((58,76,121,112,98,66,1),91,1196,1),
]
for rows, h, q2, g in fixtures:
    strong(rows)
    P = slots(rows)
    assert (P[0], 2*P[0] + P[1]) == (h, q2)
    assert gcd(*(factorial(k)*a for k,a in enumerate(P))) == g
    heights = [all_path_counts(ladder(rows, r))[-1] for r in range(3)]
    assert heights[:2] == [h,h+q2]
    print(h, q2, g, P, heights)
