"""Reconstruct the arcs and count the paths directly by their final vertex."""

CERTIFICATES = [
    (33506, (236,189,136,272,933,76,158,296,103,495), 625, 8),
    (33512, (712,821,817,374,129,465,918,270,529,184), 797, 6),
]


def apex(deletion, mask):
    n = len(deletion)
    return [row | ((1 - ((mask >> i) & 1)) << n) for i, row in enumerate(deletion)] + [mask]


def count_paths(rows):
    n = len(rows)
    for i, row in enumerate(rows):
        assert not (row >> i) & 1
        for j in range(i):
            assert ((row >> j) ^ (rows[j] >> i)) & 1
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s, ending in enumerate(dp):
        for i, ways in enumerate(ending):
            allowed = rows[i] & ~s
            while ways and allowed:
                v = allowed & -allowed
                dp[s | v][v.bit_length()-1] += ways
                allowed -= v
    return sum(dp[-1])


if __name__ == '__main__':
    for center, deletion, tags, vertex in CERTIFICATES:
        for index, mask in enumerate([tags, tags ^ (1 << vertex)]):
            rows = apex(deletion, mask)
            score = sorted(r.bit_count() for r in rows)
            assert score != sorted(len(rows)-1-d for d in score)
            H = count_paths(rows)
            assert H == center-1+2*index
            print('center', center, 'count', H, 'score sequence', score)
