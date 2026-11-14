from itertools import combinations


def tape(indices):
    p, h, word = 2, 3, '10'
    for j in indices:
        increment = (1 << j) * (p + 1)
        p, h = 2 * increment - 2, h + increment
        word += '1' * j + '0'
    return p, h, word


def count_literal(word):
    n = len(word) + 1
    full = (1 << n) - 1
    rows = [
        (full ^ ((1 << (i + 1)) - 1)) ^
        ((word[i] == '1') << (n - 1))
        for i in range(n - 1)
    ]
    rows += [sum((bit == '1') << i for i, bit in enumerate(word))]
    dp = [[0] * n for _ in range(full + 1)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, full):
        for i, v in enumerate(dp[mask]):
            if v:
                options = rows[i] & (full ^ mask)
                while options:
                    b = options & -options
                    dp[mask | b][b.bit_length() - 1] += v
                    options ^= b
    return sum(dp[full])


def decode(n):
    p, h, j, selected = 2, 3, 1, []
    while h < n:
        if ((h ^ n) >> j) & 1:
            q = (1 << j) * (p + 1)
            p, h = 2 * q - 2, h + q
            selected.append(j)
        j += 1
    return tuple(selected) if h == n else None


literal = 0
for a in range(2, 12):
    residues = set()
    for k in range(a):
        for indices in combinations(range(1, a), k):
            p, h, word = tape(indices)
            if p != int(word, 2) or decode(h) != indices:
                raise RuntimeError((a, indices, p, h))
            residues.add(h % (1 << a))
            if a == 11 and len(word) <= 11:
                literal += 1
                if h != count_literal(word):
                    raise RuntimeError((h, word))
    if residues != set(range(1, 1 << a, 2)):
        raise RuntimeError((a, len(residues)))
print(f'Bijection modulo all powers 2^a through a=11; {literal} physical DP checks.')
