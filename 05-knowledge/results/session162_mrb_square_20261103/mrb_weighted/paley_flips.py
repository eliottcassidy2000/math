from itertools import permutations
n = 7
E = [(1, 2), (2, 3), (5, 6)]
counts = []
for mask in range(8):
    out = [sum((1 << j for j in range(n) if (j - i) % 7 in (1, 2, 4))) for i in range(n)]
    for k, (i, j) in enumerate(E):
        if mask >> k & 1:
            if not out[i] >> j & 1:
                raise AssertionError()
            out[i] ^= 1 << j
            out[j] ^= 1 << i
    ps = [p for p in permutations(range(n)) if all((out[i] >> j & 1 for i, j in zip(p, p[1:])))]
    counts.append(len(ps))
    if mask == 7:
        print('masks', out)
        print('endpoint counts', [sum((p[-1] == i for p in ps)) for i in range(n)])
        cycle = (0, 1, 3, 2, 4, 6, 5, 0)
        if not all((out[i] >> j & 1 for i, j in zip(cycle, cycle[1:]))):
            raise AssertionError()
        print('spanning cycle', cycle)
if not counts == [189, 159, 159, 159, 159, 141, 141, 145]:
    raise AssertionError()
a = counts[:]
for k in range(3):
    for m in range(8):
        if m >> k & 1:
            a[m] -= a[m ^ 1 << k]
if not a == [189, -30, -30, 30, -30, 12, 12, -8]:
    raise AssertionError()
print('counts', counts, 'multilinear coefficients', a)
