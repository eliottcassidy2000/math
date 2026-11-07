"""resdp_brute.py - brute-force check of the letter-count model: enumerate literal words (exact integer states), no
automaton. For each n, every word with exactly n letters O/Z and at most two letters J (in any slots, repeats
allowed), weight 2 if it has <= 1 J and 1 if it has 2, histogram of F(t.w) = D - E + U mod q.
usage: python3 resdp_brute.py q n
"""
import sys
from itertools import product, combinations_with_replacement

q, n = int(sys.argv[1]), int(sys.argv[2])


def step(s, c):
    D, E, U, d = s
    if c == 'O':
        return (D, 0, 2 * U, d)
    if c == 'Z':
        return (D - E + U, d, 2 * U - d, d)
    return (D - E + U, D - E, D + U, D)


cls = [0] * q
for oz in product('OZ', repeat=n):
    for j in range(3):
        for slots in combinations_with_replacement(range(n + 1), j):
            w = []
            for i in range(n + 1):
                w += ['J'] * slots.count(i)
                if i < n:
                    w.append(oz[i])
            s = (1, 0, 2, 1)
            for c in w:
                s = step(s, c)
            cls[(s[0] - s[1] + s[2]) % q] += 2 if j <= 1 else 1
print('q=%d n=%d exact weighted counts (total %d): %s' % (q, n, sum(cls), ' '.join(map(str, cls))))
