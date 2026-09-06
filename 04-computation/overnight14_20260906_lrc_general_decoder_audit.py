"""Independent atlas, graph/rank and distinguished-coefficient controls."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
import sys

Q = 91**6
gates = 0


def check(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)


def atlas_sums():
    prime = [True]*357
    prime[0] = prime[1] = False
    for p in range(2, 357):
        if prime[p]:
            for n in range(p*p, 357, p):
                prime[n] = False
    products = {1}
    for p in range(2, 357):
        if prime[p] and p % 3 == 2:
            products = {x*p**e for x in products for e in range(3) if x*p**e <= 356}
    return products-{1}


SUMS = atlas_sums()


def graph(row):
    n = len(row)
    adj = [set() for _ in row]
    vectors = []
    for i, j in combinations(range(n), 2):
        d = gcd(row[i], row[j])
        a, b = row[i]//d, row[j]//d
        if a+b in SUMS:
            adj[i].add(j)
            adj[j].add(i)
            v = [0]*n
            v[i], v[j] = b, -a
            vectors.append(v)
            check(sum(x*y for x, y in zip(v, row)) == 0 and max(map(abs, v)) <= 355,
                  'literal eligible weighted edge')
    remaining, parts = set(range(n)), []
    while remaining:
        current = {min(remaining)}
        while True:
            enlarged = current | set().union(*(adj[i] for i in current))
            if enlarged == current:
                break
            current = enlarged
        remaining -= current
        parts.append(sorted(current))
    return sorted(parts, key=lambda p: (len(p), p)), vectors


def rational_rank(rows):
    if not rows:
        return 0
    a = [[F(x) for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        for r in range(rank+1, len(a)):
            multiplier = a[r][col]/a[rank][col]
            for j in range(col, len(a[0])):
                a[r][j] -= multiplier*a[rank][j]
        rank += 1
        if rank == len(a):
            break
    return rank


def domain_and_rank(row, sizes):
    check(len(row) == len(set(row)) == 13 and min(row) > 0, 'positive distinct thirteen coordinates')
    check(reduce(gcd, row) == 1 and sum(row) <= Q*Q, 'primitive physical-box domain')
    parts, vectors = graph(row)
    check(sorted(map(len, parts)) == sorted(sizes), 'actual component sizes')
    check(rational_rank(vectors) == 13-len(parts), 'independent weighted-edge rank')
    return parts, vectors


def toy_coefficients():
    cases = 0
    for q in range(1, 8):
        for a in range(1, q):
            for b in range(a+1, q+1):
                if gcd(a, b) != 1:
                    continue
                values = {a*r+b*s for r in range(-q, q+1) for s in range(-q, q+1)}
                for d in range(1, 5):
                    for y in range(1, 31):
                        delta = gcd(d, y)
                        c, x = d//delta, y//delta
                        # All distinguished coefficients, independent of the lemma.
                        brute = any(C*y % d == 0 and C*y//d in values for C in range(1, q+1))
                        check(brute == (c <= q and x in values), 'minimal distinguished coefficient iff')
                        # Opposite-variable congruence uses inverse b modulo a.
                        s0 = 0 if a == 1 else pow(b, -1, a)*x % a
                        r0 = (x-b*s0)//a
                        low = max(-((q+s0)//a), -((q-r0)//b))
                        high = min((q-s0)//a, (q+r0)//b)
                        check((low <= high) == (x in values), 'signed box with the other variable anchored')
                        cases += 1
    print('complete toy coefficient cases', cases, 'q1..7,D1..4,Y1..30')


def actual_controls():
    atlas = {(a, s-a) for s in SUMS for a in range(1, (s+1)//2) if gcd(a, s) == 1}
    check(len(atlas) == 5855 and max(b for a, b in atlas) == 355, 'multiplicatively generated atlas')
    domain_and_rank(list(range(1, 14)), [13])
    support_counts = []
    for a in range(1, 7):
        b = 13-a
        small = [3**j for j in range(a)]
        k = max(small)
        g = 2*Q*k+1
        row = small+[g*3**j for j in range(b)]
        parts, vectors = domain_and_rank(row, [a, b])
        check(g > 2*Q*k, 'independent dominance excludes every mixed relation')
        mixed = {tuple(sorted((*pair, third)))
                 for first, other in ((parts[0], parts[1]), (parts[1], parts[0]))
                 for pair in combinations(first, 2) for third in other}
        check(len(mixed) == 11*a*b//2, 'both orientations and singleton coverage')
        for part in parts:
            for i, j in combinations(part, 2):
                check(max(row[i], row[j])//gcd(row[i], row[j]) <= Q, 'all internal height gates')
        reversed_parts, _ = graph(list(reversed(row)))
        check(sorted(map(len, reversed_parts)) == [a, b], 'row-order invariance')
        support_counts.append(len(mixed))
    check(support_counts == [66, 121, 165, 198, 220, 231] and sum(support_counts) == 1001,
          'all six unordered size types')
    print('actual equality types', support_counts, 'total1001, maximum231')

    a = [3**j for j in range(10)]+[1000]
    g = 2*Q*max(a)+1
    row = a+[g, 3*g]
    _, vectors = domain_and_rank(row, [1, 2, 10])
    w = [-1000]+[0]*9+[1, 0, 0]
    check(sum(x*y for x, y in zip(w, row)) == 0 and max(map(abs, w)) <= Q, 'rank-eleven bridge relation')
    check(rational_rank(vectors+[w]) == 11 and g > 2*Q*max(a), 'literal rank11 with dominance upper bound')
    print('three components10,1,2: Vdec rank10; W rank11')

    a = [3**j for j in range(11)]
    g = Q*max(a)+1
    h = 2*g+1
    row = a+[g, h]
    _, vectors = domain_and_rank(row, [1, 1, 11])
    w1, w2 = [0]*13, [0]*13
    w1[0], w1[10], w1[11] = -1, -Q, 1
    w2[0], w2[11], w2[12] = -1, -2, 1
    for w in (w1, w2):
        check(sum(x*y for x, y in zip(w, row)) == 0 and max(map(abs, w)) <= Q,
              'rank-twelve crossing relation')
    check(rational_rank(vectors+[w1, w2]) == 12, 'two independent quotient directions')
    print('three components11,1,1: Vdec rank10; W rank12')

    u = [1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22]
    t = 3*Q+1
    row = [t*x for x in u]+[1, 3]
    _, vectors = domain_and_rank(row, [2, 11])
    w = [0]*13
    w[0], w[11], w[12] = 1, -1, -Q
    check(sum(x*y for x, y in zip(w, row)) == 0 and rational_rank(vectors+[w]) == 12,
          'opposite orientation is a real rank-twelve obstruction')

    a = [355**j for j in range(6)]+[3, 9, 27, 81, 243]
    g = 1000*max(a)+1
    row = a+[g, 3*g]
    _, vectors = domain_and_rank(row, [2, 11])
    w = [0]*13
    w[0], w[5], w[11] = -1, -1000, 1
    check(max(a) > Q and sum(x*y for x, y in zip(w, row)) == 0
          and rational_rank(vectors+[w]) == 12, 'actual internal-height failure witness')
    outside = [355**j for j in range(11)]
    g = 2*Q*max(outside)+1
    check(sum(outside)+4*g > Q*Q, 'out-of-box control cannot enter the theorem')
    print('both two-component rejection paths and physical-box hostile retained')


def main():
    sys.stdout.reconfigure(newline='\n')
    print('Independent general LRC decoder audit: actual atlas/ranks and bounded coefficients')
    toy_coefficients()
    actual_controls()
    print('PASS', gates, 'always-active gates')


if __name__ == '__main__':
    main()
