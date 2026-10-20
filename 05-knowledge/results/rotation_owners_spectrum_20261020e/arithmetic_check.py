"""Exact parity-path, arithmetic-fibre and braid-sidecar controls."""
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def words(n):
    if n == 0:
        return ((),)
    return tuple((letter,) + tail for letter, size in [('a', 1), ('b', 3), ('c', 4)]
                 if size <= n for tail in words(n - size))


def encode(word):
    i, edges = 1, []
    for letter in word:
        if letter != 'a':
            edges.append((i, i + 2))
        if letter == 'c':
            edges.append((i + 1, i + 3))
        i += {'a': 1, 'b': 3, 'c': 4}[letter]
    return tuple(edges)


def decode(n, edges):
    i, result, edges = 1, [], set(edges)
    while i <= n:
        if (i, i + 2) not in edges:
            result.append('a')
            i += 1
        elif (i + 1, i + 3) in edges:
            result.append('c')
            i += 4
        else:
            result.append('b')
            i += 3
    need(i == n + 1, 'decode overshot boundary')
    return tuple(result)


def matchings(available):
    if not available:
        yield ()
        return
    v = min(available)
    for tail in matchings(available - {v}):
        yield tail
    if v + 2 in available:
        for tail in matchings(available - {v, v + 2}):
            yield ((v, v + 2),) + tail


def polynomial_fibonacci(n):
    a, b = [0], [1]
    for _ in range(n):
        c = b[:] + [0]
        for j, coefficient in enumerate(a):
            c[j + 1] += coefficient
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        a, b = b, c
    return a


def product(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def rows_for(word):
    atoms = {'a': (0,), 'b': (2, 4, 1), 'c': (2, 4, 9, 3)}
    n = sum(len(atoms[a]) for a in word)
    result, offset = [], 0
    for a in word:
        atom = atoms[a]
        later = ((1 << n) - 1) ^ ((1 << (offset + len(atom))) - 1)
        result.extend((row << offset) | later for row in atom)
        offset += len(atom)
    return tuple(result)


def hamiltonian(rows):
    if not rows:
        return 1
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for end, value in enumerate(dp[mask]):
            candidates = rows[end] & ~mask
            while candidates:
                bit = candidates & -candidates
                candidates -= bit
                dp[mask | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def main():
    total = direct = 0
    counts = []
    for n in range(19):
        ws = words(n)
        encoded = {encode(w) for w in ws}
        independent = set(matchings(frozenset(range(1, n + 1))))
        need(len(encoded) == len(ws) and encoded == independent, 'bijection universe')
        for w in ws:
            edges = encode(w)
            need(decode(n, edges) == w, 'word round trip')
            need(len(edges) == w.count('b') + 2 * w.count('c'), 'triangle grading')
            if n <= 8:
                need(hamiltonian(rows_for(w)) == 3 ** w.count('b') * 5 ** w.count('c'),
                     'literal H versus SCC product')
                direct += 1
        counter = Counter(len(e) for e in encoded)
        expected = product(polynomial_fibonacci((n + 1) // 2 + 1),
                           polynomial_fibonacci(n // 2 + 1))
        need([counter[i] for i in range(len(expected))] == expected,
             'full parity-path matching polynomial')
        total += len(ws)
        counts.append(len(ws))
    print('Parity-path bijection: words/matchings through order18 =', total,
          '; literal Hamiltonian controls through order8 =', direct)
    print('Odd-order counts through13:', counts[1:14:2])
    e13 = [encode(w) for w in words(13)]
    need(sum((4, 6) in e for e in e13) == 63, '63 marked-edge stratum')
    need(sum((4, 6) not in e for e in e13) == 210, '210 complement stratum')
    need(sum((2, 4) in e for e in e13) == 105, 'different edge gives105')
    need(hamiltonian(rows_for(('b', 'b'))) == 9 and
         hamiltonian(rows_for(('c', 'a', 'a'))) == 5, 'equal c3 unequal H hostile')
    print('Order13: marked {4,6} gives63+210; marked {2,4} gives105+168.')

    for n in range(1, 201):
        sums = [(a, n - a) for a in range(1, n) if a <= n - a]
        factors = [(a, n // a) for a in range(1, n + 1) if n % a == 0 and a <= n // a]
        tau = sum(n % a == 0 for a in range(1, n + 1))
        need(len(sums) == n // 2, 'sum swap quotient')
        need(len(factors) == (tau + (isqrt(n) ** 2 == n)) // 2, 'product swap quotient')
    for a in range(1, 61):
        for b in range(1, 61):
            s, p = a + b, a * b
            d = isqrt(s * s - 4 * p)
            need(((s - d) // 2, (s + d) // 2) == tuple(sorted((a, b))), 'pair recovery')
            x, y, address = a, b, []
            while x != y:
                old_height = x + y
                if x > y:
                    x -= y
                    address.append('L')
                else:
                    y -= x
                    address.append('R')
                need(x + y < old_height, 'Euclidean strict rank')
            need(x == gcd(a, b), 'gcd root')
            for letter in reversed(address):
                if letter == 'L': x += y
                else: y += x
            need((x, y) == (a, b), 'retained operation address')
    need(sum((1, 6, 6)) == sum((2, 2, 9)) and 1 * 6 * 6 == 2 * 2 * 9,
         'higher-arity scalar collision')
    print('Swap-fixed fibres through200 and exact half-plane inverses for3600 pairs pass.')

    X, Y, I = ((1, 1), (0, 1)), ((1, 0), (-1, 1)), ((1, 0), (0, 1))
    need(mm(mm(X, Y), X) == mm(mm(Y, X), Y), 'braid relation')
    M = I
    strands, crossings = [0, 1, 2], Counter()
    for generator in [0, 1] * 6:
        M = mm(M, X if generator == 0 else Y)
        pair = tuple(sorted(strands[generator:generator + 2]))
        crossings[pair] += 1
        strands[generator], strands[generator + 1] = strands[generator + 1], strands[generator]
    need(M == I and strands == [0, 1, 2], 'trivial matrix and strand permutation')
    need(set(crossings.values()) == {4} and len(crossings) == 3, 'linking hostile')
    print('Braid matrix kernel: exponent12, matrixI, linking numbers(2,2,2), not unlink.')
    print('PASS: full retained maps, gradings and hostile quotient controls.')


if __name__ == '__main__':
    main()
