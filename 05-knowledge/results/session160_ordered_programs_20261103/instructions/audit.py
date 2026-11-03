from itertools import permutations, product
from math import isqrt

def edges_trans(n):
    return {(i, j) for i in range(n) for j in range(i + 1, n)}

def wrap(E, n):
    a, b = (n, n + 1)
    return (E | {(i, a) for i in range(n)} | {(a, b)} | {(b, i) for i in range(n)}, n + 2)

def hs(E, vertices):
    return sum((all(((a, b) in E for a, b in zip(p, p[1:]))) for p in permutations(vertices)))

def cover(E, n, r):
    ans = 0
    for assignment in product(range(r), repeat=n):
        value = 1
        for i in range(r):
            vertices = [j for j in range(n) if assignment[j] == i]
            value *= hs(E, vertices)
        ans += value
    return ans

def strongly_connected(E, n):
    for x in range(n):
        reachable = {x}
        while True:
            new = reachable | {v for u, v in E if u in reachable}
            if new == reachable:
                break
            reachable = new
        if len(reachable) != n:
            return False
    return True

def make_apex(n):
    word = bin(n)[2:] + '0'
    m = len(word)
    E = edges_trans(m)
    E |= {(m, i) if c == '1' else (i, m) for i, c in enumerate(word)}
    return (E, m + 1)
for n in range(1, 6):
    E, size = wrap(edges_trans(n), n)
    if not (strongly_connected(E, size) and hs(E, range(size)) == 1 + 2 ** n):
        raise AssertionError()
    if n <= 3:
        for r in range(1, 4):
            actual = cover(E, size, r)
            predicted = r ** (n + 2) + r * ((r + 1) ** n - (r - 1) ** n)
            if not actual == predicted:
                raise AssertionError((n, r, actual, predicted))
    for k in range(2 if n <= 2 else 0):
        A = E | {(i, j) for i in range(size) for j in range(size, size + k)} | {(i, j) for i in range(size, size + k) for j in range(i + 1, size + k)}
        W, v = wrap(A, size + k)
        expected = 1 + 2 ** n + 2 ** k * (2 ** (n + 2) + 2 * (3 ** n - 1))
        if not (strongly_connected(W, v) and hs(W, range(v)) == expected):
            raise AssertionError((n, k))
H = {}
for n in range(1, 16):
    E, v = make_apex(n)
    if not strongly_connected(E, v):
        raise AssertionError()
    H[n] = hs(E, range(v))
for n in range(1, 8):
    h = H[n]
    t = n % 2
    if not H[2 * n] == h + 2 * n:
        raise AssertionError()
    if not H[2 * n + 1] == h + n + 2 - t:
        raise AssertionError()
if not [H[i] for i in (2, 3, 4, 6)] == [5, 5, 9, 11]:
    raise AssertionError()
for n in range(1, 101):
    x = 1 + 2 ** n
    if not 1 + (x - 1) ** 2 == 1 + 2 ** (2 * n):
        raise AssertionError()
    if not 1 + 2 * (x - 1) ** 2 == 1 + 2 ** (2 * n + 1):
        raise AssertionError()
for x in range(1, 200):
    y = 2 * x + 1
    k = x.bit_length() - 1
    if not (x - (1 << k) >= 0 and x - (1 << k) < 1 << k):
        raise AssertionError()
    if not 1 + sum((1 << i + 1 for i in range(x.bit_length()) if x >> i & 1)) == y:
        raise AssertionError()
print('PASS independent permutation/assignment enumerator: E(n<=5), all profiles n<=3 r<=3, W(n<=2,k<=1), apex n<=15 and digit/dyadic arithmetic')
