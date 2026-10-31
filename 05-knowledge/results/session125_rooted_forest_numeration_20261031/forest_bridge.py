"""Exact rooted-path-forest / Chung-Graham / odd-tiling bridge.

Python standard library only; 0-based graph vertices and cell indices.
Words are MSB first, forest blocks are (length, 1-based root position).
"""
from fractions import Fraction
from math import isqrt


def legal(word):
    if any(x not in (0, 1, 2) for x in word):
        return False
    locked = False
    for digit in word:
        if digit == 0:
            locked = False
        elif digit == 2:
            if locked:
                return False
            locked = True
    return True


def forest_to_tiling(blocks):
    if not blocks or any(k < 1 or r < 1 or r > k for k, r in blocks):
        raise ValueError('nonempty rooted path forest required')
    tiles = []
    for j, (k, r) in enumerate(blocks):
        if j:
            tiles.append(1)
        tiles += [2] * (r - 1) + [1] + [2] * (k - r)
    return tuple(tiles)


def tiling_to_forest(tiles):
    if not tiles or any(t not in (1, 2) for t in tiles) or sum(tiles) % 2 != 1:
        raise ValueError('positive odd-length monomer/domino tiling required')
    blocks, left_d, right_d, rooted = [], 0, 0, False
    for tile in tiles:
        if tile == 2:
            if rooted:
                right_d += 1
            else:
                left_d += 1
        elif not rooted:
            rooted = True
        else:
            blocks.append((left_d + right_d + 1, left_d + 1))
            left_d, right_d, rooted = 0, 0, False
    if not rooted:
        raise ValueError('odd monomer parity violated')
    blocks.append((left_d + right_d + 1, left_d + 1))
    return tuple(blocks)


def word_to_tiling(word):
    word = tuple(word)
    if not legal(word):
        raise ValueError('illegal CG word')
    n = 2 * len(word) + 1
    placed = {}
    locked = False
    for i, digit in enumerate(reversed(word)):
        a, b = n - 2 * i - 2, n - 2 * i - 1
        if not locked:
            if digit == 0:
                placed[a] = placed[b] = 1
            elif digit == 1:
                placed[a] = 2
            else:
                placed[b], placed[a - 1] = 1, 2
                locked = True
        elif digit == 0:
            placed[a] = 1
            locked = False
        elif digit == 1:
            placed[a - 1] = 2
        else:
            raise ValueError('illegal CG transition')
    if not locked:
        placed[0] = 1
    tiles, pos = [], 0
    while pos < n:
        if pos not in placed:
            raise RuntimeError('uncovered cell')
        tiles.append(placed[pos])
        pos += placed[pos]
    if pos != n or len(tiles) != len(placed):
        raise RuntimeError('overlapping or overrun tiles')
    return tuple(tiles)


def tiling_to_word(tiles):
    tiles = tuple(tiles)
    if not tiles or any(t not in (1, 2) for t in tiles) or sum(tiles) % 2 != 1:
        raise ValueError('positive odd tiling required')
    n = sum(tiles)
    partner = {}
    pos = 0
    for tile in tiles:
        partner[pos] = pos if tile == 1 else pos + 1
        if tile == 2:
            partner[pos + 1] = pos
        pos += tile
    digits = []
    locked = False
    for i in range((n - 1) // 2):
        a, b = n - 2 * i - 2, n - 2 * i - 1
        if not locked:
            if partner[b] == a:
                digit = 1
            elif partner[a] == a:
                digit = 0
            else:
                digit = 2
                locked = True
        elif partner[a] == a:
            digit = 0
            locked = False
        else:
            digit = 1
        digits.append(digit)
    return tuple(reversed(digits))


def word_to_forest(word):
    return tiling_to_forest(word_to_tiling(word))


def forest_to_word(blocks):
    return tiling_to_word(forest_to_tiling(blocks))


def det_int(matrix):
    """Fraction-free Bareiss determinant, including determinant([])=1."""
    a = [list(row) for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError('square matrix required')
    if n == 0:
        return 1
    denom, sign = 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                if numerator % denom:
                    raise RuntimeError('nonexact Bareiss division')
                a[i][j] = numerator // denom
            a[i][k] = 0
        denom = pivot
    return sign * a[-1][-1]


def identity_laplacian(n, edges):
    a = [[int(i == j) for j in range(n)] for i in range(n)]
    for u, v in edges:
        if not 0 <= u < v < n:
            raise ValueError('distinct sorted vertex pairs required')
        a[u][u] += 1
        a[v][v] += 1
        a[u][v] -= 1
        a[v][u] -= 1
    return a


def solve(matrix, rhs):
    n = len(matrix)
    a = [[Fraction(x) for x in row] + [Fraction(rhs[i])]
         for i, row in enumerate(matrix)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            raise ValueError('singular')
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][k]
        a[k] = [x / p for x in a[k]]
        for i in range(n):
            if i != k:
                q = a[i][k]
                a[i] = [x - q*y for x, y in zip(a[i], a[k])]
    return [row[-1] for row in a]


def arrival_count(n, edges, ports):
    ports = tuple(ports)
    if len(set(ports)) != len(ports) or any(not 0 <= v < n for v in ports):
        raise ValueError('ports must be distinct old vertices')
    m = identity_laplacian(n, edges)
    b = [int(i in ports) for i in range(n)]
    for v in ports:
        m[v][v] += 1
    x = solve(m, b)
    result = det_int(m) * (1 + len(ports) - sum(bi*xi for bi, xi in zip(b, x)))
    if result.denominator != 1:
        raise RuntimeError('forest count not integral')
    return int(result)


def square_sum_edges(n):
    squares = {i*i for i in range(2, isqrt(2*n) + 1)}
    return tuple((u-1, v-1) for u in range(1, n+1) for v in range(u+1, n+1)
                 if u+v in squares)
