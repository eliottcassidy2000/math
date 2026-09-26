"""Independent exact controls for KSBFT_v7 and escaping Fibonacci boundaries.

Standard library only. All checks remain active under python -O.
The finite universe does not certify the cited external theorems.
"""
from fractions import Fraction as Q
from itertools import combinations, product


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def closure(pred):
    pred = list(pred)
    for j in range(len(pred)):
        for i in range(j):
            if pred[j] >> i & 1:
                pred[j] |= pred[i]
    return tuple(pred)


def natural_posets(n):
    pairs = list(combinations(range(n), 2))
    found = set()
    for mask in range(1 << len(pairs)):
        pred = [0] * n
        for bit, (i, j) in enumerate(pairs):
            if mask >> bit & 1:
                pred[j] |= 1 << i
        found.add(closure(pred))
    return sorted(found)


def extensions(pred):
    n = len(pred)
    result = []

    def visit(word, used):
        if len(word) == n:
            ranks = [0] * n
            for rank, v in enumerate(word, 1):
                ranks[v] = rank
            result.append(tuple(ranks))
            return
        for v in range(n):
            if not (used >> v & 1) and pred[v] & ~used == 0:
                visit(word + [v], used | (1 << v))

    visit([], 0)
    return result


def comparable(pred, a, b):
    return bool((pred[a] >> b & 1) or (pred[b] >> a & 1))


def connected_incomparability(pred):
    n = len(pred)
    reached, queue = {0}, [0]
    while queue:
        v = queue.pop()
        for u in range(n):
            if u not in reached and not comparable(pred, u, v):
                reached.add(u)
                queue.append(u)
    return len(reached) == n


def small_posets():
    counts = dict(posets=0, extensions=0, maximal_triples=0,
                  variance_triangles=0, ideal_cuts=0, antichains=0)
    for n in range(2, 7):
        allp = natural_posets(n)
        for pred in allp:
            le = extensions(pred)
            size = len(le)
            counts['posets'] += 1
            counts['extensions'] += size
            h = [sum(f[v] for f in le) for v in range(n)]
            window = [sum(f[v] - max([f[u] for u in range(n) if pred[v] >> u & 1], default=0)
                          for f in le) for v in range(n)]
            # Denominator of every scaled-order-polytope variance is (n+2)*size^2.
            variance = [[0] * n for _ in range(n)]
            for a, b in combinations(range(n), 2):
                gap_moment = sum((f[a] - f[b]) ** 2 + abs(f[a] - f[b]) for f in le)
                value = (n + 1) * size * gap_moment - (n + 2) * (h[a] - h[b]) ** 2
                check(value > 0, 'nondegenerate continuous difference')
                variance[a][b] = variance[b][a] = value
                check(3 * value >= (n + 2) * window[a] ** 2, 'window variance a')
                check(3 * value >= (n + 2) * window[b] ** 2, 'window variance b')
            for a, b, c in product(range(n), repeat=3):
                check(variance[a][c] <= variance[a][b] + variance[b][c], 'squared triangle')
                counts['variance_triangles'] += 1
            if n >= 3 and connected_incomparability(pred):
                order = sorted(range(n), key=lambda v: (h[v], v))
                displacements = [h[v] - size * (i + 1) for i, v in enumerate(order)]
                maximum = max(map(abs, displacements))
                eligible = []
                for k in range(n - 2):
                    triple = order[k:k + 3]
                    chain = all(comparable(pred, a, b) for a, b in combinations(triple, 2))
                    if (not chain and h[triple[2]] - h[triple[0]] <= 2 * size
                            and max(abs(displacements[k]), abs(displacements[k + 2])) == maximum):
                        eligible.append(k)
                check(eligible, 'maximal displacement BFT triple')
                counts['maximal_triples'] += 1
            for mask in range(1, 1 << n):
                vertices = [v for v in range(n) if mask >> v & 1]
                if all(not comparable(pred, a, b) for a, b in combinations(vertices, 2)):
                    check(2 * sum(window[v] for v in vertices) >= size * len(vertices) ** 2,
                          'antichain window bound')
                    counts['antichains'] += 1
                if mask == (1 << n) - 1 or any(pred[v] & ~mask for v in vertices):
                    continue
                outside = [v for v in range(n) if not (mask >> v & 1)]
                maxima = [v for v in vertices if not any(pred[u] >> v & 1 for u in vertices)]
                minima = [v for v in outside if not any(pred[v] >> u & 1 for u in outside)]
                check(min(h[v] for v in outside) - max(h[v] for v in vertices)
                      <= size * (len(maxima) + len(minima) - 1), 'ideal cut inequality')
                counts['ideal_cuts'] += 1
        print('naturally labelled posets', n, len(allp))
    print('small-poset exact controls', counts)


def fibonacci_controls():
    fib = [0, 1]
    for _ in range(505):
        fib.append(fib[-1] + fib[-2])
    direct = 0
    for n in range(2, 14):
        pred = tuple(sum(1 << j for j in range(i - 1)) for i in range(n))
        le = extensions(pred)
        check(len(le) == fib[n + 1], 'Fibonacci extension count')
        edge_p = [Q(0)] + [Q(fib[j] * fib[n - j], fib[n + 1]) for j in range(1, n)] + [Q(0)]
        for j in range(1, n):
            actual = Q(sum(f[j] < f[j - 1] for f in le), len(le))
            check(actual == edge_p[j], 'Fibonacci pair probability')
            direct += 1
        for j in range(1, n + 1):
            actual = Q(sum(f[j - 1] for f in le), len(le)) - j
            check(actual == edge_p[j] - edge_p[j - 1], 'Fibonacci displacement')
            direct += 1
    boundary_cases = 0
    for m in range(1, 251):
        denominator = fib[2 * m + 2]
        displacements = [Q(fib[2 * abs(i)], denominator) for i in range(-m, m + 1)]
        check(sum(displacements) == Q(2 * (fib[2 * m + 1] - 1), denominator), 'exact L1 mass')
        check(max(displacements) == Q(fib[2 * m], denominator), 'boundary maximum')
        for i, displacement in zip(range(-m, m + 1), displacements):
            distance = m - abs(i)
            check(displacement <= Q(1, 2 ** (distance + 1)), 'boundary localization')
            boundary_cases += 1
        if m in (2, 5, 20, 100, 250):
            center_pair = Q(fib[m + 1] * fib[m], denominator)
            print('Fibonacci m', m, 'center pair', float(center_pair),
                  'max displacement', float(max(displacements)),
                  'normalized L1', float(sum(displacements) / (2 * m + 1)))
    print('Fibonacci direct count/moment controls', direct, 'boundary controls', boundary_cases)


def structural_gluing():
    models = []
    for n in range(1, 4):
        for pred in natural_posets(n):
            for x in range(n):
                if not any(pred[v] >> x & 1 for v in range(n)):
                    models.append((pred, x))
    cases = 0
    for left, x in models:
        nl = len(left)
        le_left = extensions(left)
        # Dualize an independent left model to supply a right model with a marked minimum.
        for original, marked in models:
            nr = len(original)
            right = tuple(sum(1 << (nr - 1 - j) for j in range(nr) if original[j] >> (nr - 1 - i) & 1)
                          for i in range(nr))
            z = nr - 1 - marked
            le_right = extensions(right)
            y = nl
            pred = list(left) + [((1 << nl) - 1) ^ (1 << x)]
            for v in range(nr):
                pred.append(((1 << nl) - 1) | (right[v] << (nl + 1)) | (0 if v == z else 1 << y))
            le = extensions(closure(pred))
            r = Q(sum(f[x] == nl for f in le_left), len(le_left))
            s = Q(sum(f[z] == 1 for f in le_right), len(le_right))
            Z = 1 + r + s
            check(len(le) == len(le_left) * len(le_right) * Z, 'gluing normalizer')
            check(Q(sum(f[y] < f[x] for f in le), len(le)) == r / Z, 'left inversion')
            check(Q(sum(f[nl + 1 + z] < f[y] for f in le), len(le)) == s / Z, 'right inversion')
            t = Q(sum(nl - f[x] for f in le_left), len(le_left))
            u = Q(sum(f[z] - 1 for f in le_right), len(le_right))
            expected = [nl + r / Z - (1 + s) * t / Z,
                        nl + 1 + (s - r) / Z,
                        nl + 2 + (1 + r) * u / Z - s / Z]
            actual = [Q(sum(f[v] for f in le), len(le)) for v in [x, y, nl + 1 + z]]
            check(expected == actual, 'three height identities')
            cases += 1
    print('rigid triple gluing exact controls', cases)


def algebra_controls():
    count = 0
    for i, j, k in product(range(21), repeat=3):
        b1, b2, b3 = Q(i, 40), Q(j, 40), Q(k, 40)
        s = 2 * b1 + b2
        if s > 1 or b2 + b3 > Q(1, 2) or 4 * b1 * b3 > b2 ** 2:
            continue
        lhs = 7 * b1 + 3 * b2 + 3 * b3
        q = 4 * lhs - 6 + 7 * s
        check(q <= 0 or q ** 2 <= 245 * s ** 2, 'Case D convex chord')
        count += 1
    # The claimed appendix lower bound factors exactly; all grid checks use rationals.
    for j in range(101):
        v = Q(1, 2) + Q(j, 200)
        difference = (-3 * v ** 2 + 6 * v - 2) / (2 * v * (1 + v) * (1 + 2 * v))
        check(difference >= Q(1, 12), 'appendix A.5')
    check(5 - 11 * Q(691, 2500) + Q(1, 22) > 2, 'appendix final contradiction')
    print('Case D exact rational feasible grid', count, 'appendix grid', 101)


if __name__ == '__main__':
    small_posets()
    fibonacci_controls()
    structural_gluing()
    algebra_controls()
    print('PASS: exact finite controls only; no claimed verification of cited upstream theorems or global epsilon numeral')
