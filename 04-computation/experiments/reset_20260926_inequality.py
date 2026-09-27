"""Exact split-product inequalities and actual first-collision hostiles."""
from fractions import Fraction
from hashlib import sha256
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def v2(n):
    need(n > 0, 'positive valuation input')
    return (n & -n).bit_length() - 1


def U(n):
    return (3 * n + 1) >> v2(3 * n + 1)


def canonical(n):
    if n == 1:
        return None
    r = n.bit_length() - 1
    return (1, r, n - (1 << r))


def value(state):
    q, r, u = state
    return (q << r) + u


def product(state):
    if state is None:
        return 0
    q, r, u = state
    return (q * u) << r


def advance(state):
    q, r, u = state
    a = v2(3 * u + 1)
    v = (3 * u + 1) >> a
    if a < r:
        return 'consume', (3 * q, r - a, v), None, a
    if a > r:
        return 'swap', (v, a - r, 3 * q), None, r
    x = (3 * q + v) >> v2(3 * q + v)
    return 'collision', None, x, r + v2(3 * q + v)


def first_collision(n, cap=10000):
    state = canonical(n)
    trace = []
    for step in range(1, cap + 1):
        kind, following, x, m = advance(state)
        trace.append((kind, state, m))
        if kind == 'collision':
            return x, trace
        state = following
    raise RuntimeError(f'finite first-collision cap exceeded at {n}')


def main():
    counts = dict(consume=0, swap=0, collision=0)
    balanced = 0
    max_balance_ratio = Fraction(0)
    for q in range(1, 64, 2):
        for u in range(1, 128, 2):
            for r in range(1, 11):
                old = (q, r, u)
                y, k = value(old), product(old)
                kind, new, x, m = advance(old)
                counts[kind] += 1
                if kind != 'collision':
                    need(value(new) == U(y), 'noncollision source preservation')
                    ratio = Fraction(product(new), k)
                    need(ratio == Fraction(9 * u + 3, u * (1 << (2 * m))),
                         'exact product law across consume and swap')
                    if m >= 2:
                        need(ratio <= Fraction(3, 4), 'large-valuation contraction')
                    else:
                        need(ratio <= 3, 'one-valuation expansion bound')
                    p = Fraction(u, y)
                    p_new = Fraction(new[2], value(new))
                    p_expected = Fraction(3 * u + 1, 3 * y + 1)
                    need(p_new == (p_expected if kind == 'consume' else 1 - p_expected),
                         'projective role law')
                    need(abs(abs(2 * p_new - 1) - abs(2 * p - 1)) <= Fraction(2, 3 * y + 1),
                         'balance variation bound')
                else:
                    need(x == U(y) and x < y and m >= 2, 'collision contracts actual integer')
                    reset = canonical(x)
                    ratio = Fraction(product(reset), k)
                    threshold = Fraction(3 * (y - 1) * (y + 3), 16)
                    need(ratio <= Fraction(3 * y + 5, 4 * (y + 3)) * threshold / k,
                         'reset imbalance charge')
                    if y <= 4 * u <= 3 * y:
                        balanced += 1
                        need(k >= threshold, 'balanced integer split lower bound')
                        need(ratio <= Fraction(3 * y + 5, 4 * (y + 3)) < Fraction(3, 4),
                             'sharp balanced-collision reset inequality')
                        max_balance_ratio = max(max_balance_ratio, ratio)

    sharp_rows = []
    for k in range(1, 81):
        y = ((1 << (2 * k + 3)) - 5) // 3
        old = (3 * (y - 1) // 8, 1, (y + 3) // 4)
        need(value(old) == y, 'sharp source')
        kind, _, x, m = advance(old)
        need(kind == 'collision' and m == 2 and x == (1 << (2 * k + 1)) - 1,
             'sharp collision arithmetic')
        ratio = Fraction(product(canonical(x)), product(old))
        need(ratio == Fraction(3 * y + 5, 4 * (y + 3)), 'sharp attainable ratio')
        if k <= 5:
            sharp_rows.append((k, y, str(ratio)))

    hostile_rows = []
    for h in range(2, 202, 2):
        n = (1 << h) - 1
        x, trace = first_collision(n)
        need(len(trace) == h, 'all-ones first-collision clock')
        need([t[0] for t in trace] == ['consume'] * (h - 2) + ['swap', 'collision'],
             'all-ones exact transition word')
        need(x == (3 ** h - 1) >> (2 + v2(h)), 'all-ones exact endpoint')
        q, r, u = trace[-1][1]
        need((q, r, u) == ((3 ** (h - 1) - 1) // 2, 1, 3 ** (h - 1)),
             'all-ones balanced collision state')
        need(value(trace[-1][1]) <= 4 * u <= 3 * value(trace[-1][1]),
             'all-ones collision is balanced')
        if h in (2, 6, 10, 14, 18):
            hostile_rows.append((h, n, x, str(Fraction(x, n))))

    # Unbalanced reset creation in one exact collision, with an explicitly
    # noncanonical input representation (not a canonical-episode theorem).
    unbalanced = []
    for ell in range(3, 101):
        q = 3 * (1 << ell) + 1
        state = (q, 2, 1)
        kind, _, x, _ = advance(state)
        need(kind == 'collision' and x == 9 * (1 << (ell - 2)) + 1,
             'unbalanced collision')
        new_k = (1 << (ell + 1)) * ((1 << (ell - 2)) + 1)
        need(product(canonical(x)) == new_k, 'unbounded product creation')
        if ell in (3, 10, 30, 100):
            unbalanced.append((ell, str(Fraction(new_k, product(state)))))

    # Finite observation only: actual first-collision episodes under the
    # immediate canonical convention. No global collision-existence claim.
    episodes = []
    for n in range(3, 4096, 2):
        x, trace = first_collision(n)
        episodes.append((n, x, len(trace)))
    summary = dict(
        status='Exact finite controls; all-height formulas proved in accompanying note',
        legal_triple_count=sum(counts.values()),
        branches=counts,
        balanced_collision_count=balanced,
        maximum_finite_balanced_product_ratio=str(max_balance_ratio),
        sharp_family_controls=80,
        sharp_family_examples=sharp_rows,
        all_ones_family_controls=100,
        all_ones_examples=hostile_rows,
        unbalanced_collision_examples=unbalanced,
        finite_episode_count=len(episodes),
        growing_episode_count=sum(x > n for n, x, _ in episodes),
        maximum_finite_episode_length=max(j for _, _, j in episodes),
        source27_episode=next(row for row in episodes if row[0] == 27),
        episode_digest=sha256(json.dumps(episodes, separators=(',', ':')).encode()).hexdigest(),
    )
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
