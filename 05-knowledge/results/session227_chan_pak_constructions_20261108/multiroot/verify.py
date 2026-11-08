#!/usr/bin/env python3
"""Exact checks for set-root grafts; standard library; guards survive python -O."""
import itertools
import json
import math
from fractions import Fraction


def require(test, msg):
    if not test:
        raise RuntimeError(msg)


def valid(a):
    n = len(a)
    require(n > 0, 'nonempty tournament')
    for i, row in enumerate(a):
        require(0 <= row < 1 << n and not (row >> i & 1), 'row')
        for j in range(i):
            require(((row >> j) & 1) + ((a[j] >> i) & 1) == 1, 'tournament arc')


def tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for word in range(1 << len(pairs)):
        a = [0] * n
        for k, (i, j) in enumerate(pairs):
            if word >> k & 1:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
        yield a


def paths(a, root=None):
    """All induced Hamilton counts, plus root-start counts if requested."""
    valid(a)
    n, size = len(a), 1 << len(a)
    d = [[0] * n for _ in range(size)]
    s = [[0] * n for _ in range(size)] if root is not None else None
    for v in range(n):
        d[1 << v][v] = 1
    if s is not None:
        s[1 << root][root] = 1
    for mask in range(1, size):
        for u in range(n):
            x = d[mask][u]
            y = s[mask][u] if s is not None else 0
            if not (x or y):
                continue
            rest = a[u] & ~mask
            while rest:
                b = rest & -rest
                rest -= b
                v = b.bit_length() - 1
                d[mask | b][v] += x
                if s is not None:
                    s[mask | b][v] += y
    hp = [sum(row) for row in d]
    hp[0] = 1
    if root is None:
        return hp
    full = size - 1
    ss = [sum(row) for row in s]
    return [hp[full], d[full][root],
            sum(ss[m] * hp[full ^ m] for m in range(size)),
            hp[full ^ (1 << root)]]


def ordered_nonempty_covers(a):
    """Permutation-and-cut count, separate from the Hamilton subset DP."""
    n = len(a)
    out = [0] * (n + 1)
    for p in itertools.permutations(range(n)):
        bad = sum(not (a[u] >> v & 1) for u, v in zip(p, p[1:]))
        for k in range(bad + 1, n + 1):
            out[k] += math.comb(n - 1 - bad, k - 1 - bad)
    return out


def colours_from_subsets(a, max_j):
    """Allow empty ordered bins; subset convolution independent of covers()."""
    hp = paths(a)
    full = (1 << len(a)) - 1
    previous = [0] * (full + 1)
    previous[0] = 1
    out = [previous[full]]
    for j in range(1, max_j + 1):
        now = [0] * (full + 1)
        for m in range(full + 1):
            sub = m
            while True:
                now[m] += hp[sub] * previous[m ^ sub]
                if not sub:
                    break
                sub = (sub - 1) & m
        previous = now
        out.append(now[full])
    return out


def coefficients(a, r):
    """A-order and B-slot enumeration: backward A gaps are mandatory slots."""
    n = len(a)
    require(0 <= r < 1 << n, 'root-set mask')
    limit = min(r.bit_count(), n - r.bit_count()) + 1
    out = [0] * (n + 2)
    for p in itertools.permutations(range(n)):
        mandatory = 0
        optional = int(bool(r >> p[0] & 1)) + int(not (r >> p[-1] & 1))
        compatible = True
        for u, v in zip(p, p[1:]):
            eligible = not (r >> u & 1) and bool(r >> v & 1)
            if a[u] >> v & 1:
                optional += int(eligible)
            elif eligible:
                mandatory += 1
            else:
                compatible = False
                break
        if compatible:
            for ell in range(max(1, mandatory), mandatory + optional + 1):
                out[ell] += math.comb(optional, ell - mandatory)
    require(all(v == 0 for v in out[limit + 1:]), 'sharp interface degree')
    return out[:limit + 1]


def spectral(c):
    return [sum((-1) ** (ell - j) * math.comb(ell, j) * c[ell]
                for ell in range(j, len(c))) for j in range(len(c))]


def graft(a, r, b):
    n = len(a)
    new = ((1 << len(b)) - 1) << n
    out = [row | (new if not (r >> i & 1) else 0) for i, row in enumerate(a)]
    out.extend((row << n) | r for row in b)
    valid(out)
    return out


def join(*blocks):
    out = []
    total = sum(map(len, blocks))
    shift = 0
    for b in blocks:
        later = ((1 << total) - 1) ^ ((1 << (shift + len(b))) - 1)
        out.extend((row << shift) | later for row in b)
        shift += len(b)
    valid(out)
    return out


def transitive(n):
    return [((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n)]


def reachable(a, v):
    seen = todo = 1 << v
    while todo:
        bit = todo & -todo
        todo -= bit
        new = a[bit.bit_length() - 1] & ~seen
        seen |= new
        todo |= new
    return seen


def strong(a):
    full = (1 << len(a)) - 1
    return all(reachable(a, v) == full for v in range(len(a)))


def criterion(a, r):
    """Equivalent to first SCC meeting R and last SCC meeting its complement."""
    full = (1 << len(a)) - 1
    c = full ^ r
    reached = 0
    for v in range(len(a)):
        if r >> v & 1:
            reached |= reachable(a, v)
    return reached == full and all(reachable(a, v) & c for v in range(len(a)))


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    answer = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            answer = -answer
        p = a[i][i]
        answer *= p
        for j in range(i + 1, len(a)):
            f = a[j][i] / p
            for k in range(i + 1, len(a)):
                a[j][k] -= f * a[i][k]
    require(answer.denominator == 1, 'integral determinant')
    return int(answer)


def run():
    bank = [a for n in range(1, 5) for a in tournaments(n)]
    cover_bank = [ordered_nonempty_covers(a) for a in bank]
    q_bank = [colours_from_subsets(a, 5) for a in bank]
    for a, ff, qq in zip(bank, cover_bank, q_bank):
        for j in range(1, 6):
            require(qq[j] == sum(math.comb(j, ell) * ff[ell]
                                 for ell in range(1, min(j, len(a)) + 1)),
                    'binomial cover transform')
        require(ff[1] == paths(a)[-1], 'independent Hamilton check')
        if len(a) >= 2:
            require(ff[2] >= (1 << len(a)) - 2, 'two-bin floor')
    checked = strong_checked = 0
    coefficient_rows = []
    for a in bank:
        for r in range(1 << len(a)):
            cc = coefficients(a, r)
            ll = spectral(cc)
            guard = criterion(a, r)
            if guard:
                require(len(cc) > 2 and cc[2] >= 1 and cc[1] >= 3, 'strong source c1/c2 positivity')
            coefficient_rows.append([a, r, cc])
            for b, ff, qq in zip(bank, cover_bank, q_bank):
                if len(a) == 4 and len(b) == 4:
                    continue
                t = graft(a, r, b)
                h = paths(t)[-1]
                require(h == sum(cc[i] * ff[i] for i in range(1, min(len(cc), len(ff)))),
                        'positive graft cover law')
                require(h == sum(ll[j] * qq[j] for j in range(1, len(ll))),
                        'spectral graft cover law')
                require(strong(t) == guard, 'strongness iff')
                if guard:
                    require(h >= (1 << len(b)) + 1, 'sharp strong gadget order bound')
                    strong_checked += 1
                checked += 1
    sharp = []
    for r in range(1, 5):
        for c in range(1, 5):
            a = transitive(r + c)
            cc = coefficients(a, (1 << r) - 1)
            m, other = min(r, c), max(r, c)
            # Ordered surjections onto m bins by elementary inclusion-exclusion.
            onto = sum((-1) ** (m - j) * math.comb(m, j) * j ** other
                       for j in range(m + 1))
            require(cc[-1] >= math.factorial(m) * onto, 'sharp leading coefficient lower bound')
            if r == c:
                require(cc[-1] == math.factorial(m) * onto, 'equal-shore leading coefficient')
            sharp.append({'r': r, 'c': c, 'coefficients': cc, 'spectral': spectral(cc)})
    triangle = [2, 4, 1]
    four = [8, 1, 3, 6]
    # Exact inherited Q values verified rather than assumed.
    require(colours_from_subsets(triangle, 3) == [0, 3, 12, 33], 'triangle colour polynomial')
    require(colours_from_subsets(four, 3) == [0, 5, 32, 117], 'four-vertex colour polynomial')
    products = []
    for left in bank[:11]:
        for right in bank[:11]:
            direct = colours_from_subsets(join(left, right), 3)
            ql, qr = colours_from_subsets(left, 3), colours_from_subsets(right, 3)
            require(direct == [ql[j] * qr[j] for j in range(4)], 'ordered-join product')
            products.append(direct)
    observation_ranks = []
    for r in range(1, 5):
        contexts = [(transitive(2), 1), (transitive(3), 1)]
        contexts += [(transitive(2 * m), (1 << m) - 1) for m in range(2, r + 1)]
        response = [[paths(graft(a, roots, transitive(n)))[-1]
                     for n in range(1, r + 2)] for a, roots in contexts]
        expected = 4
        for m in range(2, r + 1):
            expected *= math.factorial(m) ** 2
        moment_matrix = [[j ** n for n in range(1, r + 2)] for j in range(1, r + 2)]
        expected *= determinant(moment_matrix)
        require(determinant(response) == expected and expected != 0,
                'physical full-rank gadget observer')
        observation_ranks.append({'r': r, 'dimension': r + 1,
                                  'response_matrix': response, 'determinant': expected})
    three_mode = []
    aa, rr = transitive(4), 5
    cc = coefficients(aa, rr)
    require(cc == [0, 9, 10, 2], 'explicit alternating-shell coefficients')
    require(spectral(cc)[1:] == [-5, 4, 2], 'three-mode coefficients')
    for m in range(1, 4):
        bb = join(*([triangle] * m))
        tt = graft(aa, rr, bb)
        value = -5 * 3 ** m + 4 * 12 ** m + 2 * 33 ** m
        require(paths(tt)[-1] == value and strong(tt), 'literal three-mode family')
        three_mode.append({'m': m, 'H': value, 'rows': tt})
    examples = []
    for k, ell in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 1)]:
        words = sorted(set(itertools.permutations('C' * k + 'F' * ell)))
        h, q = 3 ** k * 5 ** ell, 12 ** k * 32 ** ell
        expected = [3 * h + 12 * q, 3 * q, 6 * h + 14 * q - int(k + ell == 0), h + 4 * q]
        witnesses = []
        for word in words:
            blocks = [triangle if x == 'C' else four for x in word] + [triangle]
            b = join(*blocks)
            t = graft([2, 0], 1, b)
            root = len(t) - 3
            require(paths(t, root) == expected, 'physical collision family full state')
            require(strong(t), 'collision family strong')
            n = len(t)
            outdegrees = [v.bit_count() for v in t]
            require([i for i, d in enumerate(outdegrees) if d == 1] == [0], 'unique outdegree-one old root')
            require([i for i, d in enumerate(outdegrees) if d == n - 2] == [1], 'unique indegree-one old nonroot')
            witnesses.append({'word': ''.join(word) + 'C', 'rows': t, 'root': root})
        examples.append({'k': k, 'ell': ell, 'state': expected, 'witnesses': witnesses})
    return {'status': 'PASS', 'universe': {'labelled_gadgets': len(bank),
            'source_cuts': len(coefficient_rows),
            'pairs': checked, 'strong_pairs': strong_checked,
            'pair_filter': 'orders 1..4 on both sides, excluding 4-by-4',
            'product_checks': len(products), 'sharpness_checks': len(sharp)},
            'sharpness': sharp, 'observation_ranks': observation_ranks, 'three_mode_family': three_mode, 'physical_collision_families': examples,
            'scope': 'All-order claims are proved separately; no completeness or improved scalar cutoff.'}


if __name__ == '__main__':
    print(json.dumps(run(), sort_keys=True, indent=2))
