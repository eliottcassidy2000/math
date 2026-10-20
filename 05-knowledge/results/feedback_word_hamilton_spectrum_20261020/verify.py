#!/usr/bin/env python3
"""Independent exact checks; no imported repo implementation and no huge census."""
from itertools import combinations, permutations, product
from collections import defaultdict
import json


def need(p, why):
    if not p:
        raise RuntimeError(why)


def tournament(word, k=1, feedback_code=0):
    m = len(word)
    out = [0] * (m + k)
    for i in range(m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    for i, letter in enumerate(word):
        for a in range(k):
            if letter >> a & 1:  # feedback vertex beats chain vertex
                out[m + a] |= 1 << i
            else:
                out[i] |= 1 << (m + a)
    for bit, (a, b) in enumerate(combinations(range(k), 2)):
        u, v = (a, b) if feedback_code >> bit & 1 else (b, a)
        out[m + u] |= 1 << (m + v)
    return tuple(out)


def h_dp(out):
    n = len(out)
    if not n:
        return 1
    d = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            count = d[mask][v]
            if not count:
                continue
            opts = out[v] & ~mask
            while opts:
                b = opts & -opts
                opts -= b
                d[mask | b][b.bit_length() - 1] += count
    return sum(d[-1])


def h_perm(out):
    return sum(all(out[p[i]] >> p[i + 1] & 1 for i in range(len(p) - 1))
               for p in permutations(range(len(out))))


def h_word(word):
    # Integer four-coordinate linear representation: (h,t,last,1).
    h, t, last = 1, 0, 0
    for b in word:
        if b == 0:
            h += 2 * (t + last)
        t, last = 2 * t + last, b
    return h


def weight(d):
    return 1 if d == 1 else 1 << (d - 2)


def h_cut(word):
    return 1 + 2 * sum(weight(j - i) for i in range(len(word))
                       for j in range(i + 1, len(word))
                       if word[i] == 1 and word[j] == 0)


def feedback_arc(code, a, b, k):
    lo, hi = sorted((a, b))
    bit = list(combinations(range(k), 2)).index((lo, hi))
    return bool(code >> bit & 1) == (a == lo)


def h_automaton(word, k, code):
    # 0=empty; 1=nonempty with valid last vertex; 2=invalid last vertex.
    total = 0
    for p in permutations(range(k)):
        cur = {(0,) * (k + 1): 1}
        for letter in word:
            nxt = defaultdict(int)
            for state, count in cur.items():
                for gap in range(k + 1):
                    if state[gap] == 0 and gap and not (letter >> p[gap - 1] & 1):
                        continue
                    valid_last = gap == k or not (letter >> p[gap] & 1)
                    new = list(state)
                    new[gap] = 1 if valid_last else 2
                    # Distinct gap assignments can have the same state: sum them.
                    nxt[tuple(new)] += count
            cur = nxt
        for state, count in cur.items():
            if 2 in state:
                continue
            if any(state[j] == 0 and not feedback_arc(code, p[j - 1], p[j], k)
                   for j in range(1, k)):
                continue
            total += count
    return total


def strong(out):
    if not out:
        return False
    for s in range(len(out)):
        seen = 1 << s
        old = -1
        while old != seen:
            old = seen
            for u in range(len(out)):
                if seen >> u & 1:
                    seen |= out[u]
        if seen != (1 << len(out)) - 1:
            return False
    return True


def main():
    controls = 0
    for m in range(0, 9):
        for word in product(range(2), repeat=m):
            out = tournament(word)
            value = h_dp(out)
            need(value == h_cut(word) == h_word(word), ('word formula', word))
            need(value == h_automaton(word, 1, 0), ('k1 automaton', word))
            if m <= 5:
                need(value == h_perm(out), ('permutations', word))
            if m >= 2:
                need(strong(out) == (word[0] == 1 and word[-1] == 0), ('strong', word))
            controls += 1
    print('one-apex formula/DP/automaton controls:', controls)
    records = []
    swaps = 0
    for m in range(2, 13):
        values = []
        for mid in product(range(2), repeat=m - 2):
            w = (1,) + mid + (0,)
            value = h_word(w)
            values.append(value)
            for i in range(1, m - 2):
                if w[i:i + 2] != (1, 0):
                    continue
                v = list(w)
                v[i], v[i + 1] = 0, 1
                delta = value - h_word(v)
                exact = 2 * (1 + sum(w[j] * (weight(i + 1 - j) - weight(i - j))
                                     for j in range(i))
                             + sum((1 - w[j]) * (weight(j - i) - weight(j - i - 1))
                                   for j in range(i + 2, m)))
                need(delta == exact and delta >= 2, ('swap', w, i))
                if m >= 5:
                    need(delta >= 4, ('unit barrier', w, i))
                swaps += 1
        a, b = (m - 2) // 2, (m - 2 + 1) // 2
        expected_min = (1 << (m - 1)) + 1 if m < 4 else ((1 << a) + 1) * ((1 << b) + 1)
        r, s = m // 2, (m + 1) // 2
        expected_max = ((1 << r) - 1) * ((1 << s) - 1) + 2
        need(min(values) == expected_min and max(values) == expected_max, ('extrema', m))
        records.append({'n': m + 1, 'min': min(values), 'max': max(values),
                        'distinct_values': len(set(values))})
    print('sharp extrema:', json.dumps(records, sort_keys=True))
    print('strict adjacent-swap checks:', swaps)
    need({h_word((1,) + mid + (0,)) for mid in product(range(2), repeat=2)} == {9, 11}, 'H13 obstruction')
    for a in range(1, 17):
        for b in range(1, 17):
            w = (1,) + (0,) * a + (1,) * b + (0,)
            need(h_word(w) == ((1 << a) + 1) * ((1 << b) + 1), ('product family', a, b))
    print('strong product-family checks: 256')
    witness13 = (22, 12, 24, 17, 2)
    need(h_dp(witness13) == h_perm(witness13) == 13 and strong(witness13), 'explicit H13 witness')
    print('explicit strong H13 witness:', {'rows': witness13, 'DP': h_dp(witness13), 'permutations': h_perm(witness13)})
    aut_cases = 0
    for k, max_m in [(0, 5), (2, 4), (3, 2)]:
        for code in range(1 << (k * (k - 1) // 2)):
            for m in range(max_m + 1):
                for w in product(range(1 << k), repeat=m):
                    value = h_automaton(w, k, code)
                    out = tournament(w, k, code)
                    need(value == h_dp(out), ('general automaton', k, code, w))
                    if len(out) <= 5:
                        need(value == h_perm(out), ('general permutation', k, code, w))
                    aut_cases += 1
    print('k=0,2,3 automaton/DP controls:', aut_cases)
    # Historical THM079 Part R.2 hostile, not a counterexample to its dichotomy.
    out = (502, 508, 504, 497, 160, 192, 144, 256, 112)
    triangles = [set(c) for c in combinations(range(9), 3)
                 if all(sum(bool(out[v] >> u & 1) for u in c if u != v) == 1 for v in c)]
    a, b = {0, 1, 3}, {4, 5, 6}
    r = set(range(9)) - a - b
    poison = [(u, v) for u in sorted(r) for v in sorted(r - {u})
              if all(v in c for c in triangles if u in c)]
    need(poison == [(7, 8), (8, 7)], 'poisoning two-cycle')
    need(not any(x.isdisjoint(y) and x.isdisjoint(z) and y.isdisjoint(z)
                 for x, y, z in combinations(triangles, 3)), 'triangle matching number')
    print('historical poisoning-DAG hostile:', {'rows': out, 'R': sorted(r), 'arcs': poison, 'H': h_dp(out)})
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
