#!/usr/bin/env python3
"""Exact modular homogeneous-context algebra and branching resource audit.
Standard library only; all assertions are explicit and survive python -O.
No imports from a repository checkout. Graphs are adjacency bitmask tuples.
"""
from functools import lru_cache
from itertools import product, permutations
from math import factorial, gcd, comb
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def divisor(c):
    return factorial(c) * factorial(c - 1)


def cutoff(m):
    c = 1
    while divisor(c) % m:
        c += 1
    return c


def graph(n, mask):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if mask >> k & 1 else (j, i)
            out[u] |= 1 << v
            k += 1
    return tuple(out)


def linear(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def substitute(q, blocks):
    offsets = [0]
    for b in blocks:
        offsets.append(offsets[-1] + len(b))
    out = [0] * offsets[-1]
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            out[offsets[i] + u] = row << offsets[i]
            for j, b in enumerate(blocks):
                if q[i] >> j & 1:
                    out[offsets[i] + u] |= ((1 << len(b)) - 1) << offsets[j]
    return tuple(out)


ONE = (0,)
JOIN = (2, 0)
CYCLE = (2, 4, 1)


def all_h(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    out = [0] * (1 << n)
    out[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        for v, value in enumerate(dp[s]):
            if value:
                out[s] += value
                rest = a[v] & ~s
                while rest:
                    bit = rest & -rest
                    rest -= bit
                    dp[s | bit][bit.bit_length() - 1] += value
    return out


def hp(a):
    return all_h(a)[-1]


@lru_cache(None)
def profile(a):
    """Works also for digraphs with missing arcs, needed by the universal kernel."""
    n = len(a)
    h = all_h(a)
    @lru_cache(None)
    def go(s):
        if not s:
            return (1,) + (0,) * n
        bit = s & -s
        tail = s ^ bit
        ans = [0] * (n + 1)
        u = tail
        while True:
            block = u | bit
            old = go(s ^ block)
            for c in range(1, n + 1):
                ans[c] += h[block] * old[c - 1]
            if not u:
                break
            u = (u - 1) & tail
        return tuple(ans)
    return go((1 << n) - 1)


def literal_profile(a):
    n = len(a)
    out = [0] * (n + 1)
    for word in permutations(range(n)):
        bad = sum(not (a[word[i]] >> word[i + 1] & 1) for i in range(n - 1))
        for c in range(bad + 1, n + 1):
            out[c] += comb(n - 1 - bad, c - 1 - bad)
    for c in range(1, n + 1):
        require(out[c] % factorial(c) == 0, 'ordered cover division')
        out[c] //= factorial(c)
    return tuple(out)


def coefficients(q):
    """a_c with distinguished vertex 0, by labelled exterior-mask DP."""
    n = len(q)
    maxc = (n - 1) // 2 + 1
    full = (1 << n) - 2
    @lru_cache(None)
    def go(used, last, c):
        ans = [0] * (maxc + 1)
        if used == full and c:
            ans[c] += 1
        if c < maxc and last != 0 and (last < 0 or q[last] & 1):
            sub = go(used, 0, c + 1)
            ans = [a + b for a, b in zip(ans, sub)]
        for v in range(1, n):
            if not (used >> v & 1) and (last < 0 or q[last] >> v & 1):
                sub = go(used | 1 << v, v, c)
                ans = [a + b for a, b in zip(ans, sub)]
        return tuple(ans)
    return tuple(factorial(c) * value for c, value in enumerate(go(0, -1, 0)))


def detector(s):
    """x=0, a_i=1+i, b_i=1+s+i; triangular A->B incidence."""
    n = 2 * s + 1
    a = [0] * n
    def edge(i, j):
        a[i] |= 1 << j
    for i in range(s):
        edge(0, 1 + i)
        edge(1 + s + i, 0)
        for j in range(i + 1, s):
            edge(1 + i, 1 + j)
            edge(1 + s + i, 1 + s + j)
        for j in range(s):
            u, v = 1 + i, 1 + s + j
            edge(u, v) if i <= j else edge(v, u)
    return tuple(a)


def signature(p, m):
    return tuple((p[c] if c < len(p) else 0) % (m // gcd(m, divisor(c)))
                 for c in range(1, cutoff(m)))


@lru_cache(None)
def kernel(q, counts):
    """K_Q(counts,d): covers of labelled independent copies of each quotient vertex."""
    return profile(substitute(q, [(0,) * c for c in counts]))


@lru_cache(None)
def terms(q, m):
    k = cutoff(m)
    rows = []
    for counts in product(range(1, k), repeat=len(q)):
        p = kernel(q, counts)
        row = tuple(divisor(d) * (p[d] if d < len(p) else 0) for d in range(1, k))
        for value in row:
            for c in counts:
                require(value % divisor(c) == 0, 'universal coefficient divisibility')
        rows.append((counts, tuple(v % m for v in row)))
    return tuple(rows)


def transition(q, states, m):
    """Canonical p_c residues -> canonical p_c residues; no graph representative needed."""
    k = cutoff(m)
    weighted = [0] * (k - 1)
    for counts, row in terms(q, m):
        mult = 1
        for i, c in enumerate(counts):
            mult *= states[i][c - 1]
        if mult:
            for d, value in enumerate(row):
                weighted[d] = (weighted[d] + mult * value) % m
    out = []
    for d, value in enumerate(weighted, 1):
        g = gcd(m, divisor(d))
        require(value % g == 0, 'weighted output is admissible')
        reduced = m // g
        out.append((value // g * pow(divisor(d) // g, -1, reduced)) % reduced)
    return tuple(out)


def unary_transition(state, m, op):
    k = cutoff(m)
    p = [0] + list(state) + [0, 0]
    out = [0] * k
    for d in range(1, k):
        if op == 'P':
            out[d] = d * p[d] + p[d - 1]
        else:
            out[d] = ((p[d - 2] if d >= 2 else 0) + (2 * d - 1) * p[d - 1]
                      + d * (d + 2) * p[d] + d * (d + 1) * p[d + 1])
    return signature(out, m)


def unary_orbit(m):
    one = signature((0, 1), m)
    witnesses = {one: ('', ONE)}
    pending = [one]
    for state in pending:
        word, a = witnesses[state]
        for op in 'PC':
            nxt = unary_transition(state, m, op)
            if nxt not in witnesses:
                b = substitute(JOIN, [ONE, a]) if op == 'P' else substitute(CYCLE, [a, ONE, ONE])
                witnesses[nxt] = (word + op, b)
                pending.append(nxt)
    return witnesses


def grammar12(seed_branch=False):
    states = {(1, 0): ('1', ONE)}
    branch = substitute(CYCLE, [linear(2), linear(2), ONE])
    if seed_branch:
        states[(11, 3)] = ('B', branch)
    while True:
        new = {}
        for (h, p), (word, a) in states.items():
            nxt = ((3 * h + 2 * p) % 12, (3 * h + 8 * p) % 6)
            if nxt not in states:
                new.setdefault(nxt, ('C(' + word + ')', substitute(CYCLE, [a, ONE, ONE])))
        for x, y in product(states, repeat=2):
            h, p = x
            g, q = y
            nxt = (h * g % 12, (h * g + 2 * h * q + 2 * p * g + 2 * p * q) % 6)
            if nxt not in states:
                wx, a = states[x]
                wy, b = states[y]
                new.setdefault(nxt, ('J(' + wx + ',' + wy + ')', substitute(JOIN, [a, b])))
        if not new:
            return states
        states.update(new)


def content_words(q, counts):
    @lru_cache(None)
    def end(c, last):
        if not c[last]:
            return 0
        if sum(c) == 1:
            return 1
        old = list(c)
        old[last] -= 1
        old = tuple(old)
        return sum(end(old, v) for v in range(len(q)) if q[v] >> last & 1)
    return sum(end(counts, v) for v in range(len(q)))


def inactive_repeat_max(q, counts):
    repeated = sum(1 << i for i, c in enumerate(counts) if c >= 2)
    reach = [row & repeated for row in q]
    for k in range(len(q)):
        if repeated >> k & 1:
            for i in range(len(q)):
                if reach[i] >> k & 1:
                    reach[i] |= reach[k]
    return max((c for i, c in enumerate(counts) if c >= 2 and not (reach[i] >> i & 1)),
               default=1)


def repeat_support_audit():
    checks = 0
    for n in range(1, 5):
        for mask in range(1 << (n * (n - 1) // 2)):
            q = graph(n, mask)
            for counts in product(range(1, 4), repeat=n):
                w = content_words(q, counts)
                top = inactive_repeat_max(q, counts)
                require(w % factorial(top - 1) == 0, 'inactive repeat-support divisor')
                checks += 1
    sharp = []
    for c in range(2, 5):
        q = substitute(CYCLE, [ONE, ONE, linear(c - 1)])
        counts = (c, c) + (1,) * (c - 1)
        w = content_words(q, counts)
        require(w == factorial(c - 1), 'sharp repeat-support family')
        sharp.append({'c': c, 'W': w, 'counts': counts, 'quotient': q})
    hostile = content_words(CYCLE, (3, 3, 3))
    require(hostile == 3 and hostile % factorial(2) != 0, 'cyclic repeat-support hostile')
    return {'labelled_checks': checks, 'sharp_balanced_family': sharp,
            'cyclic_hostile_counts': [3, 3, 3], 'cyclic_hostile_W': hostile}


def main():
    # Complete labelled context universe through order five, and independent profile universe.
    bank = {}
    coefficient_gcd = [0] * 4
    context_count = profile_count = 0
    for n in range(1, 6):
        for mask in range(1 << (n * (n - 1) // 2)):
            a = graph(n, mask)
            p = profile(a)
            require(p == literal_profile(a), 'subset / literal permutation covers')
            bank.setdefault(p, a)
            profile_count += 1
            co = coefficients(a)
            require(co[1] == p[1], 'first context coefficient')
            for c in range(1, len(co)):
                require(co[c] % divisor(c) == 0, 'return-excursion divisor')
                coefficient_gcd[c] = gcd(coefficient_gcd[c], co[c])
            context_count += 1
    detector_rows = []
    for s in range(6):
        co = coefficients(detector(s))
        require(co[-1] == divisor(s + 1), 'sharp detector coefficient')
        require(all(co[c] % divisor(c) == 0 for c in range(1, len(co))), 'decoder integrality')
        detector_rows.append(co[1:])
    # Universal representative-free transitions against actual tournament substitution.
    small = [ONE, linear(2), linear(3), CYCLE]
    substitution_checks = 0
    for m in (3, 4, 8, 12):
        for q in (JOIN, CYCLE):
            for blocks in product(small, repeat=len(q)):
                pred = transition(q, tuple(signature(profile(a), m) for a in blocks), m)
                actual = signature(profile(substitute(q, blocks)), m)
                require(pred == actual, 'universal modular transition')
                substitution_checks += 1
    # Exact actual-context class counts: upper bounds are mathematical, positive coverage literal.
    exact_classes = {}
    for m in (3, 4, 8):
        found = unary_orbit(m)
        require(len(found) == {3: 5, 4: 3, 8: 14}[m], 'low-modulus closure count')
        witnesses = []
        for state, (word, a) in sorted(found.items()):
            p = profile(a)
            require(signature(p, m) == state, 'direct unary witness profile')
            require(p == literal_profile(a), 'literal unary witness covers')
            witnesses.append({'state': state, 'word': word, 'adjacency': a, 'H': p[1]})
        exact_classes[str(m)] = witnesses
    restricted = grammar12(False)
    complete12 = grammar12(True)
    require(len(restricted) == 7 and len(complete12) == 11, 'mod12 resource hierarchy')
    grammar_witnesses = []
    for state, (expr, a) in sorted(complete12.items()):
        require(signature(profile(a), 12) == state, 'direct mod12 representative')
        grammar_witnesses.append({'state': state, 'expression': expr, 'adjacency': a, 'H': hp(a)})
    # Every state seen in the exhaustive small universe lies in the claimed universal classes.
    for p in bank:
        for m in (3, 4, 8):
            require(signature(p, m) in unary_orbit(m), 'small universal class coverage')
        require(signature(p, 12) in complete12, 'small universal mod12 coverage')
    # Hostile: old cyclic detector loses p3 modulo eight, but the sharp detector sees it.
    old = substitute(CYCLE, [ONE, linear(2), linear(2)])
    sharp = detector(2)
    hostile = []
    for n in (3, 4):
        t = linear(n)
        vals = [hp(substitute(q, [t] + [ONE] * (len(q) - 1))) for q in (old, sharp)]
        hostile.append({'n': n, 'signature8': signature(profile(t), 8),
                        'cyclic_response': vals[0], 'sharp_response': vals[1]})
    require(hostile[0]['cyclic_response'] % 8 == hostile[1]['cyclic_response'] % 8,
            'old decoder hostile')
    require(hostile[0]['sharp_response'] % 8 != hostile[1]['sharp_response'] % 8,
            'sharp detector repairs hostile')
    # Mutation controls: erasing the highest modular coordinate must falsely merge this pair.
    require(hostile[0]['signature8'][:2] == hostile[1]['signature8'][:2], 'highest-coordinate mutant')
    # Grammar branch witnesses, with actual vertex-level Hamilton counts.
    b = substitute(CYCLE, [linear(2), linear(2), ONE])
    five = substitute(CYCLE, [linear(2), ONE, ONE])
    fiftyfive = substitute(JOIN, [five, b])
    require(hp(b) == 11 and hp(fiftyfive) == 55, 'branching escape values')
    output = {'status': 'PASS', 'repeat_support': repeat_support_audit(), 'labelled_contexts': context_count,
              'literal_profile_checks': profile_count, 'coefficient_gcds': coefficient_gcd[1:],
              'sharp_detector_rows': detector_rows, 'direct_modular_substitution_checks': substitution_checks,
              'exact_context_classes': exact_classes, 'exact_mod12_classes': grammar_witnesses,
              'restricted_G_mod12_states': sorted(restricted),
              'restricted_G_mod12_residues': sorted({h for h, p in restricted}),
              'branch_escapes': [{'H': 11, 'adjacency': b}, {'H': 55, 'adjacency': fiftyfive}],
              'old_cyclic_decoder_hostile_mod8': hostile,
              'modular_signature_precisions': {str(m): [m // gcd(m, divisor(c))
                  for c in range(1, cutoff(m))] for m in (2, 3, 4, 5, 8, 9, 12, 16, 32, 64)},
              'scope': 'All-order proofs in PROOF.md; finite tests are the stated universes, not a spectrum census.'}
    print(json.dumps(output, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
