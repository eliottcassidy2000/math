#!/usr/bin/env python3
"""Exact independent checks for the Session 197 language-entropy packet.
No repository data or third-party packages. Explicit checks survive python -O.
"""
import hashlib
import itertools
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (i, j) in enumerate(pairs):
            if mask >> bit & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield rows


def physical(apex, word):
    k = len(apex)
    rows = list(apex) + [0] * len(word)
    for i, letter in enumerate(word):
        v = k + i
        for a in range(k):
            if letter >> a & 1:
                rows[a] |= 1 << v
            else:
                rows[v] |= 1 << a
        for j in range(i + 1, len(word)):
            rows[v] |= 1 << (k + j)
    return rows


def induced(rows, vertices):
    return [sum(1 << j for j, v in enumerate(vertices) if rows[u] >> v & 1)
            for u in vertices]


def hamilton_dp(rows):
    n = len(rows)
    if n == 0:
        return 1
    # Independent graph-side terminal-removal recurrence.
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1)
                for v in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        remaining = mask
        while remaining:
            bit = remaining & -remaining
            remaining -= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit
            if not rest:
                dp[mask][v] = 1
            else:
                predecessors = rest & incoming[v]
                total = 0
                while predecessors:
                    b = predecessors & -predecessors
                    predecessors -= b
                    total += dp[rest][b.bit_length() - 1]
                dp[mask][v] = total
    return sum(dp[-1])


def components(rows):
    n = len(rows)
    reach = [rows[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    unseen = set(range(n))
    result = []
    while unseen:
        i = min(unseen)
        part = sorted(j for j in unseen if reach[i] >> j & 1 and reach[j] >> i & 1)
        unseen.difference_update(part)
        result.append(part)
    return result


class SlotTransfer:
    """Inherited general-k nonnegative slot transfer, generated from its bijection."""
    def __init__(self, apex):
        self.apex = tuple(apex)
        self.k = len(apex)
        self.orders = list(itertools.permutations(range(self.k)))
        self.states = [p + (z,) for p in itertools.product(range(3), repeat=self.k)
                       for z in range(2)]
        self.index = {s: i for i, s in enumerate(self.states)}
        self.transitions = {}
        self.accept = []
        for order in self.orders:
            accept = []
            for state in self.states:
                good = True
                for j in range(self.k):
                    if state[j] == 2:
                        good = False
                    if j > 0 and state[j] == 0 and not (apex[order[j - 1]] >> order[j] & 1):
                        good = False
                accept.append(good)
            self.accept.append(tuple(accept))
            for letter in range(1 << self.k):
                transitions = []
                for state in self.states:
                    targets = []
                    for j in range(self.k + 1):
                        if state[j] == 0 and j > 0 and not (letter >> order[j - 1] & 1):
                            continue
                        target = list(state)
                        target[j] = 1 if j == self.k or not (letter >> order[j] & 1) else 2
                        targets.append(self.index[tuple(target)])
                    transitions.append(tuple(targets))
                self.transitions[order, letter] = tuple(transitions)

    def initial(self):
        return tuple(tuple(int(i == 0) for i in range(len(self.states))) for _ in self.orders)

    def step(self, vectors, letter, cap=None):
        result = []
        for order, vector in zip(self.orders, vectors):
            output = [0] * len(self.states)
            for i, weight in enumerate(vector):
                if weight:
                    for target in self.transitions[order, letter][i]:
                        output[target] += weight
            if cap is not None:
                output = [min(cap, x) for x in output]
            result.append(tuple(output))
        return tuple(result)

    def value(self, vectors):
        return sum(x for vector, accept in zip(vectors, self.accept)
                   for x, good in zip(vector, accept) if good)

    def evaluate(self, word):
        state = self.initial()
        for letter in word:
            state = self.step(state, letter)
        return self.value(state)


def verify_scc_projection(apex, word, rows):
    k = len(apex)
    for component in components(rows):
        roots = [v for v in component if v < k]
        base = [v - k for v in component if v >= k]
        if not roots:
            require(len(component) == 1, 'nontrivial SCC without an apex')
            continue
        if base:
            require(base == list(range(base[0], base[-1] + 1)), 'SCC base is not an interval')
            factor = word[base[0]:base[-1] + 1]
        else:
            factor = ()
        projected = tuple(sum(((c >> r) & 1) << j for j, r in enumerate(roots))
                          for c in factor)
        normalized = roots + [k + i for i in base]
        require(induced(rows, normalized) == physical(induced(apex, roots), projected),
                'SCC factor/projection is not the actual induced graph')


def fibonacci_words(max_length):
    level = [()]
    for length in range(max_length + 1):
        yield length, level
        level = [word + (c,) for word in level for c in (1, 2)
                 if not (word and word[-1] == c == 2)]


def exact_strong_order_bound(B):
    # Busch: 5^(n-1) <= H^3 <= B^3, using exact integer comparisons.
    n = 1
    while 5 ** n <= B ** 3:
        n += 1
    return n


def main():
    cross_checks = 0
    scc_checks = 0
    saturation_checks = 0
    # Complete labelled apex banks and complete specified word universes.
    for k, max_length in [(0, 8), (1, 7), (2, 4), (3, 2)]:
        for apex in tournaments(k):
            transfer = SlotTransfer(apex)
            for length in range(max_length + 1):
                for word in itertools.product(range(1 << k), repeat=length):
                    rows = physical(apex, word)
                    h = transfer.evaluate(word)
                    require(h == hamilton_dp(rows), f'slot/graph mismatch {apex} {word}')
                    verify_scc_projection(apex, word, rows)
                    scc_checks += 1
                    vector = transfer.initial()
                    capped = transfer.initial()
                    for c in word:
                        vector = transfer.step(vector, c)
                        capped = transfer.step(capped, c, 8)
                        require(capped == tuple(tuple(min(8, z) for z in v) for v in vector),
                                'nonnegative saturation mismatch')
                    require(min(8, h) == min(8, transfer.value(capped)),
                            'capped final acceptance mismatch')
                    saturation_checks += 1
                    cross_checks += 1

    apex = [2, 0]
    transfer = SlotTransfer(apex)
    B = 5001
    N = exact_strong_order_bound(B)
    max_length = N - 2
    witness = {1: ()}
    fib = [0, 1]
    for _ in range(max_length + 3):
        fib.append(fib[-1] + fib[-2])
    word_count = strong_word_count = extra_graph_checks = 0
    length_inventory = []
    for length, level in fibonacci_words(max_length):
        require(len(level) == fib[length + 2], 'Fibonacci word count')
        length_inventory.append(len(level))
        for index, word in enumerate(level):
            h = transfer.evaluate(word)
            is_strong = (2 in word)
            if length <= 8 or index in {0, len(level) // 2, len(level) - 1}:
                rows = physical(apex, word)
                require(h == hamilton_dp(rows), 'Fibonacci graph cross-check')
                require(is_strong == (len(components(rows)) == 1), 'Fibonacci strongness')
                extra_graph_checks += 1
            if is_strong:
                strong_word_count += 1
                require(h ** 3 >= 5 ** (length + 1), 'Busch bound finite control')
            else:
                require(h == 1, 'all-1 word is not transitive')
            if h <= B:
                witness.setdefault(h, word)
            word_count += 1
    atoms = sorted(set(witness) - {1})
    closure = {1}
    queue = [1]
    factorization = {1: []}
    for x in queue:
        for a in atoms:
            if a * x > B:
                break
            if a * x not in closure:
                closure.add(a * x)
                factorization[a * x] = factorization[x] + [a]
                queue.append(a * x)
    require((N, max_length, word_count, strong_word_count) == (16, 14, 2582, 2567),
            'census size changed')
    require((len(atoms), len(closure)) == (157, 414), 'census result changed')
    require(11 not in closure and 13 in closure, 'golden-mean boundary failure')
    escape = physical(apex, (1, 2, 2))
    require(escape == [6, 24, 26, 17, 1] and hamilton_dp(escape) == 11,
            'physical escape H11')
    for m in range(9):
        require(transfer.evaluate((0,) * m + (2,)) == 3, 'weak padding hostile')
    # Exact exponent comparison: phi^3=2+sqrt(5)<5, so rho<1.
    require(5 < 9, 'sqrt(5)<3 exact comparison')
    require(5 * 55 ** 2 < 127 ** 2, 'rho<9/10 exact comparison')
    for middle in [(1, 2, 2), (2, 1, 2), (1, 2, 0)]:
        for i in range(4):
            for j in range(4):
                require(transfer.evaluate((0,) * i + middle + (3,) * j)
                        == transfer.evaluate(middle), 'finite-census padding hostile')
    payload = {
        'bound': B,
        'maximum_strong_order': N,
        'maximum_chain_length': max_length,
        'all_fibonacci_words_tested': word_count,
        'strong_words_tested': strong_word_count,
        'strong_distinct_values': len(atoms),
        'ordered_product_distinct_values_including_one': len(closure),
        'witnesses': {str(h): list(witness[h]) for h in sorted(witness)},
        'product_factorizations': {str(h): factorization[h] for h in sorted(closure)},
    }
    out = Path(__file__).resolve().parent
    data = json.dumps(payload, sort_keys=True, indent=2) + '\n'
    (out / 'bounded_image.json').write_text(data)
    summary = {
        'slot_vs_graph_checks': cross_checks,
        'scc_interval_projection_checks': scc_checks,
        'saturation_checks': saturation_checks,
        'fibonacci_extra_graph_checks': extra_graph_checks,
        'word_counts_by_length': length_inventory,
        'all_order_bound': B,
        'maximum_strong_order': N,
        'maximum_chain_length': max_length,
        'strong_generator_values': len(atoms),
        'ordered_product_values_including_one': len(closure),
        'first_extra_omission': 11,
        'escape_H11_rows': escape,
        'bounded_image_sha256': hashlib.sha256(data.encode()).hexdigest(),
    }
    receipt = json.dumps(summary, sort_keys=True, indent=2) + '\n'
    (out / 'receipt.json').write_text(receipt)
    print(receipt, end='')


if __name__ == '__main__':
    main()
