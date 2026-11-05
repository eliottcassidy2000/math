#!/usr/bin/env python3
"""Exact independent checks for the full opposed-binary sparsity packet.
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


def split_runs(word):
    m = len(word)
    if not m or all(c == word[0] for c in word):
        return None
    a = 1
    while word[a] == word[0]:
        a += 1
    b = 1
    while word[m-b-1] == word[-1]:
        b += 1
    return a, b, m-a-b


def literal_binary_count(word):
    """Separate direct three-slot implementation, without generated matrices."""
    total = 0
    for reverse in (False, True):
        values = {(0, 0, 0): 1}
        for c in word:
            P = int(c == (2 if reverse else 1))
            Q = 1-P
            new = {}
            for (x, y, z), count in values.items():
                choices = [(2 if P else 1, y, z)]
                if y or P:
                    choices.append((x, 2 if Q else 1, z))
                if z or Q:
                    choices.append((x, y, 1))
                for target in choices:
                    new[target] = new.get(target, 0) + count
            values = new
        total += sum(count for (x, y, z), count in values.items()
                     if x != 2 and (y == 1 or (y == 0 and not reverse)))
    return total


def bulk_paths(word):
    a, b, c = split_runs(word)
    require(c >= 2, 'bulk path requires disjoint anchors')
    m = len(word)
    prefix = list(range(a+1))
    suffix = list(range(m-b-1, m))
    middle = list(range(a+1, m-b-1))
    p1 = next(i for i in prefix if word[i] == 1)
    p2 = next(i for i in prefix if word[i] == 2)
    q1 = next(i for i in suffix if word[i] == 1)
    q2 = next(i for i in suffix if word[i] == 2)
    for allocation in itertools.product(range(3), repeat=len(middle)):
        slots = [[i for i in prefix if i not in (p1, p2)] + [q2],
                 [p1, q1],
                 [p2] + [i for i in suffix if i not in (q1, q2)]]
        for i, slot in zip(middle, allocation):
            slots[slot].append(i)
        for slot in slots:
            slot.sort()
        yield tuple([i+2 for i in slots[0]] + [0] + [i+2 for i in slots[1]]
                    + [1] + [i+2 for i in slots[2]])


def boundary_triangles(word):
    a, b, c = split_runs(word)
    m = len(word)
    prefix_chain = list(range(2, a+2))
    suffix_chain = list(range(m-b+2, m+2))
    prefix = ([0, 1] if word[0] == 2 else [0, a+2]) + prefix_chain
    suffix = ([0, 1] if word[-1] == 2 else [1, m-b+1]) + suffix_chain
    return prefix, suffix


def bounds(B):
    L = 0
    while 2 ** (L+1)+1 <= B:
        L += 1
    R = 0
    while 3 ** (R+1) <= B:
        R += 1
    return L, R


def candidates(B):
    L, R = bounds(B)
    for a in range(1, L+1):
        yield (2,) * a
    for a in range(1, L+1):
        for b in range(1, L+1):
            for sigma in (1, 2):
                yield (sigma,) * a + (3-sigma,) * b
                yield (sigma,) * a + (3-sigma,) + (sigma,) * b
                for tau in (1, 2):
                    for r in range(R+1):
                        for middle in itertools.product((1, 2), repeat=r):
                            yield ((sigma,) * a + (3-sigma,) + middle
                                   + (3-tau,) + (tau,) * b)


def exact_image(B, transfer):
    witness = {1: ()}
    seen_words = set()
    for word in candidates(B):
        require(word not in seen_words, 'boundary presentation duplicated')
        seen_words.add(word)
        h = transfer.evaluate(word)
        require(h == literal_binary_count(word), 'independent full candidate transfer')
        if h <= B:
            witness.setdefault(h, word)
    L, R = bounds(B)
    require(len(seen_words) == L + 8 * L * L * 2 ** R, 'candidate inventory')
    atoms = sorted(set(witness) - {1})
    closure = {1}
    queue = [1]
    factorization = {1: []}
    for x in queue:
        for a in atoms:
            if a*x > B:
                break
            if a*x not in closure:
                closure.add(a*x)
                queue.append(a*x)
                factorization[a*x] = factorization[x] + [a]
    return witness, factorization, len(seen_words)


def main():
    transfer = SlotTransfer([2, 0])
    binary_word_checks = bulk_word_checks = bulk_path_checks = boundary_checks = 0
    for m in range(9):
        for word in itertools.product((1, 2), repeat=m):
            rows = physical([2, 0], word)
            h = hamilton_dp(rows)
            require(h == transfer.evaluate(word) == literal_binary_count(word), 'small graph count')
            require((2 in word) == (len(components(rows)) == 1), 'binary strongness')
            binary_word_checks += 1
            split = split_runs(word)
            if split is None:
                require(h == (2 ** m+1 if 2 in word else 1), 'constant word value')
                continue
            a, b, c = split
            left, right = boundary_triangles(word)
            require(len(set(left)) == a+2 and len(set(right)) == b+2, 'boundary triangle distinctness')
            require(hamilton_dp(induced(rows, left)) == 2 ** a+1, 'prefix triangle clone')
            require(hamilton_dp(induced(rows, right)) == 2 ** b+1, 'suffix triangle clone')
            boundary_checks += 2
            if c >= 2:
                paths = set()
                for path in bulk_paths(word):
                    require(len(path) == len(rows) and len(set(path)) == len(rows), 'spanning path')
                    require(all(rows[u] >> v & 1 for u, v in zip(path, path[1:])), 'actual path edge')
                    require(path not in paths, 'bulk assignment collision')
                    paths.add(path)
                require(len(paths) == 3 ** (c-2), 'bulk path injection size')
                require(h >= len(paths), 'bulk lower bound')
                bulk_word_checks += 1
                bulk_path_checks += len(paths)

    witnesses, factors, candidate_count = exact_image(501, transfer)
    require(len(witnesses)-1 == 96 and len(factors) == 139, 'bounded image changed')
    positive_graph_checks = 0
    for h, word in witnesses.items():
        require(hamilton_dp(physical([2, 0], word)) == h, 'positive physical witness')
        positive_graph_checks += 1
    small, small_factors, small_candidates = exact_image(19, transfer)
    require(set(small) == {1, 3, 5, 9, 11, 13, 15, 17}, 'all-order H19 omission')
    require(19 not in factors and 67 not in factors and 173 not in factors, 'specified prime omissions')
    escape_word = (2, 0, 2, 0)
    escape_rows = physical([2, 0], escape_word)
    require(hamilton_dp(escape_rows) == transfer.evaluate(escape_word) == 19, 'aligned-letter H19 escape')
    require(2 ** 3 < 3 ** 2, 'log_3(2)<2/3 exact inequality')
    payload = {
        'bound': 501,
        'candidate_word_count': candidate_count,
        'generator_witnesses_including_one': {str(h): list(witnesses[h]) for h in sorted(witnesses)},
        'product_factorizations': {str(h): factors[h] for h in sorted(factors)},
    }
    root = Path(__file__).resolve().parent
    output = json.dumps(payload, sort_keys=True, indent=2)+'\n'
    (root/'bounded_image.json').write_text(output)
    receipt = {
        'small_binary_graph_checks': binary_word_checks,
        'outer_induced_triangle_checks': boundary_checks,
        'bulk_injection_word_checks': bulk_word_checks,
        'literal_distinct_bulk_paths_checked': bulk_path_checks,
        'ceiling501_candidates_checked_two_ways': candidate_count,
        'ceiling501_strong_generator_values': len(witnesses)-1,
        'ceiling501_product_values': len(factors),
        'positive_witness_graph_checks': positive_graph_checks,
        'ceiling19_candidates': small_candidates,
        'ceiling19_scalar_image': sorted(small),
        'H19_escape_rows': escape_rows,
        'bounded_image_sha256': hashlib.sha256(output.encode()).hexdigest(),
    }
    text = json.dumps(receipt, sort_keys=True, indent=2)+'\n'
    (root/'receipt.json').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
