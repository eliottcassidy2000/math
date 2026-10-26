"""Independent small controls for incoming S65; imports no repository code."""
import itertools
import json
import math


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


IRREGULAR = {
    ('A', 0): (1, 1, 'B'), ('A', 1): (0, 1, 'C'),
    ('B', 0): (1, -1, 'C'), ('B', 1): (1, 1, 'A'),
    ('C', 0): (1, 1, 'A'), ('C', 1): (1, -1, 'D'),
    ('D', 0): (0, -1, 'C'), ('D', 1): (0, -1, 'A'),
}


def sweep(word):
    # Dense tape with three guard cells on each side; native stop event.
    tape = [0, 0, 0] + list(word) + [0, 0, 0]
    before = tape[:]
    head, state, time = len(word) - 1, 'C', 0
    while state != 'A':
        need(-3 <= head < len(word) + 3, 'native window')
        pos = head + 3
        write, move, state = IRREGULAR[state, tape[pos]]
        tape[pos] = write
        head += move
        time += 1
        need(time <= 2*len(word) + 4, 'sweep stopping event')
    changed = [i - 3 for i, (a, b) in enumerate(zip(before, tape)) if a != b]
    return tape, head, state, time, changed


def word_graph(word):
    n = 3*len(word)
    graph = [[int(i < j) for j in range(n)] for i in range(n)]
    for cell, bit in enumerate(word):
        if bit:
            a, b = 3*cell, 3*cell+2
            graph[a][b], graph[b][a] = 0, 1
    return graph


def decode_word(graph):
    n = len(graph)
    reachable = [row[:] for row in graph]
    for i in range(n):
        reachable[i][i] = 1
    for k in range(n):
        for i in range(n):
            if reachable[i][k]:
                for j in range(n):
                    reachable[i][j] |= reachable[k][j]
    unseen, components = set(range(n)), []
    while unseen:
        a = min(unseen)
        block = {b for b in unseen if reachable[a][b] and reachable[b][a]}
        unseen -= block
        components.append(block)
    # In a tournament the SCC quotient is a transitive tournament.
    ordered = sorted(components, key=lambda block: -sum(
        graph[min(block)][min(other)] for other in components if other is not block))
    result, zeros = [], 0
    for block in ordered + [set()]:
        if len(block) == 1:
            zeros += 1
            continue
        need(zeros % 3 == 0, 'zero run remains divisible by three')
        result.extend([0] * (zeros // 3))
        zeros = 0
        if block:
            need(len(block) == 3, 'nontrivial component is C3')
            result.append(1)
    return tuple(result)


def permutation(cycle_lengths):
    out, start = [], 0
    for length in cycle_lengths:
        out.extend(start + (i+1) % length for i in range(length))
        start += length
    return tuple(out)


def signed_orbits(perm):
    remaining = set(itertools.combinations(range(len(perm)), 2))
    answer = []
    while remaining:
        first = min(remaining)
        a, b = first
        length = 0
        while True:
            remaining.remove(tuple(sorted((a, b))))
            a, b = perm[a], perm[b]
            length += 1
            if tuple(sorted((a, b))) == first:
                break
        answer.append((length, 1 if (a, b) == first else -1))
    return answer


def product_polynomial(orbits):
    poly = [1]
    for length, sign in orbits:
        factor = [math.comb(length, d) * (1 + sign * (-1)**d)
                  for d in range(length+1)]
        new = [0] * (len(poly) + length)
        for i, a in enumerate(poly):
            for j, b in enumerate(factor):
                new[i+j] += a*b
        poly = new
    return {d: c for d, c in enumerate(poly) if c}


def direct_histogram(perm):
    pairs = list(itertools.combinations(range(len(perm)), 2))
    indices = {pair: i for i, pair in enumerate(pairs)}
    comparisons = []
    for i, (a, b) in enumerate(pairs):
        x, y = perm[a], perm[b]
        comparisons.append((i, indices[tuple(sorted((x, y)))], int(x > y)))
    histogram = {}
    for mask in range(1 << len(pairs)):
        defect = sum(((mask >> i) ^ (mask >> j) ^ flip) & 1
                     for i, j, flip in comparisons)
        histogram[defect] = histogram.get(defect, 0) + 1
    return histogram


def main():
    sweeps = 0
    for k in range(51):
        for terminal in (0, 1):
            word = ([0] if terminal == 0 else [1, 1]) + [0, 1]*k
            tape, head, state, time, changed = sweep(word)
            expected = [0, 0, 0] + word + [0, 0, 0]
            expected[3] ^= 1
            need(tape == expected and changed == [0], 'one changed tape cell')
            need(state == 'A' and head == (1 if terminal == 0 else -1), 'native exit')
            need(time == 2*k + 1 + terminal, 'native clock')
            # Same fixed window before/after: exactly one internal arc changes.
            new_word = word[:]
            new_word[0] ^= 1
            # Pad both ends so both initial and final head lie in the window.
            g = word_graph([0] + word + [0])
            h = word_graph([0] + new_word + [0])
            edits = sum(g[i][j] != h[i][j] for i in range(len(g))
                        for j in range(i+1, len(g)))
            need(edits == 1, 'one literal arc edit')
            # Cyclic-triangle counts differ by one, so graphs are nonisomorphic.
            need(abs(sum(word)-sum(new_word)) == 1, 'triangle obstruction')
            sweeps += 1
    words = 0
    for n in range(6):
        for word in itertools.product((0, 1), repeat=n):
            graph = word_graph(word)
            relabelled = [row[::-1] for row in graph[::-1]]
            need(decode_word(relabelled) == word, 'unmarked SCC decoder')
            words += 1
    left, right = word_graph((0, 1)), word_graph((1, 0))
    paths = []
    for graph in (left, right):
        paths.append(sum(all(graph[a][b] for a, b in zip(p, p[1:]))
                         for p in itertools.permutations(range(6))))
    distance = min(sum(left[i][j] != right[p[i]][p[j]]
                       for i in range(6) for j in range(i+1, 6))
                   for p in itertools.permutations(range(6)))
    need(paths == [3, 3] and distance == 2, 'same H but different word')
    # Native first returns to A: only these literal rules are encountered.
    bounce = {('A', 0): (1, 1, 'B'), ('A', 1): (0, -1, 'A'),
              ('B', 1): (1, -1, 'A')}
    scalar_outputs = []
    for initial in ((0, 1), (1, 0)):
        tape, head, state, time = list(initial), 0, 'A', 0
        while time == 0 or state != 'A':
            write, move, state = bounce[state, tape[head]]
            tape[head] = write
            head += move
            time += 1
        scalar_outputs.append([tape, head, time])
    need(scalar_outputs == [[[1, 1], 0, 2], [[0, 0], -1, 1]],
         'equal scalar inputs different next macro')
    zero = word_graph((0, 0))
    outside = (0, 1, 4, 5)
    need(all(zero[2][v] == zero[3][v] for v in outside), 'source twin pair')
    need(any(left[2][v] != left[3][v] for v in outside), 'write splits twins')
    types = [(), (1,1,1,1), (3,), (4,), (2,2), (2,3), (6,), (1,2)]
    histograms, combinations = {}, 0
    for cycle_type in types:
        perm = permutation(cycle_type)
        orbits = signed_orbits(perm)
        actual = direct_histogram(perm)
        need(actual == product_polynomial(orbits), 'full signed-pair polynomial')
        even = sum(length % 2 == 0 for length in cycle_type)
        need(sum(sign == -1 for _, sign in orbits) == even, 'negative-orbit count')
        need(min(actual) == even, 'minimum defect')
        minima = 2**len(orbits) * math.prod(length//2 for length in cycle_type
                                            if length % 2 == 0)
        need(actual[even] == minima, 'minimizer count')
        need(all(d % 2 == even % 2 for d in actual), 'defect parity')
        combinations += sum(actual.values())
        histograms[str(cycle_type)] = actual
    need(histograms[str((1,2))] == {1: 4, 3: 4}, 'conditional p37 pair action')
    need(histograms[str((2,2))] == {2: 16, 4: 32, 6: 16}, 'even permutation hostile')
    def valuation(n, p):
        need(n > 0, 'positive valuation input')
        count = 0
        while n % p == 0:
            count += 1
            n //= p
        return count
    odd, credit, returned_sum = 27, 6, 0
    initial_wrong = 6 + 2*valuation(odd,3)-1
    need(initial_wrong == 11 and credit == 6, 'incoming source-index hostile')
    first_returns = []
    for index in range(257):
        need(credit == 6+returned_sum, 'sum starts at first returned odd value')
        need(credit != initial_wrong+returned_sum, 'including initial source adds five')
        if index < 5:
            first_returns.append([index, odd, credit])
        if index == 256:
            break
        k = valuation(3*odd-1,2)
        odd = 3**(k-1)*(3*odd-1)//2**k
        need(valuation(odd,3) == k-1, 'output valuation shifts the index')
        credit += 2*k-3
        returned_sum += 2*valuation(odd,3)-1
    print(json.dumps({'status': 'PASS', 'irregular_sweeps': sweeps,
                      'counter_index_controls': 257,
                      'first_odd_index_value_credit': first_returns,
                      'SCC_word_decodes': words, 'permutation_types': len(types),
                      'scalar_hostile': {'H': paths, 'distance': distance,
                                          'next_tape_head_time': scalar_outputs},
                      'clone_split_control': '00 -> 01 splits pair {2,3}',
                      'tournament_permutation_combinations': combinations,
                      'histograms': histograms}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
