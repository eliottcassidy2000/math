"""Independent phase-first controls for the coupled-move weight quotient.

Enumerate cyclic phase assignments directly, then filter native binary
differences; do not generate assignments from the author's binary words.
The all-size argument is audited separately in AUDIT_COUPLED_MOVES.md.
"""
from itertools import combinations, product


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    checks = 0
    cases = [(3, 6, {0: 0, 3: 0, 6: 0}),
             (3, 8, {0: 0, 4: 0, 8: 0}),
             (3, 8, {0: 0, 4: 1, 8: 2}),
             (5, 6, {0: 0, 6: 0})]
    reports = []
    for q, r, pins in cases:
        groups = {}
        indices = sorted(pins)
        free = [j for j in range(r+1) if j not in pins]
        for values in product(range(q), repeat=len(free)):
            phases = [pins.get(j) for j in range(r+1)]
            for j, value in zip(free, values):
                phases[j] = value
            differences = [(phases[j+1]-phases[j]) % q for j in range(r)]
            if any(value not in (0, 1) for value in differences):
                continue
            weights = tuple(sum(differences[lo:hi])
                            for lo, hi in zip(indices, indices[1:]))
            groups.setdefault(weights, []).append(tuple(phases))
        for left, right in combinations(groups, 2):
            actual = min(distance(a, b) for a in groups[left] for b in groups[right])
            predicted = sum(abs(x-y) for x, y in zip(left, right))*(q-1)//q
            if actual != predicted:
                raise RuntimeError(('layer distance', q, r, left, right, actual, predicted))
            checks += 1
        for weights, states in groups.items():
            seen = {states[0]}
            todo = list(seen)
            while todo:
                a = todo.pop()
                for b in states:
                    if b not in seen and distance(a, b) == 1:
                        seen.add(b)
                        todo.append(b)
            if len(seen) != len(states):
                raise RuntimeError(('within-layer connectivity', q, r, weights))
            checks += 1
        reports.append((q, r, len(groups), sum(map(len, groups.values()))))
    print({'status': 'PASS', 'checks': checks, 'q_r_layers_states': reports})


if __name__ == '__main__':
    main()
