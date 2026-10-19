"""Tiny permutation witnesses for packet138, not a monodromy census."""
import json
from collections import deque

checks = []


def check(name, value):
    if not value:
        raise RuntimeError(name)
    checks.append(name)


def cycle(n, entries):
    out = list(range(n))
    for left, right in zip(entries, entries[1:]+entries[:1]):
        out[left] = right
    return tuple(out)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0]*len(p)
    for i, value in enumerate(p):
        out[value] = i
    return tuple(out)


def commutator(p, q):
    return compose(compose(compose(p, q), inverse(p)), inverse(q))


def generated(generators):
    n = len(generators[0])
    identity = tuple(range(n))
    seen, queue = {identity}, deque([identity])
    while queue:
        g = queue.popleft()
        for h in generators:
            new = compose(g, h)
            if new not in seen:
                seen.add(new)
                queue.append(new)
    return seen


for left, right in [(2, 2), (2, 5), (3, 4), (4, 4), (6, 7)]:
    n = left+right-1
    sigma = cycle(n, list(range(left)))
    tau = cycle(n, [0]+list(range(left, n)))
    expected = cycle(n, [0, sigma[0], tau[0]])
    check(f'single_overlap_commutator_{left}_{right}', commutator(sigma, tau) == expected)

left = cycle(4, [0, 2, 3])
right = cycle(4, [1, 2, 3])
check('three_cycle_attachment_product', compose(left, inverse(right)) == cycle(4, [0, 2, 1]))
check('two_point_overlap_generates_A4', len(generated([cycle(4, [0, 1, 2]), cycle(4, [0, 1, 3])])) == 12)
check('one_point_overlap_generates_A5', len(generated([cycle(5, [0, 1, 2]), cycle(5, [0, 3, 4])])) == 60)
check('attachment_to_A4_generates_A6', len(generated([cycle(6, [0, 1, 2]), cycle(6, [0, 1, 3]), cycle(6, [0, 4, 5])])) == 360)

supports = [{0, 1, 2}, {2, 3, 4, 5}, {5, 6}]
n, q = 7, len(supports)
adj = {i: set() for i in range(n+q)}
for j, support in enumerate(supports):
    for vertex in support:
        adj[vertex].add(n+j)
        adj[n+j].add(vertex)
seen, queue = {0}, deque([0])
while queue:
    for neighbor in adj[queue.popleft()]:
        if neighbor not in seen:
            seen.add(neighbor)
            queue.append(neighbor)
check('support_incidence_connected', len(seen) == n+q)
check('support_incidence_tree_count', sum(map(len, supports)) == n+q-1)
check('tree_supports_intersect_at_most_once', all(len(supports[i]&supports[j]) <= 1 for i in range(q) for j in range(i)))
check('finite_branch_defect_sum', sum(len(support)-1 for support in supports) == n-1)

print(json.dumps({'packet': 138, 'status': 'PASS', 'checks': len(checks), 'verified': checks,
                  'scope': 'Tiny exact permutation and graph witnesses only; the universal incidence-tree and group arguments are proved in the note.'},
                 indent=2, sort_keys=True))
