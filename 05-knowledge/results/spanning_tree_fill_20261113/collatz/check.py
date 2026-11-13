"""Independent exact edge, boundary, and bounded-fill controls."""
from collections import deque
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def collatz(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def a_edge(x, y):
    return x > 0 and y > 0 and (collatz(x) == y or collatz(y) == x)


def extra_edge(x, y):
    return x > 1 and x % 2 == 1 and x % 3 == 1 and 3 * y == x - 1 and y > 0


def b_edge(x, y):
    return a_edge(x, y) or extra_edge(x, y)


def b_successors(n):
    out = {collatz(n), 2 * n}
    if n % 6 == 4:
        out.add((n - 1) // 3)
    if n > 1 and n % 6 == 1:
        out.add((n - 1) // 3)
    return out


def macro(n):
    need(n > 1, 'macro root is terminal')
    if n % 2 == 0:
        return [n, n // 2]
    if n % 3 == 1:
        return [n, (n - 1) // 3]
    if n % 3 == 2:
        return [n, 2 * n, (2 * n - 1) // 3]
    return [n, 3*n+1, 6*n+2, 12*n+4, 4*n+1, 4*n//3, 2*n//3, n//3]


def route(n):
    out = [n]
    steps = 0
    while out[-1] != 1:
        word = macro(out[-1])
        need(3 * word[-1] <= 2 * word[0], 'macro contraction')
        out.extend(word[1:])
        steps += 1
    return out, steps


def ceiling_log(n):
    k, two, three = 0, 1, 1
    while three < n * two:
        k += 1
        two *= 2
        three *= 3
    return k


def boundary(vertices):
    vertices = set(vertices)
    parent = {n: n for n in vertices}

    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    edges = 0
    for n in vertices:
        successor = collatz(n)
        if successor in vertices:
            edges += 1
            parent[find(n)] = find(successor)
    components = len({find(n) for n in vertices})
    sinks = sorted(n for n in vertices if collatz(n) not in vertices)
    cycles, seen = [], set()
    for n in sorted(vertices):
        path, position = [], {}
        current = n
        while current in vertices and current not in seen and current not in position:
            position[current] = len(path)
            path.append(current)
            current = collatz(current)
        if current in position:
            cycles.append(path[position[current]:])
        seen.update(path)
    need(components == len(sinks) + len(cycles), 'component-terminal identity')
    need(edges == len(vertices) - len(sinks), 'internal arc count')
    need(edges - len(vertices) + components == len(cycles), 'cycle rank')
    return dict(vertices=len(vertices), components=components, sinks=sinks,
                cycles=cycles, edges=edges)


def reverse_reachable(cap):
    incoming = [[] for _ in range(cap + 1)]
    for n in range(1, cap + 1):
        for target in b_successors(n):
            if target <= cap:
                incoming[target].append(n)
    reached = {1}
    queue = deque([1])
    while queue:
        target = queue.popleft()
        for source in incoming[target]:
            if source not in reached:
                reached.add(source)
                queue.append(source)
    return reached


def main():
    cap = 20000
    maximum_length = (0, 1)
    primitive_edges = 0
    for n in range(1, cap + 1):
        word, macro_steps = route(n)
        need(word[-1] == 1 and all(b_edge(x, y) for x, y in zip(word, word[1:])),
             'literal directed B route')
        need(max(word) <= 12*n+4, 'linear height buffer')
        need(macro_steps <= ceiling_log(n), 'logarithmic macro count')
        need(len(word)-1 <= 7*ceiling_log(n), 'logarithmic native length')
        primitive_edges += len(word)-1
        maximum_length = max(maximum_length, (len(word)-1, n))
        if n > 1:
            need(macro(n)[-1] < n, 'one earlier parent per integer')
    need(extra_edge(7, 2) and not a_edge(7, 2), 'extra edge is not genuine reverse')
    need(not extra_edge(1, 0), 'positive domain boundary')
    need([macro(n)[-1] for n in (5, 6, 9)] == [3, 3, 3], 'macro tree is not binary')
    for mask in range(1 << 10):
        boundary(n+1 for n in range(10) if mask >> n & 1)
    prefixes = []
    for n in range(1, 501):
        result = boundary(range(1, n+1))
        need(len(result['sinks']) == (n+1)//2 - (n+2)//6, 'closed boundary count')
        if n <= 12:
            prefixes.append(result)
    fill = []
    for n in (1, 3, 10, 30, 100, 1000):
        budget = 12*n+4
        reached = reverse_reachable(budget)
        need(all(i in reached for i in range(1, n+1)), 'independent finite-buffer fill')
        fill.append(dict(prefix=n, buffer=budget, all_prefix_reaches_root=True))
    need(reverse_reachable(3) == {1, 2}, 'unbuffered native prefix hostile')
    examples = {str(n): route(n)[0] for n in (1, 2, 3, 5, 7, 27, 703)}
    receipt = dict(status='PASS', route_universe=[1, cap], primitive_edges=primitive_edges,
                   largest_observed_length=maximum_length,
                   boundary_universe='all1024 subsets of[1,10], plus prefixes through500',
                   prefixes=prefixes, bounded_fill=fill, routes=examples,
                   scope='Model B reachability; model A remains classical Collatz')
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    path = HERE/'receipt.json'
    if '--write' in sys.argv:
        path.write_text(encoded)
    else:
        need(path.read_text() == encoded, 'saved receipt mismatch')
    print('PASS:20000 routes,1024 subsets,500 prefixes,6 independent buffer graphs')


if __name__ == '__main__':
    main()
