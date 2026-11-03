"""SAT search for equal bounded legality observations in the named path chart.

No Hamilton count is used. Requires python-sat; a SAT model is checked directly.
"""
from itertools import combinations
import json
import sys
from pysat.solvers import Solver


def commands(n):
    result = []
    for i in range(n):
        for j in range(i + 2, n):
            for a, b in [(i, j), (j, i)]:
                p = list(range(n))
                p.insert(b, p.pop(a))
                result.append(tuple(p))
    return result


def literals(n, word, index):
    answer = set()
    for a, b in zip(word, word[1:]):
        if abs(a - b) == 1:
            if a > b:
                return None
        else:
            v = index[min(a, b), max(a, b)]
            answer.add(v if a < b else -v)
    return frozenset(answer)


def tests(n, depth):
    free = [(i, j) for i, j in combinations(range(n), 2) if j > i + 1]
    index = {e: k + 1 for k, e in enumerate(free)}
    cmds = commands(n)
    guards = [literals(n, p, index) for p in cmds]
    result = set(g for g in guards if g is not None)
    if depth >= 2:
        # C followed by a relocation.
        result |= {frozenset(-v for v in g) for g in guards if g is not None}
        # Two relocations, retaining legality of the intermediate word.
        for p, g in zip(cmds, guards):
            if g is None:
                continue
            for q in cmds:
                h = literals(n, tuple(p[v] for v in q), index)
                if h is None:
                    continue
                both = g | h
                if not any(-v in both for v in both):
                    result.add(both)
    return free, sorted(result, key=lambda t: (len(t), sorted(t)))


def search(n, depth=2):
    free, obs = tests(n, depth)
    m = len(free)
    top = 2 * m
    clauses = []

    def conjunction(ls):
        nonlocal top
        top += 1
        a = top
        clauses.extend([[-a, v] for v in ls])
        clauses.append([a] + [-v for v in ls])
        return a

    for t in obs:
        a = conjunction(sorted(t))
        b = conjunction(sorted(v + m if v > 0 else v - m for v in t))
        clauses += [[-a, b], [a, -b]]
    differences = []
    for v in range(1, m + 1):
        top += 1
        d, w = top, v + m
        clauses += [[v, w, -d], [-v, -w, -d], [v, -w, d], [-v, w, d]]
        differences.append(d)
    clauses.append(differences)
    with Solver(name='g4', bootstrap_with=clauses) as solver:
        sat = solver.solve()
        result = dict(n=n, depth=depth, tests=len(obs), variables=top,
                      clauses=len(clauses), collision=sat)
        if sat:
            model = set(solver.get_model())
            x = sum(1 << (v - 1) for v in range(1, m + 1) if v in model)
            y = sum(1 << (v - 1) for v in range(1, m + 1) if v + m in model)
            def evaluate(z, t):
                return all(bool(z & (1 << (abs(v)-1))) == (v > 0) for v in t)
            if x == y or any(evaluate(x, t) != evaluate(y, t) for t in obs):
                raise RuntimeError('invalid SAT witness')
            result.update(masks=[x, y], differing_edges=[e for k, e in enumerate(free) if (x^y) >> k & 1])
        return result


if __name__ == '__main__':
    for n in map(int, sys.argv[1:] or [6, 7, 8, 9, 10]):
        print(json.dumps(search(n)), flush=True)
