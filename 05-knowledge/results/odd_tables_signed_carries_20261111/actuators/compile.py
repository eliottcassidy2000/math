#!/usr/bin/env python3
"""Complete finite return graph, with the unbounded translation retained."""
from collections import defaultdict
from functools import lru_cache
import json
from pathlib import Path


@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    q, a = divmod(n, 2)
    return g(q) if a else g(q) + q + q % 2


def step(state, a, multiplier=3):
    j, u, v = state
    J, b = divmod(multiplier * a + j, 2)
    nxt = (J, 2 * u + 13 * (a == 0) - 3 * multiplier * (b == 0),
           13 * (a == 0) - 3 * (1 - 2 * (J % 2)) * (b == 0))
    drift = (u + v) * a - 3 * (J + J % 2) * (b == 0)
    return nxt, drift


def graph(start, multiplier=3):
    seen, todo, edges = {start}, [start], []
    for state in todo:
        for a in (0, 1):
            nxt, drift = step(state, a, multiplier)
            if -13 <= nxt[1] <= 3 * multiplier:
                edges.append((state, nxt, a, drift))
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
    adj = defaultdict(list)
    for s, t, a, d in edges:
        adj[s].append((t, a, d))
    cycles = []
    # Every elementary cycle appears once, at its lexicographically first node.
    for root in sorted(seen):
        stack = [(root, [root], [], 0)]
        while stack:
            state, path, bits, drift = stack.pop()
            for nxt, a, d in adj[state]:
                if nxt == root:
                    cycles.append(dict(state=root,
                                       low_to_high="".join(map(str, bits + [a])),
                                       drift=drift + d))
                elif nxt > root and nxt not in path:
                    stack.append((nxt, path + [nxt], bits + [a], drift + d))
    return dict(start=start, states=sorted(seen), edges=edges, cycles=cycles)


def value(s, n, w=0, multiplier=3):
    j, u, v = s
    return 13 * g(n) - 3 * g(multiplier * n + j) + u*n + v*(n % 2) + w


def main():
    graphs = [graph((j, -4, 6)) for j in range(3)]
    checks = 0
    for chart in graphs:
        for s in chart['states']:
            for a in (0, 1):
                t, d = step(s, a)
                for n in range(128):
                    assert value(s, 2*n+a) == value(t, n, d)
                    checks += 1
    assert [len(c['states']) for c in graphs] == [20, 2, 6]
    assert graphs[0]['cycles'] == [
        dict(state=(0, -4, 10), low_to_high='0', drift=0),
        dict(state=(0, 1, 13), low_to_high='1100', drift=4)]
    assert not graphs[1]['cycles'] and not graphs[2]['cycles']
    # Independence controls justify identifying return of functions with
    # return of coefficient states, rather than merely equating symbols.
    samples = {}
    for j in range(3):
        for k in range(j+1, 3):
            # On even n, an affine n/parity/constant function has zero second
            # difference on every equally spaced triple.
            values = [g(3*n+j)-g(3*n+k) for n in (0, 2, 4, 6)]
            second = [values[i]-2*values[i+1]+values[i+2] for i in (0,1)]
            assert any(second)
            samples[f'{j},{k}'] = dict(values=values, second_differences=second)
    phase_checks = 0
    for m in range(1, 102, 2):
        if m % 5 == 0:
            continue
        for i in range(1, 5):
            ip = m*i % 5
            j = m*i // 5
            ai, ap = pow(i, -1, 5), pow(ip, -1, 5)
            for a, b in ((13,3), (25,5), (3,13), (-2,8)):
                assert (b*m*ap-a*ai) % 5 == 0
                u = (b*m*ap-a*ai)//5
                v = a*(i<=2)-b*(ip<=2)*(-1)**j
                def f(n):
                    return a*g(n)-b*g(m*n+j)+u*n+v*(n%2)
                for n in range(64):
                    assert f(16*n+3*i)-f(n) == 2*(a-b)//5
                    phase_checks += 1
    result = dict(status='FINITE COMPLETE RETURN GRAPH; algebraic proof in PROOF.md',
                  transition_controls=checks, graphs=graphs,
                  universal_phase_controls=phase_checks,
                  phase_independence=samples,
                  scope='All additive return loops of the specified ternary defect; not all target witnesses')
    path = Path(__file__).with_name('receipt.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(states=[len(c['states']) for c in graphs],
                          cycles=[c['cycles'] for c in graphs],
                          transition_controls=checks, universal_phase_controls=phase_checks)))


if __name__ == '__main__':
    main()
