#!/usr/bin/env python3
"""Discover finite cycle witnesses. Search success is not a proof obligation.

Replay is in verify_square_bridge.py and does not import this module or OR-Tools.
The insertion/2-opt strategy follows the credited Gerbicz/repo implementation.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def square(n):
    return math.isqrt(n) ** 2 == n


def solve(n):
    from ortools.sat.python import cp_model
    model = cp_model.CpModel()
    arcs = []
    for a in range(1, n + 1):
        for r in range(math.isqrt(a) + 1, math.isqrt(a + n) + 1):
            b = r * r - a
            if b != a:
                arcs.append((a - 1, b - 1, model.NewBoolVar(f'e{a}_{b}')))
    model.AddCircuit(arcs)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 3102026
    solver.parameters.max_time_in_seconds = 30
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(('finite search incomplete', n, solver.StatusName(status)))
    nxt = {a + 1: b + 1 for a, b, x in arcs if solver.Value(x)}
    result = [1]
    for _ in range(n - 1):
        result.append(nxt[result[-1]])
    return result


def extend(old, rng):
    n = len(old) + 1
    squares = {r * r for r in range(2, math.isqrt(2 * n) + 1)}
    insertions = [i for i, a in enumerate(old) if a + n in squares]
    direct = [i for i in insertions if old[(i + 1) % (n - 1)] + n in squares]
    if direct:
        i = direct[0]
        return old[:i + 1] + [n] + old[i + 1:], [i, []]
    for _ in range(8):
        i = rng.choice(insertions)
        a = old[:i + 1] + [n] + old[i + 1:]
        pos = {v: j for j, v in enumerate(a)}
        bad = (n, a[(i + 2) % n])
        moves = []
        for _ in range(1000):
            x, y = bad
            i0 = pos[x]
            if a[(i0 + 1) % n] != y:
                x, y = y, x
                i0 = pos[x]
            candidates = set()
            for s in squares:
                if s - x in pos:
                    candidates.add(pos[s - x])
                if s - y in pos:
                    candidates.add((pos[s - y] - 1) % n)
            choices = []
            best = 1
            for j in sorted(candidates):
                if j in (i0, (i0 - 1) % n, (i0 + 1) % n):
                    continue
                c, d = a[j], a[(j + 1) % n]
                gain = int(x + c in squares) + int(y + d in squares)
                if gain > best:
                    best, choices = gain, []
                if gain == best:
                    choices.append(j)
            if not choices:
                break
            j = rng.choice(choices)
            c, d = a[j], a[(j + 1) % n]
            l, r = sorted((i0, j))
            moves.append([l, r])
            a[l + 1:r + 1] = reversed(a[l + 1:r + 1])
            for z in range(l + 1, r + 1):
                pos[a[z]] = z
            if best == 2:
                return a, [i, moves]
            bad = (x, c) if x + c not in squares else (y, d)
    return None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--last', type=int, default=6242)
    args = parser.parse_args()
    rng = random.Random(20261031)
    path = solve(32)
    data = {'first': 32, 'last': args.last, 'initial': path, 'steps': []}
    resets = []
    for n in range(33, args.last + 1):
        new, trace = extend(path, rng)
        if new is None:
            new = solve(n)
            data['steps'].append({'cycle': new})
            resets.append(n)
        else:
            data['steps'].append({'insert': trace[0], 'reverse': trace[1]})
        if sorted(new) != list(range(1, n + 1)) or any(
                not square(x + y) for x, y in zip(new, new[1:] + new[:1])):
            raise RuntimeError(('invalid discovery', n))
        path = new
        if n % 1000 == 0:
            print(json.dumps({'discovered_through': n, 'reset_sizes': resets}), flush=True)
    raw = (json.dumps(data, separators=(',', ':')) + '\n').encode()
    (HERE / 'square_bridge.json.gz').write_bytes(gzip.compress(raw, mtime=0))
    print(json.dumps({'first': 32, 'last': args.last, 'reset_sizes': resets,
                      'trace_bytes': len(raw), 'status': 'discovery complete'}))


if __name__ == '__main__':
    main()
