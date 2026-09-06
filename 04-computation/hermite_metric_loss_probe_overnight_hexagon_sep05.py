#!/usr/bin/env python3
"""Exploratory exact all-dilation two-jet precision versus unlabelled metric.

No all-node metric theorem is assumed.  Universe: all translated-to-zero
integer node sets with 4..max_nodes points and diameter <=height, for every
requested prime.  H1 from THM-4435 yields L(e)=max(2S+2(n-1)e,
B+(2n-1)e), where S=max_i v_p(F'(x_i)), B=max_i(3v_p(F')-v_p(F'')).
Comparing B exposes unit differences hidden at shallow dilation depths.
"""

import argparse
import hashlib
import itertools
import math
import random


def valuation(x, p):
    if not x:
        return math.inf
    out = 0
    while x % p == 0:
        x //= p
        out += 1
    return out


def signature(nodes, p, vals):
    if len(nodes) == 1:
        return (-1, ())
    depth = min(vals[abs(x - nodes[0])] for x in nodes[1:])
    modulus = p ** (depth + 1)
    groups = {}
    for x in nodes:
        groups.setdefault(x % modulus, []).append(x)
    return (depth, tuple(sorted(signature(g, p, vals) for g in groups.values())))


def precision_constants(nodes, p):
    local = []
    for i, x in enumerate(nodes):
        differences = [x - y for j, y in enumerate(nodes) if i != j]
        derivative = math.prod(differences)
        second = 2 * sum(derivative // d for d in differences)
        s = valuation(derivative, p)
        b = 3 * s - valuation(second, p)
        local.append((s, b))
    return max(s for s, b in local), max(b for s, b in local)


def terminal_prediction(nodes, p):
    """Candidate tree formula, independently of F'' evaluation."""
    candidates = []

    def descend(cluster):
        if len(cluster) == 1:
            return
        depth = min(valuation(x - cluster[0], p) for x in cluster[1:])
        modulus = p ** (depth + 1)
        groups = {}
        for x in cluster:
            groups.setdefault(x % modulus, []).append(x)
        if all(len(g) == 1 for g in groups.values()):
            x = cluster[0]
            potential = sum(valuation(x - y, p) for y in nodes if y != x)
            candidates.append(2 * potential + depth - int(len(cluster) == p))
        else:
            for group in groups.values():
                descend(group)

    descend(nodes)
    return max(candidates)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--height', type=int, default=24)
    ap.add_argument('--max-nodes', type=int, default=6)
    ap.add_argument('--primes', type=int, nargs='+', default=[2, 3, 5, 7])
    ap.add_argument('--random', type=int, default=0,
                    help='seeded isometric high-unit-lift trials, instead of complete head')
    args = ap.parse_args()
    dig = hashlib.sha256()
    gates = 0
    if args.random:
        rng = random.Random(4435)
        for p in args.primes:
            for trial in range(args.random):
                n = rng.randrange(4, args.max_nodes + 1)
                depth = 1
                while p ** depth < 2 * n:
                    depth += 1
                depth += rng.randrange(3)
                modulus = p ** depth
                nodes = tuple(sorted(rng.sample(range(modulus), n)))
                lifted = tuple(x + modulus * rng.randrange(-100, 101) for x in nodes)
                left = precision_constants(nodes, p)
                right = precision_constants(lifted, p)
                for row, constants in ((nodes, left), (lifted, right)):
                    prediction = terminal_prediction(row, p)
                    if prediction != constants[1]:
                        print('TERMINAL FORMULA FAIL', p, row, constants, prediction, flush=True)
                        return
                gates += 1
                if left != right:
                    print('LIFT COUNTEREXAMPLE', p, nodes, lifted, left, right, flush=True)
                    return
                dig.update(repr((p, nodes, lifted, left)).encode())
            print('RANDOM', p, args.random, 'isometric high-unit lifts PASS', flush=True)
        print('GATES', gates, 'SEMANTIC', dig.hexdigest())
        print('FINITE-EXACT sampled lifts only; no unbounded theorem.')
        return
    # Canonical hostile: same metric and largest loss, different full lists.
    for p in args.primes:
        vals = [math.inf] + [valuation(d, p) for d in range(1, args.height + 1)]
        for n in range(4, args.max_nodes + 1):
            seen = {}
            count = 0
            for rest in itertools.combinations(range(1, args.height + 1), n - 1):
                nodes = (0,) + rest
                key = signature(nodes, p, vals)
                constants = precision_constants(nodes, p)
                prediction = terminal_prediction(nodes, p)
                if prediction != constants[1]:
                    print('TERMINAL FORMULA FAIL', p, nodes, constants, prediction, flush=True)
                    return
                count += 1
                gates += 1
                if key in seen:
                    old_constants, old_nodes = seen[key]
                    if constants != old_constants:
                        print('COUNTEREXAMPLE', p, n, old_nodes, nodes,
                              old_constants, constants, 'TREE', key, flush=True)
                        print('Status FINITE-EXACT; reconstruct literal inverses before promotion.')
                        return
                else:
                    seen[key] = (constants, nodes)
                dig.update(repr((p, nodes, key, constants)).encode())
                dig.update(b'\n')
            print('HEAD', p, n, 'rows', count, 'metric_types', len(seen), flush=True)
    print('GATES', gates, 'SEMANTIC', dig.hexdigest())
    print('NO COLLISION in the stated universe; no unbounded metric theorem follows.')


if __name__ == '__main__':
    main()
