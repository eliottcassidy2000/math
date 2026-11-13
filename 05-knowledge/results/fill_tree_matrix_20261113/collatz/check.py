"""Exact controls for guarded inverse trees and the literal relaxed graph."""
from collections import deque
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def oddpart(n):
    need(n > 0, 'positive odd-part input')
    return n // (n & -n)


def native_successors(n, doubled=True):
    if n % 2 == 0:
        return [n // 2, 2 * n] if doubled else [n // 2]
    out = [3 * n + 1]
    if n > 1 and n % 3 == 1:
        out.append((n - 1) // 3)
    return out


def reverse_reach(cap, mode):
    reached = {1}
    frontier = deque([1])
    while frontier:
        value = frontier.popleft()
        parents = [2 * value]
        if value % 6 == 4:
            parents.append((value - 1) // 3)
        if mode != 'original' and value % 2 == 0:
            parents.append(3 * value + 1)
        if mode == 'literal' and value % 4 == 0:
            parents.append(value // 2)
        for parent in parents:
            if 0 < parent <= cap and parent not in reached:
                reached.add(parent)
                frontier.append(parent)
    return reached


def odd_reverse_reach(cap):
    """Independent odd-word computation, retaining each macro edge's peak."""
    reached = {1}
    frontier = deque([1])
    while frontier:
        value = frontier.popleft()
        parents = []
        power = 2
        while power * value <= cap:
            numerator = power * value - 1
            if numerator % 3 == 0:
                parent = numerator // 3
                if parent > 0 and parent % 2:
                    parents.append(parent)
            power *= 2
        power = 2
        while 3 * power * value + 1 <= cap:
            parents.append(3 * power * value + 1)
            power *= 2
        for parent in parents:
            if parent not in reached:
                reached.add(parent)
                frontier.append(parent)
    return reached


def policy(n):
    if n % 2 == 0:
        return n // 2
    if n > 1 and n % 3 == 1:
        return (n - 1) // 3
    return 3 * n + 1


def path_to_one(n, limit=10000):
    path = [n]
    seen = {n}
    while path[-1] != 1:
        nxt = policy(path[-1])
        need(nxt in native_successors(path[-1]), 'illegal policy edge')
        need(nxt not in seen and len(path) < limit, 'policy cycle or exhausted cap')
        path.append(nxt)
        seen.add(nxt)
    return path


# Three equivalent finite presentations, with the original graph as control.
boxes = []
for cap in (100, 1000, 10000, 100000):
    ordinary = reverse_reach(cap, 'original')
    literal = reverse_reach(cap, 'literal')
    reduced = reverse_reach(cap, 'no_doubling')
    odd = odd_reverse_reach(cap)
    need(literal == reduced, 'even doubling changed finite reachability')
    need(ordinary <= literal, 'relaxation lost an original path')
    need(all((n in literal) == (oddpart(n) in odd)
             for n in range(1, cap + 1)), 'odd reduction mismatch')
    missing = [n for n in range(1, cap + 1) if n not in literal]
    boxes.append(dict(cap=cap, ordinary_reached=len(ordinary),
                      literal_reached=len(literal),
                      first_ten_literal_holes=missing[:10]))

# Symbolic-family controls; the proof, not this finite range, gives all K.
forced = []
for exponent in range(2, 202, 2):
    states = [3 ** j * 2 ** (exponent-j) - 1 for j in range(exponent)]
    for j, value in enumerate(states):
        need(value > 1 and value % 2 and value % 3 != 1,
             'inverse guard unexpectedly available')
        if j + 1 < exponent:
            need(3 * value + 1 == 2 * states[j+1], 'forced branch identity')
            need(states[j+1] > value, 'forced increase')
    floor = 2 * 3 ** exponent - 2
    need(3 * states[-1] + 1 == floor, 'native peak floor')
    if exponent <= 40:
        path = path_to_one(states[0])
        need(max(path) == floor, 'finite sharpness control')
        forced.append(dict(exponent=exponent, start=states[0], minimum_peak=floor,
                           native_steps=len(path)-1, path=path))

# Finite reachability certificates, not a convergence theorem.
known = {1}
max_visited = 1
for start in range(1, 100001):
    n = start
    chain = []
    local = set()
    while n not in known:
        need(n not in local and len(chain) < 10000, 'finite policy control failed')
        local.add(n)
        chain.append(n)
        max_visited = max(max_visited, n)
        n = policy(n)
    known.update(chain)


def plus_odd(n):
    return oddpart(3*n + 1)


def minus_odd(n):
    return oddpart(3*n - 1)


plus_paths = [[7, 11, 17, 13, 5, 1], [21, 1]]
minus_paths = [[5, 7, 5], [17, 25, 37, 55, 41, 61, 91, 17],
               [21, 31, 23, 17]]
for paths, transition in ((plus_paths, plus_odd), (minus_paths, minus_odd)):
    for path in paths:
        need(all(transition(a) == b for a, b in zip(path, path[1:])), 'named orbit')

sources = [
    '01-canon/theorems/THM-4514-collatz-value-time-poset-spine-blocks-and-two-place-descent-tree.md',
    '01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md',
    '05-knowledge/results/collatz_mod6_20260917_reverse_tree_pieces.md',
    '05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md',
    '05-knowledge/results/clone_clocks_negative_collatz_20261113/COLLATZ.md',
]
pins = [dict(path=path, sha256_lf=hashlib.sha256(
    (ROOT/path).read_bytes().replace(b'\r\n', b'\n')).hexdigest()) for path in sources]
report = dict(status='FINITE-EXACT controls; all-order statements proved in NOTE.md',
              boxes=boxes, family_identity_exponents=list(range(2, 202, 2)),
              sharp_minimum_peak_certificates=forced,
              policy_control=dict(all_starts_through=100000, step_cap_per_unresolved_chain=10000,
                                  distinct_certified_nodes=len(known), max_visited=max_visited),
              plus_odd_paths=plus_paths, minus_odd_paths=minus_paths,
              source_pins=pins)
suffix = '.optimized' if sys.flags.optimize else ''
(HERE/('receipt'+suffix+'.json')).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
print(json.dumps(dict(status='PASS', boxes=boxes, sharp_examples=len(forced),
                     policy_control=report['policy_control'])))
