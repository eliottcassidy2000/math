"""Independent bounded-profile audit; no producer algorithm imports.

The abstract transfer uses endpoint TYPE multiplicities, not component-slot
loops. Concrete checks enumerate edge subsets and directed paths directly.
"""
from collections import Counter, deque
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import math
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'frozen_source'
BULK = -1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def shape(counts):
    return tuple(sorted(p for p, n in counts.items() for _ in range(n)))


def independent_transfer(vector, label, chosen=None, cap=None):
    """chosen=None means forget; otherwise use aggregated type multiplicity."""
    result = Counter()
    for pattern, weight in vector:
        if chosen is None:
            target = tuple(sorted((BULK if a == label else a,
                                   BULK if b == label else b)
                                  for a, b in pattern))
            if sum(a == BULK for a, _ in target) <= 2:
                result[target] += weight
            continue
        types = Counter(pattern)
        singleton = types.copy()
        singleton[(label, label)] += 1
        result[shape(singleton)] += weight
        incoming = [t for t in types if t[1] == BULK or t[1] not in chosen]
        outgoing = [t for t in types if t[0] in chosen]
        for t in incoming:
            child = types.copy()
            child[t] -= 1
            child[(t[0], label)] += 1
            result[shape(child)] += weight * types[t]
        for t in outgoing:
            child = types.copy()
            child[t] -= 1
            child[(label, t[1])] += 1
            result[shape(child)] += weight * types[t]
        for left in incoming:
            for right in outgoing:
                choices = types[left] * (types[right] - (left == right))
                if not choices:
                    continue
                child = types.copy()
                child[left] -= 1
                child[right] -= 1
                child[(left[0], right[1])] += 1
                result[shape(child)] += weight * choices
    return tuple(sorted((p, min(n, cap) if cap is not None else n)
                        for p, n in result.items() if n))


def observe(vector, root, cap=None):
    # Direct nonnegative coefficients, including the coefficient-two overlap.
    b = s = 0
    for p, n in vector:
        if len(p) == 1:
            s += n * (p[0][0] == root)
            b += n * (p[0][0] != root)
        b += n * (len(p) in (1, 2) and any(e == root for _, e in p))
    return (min(b, cap), min(s, cap)) if cap is not None else (b, s)


def width(k):
    # Integer test rather than the producer's isqrt formula.
    def m(t):
        return max([0] + [n for n in range(1, 2 * t + 1) if n * n < 4 * t])
    return max(1 + m(a) + m(k - a) for a in range(k + 1))


def independent_bfs(k, x, limit=1000000):
    w, c = width(k), max(k, x) + 1
    start = (frozenset(), (((), 1),))
    queue = deque([start])
    visited = {start}
    answers = set()
    edges = 0
    max_patterns = 0
    status = 'COMPLETE'
    while queue:
        active, vector = queue.popleft()
        max_patterns = max(max_patterns, len(vector))
        for r in active:
            b, s = observe(vector, r, c)
            if 1 <= b <= x and 1 <= s <= k:
                answers.add((b, s))
        for label in range(w):
            if label in active:
                actions = [None]
            else:
                a = sorted(active)
                actions = [frozenset(a[j] for j in range(len(a)) if mask & (1 << j))
                           for mask in range(1 << len(a))]
            for chosen in actions:
                edges += 1
                newactive = active - {label} if chosen is None else active | {label}
                successor = (newactive, independent_transfer(vector, label, chosen, c))
                if successor in visited:
                    continue
                if len(visited) >= limit:
                    status = 'INCOMPLETE'
                    queue.clear()
                    break
                visited.add(successor)
                queue.append(successor)
            if status == 'INCOMPLETE':
                break
    return dict(status=status, K=k, X=x, W=w, cutoff=c, states=len(visited),
                transitions=edges, pairs=sorted(answers), max_patterns=max_patterns)


def concrete_forests(rows, active):
    """Enumerate actual directed edge subsets, never transfer coordinates."""
    n = len(rows)
    arcs = [(a, b) for a in range(n) for b in range(n) if rows[a] >> b & 1]
    result = Counter()
    for bits in range(1 << len(arcs)):
        succ, pred = {}, {}
        okay = True
        for j, (a, b) in enumerate(arcs):
            if bits >> j & 1:
                if a in succ or b in pred:
                    okay = False
                    break
                succ[a], pred[b] = b, a
        if not okay:
            continue
        paths, visited = [], set()
        for first in range(n):
            if first in pred:
                continue
            v = first
            while True:
                if v in visited:
                    okay = False
                    break
                visited.add(v)
                if v not in succ:
                    break
                v = succ[v]
            paths.append((active.get(first, BULK), active.get(v, BULK)))
        if okay and len(visited) == n and sum(a == BULK for a, _ in paths) <= 2:
            result[tuple(sorted(paths))] += 1
    return tuple(sorted(result.items()))


def direct_profile(rows, root):
    """DFS of all directed simple paths; uses no subset DP recurrence."""
    n, whole = len(rows), (1 << len(rows)) - 1
    counts = Counter({0: 1})
    ending = Counter()
    starts_at_root = 0
    def visit(first, last, mask):
        nonlocal starts_at_root
        counts[mask] += 1
        if last == root:
            ending[mask] += 1
        if mask == whole and first == root:
            starts_at_root += 1
        new = rows[last] & (whole ^ mask)
        while new:
            bit = new & -new
            new -= bit
            visit(first, bit.bit_length() - 1, mask | bit)
    for first in range(n):
        visit(first, first, 1 << first)
    h = counts[whole]
    ell = sum(number * counts[whole ^ mask] for mask, number in ending.items())
    return h + ell - starts_at_root, starts_at_root


def replay(word, width_bound, root_vertex):
    active, rows = {}, []
    reused = 0
    labels_seen = set()
    for op in word:
        require(op[0] in ('introduce', 'forget'), 'unrecognized word instruction')
        label = op[1]
        require(0 <= label < width_bound, 'witness label outside width')
        if op[0] == 'forget':
            require(len(op) == 2 and label in active, 'illegal forget')
            del active[label]
        else:
            chosen = op[2]
            require(label not in active and len(chosen) == len(set(chosen))
                    and set(chosen) <= active.keys(), 'illegal introduction')
            require(len(active) < width_bound, 'width overflow')
            new = len(rows)
            rows.append(0)
            beaten = {active[p] for p in chosen}
            for old in range(new):
                a, b = (new, old) if old in beaten else (old, new)
                rows[a] |= 1 << b
            active[label] = new
            reused += label in labels_seen
            labels_seen.add(label)
    require(root_vertex in active.values(), 'final root is forgotten')
    return rows, reused


def main():
    started = time.monotonic()
    results = {}
    # Independent census, including the full claimed 137 rectangle.
    cases = [(1, 1, 1000000), (1, 5, 1000000), (1, 13, 1000000),
             (1, 137, 1000000), (2, 13, 1000), (1, 1, 1)]
    independent = [independent_bfs(*case) for case in cases]
    results['independent_searches'] = independent
    for check, filename in [(independent[3], 'NORMAL137.json'),
                            (independent[4], 'GUARD_K2.json'),
                            (independent[5], 'GUARD_EMPTY.json')]:
        saved = json.loads((SOURCE / filename).read_text())
        for key in ('status', 'W', 'cutoff', 'states', 'transitions'):
            require(check[key] == saved[key], f'{filename}: {key} mismatch')
        require(check['pairs'] == [tuple(p) for p in saved['pairs']], f'{filename}: pairs mismatch')
    results['saved_result_comparison'] = 'PASS'
    # Exhaust all labelled tournaments through order 5. Use independently
    # counted edge-subset forests with both full and partly forgotten ports.
    count_graphs = count_vectors = count_roots = 0
    witness_pairs = set()
    rng = random.Random(203)
    for n in range(1, 6):
        edgepairs = list(combinations(range(n), 2))
        for bits in range(1 << len(edgepairs)):
            rows = [0] * n
            for j, (a, b) in enumerate(edgepairs):
                if bits >> j & 1:
                    a, b = b, a
                rows[a] |= 1 << b
            vector = (((), 1),)
            for u in range(n):
                vector = independent_transfer(vector, u, {v for v in range(u) if rows[u] >> v & 1})
            require(vector == concrete_forests(rows, {v: v for v in range(n)}), 'full literal forest mismatch')
            count_vectors += 1
            for r in range(n):
                actual = direct_profile(rows, r)
                require(observe(vector, r) == actual, 'active-port readout mismatch')
                count_roots += 1
                if actual[1] == 1 and actual[0] <= 137:
                    witness_pairs.add(actual)
            # Root-preserving arbitrary forget order; check after every step.
            r = rng.randrange(n)
            forgotten = [v for v in range(n) if v != r]
            rng.shuffle(forgotten)
            active = {v: v for v in range(n)}
            for u in forgotten:
                vector = independent_transfer(vector, u)
                del active[u]
                require(vector == concrete_forests(rows, active), 'forget literal forest mismatch')
                require(observe(vector, r) == direct_profile(rows, r), 'forget changed readout')
                count_vectors += 1
            count_graphs += 1
    require(witness_pairs <= set(independent[3]['pairs']), 'small tournament profile absent from completed BFS')
    results['literal_forest_checks'] = dict(tournaments=count_graphs, vectors=count_vectors, roots=count_roots)
    # Verify every saved physical witness without producer DP/runtime/replay.
    checked, reuses, max_order = 0, 0, 0
    for filename in ('NORMAL137.json', 'OPTIMIZED137.json', 'GUARD_K2.json', 'GUARD_EMPTY.json'):
        saved = json.loads((SOURCE / filename).read_text())
        require(saved['literal_witnesses_checked'] == len(saved['witnesses']), 'witness count mismatch')
        require(saved['pairs'] == [[w['B'], w['s']] for w in saved['witnesses']], 'pair/witness mismatch')
        for w in saved['witnesses']:
            rows, reuse = replay(w['word'], saved['W'], w['root'])
            require(rows == w['rows'], 'word/graph mismatch')
            require(direct_profile(rows, w['root']) == (w['B'], w['s']), 'direct witness profile mismatch')
            checked += 1
            reuses += reuse
            max_order = max(max_order, len(rows))
    results['literal_saved_witnesses'] = dict(checked=checked, label_reuse_events=reuses, max_order=max_order)
    # Abstract cap property and signed hostile control. A zero coefficient can
    # erase a huge state, but nonzero coefficients cannot expose a false small.
    cap_cases = 0
    for cap in range(2, 8):
        for values in product(range(2 * cap + 1), repeat=2):
            for coeffs in product(range(3), repeat=2):
                lhs = min(cap, sum(a * x for a, x in zip(coeffs, values)))
                rhs = min(cap, sum(a * min(cap, x) for a, x in zip(coeffs, values)))
                require(lhs == rhs, 'cap scalar congruence')
                cap_cases += 1
    require(min(5, 7 - 6) != min(5, min(5, 7) - min(5, 6)), 'signed hostile failed')
    results['cap_arithmetic_cases'] = cap_cases
    results['signed_readout_hostile'] = {'cap': 5, 'true_difference': 1, 'difference_of_caps': 0}
    # Runtime CLI guards, including -O. No generated bytecode in producer dir.
    cli_checks = []
    for optimized in (False, True):
        for k, x, limit, expected in [(1, 1, 1, 'INCOMPLETE'), (2, 13, 1000, 'INCOMPLETE'),
                                      (1, 1, 144, 'INCOMPLETE'), (1, 1, 145, 'COMPLETE'),
                                      (1, 5, 10000, 'COMPLETE'), (1, 137, 100000, 'COMPLETE')]:
            command = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(SOURCE / 'check_membership.py'),
                       '--K', str(k), '--X', str(x), '--max-states', str(limit)]
            saved = json.loads(subprocess.check_output(command, text=True))
            require(saved['status'] == expected, 'runtime guard status')
            if (k, x) == (1, 137):
                frozen = 'OPTIMIZED137.json' if optimized else 'NORMAL137.json'
                require(saved == json.loads((SOURCE / frozen).read_text()), 'full producer replay mismatch')
            for witness in saved['witnesses']:
                require(direct_profile(witness['rows'], witness['root']) == (witness['B'], witness['s']),
                        'guard witness mismatch')
            cli_checks.append(dict(optimized=optimized, K=k, X=x, limit=limit,
                                   status=saved['status'], states=saved['states'], pairs=saved['pairs']))
    for flag in ('--K', '--X', '--max-states'):
        result = subprocess.run([sys.executable, '-B', str(SOURCE / 'check_membership.py'), flag, '0'],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(result.returncode != 0 and 'positive integer' in result.stderr, 'invalid CLI input not rejected')
    results['producer_cli_checks'] = cli_checks
    results['invalid_cli_arguments_rejected'] = 3
    results['input_hashes'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in SOURCE.iterdir() if p.is_file()}
    results['status'] = 'PASS'
    results['elapsed_seconds'] = round(time.monotonic() - started, 3)
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
