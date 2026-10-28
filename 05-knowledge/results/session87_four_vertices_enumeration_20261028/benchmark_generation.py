#!/usr/bin/env python3
"""Controlled stdlib comparison: same canonicalizer/cache, different generation."""
from collections import Counter
from statistics import median
from time import perf_counter
import json
from check_generation import ParentCanonicalizer, delete_root, extend, matrix, need, produce


def full_q(parents, n):
    engine = ParentCanonicalizer(max_order=n)
    states = {}
    for raw in parents:
        code, _ = engine.maps(raw)
        rows = matrix(n, code)
        for mask in range(1 << n):
            key = engine.rooted(rows, mask)
            states.setdefault(key, extend(matrix(n, key[0]), key[1]))
    representatives = {}
    for key, child in sorted(states.items()):
        row = Counter(engine.rooted(*delete_root(child, v)) for v in range(n + 1))
        signature = tuple(sorted(row.items()))
        representatives.setdefault(signature, child)
    return list(representatives.values())


def main():
    parents = [()]
    results = []
    for n in range(8):
        if n < 6:
            parents, _ = produce(parents, n, audit=False, selector='neighbor_degrees')
            continue
        timings = {'full_Q_same_cache': [], 'pruned_neighbor_degrees': []}
        outputs = {}
        for repeat in range(3):
            methods = ['full_Q_same_cache', 'pruned_neighbor_degrees']
            if repeat % 2:
                methods.reverse()
            for method in methods:
                start = perf_counter()
                if method == 'full_Q_same_cache':
                    children = full_q(parents, n)
                else:
                    children, _ = produce(parents, n, audit=False, selector='neighbor_degrees')
                timings[method].append(perf_counter() - start)
                outputs[method] = children
        oracle = ParentCanonicalizer()
        original = {oracle.maps(t)[0] for t in outputs['full_Q_same_cache']}
        pruned = {oracle.maps(t)[0] for t in outputs['pruned_neighbor_degrees']}
        need(original == pruned, 'timed algorithms have unequal class sets')
        need(len(pruned) == len(outputs['pruned_neighbor_degrees']), 'pruned duplicate')
        medians = {k: median(v) for k, v in timings.items()}
        results.append({'child_order': n + 1, 'parents': len(parents), 'classes': len(pruned),
                        'seconds_per_run': timings, 'median_seconds': medians,
                        'ratio_of_medians': medians['full_Q_same_cache'] / medians['pruned_neighbor_degrees']})
        parents = outputs['pruned_neighbor_degrees']
    return {'status': 'MEASURED_LOCAL_NOT_ASYMPTOTIC', 'repeat_count': 3,
            'comparison': 'Both algorithms use the same exact parent canonicalizer and literal-matrix cache. '
                          'Timing excludes test oracle and lower-level generation; includes parent preparation. '
                          'Reference has full rooted bank and Q rows, not the uncached historical executable.',
            'results': results}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
