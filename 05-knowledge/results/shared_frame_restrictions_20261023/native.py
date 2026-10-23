#!/usr/bin/env python3
"""Actual native source/partner transport, with retained kernel lifts.

The clock guard is a declared synthetic port observation, NOT the repository's
electrical mod-8 invariant. An application supplies its own exact guard table.
"""
from collections import Counter, defaultdict
from itertools import combinations, product
from math import prod
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def arc(T, u, v):
    return (u, v) in T


def path_ok(T, path):
    return len(set(path)) == len(path) and all(arc(T, u, v) for u, v in zip(path, path[1:]))


def partner(T, path):
    if not path_ok(T, path):
        raise ValueError("input path is not authenticated")
    keep = set(zip(path, path[1:]))
    return {(u, v) if (u, v) in keep else (v, u) for u, v in T}


def transport(T, f):
    return {(f[u], f[v]) for u, v in T}


def converse(T):
    return {(v, u) for u, v in T}


def compile_lifts(records, guard):
    """Evaluate the predicate on each FULL lift; never one quotient representative.

    Return source-indexed records and exact multiplicities. The explicit finite
    interface makes no small-description or polynomial-complexity claim.
    """
    accepted = [r for r in records if guard(r)]
    quotient_weights, ports = Counter(), Counter()
    for r in accepted:
        weight = r.get('weight', 1)
        if not isinstance(weight, int) or weight <= 0:
            raise ValueError('source multiplicity must be a positive integer')
        quotient_weights[r['quotient']] += weight
        ports[tuple(r['ports'])] += weight
    return accepted, quotient_weights, ports


def main():
    m = 9
    n = 2 * m
    T = set()
    for u, v in combinations(range(n), 2):
        if u // m == v // m:
            forward = 1 <= (v - u) % m <= (m - 1) // 2
        else:
            forward = (u - v) % 3 == 0
        T.add((u, v) if forward else (v, u))
    P = list(range(m)) + [m + (2 + i) % m for i in range(m)]
    require(len(T) == n * (n - 1) // 2 and path_ok(T, P), "native source/path")
    Q = partner(T, P)
    require(path_ok(Q, P) and partner(Q, P) == T, "native partner involution")
    records = []
    rejected = 0
    for a, b in product(range(m), repeat=2):
        f = [(i + a) % m for i in range(m)] + [m + (i + b) % m for i in range(m)]
        actual = transport(T, f) == T
        predicted = (a - b) % 3 == 0
        require(actual == predicted, "actual block relation versus congruence")
        if not actual:
            rejected += 1
            continue
        path = [f[v] for v in P]
        require(path_ok(T, path), "transported source path")
        target = partner(T, path)
        require(target == transport(Q, f), "partner equivariance")
        require(partner(converse(T), list(reversed(path))) == converse(target), "converse square")
        clock = a + 2 * b  # Integer representative, deliberately not a group homomorphism.
        records.append(dict(frame=[a, b], quotient=a % 3, ports=[path[0], path[-1]],
                            integer_clock=clock, clock_mod8=clock % 8, path=path))
    require(len(records) == 27 and rejected == 54, "complete cyclic relation")
    require(len({tuple(r['path']) for r in records}) == 27, "distinct indexed paths")
    require(len({tuple(r['ports']) for r in records}) == 27, "all endpoint lifts retained")
    accepted, weights, ports = compile_lifts(records, lambda r: r['clock_mod8'] == 0)
    require(sum(weights.values()) == len(accepted) == sum(ports.values()), "guard source mass")
    # Direct independent guard join of the two congruences and least representatives.
    direct = {(a, b) for a in range(9) for b in range(9)
              if (a - b) % 3 == 0 and (a + 2 * b) % 8 == 0}
    require({tuple(r['frame']) for r in accepted} == direct, "guard join")
    all_hist = defaultdict(Counter)
    for r in records:
        all_hist[r['quotient']][r['clock_mod8']] += 1
    require(any(len(h) > 1 for h in all_hist.values()), "quotient-only guard hostile")
    # Each quotient has kernel C3 x C3. The chosen section does NOT preserve ports/guard.
    representatives = [r for r in records if r['frame'] == [r['quotient'], r['quotient']]]
    fake, fake_weights, _ = compile_lifts(representatives, lambda r: r['clock_mod8'] == 0)
    require(len(fake) * 9 != len(accepted), "one representative times kernel size fails")
    weighted = [dict(r, weight=1+sum(r['frame'])) for r in records]
    _, weighted_q, weighted_ports = compile_lifts(weighted, lambda r: r['clock_mod8'] == 0)
    replicated = [dict(r, weight=1) for r in weighted for _ in range(r['weight'])]
    _, replicated_q, replicated_ports = compile_lifts(replicated, lambda r: r['clock_mod8'] == 0)
    require(weighted_q == replicated_q and weighted_ports == replicated_ports,
            "source weights versus literal repeated copies")
    require(sum(weighted_q.values()) == 18, "weighted guarded mass")
    # Marking the actual first port removes one factor but keeps three second-port lifts.
    pinned, _, pinned_ports = compile_lifts(records, lambda r: r['ports'][0] == 0)
    require(len(pinned) == 3 and len(pinned_ports) == 3, "one pinned endpoint")
    require(len([r for r in pinned if r['ports'][1] == P[-1]]) == 1, "both pinned endpoints")
    # A quotient return leaves a nontrivial kernel displacement.
    return_lift = next(r for r in records if r['frame'] == [3, 3])
    require(return_lift['quotient'] == 0 and return_lift['ports'] != [P[0], P[-1]],
            "kernel carry on quotient return")
    # Strong connectivity is tested by reachability, not inferred from degree.
    for start in range(n):
        seen, todo = {start}, [start]
        while todo:
            u = todo.pop()
            for v in range(n):
                if v not in seen and (u, v) in T:
                    seen.add(v)
                    todo.append(v)
        require(len(seen) == n, "strong native source")
    out = dict(status="FINITE-EXACT PASS", vertices=n, source_arcs=len(T),
               cyclic_candidates=81, compatible_lifts=27, quotient_states=3,
               kernel_size=9, literal_paths=27, distinct_endpoint_pairs=27,
               converse_squares=27, source_strong=True,
               synthetic_guard_description="(least(a)+2 least(b)) mod8 = 0; NOT electrical mod8",
               accepted_frames=[r['frame'] for r in accepted],
               accepted_quotient_weights=dict(weights),
               guard_histogram_by_quotient={q: dict(sorted(h.items())) for q, h in all_hist.items()},
               false_section_scaled_mass=9 * len(fake), actual_guarded_mass=len(accepted),
               weighted_guarded_mass=sum(weighted_q.values()),
               one_pinned_endpoint_lifts=len(pinned), quotient_return_ports=return_lift['ports'])
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
