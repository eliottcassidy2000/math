#!/usr/bin/env python3
"""Finite exact controls for the recovered three-cycle boundary reduction.

Standard library; left permutation actions. No claimed Keller realization.
The all-size group and geometric arguments are in the companion note.
"""
from collections import Counter, deque
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod
import hashlib
import json

GATES = 0


def check(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)


def compose(a, b):
    return tuple(a[x] for x in b)


def inverse(a):
    q = [0] * len(a)
    for i, x in enumerate(a):
        q[x] = i
    return tuple(q)


def cycle(n, support):
    p = list(range(n))
    for a, b in zip(support, support[1:] + support[:1]):
        p[a] = b
    return tuple(p)


def closure(n, generators):
    identity = tuple(range(n))
    found = {identity}
    todo = deque([identity])
    while todo:
        p = todo.popleft()
        for g in generators:
            q = compose(g, p)
            if q not in found:
                found.add(q)
                todo.append(q)
    return found


def components(n, edges):
    out = []
    todo = set(range(n))
    while todo:
        part = {min(todo)}
        previous = None
        while previous != part:
            previous = set(part)
            for e in edges:
                if part.intersection(e):
                    part.update(e)
        todo.difference_update(part)
        out.append(tuple(sorted(part)))
    return out


def sign_on(p, part):
    return (-1) ** sum(p[a] > p[b] for a, b in combinations(part, 2))


def fixed(p):
    return {i for i, j in enumerate(p) if i == j}


def main():
    # Complete edge universe on five labels. Orientations need not be
    # enumerated: reversing a three-cycle replaces a generator by its inverse.
    n = 5
    triples = list(combinations(range(n), 3))
    all_perms = list(permutations(range(n)))
    histogram = Counter()
    for mask in range(1 << len(triples)):
        edges = [e for i, e in enumerate(triples) if mask >> i & 1]
        parts = components(n, edges)
        actual = closure(n, [cycle(n, e) for e in edges])
        expected = {p for p in all_perms
                    if all({p[i] for i in c} == set(c) and sign_on(p, c) == 1
                           for c in parts)}
        check(actual == expected, "all-five-label support-hypergraph theorem")
        expected_order = prod(factorial(len(c)) // 2 if len(c) >= 3 else 1
                              for c in parts)
        check(len(actual) == expected_order, "component product order")
        histogram[len(actual)] += 1
    check(sum(histogram.values()) == 1024, "complete hypergraph universe")

    # A connected sparse seven-label control tests repeated one-point gluing;
    # the disconnected six-label hostile shows why transitivity is necessary.
    check(len(closure(7, [cycle(7, e) for e in
                         ((0, 1, 2), (2, 3, 4), (4, 5, 6))])) == 2520,
          "one-point attachment chain gives A7")
    check(len(closure(6, [cycle(6, (0, 1, 2)), cycle(6, (3, 4, 5))])) == 9,
          "disconnected supports give C3 times C3")
    s4 = closure(4, [cycle(4, (0, 1, 2)), cycle(4, (0, 1, 2, 3))])
    check(len(s4) == 24, "containing a three-cycle is insufficient")

    sigma = cycle(6, (0, 3, 4))
    tau = cycle(6, (0, 1, 2))
    eta = cycle(6, (1, 2, 5))
    w = compose(sigma, tau)
    g = compose(w, w)
    a, b = fixed(sigma), fixed(tau)
    local = closure(6, [sigma, tau])
    global_group = closure(6, [sigma, tau, eta])
    check(len(local) == 60, "cusp A5")
    check({p[5] for p in local} == {5}, "cusp singleton retained")
    check({g[i] for i in a} == b, "actual-access control")
    check(compose(g, sigma) == compose(tau, g), "odd Artin relation")
    check(len(a & b) == 1, "one actual cusp sheet")
    check(len(global_group) == 360, "global A6 extension")
    check(compose(sigma, eta) == compose(eta, sigma), "node commuting meridians")
    check(not (fixed(sigma) & fixed(eta)), "node actual fibre zero")
    check(not ((set(range(6)) - fixed(sigma)) &
               (set(range(6)) - fixed(eta))), "node overlap zero")
    conjugates = {compose(compose(p, sigma), inverse(p)) for p in global_group}
    check(len(conjugates) == 40, "single-three-cycle conjugacy class")
    check(closure(6, conjugates) == global_group, "meridian normal generation")
    stabilizer = {p for p in global_group if p[0] == 0}
    check(len(stabilizer) == 60, "actual point stabilizer A5")
    for size in (2, 3):
        for tail in combinations(range(1, 6), size - 1):
            block = {0, *tail}
            check(any(({p[i] for i in block} & block) and
                      {p[i] for i in block} != block for p in global_group),
                  "no proper block/intermediate-degree shortcut")

    # Complete generic prime-data relaxation: each tuple is (e,f,actual).
    # A prime contributes f inertia cycles of length e. Actual primes have e=1.
    prime_types = [(e, f, actual) for e in range(1, 7)
                   for f in range(1, 7 // e + 1)
                   for actual in (False, True) if not actual or e == 1]
    prime_rows = []
    for count in range(1, 7):
        for row in combinations_with_replacement(prime_types, count):
            if sum(e * f for e, f, actual in row) != 6:
                continue
            if sorted(e for e, f, actual in row for _ in range(f)) != [1, 1, 1, 3]:
                continue
            if sum(f for e, f, actual in row if actual) != 3:
                continue
            deleted = [(e, f) for e, f, actual in row if not actual]
            check(deleted == [(3, 1)], "unique deleted prime from actual count")
            prime_rows.append(row)
    check(len(prime_rows) == 3, "three possible retained residue-degree partitions")
    hostile_primes = ((1, 1, False), (1, 2, True), (3, 1, False))
    check(sorted(e for e, f, _ in hostile_primes for _ in range(f)) == [1, 1, 1, 3],
          "same inertia with extra deleted fixed sheet")
    check(sum(f for e, f, actual in hostile_primes if actual) == 2,
          "lost actual-count sidecar identified")

    payload = {"hypergraph_histogram": sorted(histogram.items()),
               "sigma": sigma, "tau": tau, "eta": eta,
               "local_order": len(local), "global_order": len(global_group),
               "prime_rows": prime_rows, "hostile_primes": hostile_primes}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True,
                                      separators=(",", ":")).encode()).hexdigest()
    print("STATUS: FINITE-EXACT controls; all-size and geometric proofs in note")
    print("hypergraph_universe: all 1024 three-uniform hypergraphs on five labels")
    print("hypergraph_order_histogram=" + json.dumps(sorted(histogram.items())))
    print("named_groups: connected seven-label=2520; disconnected six-label=9; containing-only hostile=24")
    print("semilocal_fixture=" + json.dumps({k: payload[k] for k in
          ("sigma", "tau", "eta", "local_order", "global_order")}, sort_keys=True))
    print("prime_universe: all unordered (e,f,actual) rows of total degree six, inertia 3111, actual count three")
    print("prime_rows=" + json.dumps(prime_rows))
    print("scope: no full curve-complement representation or Keller realization asserted")
    print("semantic_sha256=" + digest)
    print("PASS gates=" + str(GATES))


if __name__ == "__main__":
    main()
