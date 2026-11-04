#!/usr/bin/env python3
"""
deck_parity_source_lift_check_20261104.py  (opus-2026-11-04-S72)

Checks, with brute-force canonical forms (minimum over ALL vertex permutations, independent of the C programs'
refined canonical form), the source-lift explanation of the odd-order collisions of the converse-merged deck:

  every collision class of the unweighted converse-merged deck at odd order n in {5, 7, 9} whose members have no
  Hamiltonian cycle consists of {1 -> S, 1 -> S^op} (up to converse, i.e. possibly written as sink extensions),
  where {S, S^op} is one of the same-full-deck converse pairs of order n - 1 found by the census
  (deletion_tower_order10 10 full).

Codes are the census encoding: bits ordered by (q = 1..n-1, p = 0..q-1), most significant first, bit 1 iff p -> q.
The codes below are copied from deletion_tower_order10_20261104.out (collision lines).
"""
import itertools

EVEN_PAIRS = {   # order: list of (code, code) same full deck, converse of each other
    4: [(2, 16)],
    6: [(4897, 8324), (8280, 29797), (4417, 12356), (6241, 20612)],
    8: [(122729113, 123751452), (96520905, 226797752)],
}
ODD_COLLISIONS = {   # order: list of (code, code) equal unweighted merged decks, both without Hamiltonian cycles
    5: [(2, 32)],
    7: [(35393, 65796), (65688, 364741), (33409, 98436), (41153, 295172)],
    9: [(3565331737, 3573618716), (2617447753, 11142336824)],
}


def decode(n, code):
    tot = n * (n - 1) // 2
    b = tot - 1
    adj = [[0] * n for _ in range(n)]
    for q in range(1, n):
        for p in range(q):
            bit = (code >> b) & 1
            b -= 1
            if bit:
                adj[p][q] = 1
            else:
                adj[q][p] = 1
    return adj


def canon(adj):
    """minimum adjacency word over all vertex orders (exact, slow: n! permutations)"""
    n = len(adj)
    best = None
    # prune by score order: only permutations sorting vertices by score (scores are invariant), all ties tried
    scores = [sum(r) for r in adj]
    groups = {}
    for v in range(n):
        groups.setdefault(scores[v], []).append(v)
    keys = sorted(groups)
    for parts in itertools.product(*[itertools.permutations(groups[k]) for k in keys]):
        order = [v for part in parts for v in part]
        word = tuple(adj[order[p]][order[q]] for q in range(1, n) for p in range(q))
        if best is None or word < best:
            best = word
    return (n, best)


def converse(adj):
    n = len(adj)
    return [[adj[j][i] for j in range(n)] for i in range(n)]


def strip_source_or_sink(adj):
    """returns (kind, remainder) for a tournament with a source or a sink; prefers the source"""
    n = len(adj)
    for v in range(n):
        if sum(adj[v]) == n - 1:
            return 'source', [[adj[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    for v in range(n):
        if sum(adj[v]) == 0:
            return 'sink', [[adj[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    return None, None


ok_all = True
for n, pairs in ODD_COLLISIONS.items():
    even = {}
    for a, b in EVEN_PAIRS[n - 1]:
        A, B = decode(n - 1, a), decode(n - 1, b)
        assert canon(converse(A)) == canon(B), "census even pair is not a converse pair"
        even[frozenset([canon(A), canon(B)])] = (a, b)
    for x, y in pairs:
        X, Y = decode(n, x), decode(n, y)
        kx, Rx = strip_source_or_sink(X)
        ky, Ry = strip_source_or_sink(Y)
        # bring both to source form: a sink extension S -> 1 is the converse of the source extension 1 -> S^op
        if kx == 'sink':
            Rx = converse(Rx)
        if ky == 'sink':
            Ry = converse(Ry)
        key = frozenset([canon(Rx), canon(Ry)])
        hit = even.get(key)
        ok = hit is not None and canon(Rx) != canon(Ry)
        ok_all &= ok
        print(f"order {n}: collision ({x}, {y}): extensions ({kx}, {ky}); remainders form the order-{n - 1} converse same-deck pair {hit}: {'YES' if ok else 'NO'}")
print("ALL ODD-ORDER HAMILTONIAN-CYCLE-FREE COLLISIONS ARE SOURCE LIFTS OF EVEN-ORDER CONVERSE PAIRS:", ok_all)
