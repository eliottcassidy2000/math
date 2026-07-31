#!/usr/bin/env python3
"""Exact small-p audit of the upward-event transport/cover identity.

For every monotone Boolean event on p <= 4 needles and every marginal vector
with coordinates in {0, 1/2, 1}, compare

  max P(X in A), over all laws with the prescribed one-marginals,

against

  min(1, min_w r.w),  w >= 0 and sum_{i in H} w_i >= 1

for every inclusion-minimal true set H of A.

Both LPs are solved by exhaustive vertex enumeration over Fraction arithmetic.
This is only a hostile finite control; the accompanying blocker/packing proof
establishes the identity for arbitrary p and arbitrary real marginals.
"""

from fractions import Fraction
from itertools import combinations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solve_square(rows, rhs):
    """Solve a square rational system, or return None when singular."""
    n = len(rows)
    a = [[Fraction(v) for v in row] + [Fraction(rhs[i])]
         for i, row in enumerate(rows)]
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [v / scale for v in a[col]]
        for row in range(n):
            if row == col or not a[row][col]:
                continue
            scale = a[row][col]
            a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
    return tuple(a[i][-1] for i in range(n))


def monotone_events(p):
    states = range(1 << p)
    ans = []
    for event_bits in range(1 << (1 << p)):
        good = True
        for s in states:
            if not (event_bits >> s) & 1:
                continue
            for i in range(p):
                if not (s >> i) & 1 and not (event_bits >> (s | (1 << i))) & 1:
                    good = False
                    break
            if not good:
                break
        if good:
            ans.append(event_bits)
    return ans


def minimal_true_sets(p, event_bits):
    out = []
    for s in range(1 << p):
        if not (event_bits >> s) & 1:
            continue
        if all(not ((event_bits >> (s ^ (1 << i))) & 1)
               for i in range(p) if (s >> i) & 1):
            out.append(s)
    return out


def transport_vertices(p, r):
    """Distinct basic feasible laws with total mass 1 and marginals r."""
    states = tuple(range(1 << p))
    columns = {
        s: (1,) + tuple((s >> i) & 1 for i in range(p))
        for s in states
    }
    rhs = (Fraction(1),) + tuple(r)
    vertices = set()
    for basis in combinations(states, p + 1):
        rows = [[columns[s][row] for s in basis] for row in range(p + 1)]
        weights = solve_square(rows, rhs)
        if weights is None or any(weight < 0 for weight in weights):
            continue
        vertices.add(
            tuple(
                (state, weight)
                for state, weight in zip(basis, weights)
                if weight
            )
        )
    require(vertices, ("no transport vertex", p, r))
    return tuple(sorted(vertices))


def transport_value(event_bits, vertices):
    return max(
        sum((weight for state, weight in vertex
             if (event_bits >> state) & 1), Fraction(0))
        for vertex in vertices
    )


def cover_vertices(p, hyperedges):
    """Vertices of {w >= 0 : w(H) >= 1 for H in hyperedges}."""
    constraints = []
    # A tag of 0 denotes w_i=0; a tag of 1 denotes w(H)=1.
    constraints.extend((0, i) for i in range(p))
    constraints.extend((1, h) for h in hyperedges)
    vertices = set()
    for tight in combinations(constraints, p):
        rows = []
        rhs = []
        for kind, datum in tight:
            if kind == 0:
                rows.append([int(i == datum) for i in range(p)])
                rhs.append(0)
            else:
                rows.append([(datum >> i) & 1 for i in range(p)])
                rhs.append(1)
        w = solve_square(rows, rhs)
        if w is None or any(x < 0 for x in w):
            continue
        if any(sum(w[i] for i in range(p) if (h >> i) & 1) < 1
               for h in hyperedges):
            continue
        vertices.add(w)
    require(vertices, ("no cover vertex", p, hyperedges))
    return tuple(vertices)


def cover_value(r, vertices):
    return min(sum(ri * wi for ri, wi in zip(r, w))
               for w in vertices)


def main():
    total_pairs = 0
    event_counts = []
    vertex_counts = []
    for p in range(1, 5):
        events = monotone_events(p)
        event_counts.append(len(events))
        grid = tuple(product((Fraction(0), Fraction(1, 2), Fraction(1)),
                             repeat=p))
        laws = {r: transport_vertices(p, r) for r in grid}
        vertex_counts.append(sum(len(v) for v in laws.values()))
        for event_bits in events:
            hyperedges = minimal_true_sets(p, event_bits)
            if not hyperedges:
                # Empty event.
                expected = Fraction(0)
                for r in grid:
                    got = transport_value(event_bits, laws[r])
                    require(
                        got == expected,
                        ("empty event mismatch", p, event_bits, r, got),
                    )
                    total_pairs += 1
                continue
            if hyperedges == [0]:
                # Full event; its empty minimal true set makes the cover
                # infeasible, conventionally tau = +infinity.
                expected = Fraction(1)
                for r in grid:
                    got = transport_value(event_bits, laws[r])
                    require(
                        got == expected,
                        ("full event mismatch", p, event_bits, r, got),
                    )
                    total_pairs += 1
                continue
            covers = cover_vertices(p, hyperedges)
            for r in grid:
                got = transport_value(event_bits, laws[r])
                expected = min(Fraction(1), cover_value(r, covers))
                require(
                    got == expected,
                    ("cover identity mismatch", p, event_bits, r, got, expected),
                )
                total_pairs += 1
    require(event_counts == [3, 6, 20, 168], "monotone event universe changed")
    require(total_pairs == 14211, "event/marginal universe changed")
    require(
        vertex_counts == [3, 10, 38, 192],
        "distinct basic-feasible-law census changed",
    )
    print("PASS exact upward-event transport/cover audit")
    print(f"monotone_event_counts={event_counts}")
    print("marginal_grids=[3,9,27,81]")
    print(f"event_marginal_pairs={total_pairs}")
    print(f"basic_feasible_laws_total_by_p={vertex_counts}")


if __name__ == "__main__":
    main()
