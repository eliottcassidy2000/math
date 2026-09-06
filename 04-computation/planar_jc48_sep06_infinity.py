#!/usr/bin/env python3
"""Bounded exact controls for the nodal missing-page identity.

These are local passports and Euler-incidence controls, not realized global
Keller maps. No finite check establishes the all-degree geometric theorem.
"""
from fractions import Fraction as Q
from itertools import product
import hashlib
import json


GATES = 0


def check(value, label):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(label)


def poly_mul(p, q):
    r = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def poly_eval(p, x):
    ans = Q(0)
    for a in reversed(p):
        ans = ans * x + a
    return ans


def cyclic_missing(d, size, offset):
    return {(offset + k) % d for k in range(size)}


def main():
    rows = []
    pair_counts = {}
    for d in range(2, 6):
        omega = set(range(d))
        actual = [{k for k in range(d) if mask & (1 << k)}
                  for mask in range(1, (1 << d) - 1)]
        count = 0
        for a, b in product(actual, repeat=2):
            m, n = omega - a, omega - b
            delta_p = d - len(a & b)
            overlap = len(m & n)
            check(delta_p == len(m | n), "node complement identity")
            check(delta_p == len(m) + len(n) - overlap, "node union count")
            check(overlap >= max(0, len(m) + len(n) - d), "overlap bound")
            rows.append((d, sorted(a), sorted(b), delta_p, overlap))
            count += 1
        pair_counts[d] = count

    # Each listed edge is an ordinary node, loops are self-nodes. No claim
    # that each incidence graph is realizable by a polynomial target curve.
    graphs = [(1, []), (1, [(0, 0)]), (1, [(0, 0)] * 3),
              (2, [(0, 0), (1, 1)]), (2, [(0, 1)]),
              (2, [(0, 0), (0, 1), (1, 1)]),
              (3, [(0, 1), (1, 2), (2, 0), (0, 0)])]
    graph_rows = 0
    for d in range(2, 8):
        for components, edges in graphs:
            for debt in product(range(1, d), repeat=components):
                branches = [0] * components
                local_debt = []
                overlaps = []
                for k, (i, j) in enumerate(edges):
                    branches[i] += 1
                    branches[j] += 1
                    m = cyclic_missing(d, debt[i], 2 * k)
                    n = cyclic_missing(d, debt[j], 2 * k + 1)
                    local_debt.append(len(m | n))
                    overlaps.append(len(m & n))
                stratified = sum(z * (1 - b) for z, b in zip(debt, branches))
                stratified += sum(local_debt)
                cancelled = sum(debt) - sum(overlaps)
                check(stratified == cancelled, "Euler branch cancellation")
                graph_rows += 1

    # Necessary one-component passports: generic actual degree is >=1.
    survivors = []
    for d in range(2, 18):
        for delta in range(1, d):
            for overlap_sum in range(delta + 1):
                if d - 1 == delta - overlap_sum:
                    check(delta == d - 1 and overlap_sum == 0,
                          "one-component forcing")
                    if max(0, 2 * delta - d) == 0:
                        survivors.append((d, delta, overlap_sum))
    check(survivors == [(2, 1, 0)], "only degree-two formal passport")
    check([p for p in product(range(1, 4), repeat=2) if p[0] * p[1] == 1]
          == [(1, 1)], "positive missing-length factors")

    # The inertia-fixed hostile: neither actual subset is the full fixed set.
    universe, a, b = {0, 1}, {0}, {1}
    check(len(a & b) == 0 and len(universe & universe) == 2,
          "deleted fixed sheets cannot be counted as actual")
    for nodes in range(1, 5):
        check((2 - 1) * (1 - 2 * nodes) + nodes * 2 == 1,
              "degree-two Euler hostile survives every tested node count")

    # Intrinsic (2,3) normal form, with raw polynomial coefficient laws.
    cubic_rows = []
    for c in map(Q, [-3, -1, 0, 1, 2, 5]):
        u = [Q(0), Q(0), Q(1)]
        v = [Q(0), c, Q(0), Q(1)]
        u_plus_c = [c, Q(0), Q(1)]
        check(poly_mul(v, v) == poly_mul(u, poly_mul(u_plus_c, u_plus_c)),
              "V squared equals U times (U+c) squared")
        cubic_rows.append((str(c), [str(z) for z in v]))
    for t in map(Q, [1, 2, 3, Q(1, 2)]):
        c = -t * t
        plus, minus = (t * t, t**3 + c * t), (t * t, -t**3 - c * t)
        check(plus == minus == (-c, Q(0)), "two normalization preimages")
        determinant = (2 * t) * (-2 * c) - (-2 * t) * (-2 * c)
        check(determinant == -8 * c * t and determinant != 0,
              "ordinary-node tangent determinant")
    check(poly_eval([Q(0), Q(0), Q(0), Q(1)], 0) == 0,
          "cusp parameter endpoint retained")

    digest = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
    print("STATUS: FINITE-EXACT LOCAL PASSPORT CONTROLS; geometric proof is analytic")
    print("actual-subset pairs by d:", pair_counts)
    print("raw-labeled pair digest:", digest)
    print("declared Euler graphs:", len(graphs), "exact debt rows:", graph_rows)
    print("one-component finite passport survivors:", survivors)
    print("hostiles: deleted fixed sheets; degree-two Euler survivor")
    print("intrinsic cubic parameter rows:", cubic_rows)
    print("ALL GATES PASS:", GATES)


if __name__ == "__main__":
    main()
