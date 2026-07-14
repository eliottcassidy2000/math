#!/usr/bin/env python3
"""Exact audit of HYP-6815's affine-slope threshold suspension.

For V(c)=cP+R, the identity

    min_i ||(c p_i+r_i)t|| = min_i ||p_i u+r_i t||,  u={ct},

turns a scale family into integer-slope slices of one two-torus field.  This
script verifies the identity at every exact peak in a small adversarial bank,
checks the pure-dilation cylinder and the HYP-6780 near-dilate ray, and audits
which quotients preserve covering, M, safe measure, components, cap status,
and low-support relation data.

Tournament Analysis uses representations, not runners, as vertices.  Its
pairwise observable is the fraction of proof-critical row pairs separated by
the representation.  Two gauges trade predicate preservation against
compression; the declared tie path is printed with the fingerprints.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations
from math import gcd

from lrc14_certificates import (
    LAM,
    L_exact,
    M_exact,
    PI_LO,
    good_intervals,
    is_covering,
)


def circle_dist(x):
    x %= 1
    return min(x, 1 - x)


def affine_body(P, R, c):
    return tuple(c * p + r for p, r in zip(P, R))


def affine_phi(P, R, u, t):
    return min(circle_dist(p * u + r * t) for p, r in zip(P, R))


def exact_peak(body):
    """Return the earliest exact peak and its clearance."""
    values = sorted(set(map(abs, body)))
    denoms = {2 * v for v in values}
    for i, v in enumerate(values):
        for w in values[i + 1 :]:
            denoms.add(v + w)
            denoms.add(w - v)
    best_t, best = F(0), F(0)
    for q in sorted(denoms):
        for m in range(1, q):
            t = F(m, q)
            value = min(circle_dist(v * t) for v in values)
            if value > best or (value == best and t < best_t):
                best_t, best = t, value
    assert best == M_exact(body)
    return best_t, best


def primitive(body):
    g = 0
    for v in body:
        g = gcd(g, v)
    return g == 1


def cap_data(body):
    top = max(body)
    core = list(body)
    core.remove(top)
    ivs = good_intervals(core)
    G = sum((b - a for a, b in ivs), F(0))
    r = len(ivs)
    cap = bool(G and PI_LO * top * G > r)
    ratio = F(top) * G / r if r else F(0)
    return cap, ratio


def sum_relation_signature(P, R, c=None):
    """Projected support-3 relations p_i+p_j-p_k as (a.P,a.R,owners)."""
    out = []
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            for k in range(len(P)):
                if k == i or k == j:
                    continue
                m = P[i] + P[j] - P[k]
                n = R[i] + R[j] - R[k]
                if c is None or c * m + n == 0:
                    out.append((m, n, i + 1, j + 1, k + 1))
    return tuple(out)


def fiber_mixing(rows, key, field):
    fibers = defaultdict(list)
    for row in rows:
        fibers[key(row)].append(row)
    mixed = sum(len({row[field] for row in fiber}) > 1 for fiber in fibers.values())
    return len(fibers), mixed


def separation_fraction(rows, key, field):
    total = separated = 0
    for i, a in enumerate(rows):
        for b in rows[i + 1 :]:
            if a[field] == b[field]:
                continue
            total += 1
            separated += key(a) != key(b)
    return F(separated, total) if total else F(1)


def orient(scores, tie_path):
    rank = {v: i for i, v in enumerate(tie_path)}
    vertices = list(scores)
    edge = {}
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            if scores[a] > scores[b] or (scores[a] == scores[b] and rank[a] < rank[b]):
                edge[a, b], edge[b, a] = 1, 0
            else:
                edge[a, b], edge[b, a] = 0, 1
    return vertices, edge


def tournament_fingerprint(scores, tie_path):
    vertices, edge = orient(scores, tie_path)
    wins = {a: sum(edge[a, b] for b in vertices if b != a) for a in vertices}
    cycles = 0
    for i, a in enumerate(vertices):
        for j in range(i + 1, len(vertices)):
            b = vertices[j]
            for k in range(j + 1, len(vertices)):
                c = vertices[k]
                cycles += (edge[a, b] and edge[b, c] and edge[c, a]) or (
                    edge[b, a] and edge[c, b] and edge[a, c]
                )

    # Mutual reachability gives SCCs; n=7 keeps the direct closure simplest.
    reach = {(a, b): (a == b or bool(edge.get((a, b), 0))) for a in vertices for b in vertices}
    for k in vertices:
        for a in vertices:
            for b in vertices:
                reach[a, b] = reach[a, b] or (reach[a, k] and reach[k, b])
    unseen = set(vertices)
    sccs = []
    while unseen:
        a = min(unseen)
        comp = sorted(b for b in unseen if reach[a, b] and reach[b, a])
        sccs.append(comp)
        unseen.difference_update(comp)

    hpaths = sum(
        all(edge[path[i], path[i + 1]] for i in range(len(path) - 1))
        for path in permutations(vertices)
    )
    return edge, {
        "score_hist": dict(sorted(Counter(wins.values()).items())),
        "directed_3cycles": int(cycles),
        "SCCs": sccs,
        "Hamiltonian_paths": hpaths,
    }


def main():
    P13 = tuple(range(1, 14))
    zero13 = (0,) * 13

    print("HYP-6815 AFFINE-SLOPE THRESHOLD SUSPENSION")
    print("\nExact slice identity on the pure cylinder:")
    for c in (1, 2, 3, 5):
        body = affine_body(P13, zero13, c)
        t, M = exact_peak(body)
        phi = affine_phi(P13, zero13, (c * t) % 1, t)
        assert phi == M == F(1, 14)
        print(f"  c={c}: t*={t}, Phi(ct*,t*)=M={M}, L={L_exact(body)}")

    print("\nHYP-6780 positive-measure core cylinder:")
    P12 = tuple(range(1, 13))
    base_ivs = good_intervals(P12)
    base_G = sum((b - a for a, b in base_ivs), F(0))
    base_r = len(base_ivs)
    for c in (1, 2, 3, 5):
        ivs = good_intervals(tuple(c * p for p in P12))
        G = sum((b - a for a, b in ivs), F(0))
        assert G == base_G and len(ivs) == c * base_r
        print(f"  c={c}: |G|={G}, components={len(ivs)}={c}*{base_r}")

    print("\nTransverse-shear ray R=e_13:")
    R13 = (0,) * 12 + (1,)
    for c in (26, 91):
        body = affine_body(P13, R13, c)
        t, M = exact_peak(body)
        phi = affine_phi(P13, R13, (c * t) % 1, t)
        assert primitive(body) and is_covering(body) and phi == M == F(1, 13)
        print(f"  c={c}: max={max(body)}, t*={t}, slice M={M}, covering=True")

    print("\nDual projected relation audit (a=e_i+e_j-e_k):")
    shape_counts = set()
    killed = []
    for owner in (0, 6, 12):
        R = tuple(1 if i == owner else 0 for i in range(13))
        projected = sum_relation_signature(P13, R)
        shape = [x for x in projected if x[0] == 0]
        persistent = [x for x in shape if x[1] == 0]
        actual = sum_relation_signature(P13, R, c=2)
        assert actual == tuple(persistent)
        shape_counts.add(len(shape))
        killed.append(len(shape) - len(persistent))
        print(
            f"  offset owner={owner + 1:2d}: shape-relations={len(shape)}, "
            f"survive shear={len(persistent)}, killed={len(shape)-len(persistent)}"
        )
    assert shape_counts == {36} and min(killed) == 6 and max(killed) == 11

    rows = []
    for owner in (0, 6, 12):
        R = tuple(1 if i == owner else 0 for i in range(13))
        for c in (2, 16):  # identical c mod 14, very different endpoint height.
            body = affine_body(P13, R, c)
            assert primitive(body) and len(set(body)) == 13
            t, M = exact_peak(body)
            u = (c * t) % 1
            assert affine_phi(P13, R, u, t) == M
            ivs = good_intervals(body)
            cap, cap_ratio = cap_data(body)
            rows.append(
                {
                    "P": P13,
                    "R": R,
                    "c": c,
                    "covering": is_covering(body),
                    "M": M,
                    "peak": (t, M),
                    "L": sum((b - a for a, b in ivs), F(0)),
                    "components": len(ivs),
                    "cap": cap,
                    "cap_ratio": cap_ratio,
                    "relations": sum_relation_signature(P13, R, c),
                }
            )

    carriers = {
        "shape_only": lambda x: x["P"],
        "offset_multiset": lambda x: (x["P"], tuple(sorted(x["R"]))),
        "owner_resolved_offset": lambda x: (x["P"], x["R"]),
        "dual_relation_signature": lambda x: x["relations"],
        "owner_offset_cmod14": lambda x: (x["P"], x["R"], x["c"] % 14),
        "exact_peak_certificate": lambda x: x["peak"],
        "full_affine_slope_packet": lambda x: (x["P"], x["R"], x["c"]),
    }
    fields = ("covering", "M", "L", "components", "cap", "relations")

    print("\nExact mixed-fiber audit (6 adversarial rows; c=2 and c=16):")
    print("  carrier                    fibers " + " ".join(f"{f:>10s}" for f in fields))
    for name, key in carriers.items():
        fiber_count = len({key(row) for row in rows})
        mixed = [fiber_mixing(rows, key, field)[1] for field in fields]
        print(f"  {name:26s} {fiber_count:6d} " + " ".join(f"{x:10d}" for x in mixed))

    first = rows[0]
    mate = rows[1]
    assert first["R"] == mate["R"] and first["c"] % 14 == mate["c"] % 14
    assert (first["L"], first["components"]) != (mate["L"], mate["components"])
    print("\nSame shape + owned offset + c mod 14, different metric fiber:")
    print(
        f"  owner=1,c={first['c']}: M={first['M']}, L={first['L']}, "
        f"components={first['components']}, covering={first['covering']}"
    )
    print(
        f"  owner=1,c={mate['c']}: M={mate['M']}, L={mate['L']}, "
        f"components={mate['components']}, covering={mate['covering']}"
    )

    metrics = {
        name: {field: separation_fraction(rows, key, field) for field in fields}
        for name, key in carriers.items()
    }
    compression = {
        name: F(len(rows) - len({key(row) for row in rows}), len(rows) - 1)
        for name, key in carriers.items()
    }
    predicate_scores = {
        name: 8 * sum(axis.values(), F(0)) + compression[name]
        for name, axis in metrics.items()
    }
    compression_scores = {
        name: sum(axis.values(), F(0)) + 8 * compression[name]
        for name, axis in metrics.items()
    }
    tie_path = (
        "full_affine_slope_packet",
        "exact_peak_certificate",
        "owner_offset_cmod14",
        "dual_relation_signature",
        "owner_resolved_offset",
        "offset_multiset",
        "shape_only",
    )
    edge_a, fp_a = tournament_fingerprint(predicate_scores, tie_path)
    edge_b, fp_b = tournament_fingerprint(compression_scores, tie_path)
    flips = sum(
        edge_a[a, b] != edge_b[a, b]
        for i, a in enumerate(carriers)
        for b in list(carriers)[i + 1 :]
    )

    print("\nTournament Analysis (vertices are representations, not runners):")
    print("  observable=proof-critical row pairs separated")
    print("  gauges=predicate-first and compression-first")
    print(f"  declared tie Hamiltonian path={tie_path}")
    print(f"  predicate-first: {fp_a}")
    print(f"  compression-first: {fp_b}")
    print(f"  edge_flips_between_gauges={flips}")
    print("\nVerdict: the slope packet is exact; mod-14 residue and relation shadows are not metric-complete.")


if __name__ == "__main__":
    main()
