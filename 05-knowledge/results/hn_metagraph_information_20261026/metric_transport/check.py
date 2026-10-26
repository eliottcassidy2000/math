"""Exact metric transport tests; standard library only, no floating point."""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations, product, combinations
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


PERMS = tuple(permutations(range(3)))
W = tuple((p, a) for p in PERMS for a in range(8))


def act(g, x):
    p, a = g
    return a ^ sum(((x >> i) & 1) << p[i] for i in range(3))


def parity(x):
    return x.bit_count() % 2


def vscale(c, u):
    return tuple(c * x for x in u)


def vsub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def dot(u, v, cross):
    return u[0] * v[0] + u[1] * v[1] + cross * (u[0] * v[1] + u[1] * v[0])


def marked(u, mask):
    return tuple(vscale(-1 if (mask >> i) & 1 else 1, z) for i, z in enumerate(u))


def kernel_pair(u):
    c = (det(u[1], u[2]), det(u[2], u[0]), det(u[0], u[1]))
    need(all(c), 'six distinct spokes / no parallel pair')
    need(all(sum(c[i] * u[i][j] for i in range(3)) == 0 for j in range(2)), 'kernel')
    b = sum((int(z < 0) << i) for i, z in enumerate(c))
    return frozenset((b, b ^ 7))


def triangle_pairs(u):
    # Independent criterion: the three oriented edge determinants have one
    # strict sign exactly when the origin lies inside the chosen triangle.
    inside = []
    for x in range(8):
        z = marked(u, x)
        areas = (det(z[0], z[1]), det(z[1], z[2]), det(z[2], z[0]))
        if all(t > 0 for t in areas) or all(t < 0 for t in areas):
            inside.append(x)
    return frozenset(inside)


def gram_transport(u, v, cross, g):
    p, a = g
    return all(dot(u[i], u[j], cross) ==
               (-1 if ((a >> p[i]) ^ (a >> p[j])) & 1 else 1) *
               dot(v[p[i]], v[p[j]], cross)
               for i in range(3) for j in range(3))


def distance_transport(u, v, cross, g):
    # Second path uses all six points, without the Gram-preservation formula.
    p, a = g
    src = [(i, b) for i in range(3) for b in range(2)]
    def point(frame, i, b):
        return vscale((-1)**b, frame[i])
    def target(i, b):
        return point(v, p[i], b ^ ((a >> p[i]) & 1))
    for (i, b), (j, c) in combinations(src, 2):
        ds = vsub(point(u, i, b), point(u, j, c))
        dt = vsub(target(i, b), target(j, c))
        if dot(ds, ds, cross) != dot(dt, dt, cross):
            return False
    return True


def fixed_composite(gs):
    ans = []
    for x in range(8):
        y = x
        for g in gs:
            y = act(g, y)
        if x == y:
            ans.append(x)
    return ans


def main():
    regular = ((Q(1), Q(0)), (Q(-1), Q(1)), (Q(0), Q(-1)))
    generic = ((Q(1), Q(0)), (Q(3, 5), Q(4, 5)), (Q(5, 13), Q(12, 13)))
    reflection = ((Q(1), Q(0)), (Q(3, 5), Q(4, 5)), (Q(3, 5), Q(-4, 5)))
    controls = {}
    for name, u, cross, expected in [('regular', regular, Q(1, 2), 12),
                                     ('generic', generic, Q(0), 2),
                                     ('reflection', reflection, Q(0), 4)]:
        need(all(dot(z, z, cross) == 1 for z in u), 'unit spokes')
        frames = [marked(u, a) for a in range(8)]
        accepted = 0
        for src in frames:
            need(kernel_pair(src) == triangle_pairs(src), 'two independent section constructions')
            for dst in frames:
                admissible = []
                for g in W:
                    metric = gram_transport(src, dst, cross, g)
                    need(metric == distance_transport(src, dst, cross, g), 'Gram versus all distances')
                    if metric:
                        admissible.append(g)
                        need(frozenset(act(g, x) for x in kernel_pair(src)) == kernel_pair(dst),
                             'canonical pair covariance')
                need(len(admissible) == expected, 'metric transporter size')
                accepted += len(admissible)
        sums = {tuple(sum(((x >> i) & 1) * u[i][j] for i in range(3)) for j in range(2))
                for x in range(8)}
        controls[name] = {'self_isometries': expected, 'signed_markings': 8,
                          'ordered_marking_pairs': 64, 'W_tests': 3072,
                          'accepted_transports': accepted,
                          'balanced_frames': sorted(kernel_pair(u)),
                          'subset_sum_images': len(sums)}
    need(controls['regular']['subset_sum_images'] == 7, 'flower positive control')
    need(controls['generic']['subset_sum_images'] == 8, 'theorem extends beyond collision')

    # Exhaust all regular-flower collision types and compatible seam maps on
    # one triangle: 4^3 * 12^3 = 110592 systems. Full section counts are read
    # from the composed eight-state maps, independently of parity prediction.
    labels = (0, 3, 5, 6)
    frames = {a: marked(regular, a) for a in labels}
    seams = {(a, b): tuple(g for g in W if gram_transport(frames[a], frames[b], Q(1, 2), g))
             for a in labels for b in labels}
    histogram = Counter()
    for a, b, c in product(labels, repeat=3):
        for gs in product(seams[a, b], seams[b, c], seams[c, a]):
            full = fixed_composite(gs)
            parity_pass = sum(parity(g[1]) for g in gs) % 2 == 0
            need(bool(full) == parity_pass, 'parity iff complete section')
            if parity_pass:
                need(kernel_pair(frames[a]) <= set(full), 'distinguished pair survives every loop')
            histogram[len(full)] += 1
    need(dict(histogram) == {8: 9216, 0: 55296, 4: 27648, 2: 18432}, 'triangle census')

    cycle = (1, 2, 0)
    g0, g5 = (cycle, 0), (cycle, 5)
    fix0, fix5 = fixed_composite([g0]), fixed_composite([g5])
    need(fix0 == [0, 7] and fix5 == [3, 4] and not set(fix0) & set(fix5), 'abstract hostile')
    common_metric = []
    for a in labels:
        if all(gram_transport(frames[a], frames[a], Q(1, 2), g) for g in (g0, g5)):
            common_metric.append(a)
    need(not common_metric, 'no one flower realizes both hostile holonomies')
    need(gram_transport(frames[0], frames[5], Q(1, 2), g5),
         'nontrivial marking transport is possible between different flowers')
    # In any rank-two frame, a common metric realization would have to preserve
    # its balanced pair. The hostile pair has no common invariant antipodal pair.
    antipodes = [frozenset((a, a ^ 7)) for a in range(4)]
    need(not any(all(frozenset(act(g, x) for x in pair) == pair for g in (g0, g5))
                 for pair in antipodes), 'all planar three-direction frames excluded')
    receipt = {'status': 'PASS', 'arithmetic': 'exact Fraction', 'controls': controls,
               'triangle_systems': sum(histogram.values()),
               'triangle_full_section_histogram': dict(sorted(histogram.items())),
               'hostile': {'holonomy_0_fixed': fix0, 'holonomy_5_fixed': fix5,
                           'common_realizations': common_metric},
               'scope': 'isometric transport of all six paired spokes; not arbitrary unit-distance gadgets'}
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
