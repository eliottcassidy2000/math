"""Exact marked fivefold and Petersen boundary interfaces; no discovery search.

Run: python -B -X utf8 04-computation/experiments/fivefold_boundary_20261018g.py
JSON goes to stdout; this script never modifies the repository.
"""
from fractions import Fraction as F
from itertools import product, permutations
import json


CHECKS = 0


def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ValueError(message)


def exact_tuple(x, size):
    if type(x) is not tuple or len(x) != size or any(type(v) not in (int, F) for v in x):
        raise TypeError("an exact, source-sized rational tuple is required")


def character(d, g):
    if type(d) is not int or type(g) is not int or not 0 <= d < 8 or not 0 <= g < 8:
        raise TypeError("marked three-bit character and group address required")
    return 1 if (d & g).bit_count() % 2 == 0 else -1


def golden_observer(weights):
    exact_tuple(weights, 8)
    return tuple(sum(weights[d] * character(d, g) for d in range(8)) for g in range(8))


def golden_defects(weights):
    exact_tuple(weights, 8)
    return tuple(weights[d] - weights[d ^ 2] for d in (0, 1, 4, 5))


def edge(a, b):
    return (a, b) if a < b else (b, a)


P_EDGES = tuple(sorted({edge(i, (i + 1) % 5) for i in range(5)} |
                      {edge(5 + i, 5 + (i + 2) % 5) for i in range(5)} |
                      {(i, i + 5) for i in range(5)}))
MATCHING = tuple((i, i + 5) for i in range(5))
Q_EDGES = tuple((i, j) for i in range(5) for j in range(i + 1, 5))


def divergence(vertices, edges, flow):
    out = [0] * vertices
    for (a, b), x in zip(edges, flow):
        out[a] += x
        out[b] -= x
    return tuple(out)


def external_lift(flow):
    """Orientation-preserving pullback to the ten nonmatching source edges."""
    exact_tuple(flow, 10)
    out = []
    for a, b in P_EDGES:
        if (a, b) in MATCHING:
            out.append(0)
        else:
            i, j = a % 5, b % 5
            out.append((1 if i < j else -1) * flow[Q_EDGES.index(edge(i, j))])
    return tuple(out)


def lift_cycle(flow):
    """An exact K5 circulation -> its unique integral/rational Petersen lift."""
    exact_tuple(flow, 10)
    if any(divergence(5, Q_EDGES, flow)):
        raise ValueError("quotient flow is not a circulation")
    out = list(external_lift(flow))
    div = divergence(10, P_EDGES, out)
    for i, e in enumerate(MATCHING):
        out[P_EDGES.index(e)] = -div[i]
    if any(divergence(10, P_EDGES, out)):
        raise ArithmeticError("source conservation failed")
    return tuple(out)


def restrict_cycle(flow):
    exact_tuple(flow, 15)
    if any(divergence(10, P_EDGES, flow)):
        raise ValueError("source flow is not a circulation")
    out = [0] * 10
    for e, x in zip(P_EDGES, flow):
        if e not in MATCHING:
            a, b = e[0] % 5, e[1] % 5
            out[Q_EDGES.index(edge(a, b))] = (1 if a < b else -1) * x
    return tuple(out)


def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for i in range(n):
        p = next(j for j in range(i, n) if a[j][i])
        a[i], a[p] = a[p], a[i]
        t = a[i][i]
        a[i] = [x / t for x in a[i]]
        for j in range(n):
            if j != i:
                t = a[j][i]
                a[j] = [x - t * y for x, y in zip(a[j], a[i])]
    return [row[n:] for row in a]


def laplacian(n, edges):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for i, j in edges:
        a[i][i] += 1
        a[j][j] += 1
        a[i][j] -= 1
        a[j][i] -= 1
    return a


def effective_resistance(n, edges, a, b):
    lap = laplacian(n, edges)
    retained = [v for v in range(n) if v != b]
    inv = inverse([[lap[i][j] for j in retained] for i in retained])
    return inv[retained.index(a)][retained.index(a)]


def pole(cycle):
    out = {}
    for colors in product(range(3), repeat=5):
        boundary = [None] * 5
        for i, v in enumerate(cycle):
            x, y = colors[(i - 1) % 5], colors[i]
            if x == y:
                break
            boundary[v] = 3 - x - y
        else:
            key = tuple(boundary)
            out[key] = out.get(key, 0) + 1
    return out


# Q[zeta5], exact coefficient order 1,z,z^2,z^3.
ONE = (1, 0, 0, 0)
ZETA = (0, 1, 0, 0)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, c):
    return tuple(c * x for x in a)


def mul(a, b):
    c = [0] * 7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    for k in range(6, 3, -1):
        x = c[k]
        for j in range(4):
            c[k - 4 + j] -= x
    return tuple(c[:4])


def power(a, n):
    out = ONE
    for _ in range(n):
        out = mul(out, a)
    return out


def conjugate_at(a, k):
    out = (0, 0, 0, 0)
    for i, x in enumerate(a):
        out = add(out, scale(power(ZETA, k * i), x))
    return out


def norm(a):
    return mul(a, conjugate_at(a, 4))


def a4_image(a):
    exact_tuple(a, 5)
    if any(type(x) is not int for x in a) or sum(a) != 0:
        raise ValueError("source must be an integral A4 vector")
    out = (0, 0, 0, 0)
    for i, x in enumerate(a):
        out = add(out, scale(power(ZETA, i), x))
    return out


def a4_decode(a):
    exact_tuple(a, 4)
    if any(type(x) is not int for x in a) or sum(a) % 5:
        raise ValueError("not in the A4 image ideal")
    k = -sum(a) // 5
    return tuple(x + k for x in a) + (k,)


def serial(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): serial(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [serial(v) for v in x]
    return x


def verify():
    global CHECKS
    CHECKS = 0
    # Exhaust the entire {0,1,2}^8 weighted universe, not just equal weights.
    for w in product(range(3), repeat=8):
        obs = golden_observer(w)
        for g in range(8):
            require(obs[g] == sum(character(d, g) * (w[d] + character(2, g) * w[d ^ 2])
                                 for d in (0, 1, 4, 5)), "paired weight factorization")
        require(all(obs[g] == 0 for g in range(8) if g & 2) == (not any(golden_defects(w))),
                "golden-half cancellation iff")
    hostile = (2, 0, 1, 0, 0, 0, 0, 0)
    require(golden_observer(hostile)[2] == 1, "unequal pair survives")
    for bad in ((True,) + (0,) * 7, (0.0,) + (0,) * 7, [0] * 8):
        try:
            golden_observer(bad)
        except TypeError:
            require(True, "reject aliased weights")
        else:
            require(False, "weight type alias accepted")

    quotient = {edge(a % 5, b % 5) for a, b in P_EDGES if (a, b) not in MATCHING}
    require(quotient == set(Q_EDGES), "actual matching quotient K5")
    basis = []
    for i in range(1, 5):
        for j in range(i + 1, 5):
            x = [0] * 10
            x[Q_EDGES.index((0, i))] = 1
            x[Q_EDGES.index((i, j))] = 1
            x[Q_EDGES.index((0, j))] = -1
            basis.append(tuple(x))
    lifts = [lift_cycle(x) for x in basis]
    for coeff in product((-1, 0, 1), repeat=6):
        x = tuple(sum(coeff[i] * basis[i][j] for i in range(6)) for j in range(10))
        y = lift_cycle(x)
        require(restrict_cycle(y) == x, "integral cycle-lattice inverse")
        require(y == tuple(sum(coeff[i] * lifts[i][j] for i in range(6)) for j in range(15)), "linear lift")
        require(sum(t * t for t in y) == sum(t * t for t in x) +
                sum(y[P_EDGES.index(e)] ** 2 for e in MATCHING), "removed-edge energy bill")
    gram = [[sum(a * b for a, b in zip(x, y)) for y in lifts] for x in lifts]
    qgram = [[sum(a * b for a, b in zip(x, y)) for y in basis] for x in basis]
    require(gram != qgram, "contraction not an energy isometry")
    try:
        lift_cycle((1,) + (0,) * 9)
    except ValueError:
        require(True, "reject unbalanced quotient flow")
    else:
        require(False, "unbalanced flow accepted")

    lap = laplacian(10, P_EDGES)
    innerinv = inverse([[lap[i][j] for j in range(5, 10)] for i in range(5, 10)])
    kron = [[lap[i][j] - innerinv[i][j] for j in range(5)] for i in range(5)]
    for v in product((-1, 0, 1), repeat=5):
        inner = [sum(innerinv[i][j] * v[j] for j in range(5)) for i in range(5)]
        whole = tuple(v) + tuple(inner)
        physical = sum((whole[i] - whole[j]) ** 2 for i, j in P_EDGES)
        reduced = sum(v[i] * kron[i][j] * v[j] for i in range(5) for j in range(5))
        require(physical == reduced, "Schur response preserves minimized energy")
    require(effective_resistance(10, P_EDGES, 0, 1) == F(3, 5), "Petersen unit resistance")
    require(effective_resistance(5, Q_EDGES, 0, 1) == F(2, 5), "contracted K5 resistance differs")

    # Suppress the three degree-two vertices, retaining each original path.
    paths = ((2, 1, 6), (9, 4, 3), (8, 5, 7))
    branches = (2, 3, 6, 7, 8, 9)
    rem = {e for e in P_EDGES if 0 not in e}
    direct = {e for e in rem if e[0] in branches and e[1] in branches}
    endpoints = {edge(p[0], p[-1]) for p in paths}
    require(direct | endpoints == {edge(i, j) for i in (2, 8, 9) for j in (3, 6, 7)}, "K33 path quotient")
    counts = {}
    for k in (2, 3):
        support = quotient_support = 0
        for row in product(range(k), repeat=6):
            c = dict(zip(branches, row))
            if not all(c[i] != c[j] for i, j in direct):
                continue
            path_ok = all(any(t != c[p[0]] and t != c[p[-1]] for t in range(k)) for p in paths)
            support += path_ok
            quotient_support += all(c[i] != c[j] for i, j in endpoints)
        counts[k] = (support, quotient_support)
    require(counts == {2: (0, 2), 3: (66, 42)}, "suppression changes exact boundary support")

    f = pole((0, 1, 2, 3, 4))
    histogram = {}
    for p in permutations((1, 2, 3, 4)):
        if p[0] > p[-1]:
            continue
        g = pole((0,) + p)
        dot = sum(v * g.get(k, 0) for k, v in f.items())
        histogram[dot] = histogram.get(dot, 0) + 1
    require(len(f) == 30 and set(f.values()) == {1}, "complete pentagon pole")
    require(histogram == {0: 1, 12: 5, 18: 5, 30: 1}, "all twelve cyclic gluing supports")

    require(power(ZETA, 5) == ONE and ZETA != ONE, "faithful order-five carrier")
    for xs in product(range(-2, 3), repeat=4):
        a = xs + (-sum(xs),)
        z = a4_image(a)
        require(a4_decode(z) == a, "literal A4 source decode")
        require(a4_image((a[-1],) + a[:-1]) == mul(ZETA, z), "cyclic rotation intertwines")
        require(add(norm(z), norm(conjugate_at(z, 2))) == scale(ONE, F(5, 2) * sum(x * x for x in a)),
                "four-dimensional Fourier norm identity")
    phi = add(add(ONE, ZETA), conjugate_at(ZETA, 4))
    phiinv = add(phi, scale(ONE, -1))
    require(mul(phi, phiinv) == ONE, "golden integral unit")
    require(add(phi, conjugate_at(phi, 2)) == ONE, "physical/internal golden conjugates")
    small = add(ONE, scale(ZETA, -1))
    height_rows = []
    for n in range(13):
        z = mul(small, power(phiinv, n))
        a = a4_decode(z)
        require(a4_image(a) == z and z != (0, 0, 0, 0), "nonzero shrinking projected lattice sequence")
        height_rows.append({"n": n, "a4_source": a, "physical_norm_squared": norm(z),
                            "internal_norm_squared": norm(conjugate_at(z, 2))})
    require(sum((F(p - 2, 2 * p) for p in (3, 8, 24)), F(0)) == 1, "flat hostile star")
    require(not any(all(x[i] != x[(i + 1) % 3] for i in range(3)) for x in product((8, 24), repeat=3)),
            "alternating triangle boundary cannot close")
    return serial({
        "status": "PROVED elementary maps; FINITE-EXACT declared checks; no new global conjecture result",
        "checks": CHECKS,
        "golden_character_order": [1, 3, 5, 15, 7, 21, 35, 105],
        "golden_pairs": [[1, 5], [3, 15], [7, 35], [21, 105]],
        "golden_weight_universe": "{0,1,2}^8", "unequal_weight_hostile": hostile,
        "petersen": {"edges": P_EDGES, "contracted_forest": MATCHING, "quotient_edges": Q_EDGES,
                     "quotient_cycle_basis": basis, "source_cycle_basis": lifts,
                     "source_energy_gram": gram, "quotient_energy_gram": qgram,
                     "kron_outer_laplacian": kron, "effective_resistances": {"Petersen_0_1": F(3, 5), "K5_0_1": F(2, 5)},
                     "suppressed_paths": paths, "path_support_vs_direct_support": counts,
                     "pentagon_gluing_Tait_histogram": histogram},
        "a4": {"coefficient_order": "1,zeta5,zeta5^2,zeta5^3", "shrinking_projection_sequence": height_rows},
        "scope": ["Petersen contractions preserve circulation only with a lift; not colouring support or unit energy",
                  "planar fivefold projection is injective but not discrete; this is not a Penrose tiling classifier",
                  "golden cancellation is weight-pair equality, not positivity or an unmarked count"]})


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
