"""Exact deficiency-two grammar, with complete independent controls at k<=2."""
from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import factorial
from pathlib import Path


def need(value, message):
    if not value:
        raise RuntimeError(message)


spec = spec_from_file_location(
    "strata_sparse",
    Path(__file__).resolve().parents[1] / "defect_boundary_20261020/sparse_count.py",
)
S = module_from_spec(spec)
spec.loader.exec_module(S)


def template(defect, endpoint_weights, gaps, k):
    t = len(endpoint_weights)
    slots, weights = {}, []
    for j in range(t + 1):
        if gaps[j]:
            weights.append(gaps[j])
        if j < t:
            slots[j] = len(weights)
            weights.append(endpoint_weights[j])
    n = len(weights)
    rows = [(1 << n) - (1 << (j + 1)) for j in range(n)]
    for a, b in defect:
        a, b = sorted((slots[a], slots[b]))
        rows[a] ^= 1 << b
        rows[b] ^= 1 << a
    return S.normalize(tuple(rows), weights, k)


def polynomial_coefficients(states, k):
    out = Counter()
    for (shift, h), count in S.terms(states, k).items():
        if h == 0:
            continue
        poly = [Fraction(1)]
        for j in range(1, h):
            nxt = [Fraction(0)] * (len(poly) + 1)
            for degree, value in enumerate(poly):
                nxt[degree] += (j - shift) * value
                nxt[degree + 1] += value
            poly = nxt
        for degree, value in enumerate(poly):
            out[degree] += count * value / factorial(h - 1)
    return [out[degree] for degree in (2 * k, 2 * k - 1, 2 * k - 2)]


def classify(k):
    K = k + 1
    states, provenance = set(), {}
    raw, losses, motifs = Counter(), Counter(), Counter()
    for r in range(max(0, k - 1), k + 1):
        for t in range(max(0, 2 * k - 2), 2 * r + 1):
            if (r == 0) != (t == 0):
                continue
            edges = tuple(combinations(range(t), 2))
            for defect in combinations(edges, r):
                if set(v for edge in defect for v in edge) != set(range(t)):
                    continue
                delta = 2 * r - t
                gap_loss_max = 2 - 2 * (k - r) - delta
                if gap_loss_max < 0:
                    continue
                weights = [(1,) * t]
                if k - r == 1:
                    weights += [tuple(1 + int(v == j) for v in range(t))
                                for j in range(t)]
                for ew in weights:
                    cost = sum(ew[a] * ew[b] for a, b in defect)
                    if cost > k:
                        continue
                    for u in range(gap_loss_max + 1):
                        for positions in combinations(range(t + 1), u):
                            for values in product(range(K), repeat=u):
                                gaps = [K] * (t + 1)
                                for j, value in zip(positions, values):
                                    gaps[j] = value
                                source_d = 2 * (k - r) + delta + u
                                need(source_d == 2 * k + 1 - gaps.count(K),
                                     "deficiency identity")
                                state = template(defect, ew, gaps, k)
                                final_d = 2 * k + 1 - state[2].count(K)
                                raw[source_d] += 1
                                losses[source_d, final_d] += 1
                                if final_d <= 2:
                                    states.add(state)
                                    provenance.setdefault(state, set()).add(
                                        (k - r, delta, u))
                                    if delta == 2:
                                        degrees = [sum(v in edge for edge in defect)
                                                   for v in range(t)]
                                        motifs[tuple(sorted(d for d in degrees if d > 1))] += 1
    histogram = Counter(state[2].count(K) for state in states)
    if k <= 2:
        complete, _ = S.template_states(k)
        expected = {state for state in complete if state[2].count(K) >= 2 * k - 1}
        need(states == expected, "independent complete compiler equality")
    coefficients = polynomial_coefficients(states, k)
    scale = 2 ** k * factorial(k)
    need(coefficients[0] * scale == 1, "inherited leading coefficient")
    need(coefficients[1] * scale == -3 * k, "inherited second coefficient")
    M = factorial(2 * k) // (2 ** k * factorial(k))
    P = M // (2 * k - 1)
    N2 = (M * (2 * k ** 4 + 7 * k ** 3 + 2 * k ** 2 - 3 * k + 1)
          + Fraction(P * (3 * k ** 3 - 16 * k ** 2 + 13), 3))
    need(histogram[2 * k - 1] == N2, "proved third-stratum formula")
    need(coefficients[2] * scale == Fraction(k * (17 - 9 * k), 2),
         "proved third polynomial coefficient")
    by_core_order = Counter(state[0] for state in states
                            if state[2].count(K) == 2 * k - 1)
    expected_by_core = {
        4 * k + 1: k ** 3 * (2 * k + 1) * M,
        4 * k: 4 * k ** 3 * M,
        4 * k - 1: (2 * k ** 3 - k) * M - 3 * (k - 1) * P,
        4 * k - 2: 2 * k * (k - 1) * M - 2 * (k ** 2 - 1) * P,
        4 * k - 3: Fraction((k - 1) * (k - 2) * (3 * k - 1) * P, 3) + M,
    }
    need(by_core_order == Counter({q: v for q, v in expected_by_core.items() if v}),
         "each proved collision sector by canonical core order")
    multi = sum(len(p) > 1 for p in provenance.values())
    print("k=", k, "stratum_counts=", dict(sorted(histogram.items())),
          "raw_by_source_deficiency=", dict(sorted(raw.items())), flush=True)
    print("normalization_deficiency_transitions=", dict(sorted(losses.items())),
          "cross_signature_duplicates=", multi, flush=True)
    print("collision_degree_profiles=", dict(sorted(motifs.items())),
          "top_three_scaled_coefficients=", [str(c * scale) for c in coefficients],
          flush=True)
    print("deficiency_two_by_core_order=", dict(sorted(by_core_order.items())), flush=True)
    return coefficients


def main():
    # Sharp clone split/merge: a single flip in the five-vertex total order.
    tt = S.E.matrix(5, (1 << 10) - 1)
    child = list(tt)
    child[1] ^= 1 << 3
    child[3] ^= 1 << 1
    need(len(S.G.clone_chains(tt)) == 1, "transitive control")
    need(len(S.G.clone_chains(tuple(child))) == 5, "sharp four-core jump")
    print("Sharp |b(T^e)-b(T)|=4 witness: TT5, flip positions (1,3).", flush=True)
    for k in range(1, 5):
        classify(k)
    print("PASS: complete grammar controls k=1,2; deficiency<=2 generation k=1..4.")


if __name__ == "__main__":
    main()
