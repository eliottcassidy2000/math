"""Exact deficiency-three ownership by intrinsic heavy gaps.

Run from repository root. Standard library only; checks survive python -O.
No fitted formula, Burnside count, or prior core bank is used by this grammar.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path


def need(value, message):
    if not value:
        raise RuntimeError(message)


spec = spec_from_file_location(
    "ownership_sparse",
    Path(__file__).resolve().parents[1] / "defect_boundary_20261020/sparse_count.py",
)
S = module_from_spec(spec)
spec.loader.exec_module(S)


def ordered_template(defect, endpoint_weights, gaps):
    slots, weights = {}, []
    for j in range(len(endpoint_weights) + 1):
        if gaps[j]:
            weights.append(gaps[j])
        if j < len(endpoint_weights):
            slots[j] = len(weights)
            weights.append(endpoint_weights[j])
    n = len(weights)
    rows = [(1 << n) - (1 << (j + 1)) for j in range(n)]
    edges = []
    for a, b in defect:
        a, b = sorted((slots[a], slots[b]))
        rows[a] ^= 1 << b
        rows[b] ^= 1 << a
        edges.append((a, b))
    return tuple(rows), tuple(weights), tuple(edges)


def interval_chains(n, edges):
    """Homogeneous pairs from symmetric differences of defect neighborhoods."""
    neighbors = [0] * n
    for i, j in edges:
        neighbors[i] |= 1 << j
        neighbors[j] |= 1 << i
    successor, predecessor = {}, {}
    for i, j in combinations(range(n), 2):
        outside_pair = (neighbors[i] ^ neighbors[j]) & ~((1 << i) | (1 << j))
        interval = (1 << j) - (1 << (i + 1))
        if outside_pair != interval:
            continue
        a, b = (j, i) if neighbors[i] >> j & 1 else (i, j)
        need(a not in successor and b not in predecessor, "clone branching")
        successor[a], predecessor[b] = b, a
    result = []
    for start in range(n):
        if start in predecessor:
            continue
        chain = [start]
        while chain[-1] in successor:
            chain.append(successor[chain[-1]])
        result.append(tuple(chain))
    need(sum(map(len, result)) == n, "clone cycle")
    return tuple(result)


def normalize_ordered(rows, weights, edges, k):
    chains = interval_chains(len(rows), edges)
    need(chains == S.G.clone_chains(rows), "independent clone normalization")
    core = S.G.induced(rows, [c[0] for c in chains])
    cap = tuple(min(k + 1, sum(weights[v] for v in c)) for c in chains)
    need(len(S.G.clone_chains(core)) == len(core), "pair-free normalization")
    return core, cap


def heavy_gaps(rows, weights, k):
    heavy = [v for v, w in enumerate(weights) if w == k + 1]
    heavy.sort(key=lambda v: -sum(rows[v] >> u & 1 for u in heavy))
    need(all(rows[u] >> v & 1 for u, v in combinations(heavy, 2)),
         "heavy order is transitive")
    gaps = [[] for _ in range(len(heavy) + 1)]
    for v, w in enumerate(weights):
        if w == k + 1:
            continue
        bits = tuple(rows[u] >> v & 1 for u in heavy)
        j = sum(bits)
        need(bits == (1,) * j + (0,) * (len(heavy) - j), "intrinsic light gap")
        gaps[j].append(v)
    need(all(gaps[j] for j in range(1, len(heavy))), "empty internal heavy gap")
    return tuple(heavy), tuple(tuple(g) for g in gaps)


def presentation(rows, weights, order, k):
    """Return the endpoint/gap record, or None if this is not a good order."""
    backward = tuple((i, j) for i, j in combinations(range(len(order)), 2)
                     if rows[order[j]] >> order[i] & 1)
    cost = sum(weights[order[i]] * weights[order[j]] for i, j in backward)
    if cost > k:
        return None
    endpoints = sorted({v for e in backward for v in e})
    position = {v: j for j, v in enumerate(endpoints)}
    ew = tuple(weights[order[v]] for v in endpoints)
    defects = tuple((position[i], position[j]) for i, j in backward)
    gap_weights = []
    for a, b in zip([-1] + endpoints, endpoints + [len(order)]):
        need(b - a <= 2, "multiple unmarked vertices in a canonical core gap")
        gap_weights.append(weights[order[a + 1]] if b == a + 2 else 0)
    r, t = len(defects), len(endpoints)
    q, delta = k - r, 2 * r - t
    u = sum(w != k + 1 for w in gap_weights)
    need(2 * q + delta + u == 2 * k + 1 - weights.count(k + 1),
         "intrinsic deficiency identity")
    return q, delta, u, ew, tuple(gap_weights), defects


@lru_cache(None)
def feasible_gap_orders(rows, weights, allowance):
    out = []
    for order in permutations(range(len(rows))):
        cost = sum(weights[a] * weights[b] for a, b in combinations(order, 2)
                   if rows[b] >> a & 1)
        if cost <= allowance:
            out.append((order, cost))
    return tuple(out)


@lru_cache(None)
def owner(rows, weights, k):
    """Lexicographically least good local presentation: an exact orbit key."""
    heavy, gaps = heavy_gaps(rows, weights, k)
    D = 2 * k + 1 - len(heavy)
    excess = sum(max(len(g) - 1, 0) for g in gaps)
    candidates = 1
    for gap in gaps:
        candidates *= factorial(len(gap))
    need(excess <= 2 * D, "intrinsic gap excess bound")
    need(candidates <= factorial(2 * D + 1), "bounded local ownership")
    best = None
    signatures, good_records = set(), set()
    cross_cost = sum(weights[a] * weights[b]
                     for i, j in combinations(range(len(gaps)), 2)
                     for a in gaps[i] for b in gaps[j] if rows[b] >> a & 1)
    need(cross_cost <= k, "intrinsic crossing cost exceeds budget")
    allowance = k - cross_cost
    choices = [feasible_gap_orders(S.G.induced(rows, g),
                                  tuple(weights[v] for v in g), allowance)
               for g in gaps]
    for options in product(*choices):
        if sum(cost for _, cost in options) > allowance:
            continue
        order = []
        for j, (part, _) in enumerate(options):
            order.extend(gaps[j][v] for v in part)
            if j < len(heavy):
                order.append(heavy[j])
        record = presentation(rows, weights, order, k)
        if record is not None:
            signatures.add(record[:3])
            good_records.add(record)
            if best is None or record < best:
                best = record
    need(best is not None, "no admissible owner")
    return best, excess, candidates, tuple(sorted(signatures)), len(good_records)


def menu_templates(k, maximum=3):
    """The exact D<=3 menu, generated by edge sets rather than a core bank."""
    for r in range(max(0, k - 1), k + 1):
        for t in range(2 * r + 1):
            if (r == 0) != (t == 0):
                continue
            q, delta = k - r, 2 * r - t
            remaining = maximum - 2 * q - delta
            if remaining < 0:
                continue
            for defect in combinations(tuple(combinations(range(t), 2)), r):
                if {v for e in defect for v in e} != set(range(t)):
                    continue
                degrees = [sum(v in e for e in defect) for v in range(t)]
                weights = [(1,) * t]
                if q == 1:
                    weights += [tuple(1 + int(v == j) for v in range(t))
                                for j, degree in enumerate(degrees) if degree == 1]
                for ew in weights:
                    need(sum(ew[a] * ew[b] for a, b in defect) <= k,
                         "menu weighted budget")
                    for u in range(remaining + 1):
                        for positions in combinations(range(t + 1), u):
                            for values in product(range(k + 1), repeat=u):
                                gaps = [k + 1] * (t + 1)
                                for j, value in zip(positions, values):
                                    gaps[j] = value
                                rows, cap, edges = ordered_template(defect, ew, gaps)
                                if rows:
                                    yield rows, cap, edges, (q, delta, u)


def classify(k):
    states, representatives, provenance = {}, {}, defaultdict(set)
    raw, transitions, gap_hist, permutation_hist = Counter(), Counter(), Counter(), Counter()
    cross_witness, decreasing_witness = None, None
    direct_owners = Counter()
    for rows, weights, edges, signature in menu_templates(k):
        source_D = 2 * signature[0] + signature[1] + signature[2]
        core, cap = normalize_ordered(rows, weights, edges, k)
        final_D = 2 * k + 1 - cap.count(k + 1)
        raw[source_D] += 1
        transitions[source_D, final_D] += 1
        if source_D == 3 and final_D < source_D:
            candidate = (len(rows), weights, edges, cap)
            if decreasing_witness is None or candidate < decreasing_witness:
                decreasing_witness = candidate
        if final_D > 3:
            continue
        record, excess, candidates, signatures, multiplicity = owner(core, cap, k)
        states[record] = final_D
        representatives.setdefault(record, (core, cap))
        provenance[record].add(signature)
        if len(core) == len(rows):
            raw_record = presentation(rows, weights, tuple(range(len(rows))), k)
            if record == raw_record:
                direct_owners[final_D] += 1
        if final_D == 3 and len(signatures) >= 2 and cross_witness is None:
            cross_witness = (record, signatures, multiplicity)
    legacy_to_owner, owner_to_legacy = {}, {}
    for record, (rows, cap) in representatives.items():
        legacy = S.normalize(rows, cap, k)
        need(legacy_to_owner.setdefault(legacy, record) == record,
             "two local owners for one independent weighted isomorphism class")
        need(owner_to_legacy.setdefault(record, legacy) == legacy,
             "one owner collapsed distinct weighted isomorphism classes")
        info = owner(rows, cap, k)
        if states[record] == 3:
            gap_hist[info[1]] += 1
            permutation_hist[info[2]] += 1
    if k <= 2:
        complete, _ = S.template_states(k)
        expected = {s for s in complete if s[2].count(k + 1) >= 2 * k - 2}
        need(set(legacy_to_owner) == expected, "independent complete capped census")
    histogram = Counter(states.values())
    need(direct_owners == histogram, "direct pair-free owner count without a seen-set")
    owner_hist = Counter(record[:3] for record, D in states.items() if D == 3)
    core_hist = Counter(len(representatives[record][0]) for record, D in states.items() if D == 3)
    cross_hist = Counter(len(owner(*representatives[record], k)[3])
                         for record, D in states.items() if D == 3)
    signature_hist = Counter(owner(*representatives[record], k)[3]
                             for record, D in states.items() if D == 3
                             and len(owner(*representatives[record], k)[3]) > 1)
    multiplicities = Counter(owner(*representatives[record], k)[4]
                             for record, D in states.items() if D == 3)
    print("k=", k, "unique_by_deficiency=", dict(sorted(histogram.items())), flush=True)
    print("raw_by_deficiency=", dict(sorted(raw.items())),
          "normalization_transitions=", dict(sorted(transitions.items())), flush=True)
    print("D3_owner_signatures=", dict(sorted(owner_hist.items())),
          "D3_core_orders=", dict(sorted(core_hist.items())), flush=True)
    print("D3_gap_excess=", dict(sorted(gap_hist.items())),
          "D3_local_permutations=", dict(sorted(permutation_hist.items())), flush=True)
    print("D3_number_of_good_menu_signatures=", dict(sorted(cross_hist.items())), flush=True)
    print("D3_multiple_signature_patterns=", dict(sorted(signature_hist.items())),
          "D3_distinct_good_presentations=", dict(sorted(multiplicities.items())), flush=True)
    print("direct_owners_without_deduplication=", dict(sorted(direct_owners.items())),
          "D3_to_lower_D_witness=", decreasing_witness, flush=True)
    if cross_witness:
        print("D3_cross_menu_witness=", cross_witness, flush=True)
    return representatives, states


def controls():
    tested = 0
    nonadjacent = None
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            rows = S.E.matrix(n, code)
            edges = tuple((i, j) for i, j in combinations(range(n), 2)
                          if rows[j] >> i & 1)
            chains = interval_chains(n, edges)
            need(chains == S.G.clone_chains(rows), "all-labelled interval merger audit")
            if nonadjacent is None:
                for chain in chains:
                    if any(abs(i - j) > 1 for i, j in zip(chain, chain[1:])):
                        nonadjacent = (n, edges, chains)
            tested += 1
    print("Independent interval-merger audit:", tested,
          "labelled tournaments through order 5; nonadjacent witness=", nonadjacent,
          flush=True)

    # Source deficiency is not preserved by contraction of light weights.
    rows, weights, edges = ordered_template(((0, 1),), (1, 1), (0, 0, 0))
    core, cap = normalize_ordered(rows, weights, edges, 1)
    need(cap == (2,), "minimal D3 to D2 merger")

    # One good order has a shared endpoint; another spends one fewer defect.
    rows, weights, _ = ordered_template(((0, 1), (0, 2)), (1, 1, 1), (0, 0, 3, 3))
    first = presentation(rows, weights, tuple(range(5)), 2)
    second = presentation(rows, weights, (1, 0, 2, 3, 4), 2)
    need(first[:3] == (0, 1, 2) and second[:3] == (1, 0, 1),
         "shared endpoint versus smaller matching budget")
    print("Changed-gap hostile, two good records of one state:", first, second,
          flush=True)

    # Identical intrinsic gap tournaments and weights do not encode ports.
    examples = []
    for defect in (((0, 1), (2, 3)), ((0, 2), (1, 3))):
        rows, weights, edges = ordered_template(defect, (1, 1, 1, 1), (1, 3, 1, 3, 1))
        need(len(interval_chains(len(rows), edges)) == len(rows), "port hostile pair free")
        heavy, gaps = heavy_gaps(rows, weights, 2)
        local_types = tuple((S.E.canonical(S.G.induced(rows, g)),
                             tuple(sorted(weights[v] for v in g))) for g in gaps)
        examples.append((owner(rows, weights, 2)[0], local_types))
    need(examples[0][1] == examples[1][1], "same induced gap data")
    need(examples[0][0] != examples[1][0], "external ports cannot be discarded")
    print("Port hostile: k=2, gaps=(1,3,1,3,1), endpoint weights all 1;",
          "matchings {01,23} and {02,13} have identical local gap data but distinct owners.",
          flush=True)


def main():
    controls()
    for k in range(1, 4):
        classify(k)
    print("PASS: D<=3 grammar k=1..3, complete independent k<=2 census, exact weighted-class key agreement.")


if __name__ == "__main__":
    main()
