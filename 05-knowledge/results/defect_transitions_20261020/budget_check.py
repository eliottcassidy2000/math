"""Exact cap-state transitions under one actual arc reversal.

No complete next-budget core census is used by the transition. The previous
independent endpoint compiler is used only as a set-equality control.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from pathlib import Path


def need(test, message):
    if not test:
        raise RuntimeError(message)


path = Path(__file__).resolve().parents[1] / "defect_boundary_20261020/sparse_count.py"
spec = spec_from_file_location("transition_sparse", path)
S = module_from_spec(spec)
spec.loader.exec_module(S)


@lru_cache(None)
def fragments(cap_weight, marked, k):
    """Gap caps around marked singleton endpoints, with target cap L=k+2."""
    K, L = k + 1, k + 2
    result = []
    for gaps in product(range(L + 1), repeat=marked + 1):
        minimum = sum(gaps) + marked
        if (minimum == cap_weight if cap_weight < K else minimum >= K):
            result.append(gaps)
    return tuple(result)


def refine_flip(core, gaps_by_vertex, target_budget, flip=True):
    """Build only the bounded refinement, reverse singleton arc, normalize."""
    owners, weights, marked = [], [], []
    for v, gaps in enumerate(gaps_by_vertex):
        for j, gap in enumerate(gaps):
            if gap:
                owners.append(v)
                weights.append(gap)
            if j + 1 < len(gaps):
                marked.append(len(owners))
                owners.append(v)
                weights.append(1)
    rows = []
    for a, u in enumerate(owners):
        row = 0
        for b, v in enumerate(owners):
            if a == b:
                continue
            arc = (a < b) if u == v else (core[u] >> v) & 1
            row |= int(arc) << b
        rows.append(row)
    if flip:
        need(len(marked) == 2, "arc reversal needs two distinct singleton marks")
        a, b = marked
        rows[a] ^= 1 << b
        rows[b] ^= 1 << a
    else:
        need(not marked, "unmarked carry")
    need(len(rows) <= len(core) + 4, "bounded temporary refinement")
    state = S.normalize(tuple(rows), weights, target_budget)
    if flip:
        need(abs(state[0] - len(core)) <= 4, "symmetric core jump bound")
    return state


def successors(state, k):
    q, code, cap = state
    core = S.E.decode(q, code)
    result = set()
    raw = Counter()
    for gaps in product(*(fragments(c, 0, k) for c in cap)):
        result.add(refine_flip(core, gaps, k + 1, flip=False))
        raw["carry"] += 1
    for a in range(q):
        for b in range(a, q):
            marks = [int(v == a) + int(v == b) for v in range(q)]
            options = [fragments(c, marks[v], k) for v, c in enumerate(cap)]
            for gaps in product(*options):
                result.add(refine_flip(core, gaps, k + 1))
                raw["same chain" if a == b else "different chains"] += 1
    return result, raw


def flip_arc(rows, a, b):
    result = list(rows)
    result[a] ^= 1 << b
    result[b] ^= 1 << a
    return tuple(result)


def exact_marked_refinement(rows, a, b, budget):
    """Extract actual gaps from a labelled tournament, without cap guessing."""
    chains = S.G.clone_chains(rows)
    core = S.G.induced(rows, [chain[0] for chain in chains])
    gap_data = []
    for chain in chains:
        positions = [i for i, v in enumerate(chain) if v in (a, b)]
        boundaries = [-1] + positions + [len(chain)]
        gap_data.append(tuple(boundaries[j + 1] - boundaries[j] - 1
                              for j in range(len(boundaries) - 1)))
    return refine_flip(core, tuple(gap_data), budget)


def main():
    states = {(1, 0, (1,))}
    layers = [states]
    last_successor_counts = {}
    for k in range(2):
        reached, raw, incoming = set(), Counter(), defaultdict(set)
        for parent in sorted(states):
            children, tally = successors(parent, k)
            if k == 1:
                last_successor_counts[parent] = len(children)
            raw.update(tally)
            reached.update(children)
            for child in children:
                incoming[child].add(parent)
        expected, _ = S.template_states(k + 1)
        need(reached == expected, "cap transition vs independent endpoint compiler")
        for q, code, cap in reached:
            core = S.E.decode(q, code)
            need(S.fas(core, cap) <= k + 1, "intrinsic weighted feedback budget")
            need(q <= 4 * (k + 1) + 1, "target core-size bound")
        incoming_histogram = Counter(map(len, incoming.values()))
        print("budget", k, "->", k + 1, "parent cap states=", len(states),
              "child cap states=", len(reached), "raw branches=", dict(raw))
        print("distinct-parent multiplicity histogram:",
              dict(sorted(incoming_histogram.items())))
        states = reached
        layers.append(states)

    # Independent labelled control: extract exact integer gap positions, and
    # compare bounded refinement with flipping the original full tournament.
    controls = 0
    jump_histogram = Counter()
    for n in range(2, 7):
        for rows in S.E.census(n)["reps"].values():
            for a, b in combinations(range(n), 2):
                child = flip_arc(rows, a, b)
                expected = S.normalize(child, (1,) * n, n)
                observed = exact_marked_refinement(rows, a, b, n)
                need(observed == expected, "literal arc flip vs marked refinement")
                old_b = len(S.G.clone_chains(rows))
                new_b = len(S.G.clone_chains(child))
                jump_histogram[new_b - old_b] += 1
                controls += 1
    need(-4 in jump_histogram and 4 in jump_histogram, "sharp core jumps occur")
    print("Literal full-vertex arc controls through order6=", controls,
          "core-jump histogram=", dict(sorted(jump_histogram.items())))

    # Cap-state precision cannot be interpreted as exact old length.
    transitive3 = (6, 4, 0)
    triangle = S.normalize(flip_arc(transitive3, 0, 2), (1, 1, 1), 1)
    need(S.normalize(transitive3, (1, 1, 1), 0) == (1, 0, (1,)),
         "transitive length3 has old cap1")
    need(triangle[0] == 3 and triangle in layers[1], "same-chain nontrivial child")
    need((1, 0, (2,)) in layers[1], "untouched heavy cap refinement")
    print("Cap hostile: old singleton cap1 carries actual lengths1,2,3,...;"
          " length3 creates C3, and untouched length2 requires new cap2.")

    # Reversing a quotient arc reverses every cross edge, not one actual arc.
    core = (2, 4, 1)
    weights = (2, 1, 1)
    inflated = S.inflate(core, weights)
    first_offset = S.normalize(flip_arc(inflated, 0, 2), (1,) * 4, 2)
    one_arc = S.normalize(flip_arc(inflated, 1, 2), (1,) * 4, 2)
    whole_pair = S.normalize(flip_arc(core, 0, 1), weights, 2)
    need(first_offset[0] == 3 and one_arc[0] == 4 and whole_pair[0] == 1,
         "endpoint offset and whole-block reversal hostile")
    print("C3[2,1,1] hostile: first/second clone arc reversals give core3/core4;"
          " quotient-arc reversal gives transitive core1.")
    sink_core = (0, 5, 9, 3)
    heavy_sink = S.normalize(sink_core, (2, 1, 1, 1), 1)
    heavy_triangle = S.normalize(sink_core, (1, 1, 1, 2), 1)
    summary = lambda p: (p[0], p[1], p[2].count(2), sum(p[2]))
    need(summary(heavy_sink) == summary(heavy_triangle), "coarse transition summary")
    need(last_successor_counts[heavy_sink] == 66
         and last_successor_counts[heavy_triangle] == 95,
         "heavy location determines successor profile")
    print("Same C3->sink core, one heavy, cap mass5: heavy sink has66 successors;"
          " heavy triangle vertex has95.")
    mass_controls = 0
    for k in range(3):
        K, L = k + 1, k + 2
        for c in range(1, K + 1):
            for marks in range(3):
                for gaps in fragments(c, marks, k):
                    minimum = marks + sum(gaps)
                    choices = [range(L, L + 3) if d == L else (d,) for d in gaps]
                    masses = set()
                    for actual in product(*choices):
                        need(min(marks + sum(actual), K) == c, "old cap retained")
                        masses.add(marks + sum(actual))
                    for n in range(max(0, minimum - 1), minimum + 3):
                        predicted = n >= minimum and (n == minimum or L in gaps)
                        need((n in masses) == predicted, "fixed-mass refinement guard")
                        mass_controls += 1
    print("Fixed-mass refinement guard boundary controls=", mass_controls)
    print("PASS (explicit exceptions active under normal and -O Python)")


if __name__ == "__main__":
    main()
