"""Finite endpoint templates for unlabelled tournaments with feedback <= k.

No Burnside formula is used by this compiler. It canonicalizes finite capped
weighted cores from defect templates, not every tournament of the target order.
The independent control generates all tournament classes through order seven
by one-vertex augmentation and solves feedback exactly by subset DP.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path


def need(value, message):
    if not value:
        raise RuntimeError(message)


def load(name, filename):
    p = Path(__file__).resolve().parents[1] / "metagraph_growth_20261020" / filename
    spec = spec_from_file_location(name, p)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load("sparse_enum", "enumeration.py")
G = load("sparse_growth", "growth_check.py")


@lru_cache(None)
def canonical_orders(rows):
    cells = defaultdict(list)
    for v, row in enumerate(rows):
        cells[row.bit_count()].append(v)
    best, orders = None, []
    for blocks in product(*(permutations(cells[d]) for d in sorted(cells))):
        order = tuple(v for b in blocks for v in b)
        code = E.ordercode(rows, order)
        if best is None or code < best:
            best, orders = code, [order]
        elif code == best:
            orders.append(order)
    return best, tuple(orders)


def normalize(rows, weights, k):
    chains = G.clone_chains(rows)
    core = G.induced(rows, [c[0] for c in chains])
    cap = tuple(min(k + 1, sum(weights[v] for v in c)) for c in chains)
    code, orders = canonical_orders(core)
    w = min(tuple(cap[v] for v in order) for order in orders)
    return len(core), code, w


def template_states(k):
    """Complete finite endpoint/gap grammar; no core bank is an input."""
    states = {(1, 0, (v,)) for v in range(1, k + 2)}
    templates = 0
    for t in range(2, 2 * k + 1):
        for r in range(1, k + 1):
            for defect in combinations(list(combinations(range(t), 2)), r):
                if set(v for edge in defect for v in edge) != set(range(t)):
                    continue
                for endpoint_weights in product(range(1, k + 1), repeat=t):
                    if sum(endpoint_weights[a] * endpoint_weights[b]
                           for a, b in defect) > k:
                        continue
                    for gaps in product(range(k + 2), repeat=t + 1):
                        slots, weights = {}, []
                        for a in range(t + 1):
                            if gaps[a]:
                                weights.append(gaps[a])
                            if a < t:
                                slots[a] = len(weights)
                                weights.append(endpoint_weights[a])
                        n = len(weights)
                        rows = [(1 << n) - (1 << (a + 1)) for a in range(n)]
                        for a, b in defect:
                            a, b = slots[a], slots[b]
                            rows[a] ^= 1 << b
                            rows[b] ^= 1 << a
                        states.add(normalize(tuple(rows), weights, k))
                        templates += 1
    return states, templates


def fas(rows, weights=None):
    n = len(rows)
    if weights is None:
        weights = (1,) * n
    dp = [0] + [10 ** 30] * ((1 << n) - 1)
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                cost = weights[v] * sum(weights[u] for u in range(n)
                                        if rest >> u & 1 and rows[v] >> u & 1)
                dp[mask] = min(dp[mask], dp[rest] + cost)
    return dp[-1]


def terms(states, k):
    result = Counter()
    for _, _, weights in states:
        h = weights.count(k + 1)
        result[sum(weights), h] += 1
    return result


def coefficient(ts, n):
    total = int(n == 0)
    for (shift, h), count in ts.items():
        if h == 0:
            total += count * (n == shift)
        elif n >= shift:
            total += count * comb(n - shift + h - 1, h - 1)
    return total


def numerator(ts, degree):
    out = Counter()
    out.update({j: (-1) ** j * comb(degree, j) for j in range(degree + 1)})
    for (shift, h), count in ts.items():
        need(h <= degree, "common denominator")
        for j in range(degree - h + 1):
            out[shift + j] += count * (-1) ** j * comb(degree - h, j)
    return {j: c for j, c in sorted(out.items()) if c}


def inflate(rows, weights):
    block = [v for v, w in enumerate(weights) for _ in range(w)]
    return tuple(sum(int((a < b) if u == v else (rows[u] >> v & 1)) << b
                     for b, v in enumerate(block) if a != b)
                 for a, u in enumerate(block))


def main():
    all_terms = []
    for k in range(3):
        states, tested = template_states(k)
        ts = terms(states, k)
        all_terms.append(ts)
        heavy_counts = Counter()
        for q, code, weights in states:
            rows = E.decode(q, code)
            need(len(G.clone_chains(rows)) == q, "normalized core is pair free")
            need(fas(rows, weights) <= k, "cap state exceeds budget")
            need(q <= 4 * k + 1, "core-size bound")
            h = weights.count(k + 1)
            need(h <= 2 * k + 1, "heavy coordinate bound")
            heavy_counts[h] += 1
            _, orders = canonical_orders(rows)
            for order in orders:
                # The displayed canonical rows use the same degree-cell ordering.
                if tuple(weights[v] for v in order) == weights:
                    need(all(order[v] == v for v in range(q)
                             if weights[v] == k + 1), "heavy coordinate permuted")
        need(heavy_counts[2 * k + 1] ==
             factorial(2 * k) // (2 ** k * factorial(k)), "top matching count")
        print("k=", k, "templates=", tested, "cap_orbits=", len(states),
              "heavy_histogram=", dict(sorted(heavy_counts.items())))
        print("terms (shift,heavy):", sorted(ts.items()))
        print("numerator over (1-x)^" + str(2 * k + 1) + ":",
              numerator(ts, 2 * k + 1))
        print("counts n=0..20:", [coefficient(ts, n) for n in range(21)])
    for n in range(4, 1001):
        p = n ** 4 - 6 * n ** 3 - n ** 2 + 62 * n - 72
        need(p % 8 == 0 and coefficient(all_terms[2], n) == p // 8,
             "budget-two polynomial consequence")
    print("Budget-two exact polynomial consequence checked at n=4..1000")
    reps = {0: ()}
    for n in range(1, 8):
        nxt = {}
        for rows in reps.values():
            for mask in range(1 << (n - 1)):
                child = E.extension(rows, mask)
                nxt.setdefault(E.canonical(child), child)
        reps = nxt
        observed = Counter(fas(rows) for rows in reps.values())
        for k in range(3):
            need(sum(c for d, c in observed.items() if d <= k) ==
                 coefficient(all_terms[k], n), f"independent census n={n},k={k}")
        need(len(reps) == [1, 1, 2, 4, 12, 56, 456][n - 1], "A000568 census")
        print("Independent augmentation n=", n, "classes=", len(reps),
              "exact feedback histogram=", dict(sorted(observed.items())))
    # Independent all-labelled control for weighted block straightening.
    weighted = 0
    for n in range(1, 5):
        for mask in range(1 << comb(n, 2)):
            rows = E.matrix(n, mask)
            for weights in product((1, 2), repeat=n):
                need(fas(rows, weights) == fas(inflate(rows, weights)),
                     "weighted straightening mismatch")
                weighted += 1
    print("Weighted straightening independently checked:", weighted,
          "labelled base/weight cases, inflated order <=8")
    print("PASS")


if __name__ == "__main__":
    main()
