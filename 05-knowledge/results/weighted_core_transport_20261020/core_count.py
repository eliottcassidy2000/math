"""Exact positive-weight orbit census on canonical pair-free tournament cores.

Imports the frozen 76caae6f1e metagraph-growth implementation for an independent
source of all small tournament classes and maximal transitive modules.
Run: python 05-knowledge/results/weighted_core_transport_20261020/core_count.py
"""
from collections import Counter, defaultdict
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from math import comb, factorial, gcd
from pathlib import Path


def need(test, message):
    if not test:
        raise RuntimeError(message)


def load(name, filename):
    p = Path(__file__).resolve().parents[1] / "metagraph_growth_20261020" / filename
    spec = spec_from_file_location(name, p)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


e = load("core_enum", "enumeration.py")
g = load("core_growth", "growth_check.py")


def automorphisms(rows):
    cells = defaultdict(list)
    for v, r in enumerate(rows):
        cells[r.bit_count()].append(v)
    blocks = tuple(tuple(vs) for _, vs in sorted(cells.items()))
    original = tuple(v for b in blocks for v in b)
    result = []
    for selections in product(*(permutations(b) for b in blocks)):
        image = tuple(v for b in selections for v in b)
        p = [0] * len(rows)
        for v, w in zip(original, image):
            p[v] = w
        if all(((rows[v] >> w) & 1) == ((rows[p[v]] >> p[w]) & 1)
               for v in range(len(rows)) for w in range(len(rows))):
            result.append(tuple(p))
    need(len(result) == e.canonical(rows, aut=True)[1], "automorphism count")
    return tuple(result)


def antimorphisms(rows):
    n = len(rows)
    target = e.converse(rows)
    source_cells, target_cells = defaultdict(list), defaultdict(list)
    for v in range(n):
        source_cells[rows[v].bit_count()].append(v)
        target_cells[target[v].bit_count()].append(v)
    if {d: len(c) for d, c in source_cells.items()} != {
            d: len(c) for d, c in target_cells.items()}:
        return ()
    degrees = sorted(source_cells)
    original = tuple(v for d in degrees for v in source_cells[d])
    result = []
    for selections in product(*(permutations(target_cells[d]) for d in degrees)):
        image = tuple(v for b in selections for v in b)
        p = [0] * n
        for v, w in zip(original, image):
            p[v] = w
        if all(((rows[v] >> w) & 1) == ((target[p[v]] >> p[w]) & 1)
               for v in range(n) for w in range(n)):
            result.append(tuple(p))
    return tuple(result)


def cycle_type(p):
    used, sizes = set(), []
    for v in range(len(p)):
        if v in used:
            continue
        size, w = 0, v
        while w not in used:
            used.add(w)
            size += 1
            w = p[w]
        sizes.append(size)
    return tuple(sorted(sizes))


def canonical_weight(weights, group):
    return min(tuple(weights[p[v]] for v in range(len(weights))) for p in group)


def weight_layers(q, group, maximum):
    """Constructive increment recurrence; no orbit-averaging formula."""
    layers = {q: {(1,) * q}}
    for total in range(q + 1, maximum + 1):
        nxt = set()
        for old in layers[total - 1]:
            for v in range(q):
                weights = list(old)
                weights[v] += 1
                nxt.add(canonical_weight(weights, group))
        layers[total] = nxt
    return layers


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def fixed_weight_count(total, cycles):
    """Coefficient of x^q / product_cycles(1-x^length)."""
    q = sum(cycles)
    if total < q:
        return 0
    target = total - q
    dp = [1] + [0] * target
    for length in cycles:
        for j in range(length, target + 1):
            dp[j] += dp[j - length]
    return dp[target]


def cycle_formula(total, types, group_order):
    out = Fraction(sum(count * fixed_weight_count(total, cycles)
                       for cycles, count in types.items()), group_order)
    need(out.denominator == 1, "cycle formula integral")
    return out.numerator


def inflate(rows, weights):
    blocks = tuple(v for v, size in enumerate(weights) for _ in range(size))
    return tuple(sum((int(i < j) if v == w else ((rows[v] >> w) & 1)) << j
                     for j, w in enumerate(blocks)) for i, v in enumerate(blocks))


def weight_gcd(weights):
    z = 0
    for w in weights:
        z = gcd(z, w)
    return z


def mu(n):
    value, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            value = -value
            if n % p == 0:
                return 0
        p += 1
    return -value if n > 1 else value


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def jordan(n, k):
    return sum(mu(d) * (n // d) ** k for d in divisors(n))


def cycle_monomial(cycles, variables=7):
    count = Counter(cycles)
    return tuple(count[j] for j in range(1, variables + 1))


def degree(monomial):
    return sum((j + 1) * power for j, power in enumerate(monomial))


def trim(poly, bound):
    return {m: c for m, c in poly.items() if c and degree(m) <= bound}


def substitute_chain(poly, bound, scale=1, sign=1):
    """Substitute p_i -> scale*p_i/(1-sign*p_i), exact weighted truncation."""
    out = defaultdict(Fraction)
    for monomial, coefficient in poly.items():
        terms = {(0,) * len(monomial): coefficient}
        for i, power in enumerate(monomial):
            if not power:
                continue
            nxt = defaultdict(Fraction)
            for old, c in terms.items():
                for k in range(power, (bound - degree(old)) // (i + 1) + 1):
                    m = list(old)
                    m[i] = k
                    nxt[tuple(m)] += c * scale ** power * sign ** (k - power) * comb(k - 1, power - 1)
            terms = nxt
        for m, c in terms.items():
            out[m] += c
    return trim(out, bound)


def derivative_one(poly):
    out = defaultdict(Fraction)
    for m, c in poly.items():
        if m[0]:
            shifted = (m[0] - 1,) + m[1:]
            out[shifted] += m[0] * c
    return dict(out)


def print_cycle_poly(poly):
    return sorted((tuple(j + 1 for j, power in enumerate(m) for _ in range(power)),
                   str(c)) for m, c in poly.items())


def main():
    maximum = 12
    atlases, core_data = {}, {}
    weighted_by_order = defaultdict(Counter)
    generated_classes = defaultdict(set)
    primitive = Counter()
    self_converse = Counter()
    total_orbits = 0
    core_index = defaultdict(Fraction)
    core_index[(0,) * 7] = Fraction(1)
    print("FINITE-EXACT: all pair-free tournament cores through order 7.")
    print("Weight orbits generated constructively by increment/canonicalize through total 12.")
    for q in range(1, 8):
        atlas = atlases[q] = e.census(q)
        cores = {code: rows for code, rows in atlas["reps"].items()
                 if len(g.clone_chains(rows)) == q}
        hist = Counter()
        core_sc = 0
        mass = Fraction(0)
        for code, rows in sorted(cores.items()):
            group = automorphisms(rows)
            types = Counter(cycle_type(p) for p in group)
            for cycles, count in types.items():
                core_index[cycle_monomial(cycles)] += Fraction(count, len(group))
            anti = antimorphisms(rows)
            anti_types = Counter(cycle_type(p) for p in anti)
            need(bool(anti) == (e.canonical(e.converse(rows)) == code),
                 "anti-isomorphism existence")
            need(not anti or len(anti) == len(group), "anti coset cardinality")
            core_sc += bool(anti)
            mass += Fraction(1, len(group))
            hist[(len(group), tuple(sorted(types.items())))] += 1
            layers = weight_layers(q, group, maximum)
            core_data[(q, code)] = dict(rows=rows, group=group, types=types, layers=layers)
            for n, layer in layers.items():
                direct = {canonical_weight(w, group) for w in compositions(n, q)}
                need(layer == direct, "increment recurrence vs direct positive compositions")
                need(len(layer) == cycle_formula(n, types, len(group)),
                     "constructive layers vs independent cycle formula")
                sc = sum(canonical_weight(tuple(w[anti[0][v]] for v in range(q)), group) == w
                         for w in layer) if anti else 0
                need(sc == cycle_formula(n, anti_types, len(group)),
                     "self-converse weight orbits vs anti-coset formula")
                weighted_by_order[n][q] += len(layer)
                total_orbits += len(layer)
                prim = sum(weight_gcd(w) == 1 for w in layer)
                mobius = sum(mu(d) * len(layers.get(n // d, ())) for d in divisors(n))
                need(prim == mobius, "primitive gcd orbit inversion")
                if n <= 7:
                    primitive[n] += prim
                    self_converse[n] += sc
                    for w in layer:
                        expanded = inflate(rows, w)
                        need(sorted(map(len, g.clone_chains(expanded))) == sorted(w),
                             "weighted inflation maximal module lengths")
                        key = e.canonical(expanded)
                        need(key not in generated_classes[n], "weighted core duplicate class")
                        generated_classes[n].add(key)
        print(f"\ncore_order={q}: all_classes={len(atlas['reps'])}, pair_free_cores={len(cores)}")
        inverse_mass = sum((Fraction((-1) ** (q - k) * comb(q - 1, k - 1)
                                     * (1 << comb(k, 2)), factorial(k))
                            for k in range(1, q + 1)), Fraction(0))
        need(mass == inverse_mass, "formal labelled-core inverse composition")
        print("self_converse_cores=", core_sc, "automorphism_mass=", mass,
              "labelled_cores=", int(mass * factorial(q)))
        print("core automorphism cycle-index histogram:", sorted(hist.items()))

    print("\nEXACT TRIANGULAR COUNT: row n lists core_order: contribution.")
    for n in range(1, 8):
        need(generated_classes[n] == set(atlases[n]["reps"]),
             "weighted core class set vs independent full census")
        total = sum(weighted_by_order[n].values())
        need(total == len(atlases[n]["reps"]), "weighted count completeness")
        need(self_converse[n] == sum(c == d for c, d in atlases[n]["dual"].items()),
             "self-converse weighted-core reconstruction")
        need(primitive[n] == sum(mu(d) * len(atlases[n // d]["reps"])
                                 for d in divisors(n)), "global primitive inversion")
        print("n=", n, "contributions=", sorted(weighted_by_order[n].items()),
              "A(n)=", total, "gcd_primitive=", primitive[n])
        print("self_converse=", self_converse[n], "merged=",
              (total + self_converse[n]) // 2)

    print("\nSAME-ORDER, SAME-GROUP-ORDER RESPONSE HOSTILE:")
    witness = None
    for q in range(1, 8):
        keys = [key for key in core_data if key[0] == q]
        for i, key_a in enumerate(keys):
            a = core_data[key_a]
            for key_b in keys[i + 1:]:
                b = core_data[key_b]
                if len(a["group"]) != len(b["group"]):
                    continue
                responses_a = tuple(len(a["layers"][n]) for n in range(q, maximum + 1))
                responses_b = tuple(len(b["layers"][n]) for n in range(q, maximum + 1))
                if responses_a != responses_b:
                    first = next(n for n in range(q, maximum + 1)
                                 if len(a["layers"][n]) != len(b["layers"][n]))
                    witness = q, key_a, key_b, first
                    break
            if witness:
                break
        if witness:
            break
    need(witness is not None, "expected response hostile absent")
    q, key_a, key_b, first = witness
    for label, key in (("A", key_a), ("B", key_b)):
        d = core_data[key]
        print(label, "core_order", q, "rows", d["rows"], "group_order", len(d["group"]),
              "automorphisms", d["group"],
              "cycles", sorted(d["types"].items()),
              "weight_response_n=q..12", [len(d["layers"][n]) for n in range(q, maximum + 1)])
    print("first separating total=", first)

    print("\nFULL FORMAL CYCLE-INDEX CONTROLS THROUGH WEIGHTED DEGREE 7:")
    tournament_index = defaultdict(Fraction)
    tournament_index[(0,) * 7] = Fraction(1)
    for atlas in atlases.values():
        for rows in atlas["reps"].values():
            group = automorphisms(rows)
            for p in group:
                tournament_index[cycle_monomial(cycle_type(p))] += Fraction(1, len(group))
    need(substitute_chain(core_index, 7) == dict(tournament_index),
         "T(p)=C(p_i/(1-p_i))")
    need(substitute_chain(tournament_index, 7, sign=-1) == dict(core_index),
         "C(p)=T(p_i/(1+p_i))")
    left = trim(derivative_one(tournament_index), 6)
    right = trim({m: c * 2 ** sum(m) for m, c in tournament_index.items()}, 6)
    need(left == right, "rooted tournament cycle-index identity")
    core_derivative = derivative_one(core_index)
    left = defaultdict(Fraction)
    for m, c in core_derivative.items():
        for shift, factor in ((0, 1), (1, 2), (2, 1)):
            left[(m[0] + shift,) + m[1:]] += factor * c
    need(trim(left, 6) == substitute_chain(core_index, 6, scale=2),
         "core functional differential recurrence")
    boundary = {m: c for m, c in core_index.items() if m[0] == 0}
    print("C boundary p1=0 (cycle-length multiset, coefficient):", print_cycle_poly(boundary))
    print("Core and tournament substitution/inversion: passed through weighted degree7.")
    print("Rooted T derivative and core PDE: passed through weighted degree6.")

    print("\nC3 PRIMITIVE NUMBER-THEORY CONTROL:")
    triangle = next(d for (q, _), d in core_data.items() if q == 3)
    for n in range(1, 61):
        # Ordered primitive positive compositions in 3 parts, divided by
        # cyclic rotations; only (1,1,1) supplies nonidentity fixed points.
        value = Fraction(jordan(n, 2) - 3 * jordan(n, 1), 6)
        if n == 1:
            value = Fraction(0)
        if n == 3:
            value += Fraction(2, 3)
        direct = {canonical_weight(w, triangle["group"])
                  for w in compositions(n, 3) if weight_gcd(w) == 1}
        need(value.denominator == 1 and value.numerator == len(direct),
             "Jordan/phi primitive C3 formula")
    print("P_C3(n)=(J2(n)-3*phi(n))/6+(2/3)[n=3], n>1; P_C3(1)=0.")
    print("Independent primitive composition orbit audit n=1..60 passed.")

    print("\nFROZEN-CORE TRUNCATION, NOT THE FULL A(n) FOR n>7:")
    for n in range(8, maximum + 1):
        print("n=", n, "classes with pair-free core order <=7=",
              sum(weighted_by_order[n].values()),
              "by_core_order=", sorted(weighted_by_order[n].items()))
    print("Total constructed weight orbits in the checked finite universe=", total_orbits)
    print("ALL CHECKS PASSED (explicit exceptions; normal and -O active).")


if __name__ == "__main__":
    main()
