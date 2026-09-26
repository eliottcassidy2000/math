"""Finite prime-incidence carriers for Collatz rank; exact and factor-free controls.

Run with python and python -O. Standard library only; explicit checks survive -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import gcd, prod


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def rank(rows, columns):
    rows = [[Fraction(x) for x in row] for row in rows]
    r = 0
    for j in range(columns):
        pivot = next((i for i in range(r, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        v = rows[r][j]
        rows[r] = [x / v for x in rows[r]]
        for i in range(r + 1, len(rows)):
            if rows[i][j]:
                v = rows[i][j]
                rows[i] = [x - v * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


def valuation_rows(values):
    fs = [factor(a) for a in values]
    primes = sorted(set().union(*(f.keys() for f in fs)))
    return {p: [f.get(p, 0) for f in fs] for p in primes}


def coprime_rows(values):
    """Naive exact gcd-free basis; atoms need not be prime."""
    atoms, splits = set(a for a in values if a > 1), 0
    while True:
        conflict = None
        ordered = sorted(atoms)
        for i, b in enumerate(ordered):
            for c in ordered[i + 1:]:
                g = gcd(b, c)
                if g > 1:
                    conflict = b, c, g
                    break
            if conflict:
                break
        if conflict is None:
            break
        b, c, g = conflict
        old_product = prod(atoms)
        atoms.difference_update((b, c))
        atoms.update(x for x in (g, b // g, c // g) if x > 1)
        check(prod(atoms) < old_product, "gcd refinement termination measure")
        splits += 1
    rows, residual = {}, list(values)
    for atom in sorted(atoms):
        row = []
        for i, value in enumerate(residual):
            exponent = 0
            while value % atom == 0:
                value //= atom
                exponent += 1
            residual[i] = value
            row.append(exponent)
        rows[atom] = row
    check(all(x == 1 for x in residual), "coprime-atom reconstruction")
    check(all(gcd(a, b) == 1 for a, b in combinations(rows, 2)), "coprime atoms")
    return rows, splits


def graph_reduce(rows, count):
    """Exact nullity using singleton/pair rows then compressed larger rows."""
    adj = [[] for _ in range(count)]
    pins, wide = set(), []
    for row in rows.values():
        support = [i for i, x in enumerate(row) if x]
        if len(support) == 1:
            pins.add(support[0])
        elif len(support) == 2:
            i, j = support
            adj[i].append((j, -Fraction(row[i], row[j])))
            adj[j].append((i, -Fraction(row[j], row[i])))
        elif len(support) > 2:
            wide.append(row)
    seen, free_components, killed_components = set(), [], 0
    for root in range(count):
        if root in seen:
            continue
        weights, queue, killed = {root: Fraction(1)}, [root], False
        seen.add(root)
        for i in queue:
            killed |= i in pins
            for j, gain in adj[i]:
                proposed = gain * weights[i]
                if j in weights:
                    killed |= proposed != weights[j]
                else:
                    weights[j] = proposed
                    seen.add(j)
                    queue.append(j)
        if killed:
            killed_components += 1
        else:
            free_components.append(weights)
    compressed = [[sum(row[i] * weight for i, weight in comp.items())
                   for comp in free_components] for row in wide]
    b = len(free_components)
    residual_rank = rank(compressed, b)
    return count - b + residual_rank, b, residual_rank, killed_components


def strip_shared(value, others):
    for other in others:
        while (g := gcd(value, other)) > 1:
            value //= g
    return value


def exclusive_graph(values):
    """Factor-free private vertices and pair-exclusive-prime edges."""
    count = len(values)
    private = {i for i, a in enumerate(values)
               if strip_shared(a, [b for j, b in enumerate(values) if i != j]) > 1}
    edges, adj = [], [[] for _ in values]
    for i, j in combinations(range(count), 2):
        residual = strip_shared(gcd(values[i], values[j]),
                                [a for k, a in enumerate(values) if k not in (i, j)])
        if residual > 1:
            edges.append((i, j))
            adj[i].append(j)
            adj[j].append(i)
    colors, forced, odd_components = {}, set(), 0
    for root in range(count):
        if root in colors:
            continue
        colors[root], queue, odd = 0, [root], False
        for i in queue:
            for j in adj[i]:
                if j in colors:
                    odd |= colors[i] == colors[j]
                else:
                    colors[j] = 1 - colors[i]
                    queue.append(j)
        odd_components += odd
        if odd or private.intersection(queue):
            forced.update(queue)
    return private, edges, forced, odd_components


def graph_peel(values):
    """Repeatedly remove coefficients forced zero; no factorization used."""
    remaining, rounds = list(range(len(values))), 0
    while remaining:
        _, _, forced, _ = exclusive_graph([values[i] for i in remaining])
        if not forced:
            break
        remaining = [i for j, i in enumerate(remaining) if j not in forced]
        rounds += 1
    return remaining, rounds


def odd_orbit(n, count):
    nodes = [n]
    while len(nodes) < count:
        x = 3 * nodes[-1] + 1
        while x % 2 == 0:
            x //= 2
        nodes.append(x)
    return nodes


def paired_rows(nodes):
    count = len(nodes)
    atoms, _ = coprime_rows([3 * n for n in nodes] + [3 * n + 1 for n in nodes])
    # Coordinates are separate: do not add the two valuation rows together.
    return {(atom, side): row[side * count:(side + 1) * count]
            for atom, row in atoms.items() for side in (0, 1)}


def graph_values(count, edges, twists=None):
    primes, p = [], 2
    while len(primes) < len(edges):
        if len(factor(p)) == 1 and factor(p).get(p) == 1:
            primes.append(p)
        p += 1
    values = [1] * count
    for e, (i, j) in enumerate(edges):
        values[i] *= primes[e] ** (2 if twists == (e, i) else 1)
        values[j] *= primes[e] ** (2 if twists == (e, j) else 1)
    return values


def order(a, p):
    for k in range(1, p):
        if pow(a, k, p) == 1:
            return k
    raise RuntimeError("order not found")


def main():
    for p, d, log in ((223, 37, 180), (233, 29, 72)):
        check(factor(p) == {p: 1}, "prime control")
        check(order(2, p) == d and order(3, p) == p - 1, "orders")
        check(pow(3, log, p) == 2, "discrete log")
        chi = tuple(pow(x, (p - 1) // 2, p) for x in (2, 3, p - 1))
        print(f"prime {p}: ord2={d}; ord3={p-1}; log3(2)={log}; index<2>={(p-1)//d}; chi(2,3,-1)={chi}")
    check(2 ** 37 - 1 == 223 * 616318177, "223 Mersenne identity")
    check(2 ** 29 - 1 == 233 * 1103 * 2089, "233 Mersenne identity")

    cycle = [(i, (i + 1) % 10) for i in range(10)]
    petersen = [(i, (i + 1) % 5) for i in range(5)]
    petersen += [(i, i + 5) for i in range(5)]
    petersen += [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    for name, values, expected in (("C10", graph_values(10, cycle), 9),
                                   ("twisted C10", graph_values(10, cycle, (0, 0)), 10),
                                   ("Petersen", graph_values(10, petersen), 10)):
        rows = valuation_rows(values)
        direct = rank(rows.values(), 10)
        reduced = graph_reduce(rows, 10)
        private, edges, forced, odd = exclusive_graph(values)
        check(direct == reduced[0] == expected, "ten-vertex exact rank")
        check(len(private) == 0, "no private primes in graph realizations")
        check(len(forced) == (10 if name == "Petersen" else 0), "odd-cycle certificate boundary")
        print(f"ten-node {name}: rank={direct}; pair edges={len(edges)}; private=0; "
              f"sign-forced={len(forced)}; free graph components={reduced[1]}")
    base, twist = graph_values(10, cycle), graph_values(10, cycle, (0, 0))
    check(exclusive_graph(base)[1] == exclusive_graph(twist)[1], "same unweighted graph, different rank")
    check(prod(base[::2]) == prod(base[1::2]), "C10 alternating relation")
    check(prod(twist[::2]) != prod(twist[1::2]), "twisted C10 loses relation")

    orbit = odd_orbit(27, 45)
    wanted = {41, 31, 47, 71, 107, 121, 91, 137, 175, 325}
    selected = [(i, n) for i, n in enumerate(orbit) if n in wanted]
    check(len(selected) == 10, "ten actual selected nodes")
    values = [3 * n for _, n in selected]
    private, edges, forced, odd = exclusive_graph(values)
    check(len(private) == 7 and len(forced) == 10 and odd == 1, "actual triangle gain")
    check(rank(valuation_rows(values).values(), 10) == 10, "actual ten-node full rank")
    triple = [91, 175, 325]
    check(all(a == prod(gcd(a, b) for j, b in enumerate(triple) if i != j)
              for i, a in enumerate(triple)), "actual triangle gcd equality")
    check(exclusive_graph(triple)[0] == set(), "actual triangle no private prime")
    print(f"orbit27 actual ten selected (odd index,value): {selected}")
    print(f"orbit27 carrier: private={len(private)}; pair edges={edges}; odd components={odd}; rank=10")
    print("orbit27 triangle: 91=7*13, 175=5^2*7, 325=5^2*13; no private prime; every gcd dominance equality")

    tuples, certs, refined_splits, reduced_cases = 0, 0, 0, Counter()
    for count in (2, 3, 4):
        for values in combinations_with_replacement(range(2, 21), count):
            rows = valuation_rows(values)
            direct = rank(rows.values(), count)
            atom_rows, splits = coprime_rows(values)
            reduced = graph_reduce(atom_rows, count)
            refined_splits += splits
            check(direct == reduced[0], "compressed graph exact rank")
            private, edges, forced, odd = exclusive_graph(values)
            if len(forced) == count:
                check(direct == count, "factor-free graph soundness")
                certs += 1
            remaining, rounds = graph_peel(values)
            if not remaining:
                check(direct == count, "factor-free iterative carrier")
            reduced_cases[reduced[1]] += 1
            tuples += 1
    print(f"complete small universe: {tuples} multisets, graph certificate={certs}; "
          f"{refined_splits} gcd splits; atom/prime/compressed rank all pass")

    counts, examples, rank_hist = Counter(), {}, Counter()
    for source in range(3, 20001, 2):
        nodes = odd_orbit(source, 10)
        if len(set(nodes)) != 10:
            continue
        values = [3 * n for n in nodes]
        private, edges, forced, odd = exclusive_graph(values)
        remaining, rounds = graph_peel(values)
        atom_rows, splits = coprime_rows(values)
        reduced = graph_reduce(atom_rows, 10)
        full = reduced[0] == 10
        pair_rank = graph_reduce(paired_rows(nodes), 10)[0]
        rank_hist[reduced[0]] += 1
        counts["injective_prefixes"] += 1
        counts["all_private"] += len(private) == 10
        counts["graph_full"] += len(forced) == 10
        counts["full_rank"] += full
        counts["full_paired_rank"] += pair_rank == 10
        counts["odd_component"] += odd > 0
        counts["graph_gain"] += len(forced) == 10 and len(private) < 10
        counts["peel_full"] += not remaining
        counts["nodip"] += min(nodes) >= source
        counts["nodip_odd_component"] += min(nodes) >= source and odd > 0
        if len(forced) == 10:
            check(full, "actual graph certificate")
        if not remaining:
            check(full, "actual iterative graph certificate")
        if source <= 2001 or not full:
            check(rank(valuation_rows(values).values(), 10) == reduced[0], "actual direct rank audit")
            direct_pair = {(p, 0): row for p, row in valuation_rows(values).items()}
            direct_pair.update({(p, 1): row for p, row in valuation_rows([3 * n + 1 for n in nodes]).items()})
            check(rank(direct_pair.values(), 10) == pair_rank, "separate-coordinate direct prime audit")
            counts["direct_audits"] += 1
        if len(forced) == 10 and len(private) < 10 and "graph_gain" not in examples:
            examples["graph_gain"] = (source, nodes, len(private), edges)
        if odd and "odd" not in examples:
            examples["odd"] = (source, nodes, len(private), edges)
        if odd and min(nodes) >= source and "nodip_odd" not in examples:
            examples["nodip_odd"] = (source, nodes, len(private), edges)
    print(f"consecutive ten-node census, odd sources3..20000: {dict(sorted(counts.items()))}")
    print(f"first-slot rank histogram: {dict(sorted(rank_hist.items()))}")
    print(f"first actual graph-gain example: {examples.get('graph_gain')}")
    print(f"first actual odd-component example: {examples.get('odd')}")
    print(f"first actual no-dip odd-component example: {examples.get('nodip_odd')}")
    values = [117, 39]
    check(graph_peel(values)[0] == [0, 1], "unweighted graph stopping witness")
    atom_rows, splits = coprime_rows(values)
    check(graph_reduce(atom_rows, 2)[0] == 2, "parallel weighted-edge obstruction")
    print(f"unweighted miss: first slots117,39 have atom rows{atom_rows}; exact rank2 despite one unweighted edge")
    check(graph_reduce(valuation_rows((9, 15, 3)), 3)[0] == 2, "first-slot dependence does not settle paired rank")
    check(graph_reduce(paired_rows([3, 5, 1]), 3)[0] == 3, "two-coordinate repair")
    check(rank(valuation_rows((5, 75)).values(), 2) == 2, "raw-node scaling hostile")
    check(rank(valuation_rows((15, 225)).values(), 2) == 1, "common factor3 must be retained")
    print("scale hostile: (5,75) has rank2, but (15,225) has rank1; coefficient3 is load-bearing")
    print("PASS: exact graph compression and factor-free sufficient carrier; no finite-modulus Collatz proof")


if __name__ == "__main__":
    main()
