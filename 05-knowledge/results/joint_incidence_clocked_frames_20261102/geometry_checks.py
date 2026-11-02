"""Exact PG(2,3) transport, arithmetic folds, and Collatz clock controls.

Uses SymPy only to find a small generating list. A separate literal BFS
enumerates its closure and all marked-observation fibres.
"""
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, gcd, prod
from pathlib import Path
import json


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def compose(p, q):
    """Apply p and then q; tuples map old positions to new positions."""
    return tuple(q[x] for x in p)


def conway():
    from sympy.combinatorics import Permutation, PermutationGroup
    points = [v for v in product(range(3), repeat=3)
              if any(v) and next(x for x in v if x) == 1]
    lines = [tuple(i for i, v in enumerate(points) if sum(a*b for a, b in zip(v, normal)) % 3 == 0)
             for normal in points]
    require(len(points) == len(lines) == 13 and all(len(L) == 4 for L in lines), "projective plane")
    identity = tuple(range(13))
    moves = {}
    for a, b in combinations(range(13), 2):
        choices = [L for L in lines if a in L and b in L]
        require(len(choices) == 1, "unique line")
        c, d = [x for x in choices[0] if x not in (a, b)]
        p = list(identity)
        p[a], p[b], p[c], p[d] = b, a, d, c
        moves[a, b] = moves[b, a] = tuple(p)
    def word(vertices):
        p = identity
        for a, b in zip(vertices, vertices[1:]):
            if a != b:
                p = compose(p, moves[a, b])
        return p
    loops = [(0, a, b, 0) for a in range(1, 13) for b in range(1, 13) if a != b]
    # These are ALL basepoint triangular loops; hence they generate all closed paths.
    full = PermutationGroup([Permutation(word(w)) for w in loops])
    require(full.order() == 95040, "closed path group")
    gens, words = [], []
    subgroup = PermutationGroup([Permutation(identity)])
    for w in loops:
        p = word(w)
        if not subgroup.contains(Permutation(p)):
            gens.append(p)
            words.append(w)
            subgroup = PermutationGroup([Permutation(g) for g in gens])
            if subgroup.order() == 95040:
                break
    seen = {identity}
    todo = deque([identity])
    while todo:
        p = todo.popleft()
        for g in gens:
            h = compose(p, g)
            if h not in seen:
                seen.add(h)
                todo.append(h)
    require(len(seen) == full.order(), "independent closure")
    fibres = []
    for k in range(6):
        counts = Counter(p[1:k+1] for p in seen)
        expected = factorial(12)//factorial(12-k)
        require(len(counts) == expected and set(counts.values()) == {95040//expected}, "marked fibres")
        fibres.append({"tokens": k, "fixed_hole_states": expected,
                       "moving_hole_states": 13*expected, "fibre_size": 95040//expected})
    hidden4 = min(p for p in seen if p != identity and p[1:5] == identity[1:5])
    require(any(hidden4[i] != i for i in range(5, 13)), "remarking hostile")
    stab4 = [p for p in seen if p[1:5] == identity[1:5]]
    def order(p):
        h, k = p, 1
        while h != identity:
            h = compose(h, p)
            k += 1
        return k
    order_histogram = dict(sorted(Counter(order(p) for p in stab4).items()))
    require(order_histogram == {1: 1, 2: 1, 4: 6}, "quaternion four-marker stabilizer")
    noncommuting = next((p, q) for p in stab4 for q in stab4 if compose(p, q) != compose(q, p))
    collinear = next(w for w in loops if word(w) == identity)
    noncollinear = next(w for w in loops if word(w) != identity)
    require(not any(all(v in L for v in noncollinear[:3]) for L in lines), "triangle holonomy")
    # Explicitly relabel one line as the paper's {4,13,10,12}; this is a
    # local line identification, not a claim to reconstruct its whole figure.
    line = lines[0]
    require(word((line[0], line[1], line[0])) == identity, "involutory step")
    return {"points": points, "lines": lines, "generating_hole_words": words,
            "order_schreier_sims": full.order(), "order_literal_bfs": len(seen),
            "marker_fibres": fibres, "four_marker_hidden_permutation": hidden4,
            "four_marker_stabilizer_orders": order_histogram,
            "four_marker_noncommuting_pair": noncommuting,
            "identity_triangle": collinear, "nonidentity_triangle": noncollinear,
            "nonidentity_transport": word(noncollinear), "full_moving_hole_states": 13*len(seen)}


def collatz():
    # Reuse SF5's strictly increasing children, not a newly claimed inverse-tree construction.
    def children(n):
        powers = (2, 4, 6) if n % 3 == 1 else (3, 5, 7)
        result = [((2**k*n-1)//3, k) for k in powers if (2**k*n-1) % 9 != 0]
        require(len(result) == 2, "two unit children")
        for child, k in result:
            require(child > n and child % 2 == 1 and child % 3 != 0, "proper unit tree")
            require(3*child+1 == 2**k*n, "source preserving inverse edge")
            cur = child
            for _ in range(k+1):
                cur = 3*cur+1 if cur % 2 else cur//2
            require(cur == n and k+1 <= 8, "literal clock")
        return result
    levels, all_nodes = [[5]], {5}
    for depth in range(1, 13):
        new = [child for n in levels[-1] for child, _ in children(n)]
        require(len(new) == len(set(new)) == 2**depth and not (set(new) & all_nodes), "tree injectivity")
        all_nodes.update(new)
        levels.append(new)
    table = {}
    for residue in (1, 2, 4, 5, 7, 8):
        n = residue if residue % 2 else residue+9
        if n == 1:
            n += 18
        table[str(residue)] = [k for _, k in children(n)]
    return {"root": 5, "inverse_exponents_by_mod9": table,
            "level_counts": [len(L) for L in levels], "unique_nodes": len(all_nodes),
            "max_native_steps_per_inverse_edge": 8,
            "triangle_cycle": [1, 4, 2, 1], "wild_35_path": [5, 16, 8, 4, 2, 1]}


def arithmetic():
    exponents = {2: 5, 3: 7, 5: 3, 7: 4, 11: 1, 13: 2}
    N = prod(p**e for p, e in exponents.items())
    def height(b):
        k, t = 0, N
        while t % b == 0:
            t //= b
            k += 1
        return k
    cycle = [4, 2, 5, 3, 7]
    require(all(height(b) == cycle[(i+1) % 5] for i, b in enumerate(cycle)), "inherited height cycle")
    require(height(13) == height(4) == 2, "inherited forward fold")
    def phi(n):
        return sum(gcd(i, n) == 1 for i in range(1, n+1))
    curves = []
    for n in (4, 13):
        divs = [d for d in range(1, n+1) if n % d == 0]
        primes = [p for p in divs if p >= 2 and all(p % q for q in range(2, p))]
        mu = n*prod(F(p+1, p) for p in primes)
        cusps = sum(phi(gcd(d, n//d)) for d in divs)
        e2 = sum((x*x+1) % n == 0 for x in range(n))
        e3 = sum((x*x+x+1) % n == 0 for x in range(n))
        genus = 1+mu/12-F(e2, 4)-F(e3, 3)-F(cusps, 2)
        require(genus == 0, "genus zero")
        curves.append({"level": n, "index": int(mu), "cusps": cusps,
                       "elliptic2": e2, "elliptic3": e3, "genus": int(genus)})
    require(curves[0]["cusps"] != curves[1]["cusps"], "marked modular hostile")
    # No ordering of these exact successive ratios realizes the requested six-state orbit.
    fractions = [F(1, 2), F(5, 2), F(3, 5), F(7, 3), F(13, 7), F(4, 13)]
    require(fractions[0].denominator == fractions[1].denominator == 2, "priority obstruction")
    # The inherited affine clock: y=2x+1 sends 3x+1 to 3y modulo every odd m.
    clocks = []
    for m in (15, 17):
        require(all((2*(3*x+1)+1-3*(2*x+1)) % m == 0 for x in range(m)), "affine conjugacy")
        orbit, x = [], 0
        while x not in orbit:
            orbit.append(x)
            x = (3*x+1) % m
        clocks.append({"modulus": m, "orbit_from_zero": orbit, "repeat_index": orbit.index(x),
                       "image_size": len({(3*x+1) % m for x in range(m)})})
    require(clocks[1]["image_size"] == 17 and len(clocks[1]["orbit_from_zero"]) == 16, "17 clock")
    require(clocks[0]["image_size"] == 5, "15 nonunit hostile")
    return {"group_order_core": N, "inherited_height_cycle": cycle,
            "fold_values": {"4": height(4), "13": height(13)}, "modular_curves": curves,
            "affine_clock_controls": clocks,
            "conflicting_reduced_fractions": [str(fractions[0]), str(fractions[1])]}


def main():
    result = {"conway": conway(), "collatz": collatz(), "arithmetic": arithmetic()}
    Path(__file__).with_name("geometry_verification.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"closed_group": result["conway"]["order_literal_bfs"],
                      "marked_states": result["conway"]["marker_fibres"],
                      "collatz_checked_nodes": result["collatz"]["unique_nodes"]}, sort_keys=True))


if __name__ == "__main__":
    main()
