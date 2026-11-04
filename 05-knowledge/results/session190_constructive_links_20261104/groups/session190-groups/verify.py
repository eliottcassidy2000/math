#!/usr/bin/env python3
"""Exact, dependency-free checks for the Session 190 group-action packet.

All gates are explicit exceptions, so python3 -O preserves every check.
The all-order statements are proved in PROOF.md; these are finite controls.
"""
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
import hashlib
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    ans = [0] * len(p)
    for i, x in enumerate(p):
        ans[x] = i
    return tuple(ans)


def order(p):
    q, count = tuple(range(len(p))), 0
    while True:
        q = compose(p, q)
        count += 1
        if q == tuple(range(len(p))):
            return count


def orbit_partition(items, group, action):
    remaining, out = set(items), []
    while remaining:
        x = min(remaining)
        orb = {action(g, x) for g in group}
        require(orb <= remaining, "orbits must partition")
        remaining -= orb
        out.append(sorted(orb))
    return out


def tournament_tools(n):
    pairs = list(combinations(range(n), 2))
    index = {e: i for i, e in enumerate(pairs)}
    full = (1 << len(pairs)) - 1

    def arc(mask, i, j):
        if i == j:
            return False
        value = (mask >> index[tuple(sorted((i, j)))]) & 1
        return bool(value) if i < j else not value

    def relabel(mask, p):
        out = 0
        for k, (i, j) in enumerate(pairs):
            u, v = p[i], p[j]
            bit = ((mask >> k) & 1) ^ (u > v)
            out |= bit << index[tuple(sorted((u, v)))]
        return out

    def path_mask(p):
        return sum(1 << index[tuple(sorted((a, b)))] for a, b in zip(p, p[1:]))

    return pairs, full, arc, relabel, path_mask


def hamilton_dp(n, arc, mask):
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for support in range(1, 1 << n):
        for end in range(n):
            if not dp[support][end]:
                continue
            for nxt in range(n):
                if not (support >> nxt) & 1 and arc(mask, end, nxt):
                    dp[support | (1 << nxt)][nxt] += dp[support][end]
    return sum(dp[-1])


def small_metagraphs():
    rows = []
    for n in range(3, 6):
        pairs, full, arc, relabel, path_mask = tournament_tools(n)
        perms = list(permutations(range(n)))
        transforms = {p: [relabel(t, p) for t in range(full + 1)] for p in perms}
        edges = set()
        auts = {}
        degrees = Counter()
        for t in range(full + 1):
            auts[t] = {p for p in perms if transforms[p][t] == t}
            paths = [p for p in perms if all(arc(t, a, b) for a, b in zip(p, p[1:]))]
            require(len(paths) == hamilton_dp(n, arc, t), "HP permutation/DP mismatch")
            neighbours = {t ^ (full ^ path_mask(p)) for p in paths}
            require(len(paths) == len(neighbours), "path-to-neighbour injectivity")
            degrees[t] = len(neighbours)
            for u in neighbours:
                require(t != u, "n >= 3 has no labelled metagraph loops")
                edges.add(tuple(sorted((t, u))))
        identity = tuple(range(n))
        signed_stabilizer_counts = Counter()
        for t, u in edges:
            require(auts[t] & auts[u] == {identity}, "common ordinary automorphism")
            ordinary = [p for p in perms if set((transforms[p][t], transforms[p][u])) == {t, u}]
            require(ordinary == [identity], "ordinary edge stabilizer is not trivial")
            negative = [p for p in perms if set((transforms[p][t] ^ full, transforms[p][u] ^ full)) == {t, u}]
            require(len(negative) <= 1, "signed edge stabilizer too large")
            for p in negative:
                require(transforms[p][t] ^ full == t, "signed symmetry must not swap endpoints")
                require(order(p) == 2, "signed edge symmetry must be an involution")
            signed_stabilizer_counts[1 + len(negative)] += 1
        ordinary_orbits = orbit_partition(edges, perms, lambda p, e: tuple(sorted((transforms[p][e[0]], transforms[p][e[1]]))))
        signed_group = [(p, eps) for p in perms for eps in (0, 1)]
        signed_orbits = orbit_partition(edges, signed_group, lambda g, e: tuple(sorted((transforms[g[0]][e[0]] ^ (full * g[1]), transforms[g[0]][e[1]] ^ (full * g[1])))))
        q = (n - 1) * (n - 2) // 2
        f = (n - 1) // 2
        c = (q + f) // 2
        require(len(ordinary_orbits) == 2 ** (q - 1), "ordinary edge orbit formula")
        require(len(signed_orbits) == (2 ** (q - 1) + 2 ** (c - 1)) // 2, "signed edge orbit formula")
        rows.append({"n": n, "labelled_tournaments": full + 1, "labelled_edges": len(edges),
                     "ordinary_edge_orbits": len(ordinary_orbits), "signed_edge_orbits": len(signed_orbits),
                     "signed_edge_stabilizer_sizes": dict(sorted(signed_stabilizer_counts.items())),
                     "checked_all_paths_by_permutations_and_subset_dp": True})
    return rows


def normalized_edge_cubes():
    rows = []
    for n in range(3, 9):
        off = [(i, j) for i in range(n) for j in range(i + 2, n)]
        q = len(off)
        index = {e: i for i, e in enumerate(off)}
        rho = [index[(n - 1 - j, n - 1 - i)] for i, j in off]
        f = sum(i == j for i, j in enumerate(rho))
        c = (q + f) // 2
        # Exhaustive through n=7; n=8 has 2^21 states and is formula-only here.
        tested = n <= 7
        if tested:
            full = (1 << q) - 1
            def reflect(x):
                return sum(((x >> i) & 1) << rho[i] for i in range(q))
            fixed = sum(reflect(x) == x for x in range(1 << q))
            anti_fixed = sum(reflect(x) == (x ^ full) for x in range(1 << q))
            require(fixed == 1 << c, "reflection fixed cube size")
            require(anti_fixed == 0, "reflection anti-fixed class forbidden")
        require(f == (n - 1) // 2, "fixed off-path coordinate count")
        rows.append({"n": n, "q": q, "fixed_coordinates": f, "reflection_cycles": c,
                     "ordinary_edge_orbits": 1 << (q - 1),
                     "signed_edge_orbits": ((1 << (q - 1)) + (1 << (c - 1))) // 2,
                     "C2_stabilizer_ordinary_edge_orbits": 1 << (c - 1),
                     "exhaustive_bit_cube": tested})
    return rows


def paley_seven():
    n = 7
    pairs, full, arc, relabel, path_mask = tournament_tools(n)
    squares = {1, 2, 4}
    t = sum(((j - i) % 7 in squares) << k for k, (i, j) in enumerate(pairs))
    perms = list(permutations(range(n)))
    aut = [p for p in perms if relabel(t, p) == t]
    anti = [p for p in perms if relabel(t, p) == (t ^ full)]
    affine = [tuple((a * x + b) % 7 for x in range(7)) for a in range(1, 7) for b in range(7)]
    affine_positive = [tuple((a * x + b) % 7 for x in range(7)) for a in squares for b in range(7)]
    require(set(aut) == set(affine_positive), "Paley automorphism group")
    require(set(aut + anti) == set(affine), "Paley signed automorphism group")
    stabilizers = {
        "vertex": [p for p in affine if p[0] == 0],
        "unordered_pair": [p for p in affine if {p[0], p[1]} == {0, 1}],
        "ordered_pair": [p for p in affine if (p[0], p[1]) == (0, 1)],
    }
    require([len(stabilizers[k]) for k in stabilizers] == [6, 2, 1], "affine stabilizer indices")
    paths = [p for p in perms if all(arc(t, a, b) for a, b in zip(p, p[1:]))]
    require(len(paths) == 189 == hamilton_dp(n, arc, t), "Paley H")
    endpoints = [[0] * n for _ in range(n)]
    for p in paths:
        endpoints[p[0]][p[-1]] += 1
    require(all(endpoints[i][j] == (0 if i == j else (1 if arc(t, i, j) else 8)) for i in range(n) for j in range(n)), "endpoint matrix 1/8")
    neighbours = {t ^ (full ^ path_mask(p)) for p in paths}
    require(len(neighbours) == 189, "Paley distinct neighbours")
    for u in neighbours:
        require(sum(relabel(u, g) == u for g in aut) == 1, "Paley neighbour retains a symmetry")
    neighbour_orbits = orbit_partition(neighbours, aut, lambda g, u: relabel(u, g))
    require([len(o) for o in neighbour_orbits] == [21] * 9, "Paley nine free neighbour orbits")
    signed_group = [(g, 0) for g in aut] + [(g, 1) for g in anti]
    signed_path_orbits = orbit_partition(paths, signed_group, lambda ge, p: tuple(ge[0][x] for x in (p[::-1] if ge[1] else p)))
    require(sorted(map(len, signed_path_orbits)) == [21, 21, 21, 42, 42, 42], "Paley signed paths")
    involutions = [g for g in anti if order(g) == 2]
    fixed_counts = [sum(tuple(g[x] for x in p[::-1]) == p for p in paths) for g in involutions]
    require(len(involutions) == 7 and fixed_counts == [9] * 7, "anti involution fixed paths")
    # Integer identity gives the real six-dimensional complex structure.
    S = [[0 if i == j else (1 if arc(t, i, j) else -1) for j in range(n)] for i in range(n)]
    S2 = [[sum(S[i][k] * S[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    require(S2 == [[1 - 7 * (i == j) for j in range(n)] for i in range(n)], "S squared = J - 7 I")
    require(all(sum(row) == 0 for row in S), "S annihilates constants")
    for g, eps in signed_group:
        require(all(S[g[i]][g[j]] == ((-1) ** eps) * S[i][j] for i in range(n) for j in range(n)), "unitary/antiunitary symmetry")
    # Conjugacy-class controls for the minimal faithful real dimension proof.
    classes = orbit_partition(aut, aut, lambda g, h: compose(compose(g, h), inverse(g)))
    require(sorted(map(len, classes)) == [1, 3, 3, 7, 7], "F21 conjugacy classes")
    require(Counter(order(g) for g in aut) == {1: 1, 3: 14, 7: 6}, "F21 element orders")
    require(not any(order(g) == 21 for g in aut), "F21 is not C21")
    return {"mask": t, "aut_order": len(aut), "signed_aut_order": len(aut + anti),
            "affine_stabilizer_orders": {k: len(v) for k, v in stabilizers.items()},
            "affine_orbit_sizes": {k: 42 // len(v) for k, v in stabilizers.items()},
            "H": len(paths), "endpoint_matrix": endpoints, "neighbour_orbit_sizes": [len(o) for o in neighbour_orbits],
            "signed_path_orbit_sizes": sorted(map(len, signed_path_orbits)),
            "anti_involutions": len(involutions), "fixed_paths_per_anti_involution": fixed_counts,
            "conjugacy_class_sizes": sorted(map(len, classes)), "S_squared": S2,
            "neighbour_orbit_representatives": [o[0] for o in neighbour_orbits]}


def binary_representation():
    def act(v, p):
        return sum(((v >> i) & 1) << p[i] for i in range(7))
    trans = tuple((i + 1) % 7 for i in range(7))
    double = tuple(2 * i % 7 for i in range(7))
    def polynomial(v, exponents):
        out, w = 0, v
        for i in range(max(exponents) + 1):
            if i in exponents:
                out ^= w
            w = act(w, trans)
        return out
    polys = [(0, 1), (0, 1, 3), (0, 2, 3)]
    kernels = [{v for v in range(128) if polynomial(v, p) == 0} for p in polys]
    require(list(map(len, kernels)) == [2, 8, 8], "binary 1+3+3 dimensions")
    require(len({a ^ b ^ c for a in kernels[0] for b in kernels[1] for c in kernels[2]}) == 128, "binary direct sum")
    for kernel in kernels:
        require(all(act(v, g) in kernel for v in kernel for g in (trans, double)), "binary F21 invariance")
    for kernel in kernels[1:]:
        v = min(kernel - {0})
        orbit = {v}
        for _ in range(6):
            v = act(v, trans)
            orbit.add(v)
        require(orbit == kernel - {0}, "binary three-dimensional irreducibility")
    negate = tuple((-i) % 7 for i in range(7))
    require({act(v, negate) for v in kernels[1]} == kernels[2], "converse exchanges the cubic binary summands")
    def mul(x, y):
        out = 0
        while y:
            if y & 1:
                out ^= x
            y >>= 1
            x <<= 1
            if x & 8:
                x ^= 11
        return out
    powers = [1]
    for _ in range(6):
        powers.append(mul(powers[-1], 2))
    require(set(powers) == set(range(1, 8)), "GF8 primitive element")
    logarithm = {x: i for i, x in enumerate(powers)}
    semilinear = []
    for b in range(7):
        for j in range(3):
            def fn(x):
                for _ in range(j):
                    x = mul(x, x)
                return mul(powers[b], x)
            require(all(fn(x ^ y) == fn(x) ^ fn(y) for x in range(8) for y in range(8)), "GF8 semilinear map is F2-linear")
            p = tuple(logarithm[fn(powers[i])] for i in range(7))
            require(p == tuple(((2 ** j) * i + b) % 7 for i in range(7)), "GF8/F7 affine dictionary")
            semilinear.append(p)
    require(len(set(semilinear)) == 21, "GammaL1(8) order")
    gl = []
    for basis in permutations(range(1, 8), 3):
        images = [0] * 8
        for x in range(8):
            for k in range(3):
                if (x >> k) & 1:
                    images[x] ^= basis[k]
        if len(set(images)) == 8:
            gl.append(tuple(logarithm[images[powers[i]]] for i in range(7)))
    require(len(set(gl)) == 168, "GL3(2) order")
    C7 = {tuple((i + b) % 7 for i in range(7)) for b in range(7)}
    normalizer = [g for g in gl if {compose(compose(g, h), inverse(g)) for h in C7} == C7]
    require(set(normalizer) == set(semilinear), "Singer normalizer F21")
    return {"field": "F2", "kernel_dimensions": [1, 3, 3], "kernel_masks": [sorted(k) for k in kernels],
            "GF8_primitive_powers": powers, "GL3_F2_order": len(gl), "Singer_normalizer_order": len(normalizer),
            "Maschke_valid_for_F21": 21 % 2 != 0, "Maschke_not_available_for_GL3_F2": 168 % 2 == 0}


def circulant_controls():
    rows = []
    for n in (3, 5, 7):
        pairs, full, arc, relabel, path_mask = tournament_tools(n)
        perms = list(permutations(range(n)))
        for bits in product((0, 1), repeat=(n - 1) // 2):
            steps = {i + 1 if b else n - i - 1 for i, b in enumerate(bits)}
            mask = sum(((j - i) % n in steps) << k for k, (i, j) in enumerate(pairs))
            h = hamilton_dp(n, arc, mask)
            aut = sum(relabel(mask, p) == mask for p in perms)
            require(h % aut == 0, "automorphism/H divisibility")
            rows.append({"n": n, "steps": sorted(steps), "H": h, "aut_order": aut})
    require(Counter((r["H"], r["aut_order"]) for r in rows if r["n"] == 7) == {(175, 7): 6, (189, 21): 2}, "C7 cube classification")
    return rows


def main():
    result = {"status": "FINITE-EXACT PASS", "small_labelled_metagraphs": small_metagraphs(),
              "normalized_edge_cubes": normalized_edge_cubes(), "paley_seven": paley_seven(),
              "binary_representation": binary_representation(), "circulants": circulant_controls()}
    payload = json.dumps(result, sort_keys=True, indent=2) + "\n"
    target = Path(__file__).with_name("results.json")
    if "--check" in __import__("sys").argv:
        require(target.read_text() == payload, "stored result differs")
    else:
        target.write_text(payload)
    print("FINITE-EXACT PASS")
    print("results_sha256=" + hashlib.sha256(payload.encode()).hexdigest())
    print("All labelled tournament metagraph edges checked at n=3,4,5.")
    print("All normalized edge cubes checked at n=3,...,7; n=8 formula only.")
    print("Paley7: Aut=21, signed Aut=42, H=189, nine free neighbour orbits.")
    print("Paley7 signed path orbits: 21,21,21,42,42,42.")
    print("Binary decomposition 1+3+3 and Singer normalizer order21 verified.")


if __name__ == "__main__":
    main()
