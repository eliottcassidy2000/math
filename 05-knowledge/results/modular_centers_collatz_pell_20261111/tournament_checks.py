"""Literal controls for TOURNAMENTS.md; no external dependencies."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd


def require(ok, message):
    if not ok:
        raise ValueError(message)


def prime_factors(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def jacobi(a, m):
    value = 1
    for p, exponent in prime_factors(m).items():
        if a % p == 0:
            return 0
        legendre = 1 if pow(a, (p - 1) // 2, p) == 1 else -1
        value *= legendre ** exponent
    return value


def level(d, p):
    v = 0
    while d % p == 0:
        v += 1
        d //= p
    return v, d % p


def recursive_rows(p, exponent):
    n = p ** exponent
    return [sum(1 << b for b in range(n) if b != a and
                pow(level((b - a) % n, p)[1], (p - 1) // 2, p) == 1)
            for a in range(n)]


def path_tables(rows):
    # Inherited subset-DP recurrence, independently checked by permutations below.
    n = len(rows)
    end = [[0] * n for _ in range(1 << n)]
    start = [[0] * n for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        for v in range(n):
            if not mask >> v & 1:
                continue
            rest = mask ^ (1 << v)
            if rest == 0:
                end[mask][v] = start[mask][v] = 1
            else:
                end[mask][v] = sum(end[rest][u] for u in range(n) if rows[u] >> v & 1)
                start[mask][v] = sum(start[rest][u] for u in range(n) if rows[v] >> u & 1)
    full = (1 << n) - 1
    q = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        for i in range(n):
            for j in range(n):
                q[i][j] += end[mask][i] * start[full ^ mask][j]
    u = [start[full][j] + sum(q[i][j] for i in range(n)) for j in range(n)]
    w = [[q[i][j] + q[j][i] for j in range(n)] for i in range(n)]
    return sum(end[full]), u, w, q, start[full], end[full]


def literal_cover_counts(rows):
    n = len(rows)
    h, u = 0, [0] * n
    w = [[0] * n for _ in range(n)]
    for ordering in permutations(range(n)):
        bad = [i for i in range(n - 1) if not rows[ordering[i]] >> ordering[i + 1] & 1]
        if not bad:
            h += 1
            u[ordering[0]] += 1
            cuts = range(1, n)
        elif len(bad) == 1:
            cuts = [bad[0] + 1]
        else:
            continue
        for cut in cuts:
            a, b = ordering[cut - 1], ordering[cut]
            u[b] += 1
            w[a][b] += 1
            w[b][a] += 1
    return h, u, w


def automorphism_count(rows):
    n = len(rows)
    degrees = [row.bit_count() for row in rows]
    candidates = [[j for j in range(n) if degrees[i] == degrees[j]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(candidates[i]))
    mapping = [-1] * n
    def visit(depth, used):
        if depth == n:
            return 1
        i, total = order[depth], 0
        for j in candidates[i]:
            if used >> j & 1:
                continue
            if any((rows[i] >> a & 1) != (rows[j] >> mapping[a] & 1) for a in order[:depth]):
                continue
            mapping[i] = j
            total += visit(depth + 1, used | (1 << j))
        mapping[i] = -1
        return total
    return visit(0, 0)


def field27_rows_and_fixers():
    # Polynomial representatives over F3, with z^3=z+1.
    coeff = [tuple((x // 3 ** k) % 3 for k in range(3)) for x in range(27)]
    def add(a, b):
        return sum(((coeff[a][k] + coeff[b][k]) % 3) * 3 ** k for k in range(3))
    def negative(a):
        return sum((-coeff[a][k] % 3) * 3 ** k for k in range(3))
    def multiply(a, b):
        c = [0] * 5
        for i in range(3):
            for j in range(3):
                c[i + j] += coeff[a][i] * coeff[b][j]
        for degree in (4, 3):
            c[degree - 2] += c[degree]
            c[degree - 3] += c[degree]
        return sum(c[k] % 3 * 3 ** k for k in range(3))
    mult = [[multiply(a, b) for b in range(27)] for a in range(27)]
    def cube(a):
        return mult[mult[a][a]][a]
    square_set = {mult[a][a] for a in range(1, 27)}
    inverses = {a: next(b for b in range(1, 27) if mult[a][b] == 1) for a in range(1, 27)}
    rows = [sum(1 << b for b in range(27) if add(b, negative(a)) in square_set) for a in range(27)]
    checked = 0
    for a, b in combinations(range(27), 2):
        d = add(b, negative(a))
        g = [add(a, mult[d][cube(mult[add(x, negative(a))][inverses[d]])]) for x in range(27)]
        require(g[a] == a and g[b] == b and sorted(g) == list(range(27)), "field pair fixer")
        require(any(g[x] != x for x in range(27)) and all(g[g[g[x]]] == x for x in range(27)), "field order three")
        for x in range(27):
            require(sum(1 << g[y] for y in range(27) if rows[x] >> y & 1) == rows[g[x]], "field automorphism")
        checked += 1
    common = sorted({(rows[a] & rows[b]).bit_count() for a, b in combinations(range(27), 2)})
    require(common == [6], "field common outneighbors")
    return {"pair_fixers_of_order_three": checked, "common_outneighbor_counts": common,
            "polynomial": "z^3-z-1 over F3"}


def run():
    for m in range(3, 300, 4):
        chi = [jacobi(a, m) for a in range(m)]
        require(all(chi[(-a) % m] == -chi[a] for a in range(m)), "Jacobi skewness")
        units = sum(gcd(a, m) == 1 for a in range(m))
        require(Counter(chi)[1] == Counter(chi)[-1] == units // 2, "Jacobi degrees")
        require(Counter(chi)[0] == m - units, "Jacobi missing edges")
    for m in (15, 27):
        least = min(prime_factors(m))
        require(all(gcd(a - b, m) == 1 for a, b in combinations(range(least), 2)), "unit clique witness")
        require(not any(all(gcd(a - b, m) == 1 for a, b in combinations(s, 2))
                        for s in combinations(range(m), least + 1)), "no larger unit clique")
    levels = []
    for p, exponent in ((3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (7, 1), (7, 2), (7, 3), (11, 1), (11, 2)):
        n = p ** exponent
        rows = recursive_rows(p, exponent)
        counts = Counter()
        for a in range(n):
            require(rows[a].bit_count() == (n - 1) // 2, "recursive degree")
            for b in range(n):
                if a == b:
                    continue
                require(((rows[a] >> b & 1) + (rows[b] >> a & 1)) == 1, "recursive orientation")
                v, unit = level((b - a) % n, p)
                # Independently compare the first unequal base-p digits.
                aa, bb = a, b
                while aa % p == bb % p:
                    aa //= p
                    bb //= p
                require((rows[a] >> b & 1) == (pow((bb - aa) % p, (p - 1) // 2, p) == 1), "lexicographic match")
                if rows[a] >> b & 1:
                    counts[v] += 1
        expected = [n * (p - 1) * p ** (exponent - v - 1) // 2 for v in range(exponent)]
        require([counts[v] for v in range(exponent)] == expected, "valuation arc counts")
        levels.append({"p": p, "e": exponent, "n": n, "arcs_by_valuation": expected})
    seed_controls = []
    for exponent in (1, 2):
        n = 3 ** exponent
        rows = recursive_rows(3, exponent)
        h, u, w, q, starts, ends = path_tables(rows)
        require((h, u, w) == literal_cover_counts(rows), "independent path-cover counts")
        if exponent == 1:
            group = [tuple((x + c) % 3 for x in range(3)) for c in range(3)]
        else:
            group = [tuple((x % 3 + outer) % 3 + 3 * ((x // 3 + shifts[x % 3]) % 3) for x in range(9))
                     for outer in range(3) for shifts in product(range(3), repeat=3)]
        for g in group:
            require(all(sum(1 << g[j] for j in range(n) if rows[i] >> j & 1) == rows[g[i]] for i in range(n)), "wreath automorphism")
        free, by_size, offsets, extensions = 0, Counter(), {}, Counter()
        for mask in range(1 << n):
            selected = [v for v in range(n) if mask >> v & 1]
            stabilizer = sum(sum(1 << g[v] for v in selected) == mask for g in group)
            offset = h - sum(w[i][j] for i, j in combinations(selected, 2))
            if exponent >= 2:
                require(offset % 3 == 0, "all-subset offset obstruction")
            if stabilizer == 1:
                free += 1
                by_size[len(selected)] += 1
                offsets.setdefault(str(len(selected)), set()).add(offset)
            ext = sum(starts[i] if mask >> i & 1 else ends[i] for i in range(n))
            ext += sum(q[i][j] for i in range(n) for j in range(n) if not mask >> i & 1 and mask >> j & 1)
            extensions[ext] += 1
            if exponent == 2:
                require(ext % 3 == 0, "all one-vertex extensions")
        expected_group = 3 ** ((n - 1) // 2)
        require(len(group) == expected_group and free == 2 * expected_group, "free subset count")
        require(dict(by_size) == {(n - 1) // 2: expected_group, (n + 1) // 2: expected_group}, "free subset sizes")
        if exponent == 2:
            # A second end-to-end method: direct subset-DP on representative actual extensions.
            for mask in (0, 1, 7, 15, 31, 85, 170, 255, 511):
                child = [rows[i] | (0 if mask >> i & 1 else 1 << n) for i in range(n)] + [mask]
                count = path_tables(child)[0]
                expected = sum(starts[i] if mask >> i & 1 else ends[i] for i in range(n))
                expected += sum(q[i][j] for i in range(n) for j in range(n) if not mask >> i & 1 and mask >> j & 1)
                require(count == expected, "actual child replay")
        seed_controls.append({"n": n, "H": h, "uniform_U": sorted(set(u)),
                              "pair_W_values": sorted({w[i][j] for i, j in combinations(range(n), 2)}),
                              "wreath_group_order": len(group), "free_subsets": free,
                              "free_subsets_by_size": dict(by_size),
                              "free_subset_offsets": {k: sorted(v) for k, v in offsets.items()},
                              "one_vertex_extension_histogram": dict(sorted(extensions.items()))})
    ring27 = recursive_rows(3, 3)
    common = sorted({(ring27[a] & ring27[b]).bit_count() for a, b in combinations(range(27), 2)})
    field = field27_rows_and_fixers()
    require(common != field["common_outneighbor_counts"], "ring vs field nonisomorphism")
    base = recursive_rows(3, 2)
    asymmetric_child = [base[i] | (0 if 15 >> i & 1 else 1 << 9) for i in range(9)] + [15]
    witness_h = path_tables(asymmetric_child)[0]
    witness_aut = automorphism_count(asymmetric_child)
    require(witness_h == 15033 and witness_aut == 1, "asymmetric divisible child")
    return {"Jacobi_moduli": "3 mod4 through299", "unit_clique_controls": [15, 27],
            "recursive_families": levels, "seed_controls": seed_controls,
            "ring27_common_outneighbor_counts": common, "field27": field,
            "asymmetric_child": {"n": 10, "rows": asymmetric_child, "new_out_mask": 15,
                                 "H": witness_h, "automorphism_group_order": witness_aut}}
