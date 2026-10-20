"""Exact 105-carrier incidence, seven-port blindspace, and 17-spread controls.

No floating arithmetic and no assertions disabled by python -O. The finite
universe is all matchings on eight labelled points, all 3+2+2 junctions on
seven, and all 91 odd ternary colour orbits. Spread tests construct one
F16 spread and all 255 single symplectic transvections of that spread;
they are not a census of all symplectic spreads.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import json

import sympy as sp

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for k, b in enumerate(items[1:], 1):
        for tail in matchings(items[1:k] + items[k + 1:]):
            yield ((a, b),) + tail


def matching_key(pairs):
    return tuple(sorted(tuple(sorted(e)) for e in pairs))


def canonical_word(word):
    names = {}
    return tuple(names.setdefault(c, len(names)) for c in word)


def partitions(n, maximum=None):
    if not n:
        yield ()
        return
    maximum = n if maximum is None else min(n, maximum)
    for a in range(maximum, 0, -1):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def z_value(cycle_type):
    answer = 1
    for k, count in Counter(cycle_type).items():
        answer *= k ** count * factorial(count)
    return answer


@lru_cache(None)
def character(shape, cycle_type):
    """Murnaghan--Nakayama, with rim strips checked as actual cell sets."""
    if not cycle_type:
        return int(not shape)
    size = cycle_type[0]
    cells = {(i, j) for i, row in enumerate(shape) for j in range(row)}
    total = 0
    for smaller in partitions(sum(shape) - size):
        if any(row > (shape[i] if i < len(shape) else 0)
               for i, row in enumerate(smaller)):
            continue
        removed = cells - {(i, j) for i, row in enumerate(smaller)
                           for j in range(row)}
        if len(removed) != size:
            continue
        if any({(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)}
               <= removed for i, j in removed):
            continue
        seen = {next(iter(removed))}
        todo = list(seen)
        while todo:
            i, j = todo.pop()
            for other in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                if other in removed and other not in seen:
                    seen.add(other)
                    todo.append(other)
        if seen != removed:
            continue
        total += (-1) ** (len({i for i, j in removed}) - 1) * character(
            smaller, cycle_type[1:])
    return total


def permutation_of_type(cycle_type):
    answer = list(range(sum(cycle_type)))
    first = 0
    for size in cycle_type:
        for i in range(size):
            answer[first + i] = first + (i + 1) % size
        first += size
    return tuple(answer)


def moved_word(word, permutation):
    answer = [0] * len(word)
    for i, j in enumerate(permutation):
        answer[j] = word[i]
    return canonical_word(answer)


def decomposition(characters, types, n=7):
    answer = {}
    for shape in partitions(n):
        multiplicity = sum(sp.Rational(characters[t] * character(shape, t),
                                       z_value(t)) for t in types)
        need(multiplicity.is_Integer and multiplicity >= 0, "character multiplicity")
        if multiplicity:
            answer[str(shape)] = int(multiplicity)
    return answer


def modular_pivots(matrix, prime=101):
    rows = [[int(x) % prime for x in row] for row in matrix.tolist()]
    position, pivots = 0, []
    for column in range(matrix.cols):
        pivot = next((i for i in range(position, matrix.rows) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[position], rows[pivot] = rows[pivot], rows[position]
        inverse = pow(rows[position][column], -1, prime)
        rows[position] = [(x * inverse) % prime for x in rows[position]]
        for i in range(position + 1, matrix.rows):
            scale = rows[i][column]
            if scale:
                rows[i] = [(x - scale * y) % prime
                           for x, y in zip(rows[i], rows[position])]
        position += 1
        pivots.append(column)
    return pivots


def polynomial(matrix):
    return str(sp.factor(matrix.charpoly().as_expr()))


def psd_rank(matrix):
    """A real symmetric matrix is diagonalizable, so zero root multiplicity
    equals nullity. All callers pass an integer Gram matrix."""
    coefficients = matrix.charpoly().all_coeffs()
    nullity = 0
    while coefficients[-1 - nullity] == 0:
        nullity += 1
    return matrix.rows - nullity


def eight_gram(matchings8):
    rows = []
    for first in matchings8:
        row = []
        for second in matchings8:
            parent = list(range(8))

            def find(a):
                while parent[a] != a:
                    a = parent[a]
                return a
            for a, b in first + second:
                parent[find(a)] = find(b)
            row.append(3 ** len({find(i) for i in range(8)}))
        rows.append(row)
    return sp.Matrix(rows)


def binary_dot(x, y):
    return (x & y).bit_count() % 2


def field_multiply(a, b):
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19                    # X^4+X+1
    return answer


def field_trace(a):
    answer = a
    for _ in range(3):
        a = field_multiply(a, a)
        answer ^= a
    need(answer in (0, 1), "F16 trace")
    return answer


def make_spread():
    remaining = [(1 << i) | (1 << 8) for i in range(8)]
    e_basis, f_basis = [], []
    while remaining:
        e = remaining.pop(0)
        index = next(i for i, x in enumerate(remaining) if binary_dot(e, x))
        f = remaining.pop(index)
        e_basis.append(e)
        f_basis.append(f)
        remaining = [x ^ (f if binary_dot(e, x) else 0)
                     ^ (e if binary_dot(f, x) else 0) for x in remaining]

    def vector(x, y):
        answer = 0
        for i in range(4):
            if x >> i & 1:
                answer ^= e_basis[i]
            if field_trace(field_multiply(1 << i, y)):
                answer ^= f_basis[i]
        return answer

    need(len({vector(x, y) for x in range(16) for y in range(16)}) == 256,
         "symplectic coordinate bijection")
    spread = [tuple(vector(x, field_multiply(t, x)) for x in range(16))
              for t in range(16)]
    spread.append(tuple(vector(0, x) for x in range(16)))
    return spread, e_basis, f_basis


def audit_spread(spread):
    need(len(spread) == 17 and all(len(set(l)) == 16 for l in spread), "spread sizes")
    nonzero = [x for l in spread for x in l if x]
    need(len(nonzero) == len(set(nonzero)) == 255, "spread partitions nonzero vectors")
    codes = []
    matching_codes, pairs = [], []
    for l in spread:
        need(all(binary_dot(x, y) == 0 for x in l for y in l), "isotropy")
        need(all(x ^ y in l for x in l for y in l), "linear Lagrangian")
        code = sorted(set(l) | {x ^ 1023 for x in l})
        need(len(code) == 32, "lifted self-dual code size")
        duads = sorted(tuple(i for i in range(10) if x >> i & 1)
                       for x in code if x.bit_count() == 2)
        need(len(duads) in (1, 5), "inherited length-ten classification")
        if len(duads) == 5:
            matching_codes.append(duads)
        else:
            pairs.append(duads[0])
        codes.append(code)
    need(len(matching_codes) == 7 and len(pairs) == 10, "17=7+10")
    for first, second in combinations(matching_codes, 2):
        reached = {0}
        while True:
            extended = reached | {v for edge in first + second
                                  if reached.intersection(edge) for v in edge}
            if extended == reached:
                break
            reached = extended
        need(len(reached) == 10, "any two matching factors form a Hamiltonian cycle")
    adjacency = [[] for _ in range(10)]
    for a, b in pairs:
        adjacency[a].append(b)
        adjacency[b].append(a)
    need(all(len(a) == 2 for a in adjacency), "residual is a 2-factor")
    need(len(set(pairs) | {e for m in matching_codes for e in m}) == 45,
         "duad partition")
    seen, cycles = set(), []
    for i in range(10):
        if i in seen:
            continue
        cycle, previous, current = [], -1, i
        while current not in seen:
            cycle.append(current)
            seen.add(current)
            nxt = next(x for x in adjacency[current] if x != previous)
            previous, current = current, nxt
        need(current == i, "closed residual cycle")
        cycles.append(cycle)
    cubic_colourings = []
    if all(len(cycle) % 2 == 0 for cycle in cycles):
        for matching in matching_codes:
            colouring = [(a, b, 2) for a, b in matching]
            for cycle in cycles:
                colouring.extend((v, cycle[(j + 1) % len(cycle)], j % 2)
                                  for j, v in enumerate(cycle))
            need(len({tuple(sorted((a, b))) for a, b, colour in colouring}) == 15,
                 "simple cubic completion")
            need(all({colour for a, b, colour in colouring if v in (a, b)} == {0, 1, 2}
                     for v in range(10)), "explicit proper three-edge-colouring")
            cubic_colourings.append(colouring)
    return dict(codes=codes, matching_codes=matching_codes,
                matching_pairs_with_Hamiltonian_union=21,
                hamming_distinguished_pairs=pairs, residual_cycles=cycles,
                cycle_type=sorted(map(len, cycles)),
                cubic_completion_colourings=cubic_colourings)


def main():
    near = [(i, m) for i in range(7)
            for m in matchings(tuple(x for x in range(7) if x != i))]
    full = [matching_key(m + ((i, 7),)) for i, m in near]
    junctions = [(t, m) for t in combinations(range(7), 3)
                 for m in matchings(tuple(x for x in range(7) if x not in t))]
    words = sorted({canonical_word(w) for w in product(range(3), repeat=7)
                    if all(w.count(c) % 2 for c in range(3))})
    need((len(near), len(junctions), len(words)) == (105, 105, 91), "universes")
    r = sp.Matrix([[int(set(pairs2) <= set(pairs3))
                    for triple, pairs2 in junctions] for i, pairs3 in near])
    a = sp.Matrix([[int(len({w[i] for i in triple}) == 3 and
                        all(w[i] == w[j] for i, j in pairs))
                    for w in words] for triple, pairs in junctions])
    need(set(map(sum, r.tolist())) == set(map(sum, r.T.tolist())) == {3},
         "three-to-three incidence")
    full_index = {m: i for i, m in enumerate(full)}
    star_sum = sp.zeros(105)
    for point in range(7):
        perm = list(range(8))
        perm[point], perm[7] = perm[7], perm[point]
        for j, m in enumerate(full):
            moved = matching_key(tuple((perm[x], perm[y]) for x, y in m))
            star_sum[full_index[moved], j] += 1
    need(r * r.T == 2 * sp.eye(105) + star_sum, "rooted transposition identity")
    determinants = []
    for left in combinations(range(8), 4):
        if 0 not in left:
            continue
        right = tuple(i for i in range(8) if i not in left)
        relation = [0] * 105
        for p in permutations(range(4)):
            inversions = sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
            matching = matching_key(tuple((left[i], right[p[i]]) for i in range(4)))
            relation[full_index[matching]] = (-1) ** inversions
        determinants.append(relation)
    d = sp.Matrix(determinants)
    gram8 = eight_gram(full)
    need(psd_rank(d * d.T) == 14 and d * r == sp.zeros(35, 105),
         "determinant incidence kernel")
    need(d * gram8 == sp.zeros(35, 105), "determinant colour kernel")
    need(psd_rank(r * r.T) == psd_rank(gram8) == 91,
         "same complete 14-dimensional kernel")
    w = sp.Matrix([[(1 if sorted(Counter(word).values()) == [1, 1, 5] else -1)
                    if word.count(word[i]) == 1 else 0 for word in words]
                   for i in range(7)])
    need(w * w.T == 15 * sp.eye(7) + sp.ones(7), "Kempe Gram")
    need(a * w.T == sp.zeros(105, 7), "seven-port blind vectors")
    edges = list(combinations(range(7), 2))
    vertex_incidence = sp.Matrix([[int(i in edge) for edge in edges] for i in range(7)])
    edge_basis = vertex_incidence.nullspace()
    need(len(edge_basis) == 14, "zero-row-sum edge space")
    extra = []
    for coefficient in edge_basis:
        edge_value = dict(zip(edges, coefficient))
        row = []
        for word in words:
            blocks = [tuple(i for i in range(7) if word[i] == c) for c in range(3)]
            singletons = tuple(sorted(b[0] for b in blocks if len(b) == 1))
            row.append(-8 * edge_value[singletons] if len(singletons) == 2 else
                       sum(edge_value[e] for block in blocks if len(block) == 3
                           for e in combinations(block, 2)))
        extra.append(row)
    u = sp.Matrix(extra)
    ra = r * a
    need(ra * u.T == sp.zeros(105, 14), "14 additional incidence losses")
    blind = w.col_join(u)
    need(psd_rank(blind * blind.T) == 21 and psd_rank(ra * ra.T) == 70,
         "full composed blindspace")
    need(psd_rank(a * a.T) == 84, "junction rank")
    need(set(ra) == {0, 1} and set(map(sum, ra.tolist())) == {27}, "union kernels")
    root_incidence = sp.Matrix([[int(i == free) for free, pairs in near]
                                for i in range(7)])
    triple_incidence = sp.Matrix([[int(i in triple) for triple, pairs in junctions]
                                  for i in range(7)])
    need(root_incidence * r == triple_incidence, "root-to-junction triple map")
    visible_seven = triple_incidence * a
    need(psd_rank(visible_seven * visible_seven.T) == 7
         and visible_seven * w.T == sp.zeros(7),
         "visible and Kempe seven-dimensional copies are distinct")
    for i in range(7):
        for j, word in enumerate(words):
            factor = 12 if sorted(Counter(word).values()) == [1, 1, 5] else 6
            need(visible_seven[i, j] == 3 + factor * int(word.count(word[i]) == 1),
                 "visible singleton formula")
    types = tuple(partitions(7))
    need(all(sum(sp.Rational(character(l, t) * character(k, t), z_value(t))
                 for t in types) == int(l == k)
             for l in types for k in types), "independent character orthogonality")
    xchars, jchars, vchars, pchars, echars = {}, {}, {}, {}, {}
    for t in types:
        p = permutation_of_type(t)
        xchars[t] = sum((p[i], matching_key(tuple((p[x], p[y]) for x, y in m))) == (i, m)
                        for i, m in near)
        jchars[t] = sum((tuple(sorted(p[i] for i in triple)),
                         matching_key(tuple((p[x], p[y]) for x, y in m))) == (triple, m)
                        for triple, m in junctions)
        vchars[t] = sum(moved_word(word, p) == word for word in words)
        pchars[t] = sum(p[i] == i for i in range(7))
        echars[t] = sum(tuple(sorted(p[i] for i in edge)) == edge for edge in edges)
    need((xchars[(2, 1, 1, 1, 1, 1)], jchars[(2, 1, 1, 1, 1, 1)]) == (15, 25),
         "non-equivariant raw carriers")
    image_chars = {t: vchars[t] - pchars[t] for t in types}
    composed_chars = {t: vchars[t] - echars[t] for t in types}
    modules = {
        "matching105_restricted_S7": decomposition(xchars, types),
        "junction105": decomposition(jchars, types),
        "odd_boundary91": decomposition(vchars, types),
        "junction_image84": decomposition(image_chars, types),
        "junction_relations21": decomposition({t: jchars[t] - image_chars[t] for t in types}, types),
        "composed_image70": decomposition(composed_chars, types),
    }
    independent_modular_ranks = {
        label: len(modular_pivots(matrix)) for label, matrix in
        (("R", r), ("A", a), ("RA", ra), ("determinant_relations", d))
    }
    need(independent_modular_ranks == dict(R=91, A=84, RA=70, determinant_relations=14),
         "independent modular lower bounds match explicit rational upper bounds")
    # Authenticate the S8 type of the determinant kernel, including the two
    # distinct 14-dimensional constituents of the full matching module.
    determinant_basis = d[list(modular_pivots(d.T)), :]
    pivot_columns = modular_pivots(determinant_basis)
    square_inverse = determinant_basis[:, pivot_columns].inv()
    matching8_chars, kernel8_chars = {}, {}
    types8 = tuple(partitions(8))
    need(all(sum(sp.Rational(character(l, t) * character(k, t), z_value(t))
                 for t in types8) == int(l == k)
             for l in types8 for k in types8), "S8 character orthogonality")
    for t in types8:
        p = permutation_of_type(t)
        action = [full_index[matching_key(tuple((p[x], p[y]) for x, y in m))]
                  for m in full]
        matching8_chars[t] = sum(i == image for i, image in enumerate(action))
        inverse = [0] * 105
        for i, image in enumerate(action):
            inverse[image] = i
        transformed = determinant_basis[:, [inverse[j] for j in pivot_columns]]
        trace = sp.trace(transformed * square_inverse)
        need(trace.is_Integer, "determinant representation integral character")
        kernel8_chars[t] = int(trace)
    modules8 = {
        "matching105": decomposition(matching8_chars, types8, 8),
        "determinant_kernel14": decomposition(kernel8_chars, types8, 8),
        "ternary_pairing_image91": decomposition(
            {t: matching8_chars[t] - kernel8_chars[t] for t in types8}, types8, 8),
    }
    spread, e_basis, f_basis = make_spread()
    original_spread = audit_spread(spread)
    alternate_spread = audit_spread([[x ^ (23 if binary_dot(x, 23) else 0)
                                     for x in l] for l in spread])
    need(original_spread["cycle_type"] == [10], "first spread control")
    need(alternate_spread["cycle_type"] == [4, 6], "changed marking hostile")
    need(len(original_spread["cubic_completion_colourings"]) == 7,
         "seven explicitly colourable Hamiltonian cubic completions")
    transvection_histogram = Counter()
    for v in range(1, 512):
        if v.bit_count() % 2:
            continue
        moved = [[x ^ (v if binary_dot(x, v) else 0) for x in l] for l in spread]
        transvection_histogram[tuple(audit_spread(moved)["cycle_type"])] += 1
    need(sum(transvection_histogram.values()) == 255, "single-transvection universe")
    lines = sorted({tuple(sorted((i, j, i ^ j)))
                    for i in range(1, 8) for j in range(i + 1, 8)})
    fano = sp.Matrix([[int(point in line) for point in range(1, 8)] for line in lines])
    need(fano.T * fano == 2 * sp.eye(7) + sp.ones(7), "Fano incidence Gram")
    fano_inverse = sp.Rational(1, 2) * fano.T - sp.Rational(1, 6) * sp.ones(7)
    need(fano_inverse * fano == sp.eye(7), "lossless marked Fano recoding")
    receipt = {
        "universe": {"matchings8": 105, "junctions7": 105, "odd_colour_orbits7": 91,
                     "single_transvections_of_one_F16_spread": 255,
                     "all_spreads_enumerated": False},
        "near_matchings": near, "junctions": junctions,
        "incidence_neighbours": [[j for j in range(105) if r[i, j]] for i in range(105)],
        "spectra": {"R_Rt": polynomial(r * r.T), "A_At": polynomial(a * a.T),
                    "RA_RAt": polynomial(ra * ra.T), "eight_colour_Gram": polynomial(gram8)},
        "ranks": {"R": 91, "eight_colour_Gram": 91, "determinant_relations": 14,
                  "A": 84, "RA": 70, "Kempe_blindspace": 7, "extra_edge_blindspace": 14},
        "raw_transposition_fixed_points": {"matchings": 15, "junctions": 25},
        "S7_modules": modules,
        "S8_modules": modules8,
        "independent_modular_ranks_prime101": independent_modular_ranks,
        "seven_Kempe_Gram": "15 I + J; eigenvalues 22 once, 15 six times",
        "seven_visible_singleton_functions":
            "3+12[singleton i] on (5,1,1); 3+6[singleton i] on (3,3,1); rank7 and orthogonal to every w_i",
        "Fano_lines": lines,
        "spread_coordinate_basis": {"e": e_basis, "f": f_basis, "field_polynomial": 19},
        "original_spread": original_spread, "transvection23_spread": alternate_spread,
        "single_transvection_cycle_types": {str(k): v for k, v in sorted(transvection_histogram.items())},
        "source_sha256_LF": hashlib.sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
    }
    (HERE / "105_17_probe.json").write_text(json.dumps(receipt, indent=2) + "\n",
                                           encoding="utf-8", newline="\n")
    for name in ("spectra", "ranks", "raw_transposition_fixed_points", "S7_modules", "S8_modules",
                 "single_transvection_cycle_types"):
        print(name + ": " + json.dumps(receipt[name], sort_keys=True))
    print("PASS: exact incidence, Gram, kernel, character, Fano, and spread controls.")


if __name__ == "__main__":
    main()
