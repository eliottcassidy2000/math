#!/usr/bin/env python3
"""Independent exact audit for provisional THM-4114.

This file is self-contained: it imports no repository computation, candidate,
or probe.  It distinguishes THM-002's conflict graph from the weighted
inserted-footprint polynomial derived from its odd-cycle packings.
Tournament codes use THM-4097's lexicographic upper-pair, LSB-first convention:
bit (i,j)=1 means i->j for i<j.
"""

from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
import math

from sympy import Matrix, kronecker_product, sqrt
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import GF, ZZ
from sympy.polys.matrices import DomainMatrix


EXPECTED_SEMANTIC = (
    "3917956b00891bdfb012076298867a484f4b0a41950e0e18326a782067131c5e"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bit_index(n, i, j):
    require(0 <= i < j < n, ("bad pair", n, i, j))
    return i * (2 * n - i - 1) // 2 + j - i - 1


def arc(code, n, i, j):
    require(i != j, "loop queried")
    if i < j:
        return bool((code >> bit_index(n, i, j)) & 1)
    return not arc(code, n, j, i)


def submasks(mask):
    current = mask
    while True:
        yield current
        if current == 0:
            break
        current = (current - 1) & mask


def popcount(mask):
    return mask.bit_count()


@lru_cache(maxsize=None)
def induced_hamilton_count(code, n, universe):
    """Hamiltonian-path count of the tournament induced on universe."""
    if universe == 0:
        return 1
    dp = {}
    for v in range(n):
        if (universe >> v) & 1:
            dp[(1 << v, v)] = 1
    for size in range(1, popcount(universe)):
        for mask in (m for m in submasks(universe) if popcount(m) == size):
            for last in range(n):
                count = dp.get((mask, last), 0)
                if not count:
                    continue
                remaining = universe ^ mask
                for nxt in range(n):
                    if ((remaining >> nxt) & 1) and arc(code, n, last, nxt):
                        key = (mask | (1 << nxt), nxt)
                        dp[key] = dp.get(key, 0) + count
    return sum(dp.get((universe, last), 0) for last in range(n))


def directed_odd_cycles(code, n):
    """Permutation-cycle objects, quotiented by cyclic rotation only."""
    cycles = []
    vertices = tuple(range(n))
    for length in range(3, n + 1, 2):
        for support in combinations(vertices, length):
            root = min(support)
            others = tuple(v for v in support if v != root)
            for tail in permutations(others):
                cycle = (root,) + tail
                if all(arc(code, n, cycle[i], cycle[(i + 1) % length])
                       for i in range(length)):
                    cycles.append(cycle)
    return tuple(cycles)


@lru_cache(maxsize=None)
def odd_cycle_packings(code, n):
    cycles = directed_odd_cycles(code, n)
    masks = tuple(sum(1 << v for v in cycle) for cycle in cycles)
    packings = []

    def search(position, used, chosen):
        if position == len(cycles):
            packings.append(tuple(chosen))
            return
        search(position + 1, used, chosen)
        if not (used & masks[position]):
            chosen.append(position)
            search(position + 1, used | masks[position], chosen)
            chosen.pop()

    search(0, 0, [])
    return cycles, masks, tuple(packings)


def mobius_from_values(values, variable_count):
    atoms = list(values)
    for bit in range(variable_count):
        for mask in range(1 << variable_count):
            if (mask >> bit) & 1:
                atoms[mask] -= atoms[mask ^ (1 << bit)]
    return tuple(atoms)


def variable_vertices(mask, vertices):
    return sum(1 << vertex for position, vertex in enumerate(vertices)
               if (mask >> position) & 1)


def packing_atoms(code, n, inserted):
    cycles, cycle_masks, packings = odd_cycle_packings(code, n)
    del cycles
    atoms = [0] * (1 << len(inserted))
    inserted_mask = sum(1 << v for v in inserted)
    position = {vertex: index for index, vertex in enumerate(inserted)}
    for packing in packings:
        used = 0
        for cycle_index in packing:
            used |= cycle_masks[cycle_index]
        footprint_vertices = used & inserted_mask
        footprint = sum(1 << position[v] for v in inserted
                        if (footprint_vertices >> v) & 1)
        atoms[footprint] += 1 << len(packing)
    return tuple(atoms)


def matrix_rank_mod(matrix, prime):
    return DomainMatrix.from_Matrix(Matrix(matrix)).convert_to(GF(prime)).rank()


def snf_diagonal(matrix):
    smith = smith_normal_form(Matrix(matrix), domain=ZZ)
    return tuple(abs(int(smith[i, i]))
                 for i in range(min(smith.rows, smith.cols)))


def zeta_matrix(bits):
    return Matrix([
        [int(column & ~row == 0) for column in range(1 << bits)]
        for row in range(1 << bits)
    ])


def flatten(values, left_positions, right_positions):
    rows = []
    for left in range(1 << len(left_positions)):
        row = []
        for right in range(1 << len(right_positions)):
            mask = sum(((left >> j) & 1) << position
                       for j, position in enumerate(left_positions))
            mask |= sum(((right >> j) & 1) << position
                        for j, position in enumerate(right_positions))
            row.append(values[mask])
        rows.append(row)
    return Matrix(rows)


def child_hamilton_count(code, n, cut):
    """Add x=n, with x->v iff v is in cut, and count by subset DP."""
    size = n + 1

    def child_arc(i, j):
        require(i != j, "child loop")
        if i < n and j < n:
            return arc(code, n, i, j)
        if i == n:
            return bool((cut >> j) & 1)
        return not bool((cut >> i) & 1)

    dp = {(1 << v, v): 1 for v in range(size)}
    full = (1 << size) - 1
    for cardinality in range(1, size):
        for mask in range(1 << size):
            if popcount(mask) != cardinality:
                continue
            for last in range(size):
                count = dp.get((mask, last), 0)
                if not count:
                    continue
                for nxt in range(size):
                    if not ((mask >> nxt) & 1) and child_arc(last, nxt):
                        key = (mask | (1 << nxt), nxt)
                        dp[key] = dp.get(key, 0) + count
    return sum(dp.get((full, last), 0) for last in range(size))


def ear_boundary(code, n):
    starts = [0] * n
    ends = [0] * n
    q = [[0] * n for _ in range(n)]
    for word in permutations(range(n)):
        good = [arc(code, n, word[i], word[i + 1]) for i in range(n - 1)]
        if all(good):
            starts[word[0]] += 1
            ends[word[-1]] += 1
        for position in range(n - 1):
            if all(good[j] for j in range(n - 1) if j != position):
                q[word[position]][word[position + 1]] += 1
    return tuple(starts), tuple(ends), tuple(tuple(row) for row in q)


def ear_cut_value(starts, ends, q, cut):
    n = len(starts)
    return (
        sum(starts[v] for v in range(n) if (cut >> v) & 1)
        + sum(ends[v] for v in range(n) if not ((cut >> v) & 1))
        + sum(q[a][b] for a in range(n) for b in range(n) if a != b
              and not ((cut >> a) & 1) and ((cut >> b) & 1))
    )


def strongly_connected(code, n):
    reach = [[i == j or arc(code, n, i, j) for j in range(n)]
             for i in range(n)]
    for middle in range(n):
        for left in range(n):
            for right in range(n):
                reach[left][right] = (reach[left][right]
                                      or (reach[left][middle]
                                          and reach[middle][right]))
    return all(all(row) for row in reach)


def local_gap_factor(code, n, left, right, inserted):
    """Exact kappa vector for one boundary/internal gap."""
    factor = []
    for mask in range(1 << len(inserted)):
        block = tuple(inserted[position] for position in range(len(inserted))
                      if (mask >> position) & 1)
        count = 0
        for ordering in permutations(block):
            word = ((left,) if left is not None else ()) + ordering \
                + ((right,) if right is not None else ())
            if all(arc(code, n, word[j], word[j + 1])
                   for j in range(len(word) - 1)):
                count += 1
        factor.append(count)
    return tuple(factor)


def squarefree_product(left, right):
    require(len(left) == len(right), "squarefree shape")
    product = [0] * len(left)
    for a, left_value in enumerate(left):
        if not left_value:
            continue
        for b, right_value in enumerate(right):
            if right_value and not (a & b):
                product[a | b] += left_value * right_value
    return tuple(product)


def fixed_word_gap_polynomial(code, n, word, inserted):
    gaps = [(None, word[0])]
    gaps.extend((word[j], word[j + 1]) for j in range(len(word) - 1))
    gaps.append((word[-1], None))
    factors = tuple(local_gap_factor(code, n, left, right, inserted)
                    for left, right in gaps)
    product = (1,) + (0,) * ((1 << len(inserted)) - 1)
    for factor in factors:
        product = squarefree_product(product, factor)
    return factors, product


ledger = {
    "conventions": {
        "conflict_graph": (
            "Omega_conf(W): directed odd permutation-cycles modulo cyclic "
            "rotation, adjacent exactly when their vertex sets meet"
        ),
        "footprint_polynomial": (
            "Phi_W_over_B(Y)=sum_C 2^|C| Y_fp_I(C), with C a set of "
            "pairwise vertex-disjoint vertices of Omega_conf(W); A is not "
            "the conflict graph"
        ),
        "zeta_lift": (
            "Z_W_over_B(X)=sum_D mu(D) X_D product_(i notin D)(1+X_i)"
        ),
        "promotion": (
            "coefficient extraction in the unique multiaffine normal form; "
            "not a derivation or Xi=1 homomorphism of Z[X]/(Xi^2)"
        ),
        "rank_scope": (
            "Smith normal form over Z; ordinary matrix rank over a declared "
            "field; no Boolean, tropical, or nonnegative-rank claim"
        ),
        "boolean_difference": (
            "Delta_A h(S)=sum_(T subset A)(-1)^(|A|-|T|)h(S union T), "
            "with A disjoint from S"
        ),
        "ear_domain": (
            "the cut response is submodular quadratic on the full Boolean "
            "cube; nonconstant strong-ear signatures form a punctured domain"
        ),
    }
}
split_count = 0
mixed_difference_count = 0
deletion_promotion_count = 0
packing_atom_count = 0

# Exhaust every labelled tournament and every nonempty-base split through n=5.
for n in range(1, 6):
    for code in range(1 << (n * (n - 1) // 2)):
        for base in range(1, 1 << n):
            inserted = tuple(v for v in range(n) if not ((base >> v) & 1))
            variable_count = len(inserted)
            values = tuple(
                induced_hamilton_count(
                    code, n, base | variable_vertices(mask, inserted)
                )
                for mask in range(1 << variable_count)
            )
            atoms = mobius_from_values(values, variable_count)
            pack_atoms = packing_atoms(code, n, inserted)
            require(atoms == pack_atoms,
                    ("OCF footprint atom mismatch", n, code, base,
                     atoms, pack_atoms))
            packing_atom_count += len(atoms)
            require(atoms[0] % 2 == 1,
                    ("empty atom should be odd", n, code, base, atoms[0]))
            require(all(atom >= 0 for atom in atoms),
                    ("negative Mobius atom", n, code, base, atoms))
            require(all(atoms[mask] % 2 == 0
                        for mask in range(1, 1 << variable_count)),
                    ("nonempty atom should be even", n, code, base, atoms))

            universe = (1 << variable_count) - 1
            for a in range(1 << variable_count):
                available = universe ^ a
                for s in submasks(available):
                    difference = sum(
                        (-1 if (popcount(a) - popcount(t)) % 2 else 1)
                        * values[s | t]
                        for t in submasks(a)
                    )
                    require(difference >= 0,
                            ("negative mixed difference", n, code, base, a, s))
                    if a:
                        require(difference % 2 == 0,
                                ("odd nonempty mixed difference",
                                 n, code, base, a, s, difference))
                    expected = sum(
                        atoms[d] for d in range(1 << variable_count)
                        if (d & a) == a and (d & ~(s | a)) == 0
                    )
                    require(difference == expected,
                            ("mixed-difference packing mismatch",
                             n, code, base, a, s, difference, expected))
                    mixed_difference_count += 1

            for position in range(variable_count):
                rest = tuple(j for j in range(variable_count) if j != position)
                promoted_base = base | (1 << inserted[position])
                promoted_inserted = tuple(inserted[j] for j in rest)
                for rest_mask in range(1 << len(rest)):
                    old_mask = sum(((rest_mask >> j) & 1) << old
                                   for j, old in enumerate(rest))
                    promoted_mask = old_mask | (1 << position)
                    promoted_value = induced_hamilton_count(
                        code, n,
                        promoted_base
                        | variable_vertices(rest_mask, promoted_inserted)
                    )
                    require(values[promoted_mask] == promoted_value,
                            ("Z promotion mismatch", n, code, base,
                             position, rest_mask))
                    deletion_value = induced_hamilton_count(
                        code, n,
                        base | variable_vertices(rest_mask, promoted_inserted)
                    )
                    require(values[old_mask] == deletion_value,
                            ("Z deletion mismatch", n, code, base,
                             position, rest_mask))
                    deletion_promotion_count += 1
                promoted_values = tuple(
                    induced_hamilton_count(
                        code, n,
                        promoted_base
                        | variable_vertices(mask, promoted_inserted)
                    )
                    for mask in range(1 << len(rest))
                )
                promoted_atoms = mobius_from_values(promoted_values, len(rest))
                for rest_mask in range(1 << len(rest)):
                    old_mask = sum(((rest_mask >> j) & 1) << old
                                   for j, old in enumerate(rest))
                    require(promoted_atoms[rest_mask]
                            == atoms[old_mask] + atoms[old_mask | (1 << position)],
                            ("Phi promotion mismatch", n, code, base,
                             position, rest_mask))

            for left_bits in range(variable_count + 1):
                left = tuple(range(left_bits))
                right = tuple(range(left_bits, variable_count))
                value_matrix = flatten(values, left, right)
                atom_matrix = flatten(atoms, left, right)
                require(value_matrix
                        == zeta_matrix(len(left)) * atom_matrix
                           * zeta_matrix(len(right)).T,
                        ("zeta flattening mismatch", n, code, base, left_bits))
                split_count += 1

ledger["finite_footprint_audit"] = {
    "labelled_tournaments": sum(1 << (n * (n - 1) // 2) for n in range(1, 6)),
    "packing_atoms_checked": packing_atom_count,
    "mixed_differences_checked": mixed_difference_count,
    "deletion_promotion_coefficients_checked": deletion_promotion_count,
    "zeta_flattenings_checked": split_count,
}

c3_rows = []
for k in range(1, 5):
    local = Matrix([[1, 1], [1, 3]])
    matrix = local
    for _ in range(1, k):
        matrix = kronecker_product(matrix, local)
    atom = Matrix.zeros(1 << k)
    for subset in range(1 << k):
        atom[subset, subset] = 1 << popcount(subset)
    require(matrix == zeta_matrix(k) * atom * zeta_matrix(k).T,
            ("C3 zeta factorization", k))
    expected_determinant = 1 << (k * (1 << (k - 1)))
    require(matrix.det() == expected_determinant,
            ("C3 determinant", k, matrix.det(), expected_determinant))
    expected_snf = tuple(
        value for exponent in range(k + 1)
        for value in [1 << exponent] * math.comb(k, exponent)
    )
    require(snf_diagonal(matrix) == expected_snf,
            ("C3 SNF", k, snf_diagonal(matrix), expected_snf))
    value_range = sorted(set(int(entry) for entry in matrix))
    require(value_range == [3 ** exponent for exponent in range(k + 1)],
            ("C3 sparse range", k, value_range))
    c3_rows.append({
        "k": k,
        "shape": [1 << k, 1 << k],
        "range_size": len(value_range),
        "determinant": int(matrix.det()),
        "rank_Q": matrix.rank(),
        "rank_F2": matrix_rank_mod(matrix, 2),
        "rank_F3": matrix_rank_mod(matrix, 3),
        "snf": list(expected_snf),
    })
ledger["C3_ordinal_sum"] = c3_rows

code = 42
n = 4
base = 1 << 0
inserted = (1, 2, 3)
code42_values = tuple(
    induced_hamilton_count(code, n, base | variable_vertices(mask, inserted))
    for mask in range(8)
)
code42_atoms = mobius_from_values(code42_values, 3)
require(code42_values == (1, 1, 1, 1, 1, 1, 3, 5),
        ("code42 values", code42_values))
require(code42_atoms == (1, 0, 0, 0, 0, 0, 2, 2),
        ("code42 atoms", code42_atoms))
ledger["code42_base0"] = {
    "coefficient_order": "masks 000,001,010,011,100,101,110,111 on I=(1,2,3)",
    "H_values": list(code42_values),
    "Mobius_atoms": list(code42_atoms),
    "directed_triangles": [[0, 2, 3], [1, 2, 3]],
}

# The sharper tropical/defect witness uses B={0,1}, I=(2,3).
bad_factors, bad_product = fixed_word_gap_polynomial(
    42, 4, (0, 1), (2, 3)
)
good_factors, good_product = fixed_word_gap_polynomial(
    42, 4, (1, 0), (2, 3)
)
require(bad_factors == ((1, 0, 1, 1), (0, 0, 0, 1), (1, 1, 0, 1)),
        ("code42 bad factors", bad_factors))
require(bad_product == (0, 0, 0, 1),
        ("code42 delayed product", bad_product))
require(good_product == (1, 1, 1, 4),
        ("code42 good product", good_product))
require(tuple(bad_product[j] + good_product[j] for j in range(4))
        == (1, 1, 1, 5), "code42 split reconstruction")
ledger["code42_defect_delay"] = {
    "base": [0, 1],
    "inserted": [2, 3],
    "bad_word": [0, 1],
    "bad_adjacencies": 1,
    "left_internal_right_factors": [list(row) for row in bad_factors],
    "fixed_word_product": list(bad_product),
    "first_live_degree": 2,
    "singleton_repair_row": [0, 0],
    "two_vertex_repair": "0->2->3->1",
}

z_hessian = Matrix([[2, 1, 1], [1, 0, 3], [1, 3, 0]])
atom_hessian = Matrix([[2, 0, 0], [0, 0, 2], [0, 2, 0]])
z_factorial_hessian = Matrix([[1, 1, 1], [1, 0, 3], [1, 3, 0]])
atom_factorial_hessian = Matrix([[1, 0, 0], [0, 0, 2], [0, 2, 0]])
require(z_hessian.eigenvals() == {-3: 1, 1: 1, 4: 1},
        ("Z homogenization signature", z_hessian.eigenvals()))
require(atom_hessian.eigenvals() == {-2: 1, 2: 2},
        ("atom homogenization signature", atom_hessian.eigenvals()))
require(z_factorial_hessian.eigenvals()
        == {-3: 1, 2 - sqrt(3): 1, 2 + sqrt(3): 1},
        ("Z factorial homogenization signature",
         z_factorial_hessian.eigenvals()))
require(atom_factorial_hessian.eigenvals() == {-2: 1, 1: 1, 2: 1},
        ("atom factorial homogenization signature",
         atom_factorial_hessian.eigenvals()))
ledger["stability_hostile"] = {
    "Z": "1+U+V+3UV",
    "Z_upper_half_plane_zero": "U=V=(-1+i*sqrt(2))/3",
    "Z_Rayleigh_bc_minus_ad": -2,
    "Phi": "1+2UV",
    "Phi_upper_half_plane_zero": "U=V=i/sqrt(2)",
    "Phi_Rayleigh_bc_minus_ad": -2,
    "ordinary_homogenization_Hessian_signatures": {
        "Z": "(+,+,-), eigenvalues 4,1,-3",
        "Phi": "(+,+,-), eigenvalues 2,2,-2",
    },
    "factorial_normalized_Hessian_signatures": {
        "Z": "(+,+,-), eigenvalues 2+sqrt(3),2-sqrt(3),-3",
        "Phi": "(+,+,-), eigenvalues 2,1,-2",
    },
}

# The theorem's additional homogeneous-slice hostile: code 10, base {0}.
code10_inserted = (1, 2, 3, 4)
code10_hessian = Matrix([
    [
        0 if left == right else induced_hamilton_count(
            10, 5, (1 << 0) | (1 << left) | (1 << right)
        )
        for right in code10_inserted
    ]
    for left in code10_inserted
])
expected_code10_hessian = Matrix([
    [0, 3, 1, 3],
    [3, 0, 1, 1],
    [1, 1, 0, 3],
    [3, 1, 3, 0],
])
require(code10_hessian == expected_code10_hessian,
        ("code10 degree-two Hessian", code10_hessian))
require(tuple(code10_hessian.charpoly().all_coeffs()) == (1, 0, -30, -48, 13),
        ("code10 characteristic polynomial",
         code10_hessian.charpoly().all_coeffs()))
ledger["stability_hostile"]["code10_degree2_slice"] = {
    "base": [0],
    "inserted": list(code10_inserted),
    "Hessian": [[int(entry) for entry in row]
                for row in code10_hessian.tolist()],
    "characteristic_polynomial": "(t^2-4t-13)(t^2+4t-1)",
    "signature": "(+,+,-,-)",
}

ear_parents = 0
ear_cuts = 0
ear_pair_differences = 0
ear_third_differences = 0
strong_parent_counts = {3: 0, 4: 0, 5: 0}
nonsolid_strong_images = {3: 0, 4: 0, 5: 0}
for n in range(1, 6):
    for code in range(1 << (n * (n - 1) // 2)):
        starts, ends, q = ear_boundary(code, n)
        values = tuple(ear_cut_value(starts, ends, q, cut)
                       for cut in range(1 << n))
        direct = tuple(child_hamilton_count(code, n, cut)
                       for cut in range(1 << n))
        require(values == direct, ("ear cut mismatch", n, code, values, direct))
        ear_parents += 1
        ear_cuts += len(values)
        if n >= 3 and strongly_connected(code, n):
            strong_parent_counts[n] += 1
            image = sorted(set(values[1:-1]))
            solid = image == list(range(image[0], image[-1] + 1, 2))
            if not solid:
                nonsolid_strong_images[n] += 1
        for i, j in combinations(range(n), 2):
            rest = ((1 << n) - 1) ^ (1 << i) ^ (1 << j)
            for base_cut in submasks(rest):
                second = (values[base_cut | (1 << i) | (1 << j)]
                          - values[base_cut | (1 << i)]
                          - values[base_cut | (1 << j)]
                          + values[base_cut])
                require(second == -(q[i][j] + q[j][i]),
                        ("ear second difference", n, code, i, j,
                         base_cut, second, q[i][j], q[j][i]))
                require(second <= 0, ("ear non-submodular", n, code, i, j))
                ear_pair_differences += 1
        for i, j, k in combinations(range(n), 3):
            rest = ((1 << n) - 1) ^ (1 << i) ^ (1 << j) ^ (1 << k)
            toggles = (1 << i) | (1 << j) | (1 << k)
            for base_cut in submasks(rest):
                third = sum(
                    (-1 if (3 - popcount(t)) % 2 else 1)
                    * values[base_cut | t]
                    for t in submasks(toggles)
                )
                require(third == 0,
                        ("ear third difference", n, code, i, j, k,
                         base_cut, third))
                ear_third_differences += 1

require(strong_parent_counts == {3: 2, 4: 24, 5: 544},
        ("strong parent census", strong_parent_counts))
require(nonsolid_strong_images == {3: 0, 4: 0, 5: 544},
        ("nonsolid strong-ear census", nonsolid_strong_images))

ledger["ear_cut_audit"] = {
    "labelled_parents": ear_parents,
    "cut_values_checked": ear_cuts,
    "second_differences_checked": ear_pair_differences,
    "third_differences_checked": ear_third_differences,
    "formula": "Delta_i Delta_j F=-(Q_ij+Q_ji), Delta_i Delta_j Delta_k F=0",
    "strong_parent_counts": strong_parent_counts,
    "nonsolid_strong_images": nonsolid_strong_images,
}

n = 5
code = 8
require(strongly_connected(code, n), "code 8 should be strong")
starts, ends, q = ear_boundary(code, n)
ear_values = tuple(ear_cut_value(starts, ends, q, cut)
                   for cut in range(1 << n))
require(ear_values[0] == ear_values[-1] == 9,
        ("code8 constant cuts", ear_values[0], ear_values[-1]))
nonconstant_image = sorted(set(ear_values[1:-1]))
require(nonconstant_image
        == [15, 17, 19, 23, 25, 27, 29, 33, 37, 41],
        ("code8 image", nonconstant_image))
q_matrix = Matrix(q)
sym_q = q_matrix + q_matrix.T
split_left = (0, 2)
split_right = (1, 3, 4)
response_flattening = flatten(ear_values, split_left, split_right)
response_atoms = mobius_from_values(ear_values, n)
atom_flattening = flatten(response_atoms, split_left, split_right)
require(response_flattening
        == zeta_matrix(2) * atom_flattening * zeta_matrix(3).T,
        "code8 zeta flattening")
require(snf_diagonal(response_flattening) == (1, 2, 2, 8),
        ("code8 flatten SNF", snf_diagonal(response_flattening)))
ledger["strong_order5_code8"] = {
    "code_convention": "lexicographic upper pairs, LSB first, bit 1 means i->j",
    "H_parent": 9,
    "nonconstant_image": nonconstant_image,
    "missing_odds_between_extremes": [
        value for value in range(15, 42, 2) if value not in nonconstant_image
    ],
    "Q_rank_Q": q_matrix.rank(),
    "Q_det": int(q_matrix.det()),
    "Q_snf": list(snf_diagonal(q_matrix)),
    "symmetric_quadratic_rank_Q": sym_q.rank(),
    "symmetric_quadratic_det": int(sym_q.det()),
    "symmetric_quadratic_snf": list(snf_diagonal(sym_q)),
    "flattening_split": [list(split_left), list(split_right)],
    "flattening_shape": [response_flattening.rows, response_flattening.cols],
    "flattening_rank_Q": response_flattening.rank(),
    "flattening_rank_F2": matrix_rank_mod(response_flattening, 2),
    "flattening_rank_F3": matrix_rank_mod(response_flattening, 3),
    "flattening_snf": list(snf_diagonal(response_flattening)),
}

canonical = json.dumps(ledger, sort_keys=True, separators=(",", ":"))
semantic_digest = sha256(canonical.encode("utf-8")).hexdigest()
require(semantic_digest == EXPECTED_SEMANTIC,
        ("semantic drift", semantic_digest, EXPECTED_SEMANTIC))
print("THM-4114 INDEPENDENT OCF-MOBIUS/EAR-CURVATURE AUDIT")
for key in sorted(ledger):
    print(key + "=" + json.dumps(ledger[key], sort_keys=True,
                                 separators=(",", ":")))
print("semantic_sha256=" + semantic_digest)
print("all_exact_checks=PASS")
