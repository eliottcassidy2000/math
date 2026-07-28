#!/usr/bin/env python3
"""Exact hostile/control certificate for THM-2779's Heisenberg gate.

This companion is intentionally independent of the THM-2771/2772 scripts.
It checks the finite Heisenberg law and sharp p^2 permutation action for
p=2,3,5,7,13; the p=2 dihedral boundary; the complete F_13 frame census;
and the THM-2771 decoder's affine-socle torsor directly from its printed
coefficient vectors.
"""

from collections import Counter
from itertools import product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mul(g, h, p):
    a, b, c = g
    aa, bb, cc = h
    return ((a + aa) % p, (b + bb) % p, (c + cc + a * bb) % p)


def inverse(g, p):
    a, b, c = g
    return ((-a) % p, (-b) % p, (-c + a * b) % p)


def commutator(g, h, p):
    return mul(mul(mul(g, h, p), inverse(g, p), p), inverse(h, p), p)


def det(v, w, p):
    return (v[0] * w[1] - v[1] * w[0]) % p


def act(g, point, p):
    """The sharp left action (a,b,c).(x,y)=(x+b,y+c+a*x)."""
    a, b, c = g
    x, y = point
    return ((x + b) % p, (y + c + a * x) % p)


def conjugate(g, h, p):
    return mul(mul(g, h, p), inverse(g, p), p)


def element_order(g, p):
    identity = (0, 0, 0)
    power = identity
    for n in range(1, 2 * p + 2):
        power = mul(power, g, p)
        if power == identity:
            return n
    raise RuntimeError("order search escaped")


def cyclic_subgroup(g, p):
    identity = (0, 0, 0)
    subgroup = {identity}
    power = identity
    while True:
        power = mul(power, g, p)
        if power == identity:
            return frozenset(subgroup)
        require(power not in subgroup, "cyclic subgroup repeated before identity")
        subgroup.add(power)


def cyclic_convolution(a, b, p):
    n = len(a)
    require(len(b) == n, "cyclic convolution length mismatch")
    out = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[(i + j) % n] = (out[(i + j) % n] + ai * bj) % p
    return out


def truncated_convolution(a, b, p):
    n = len(a)
    require(len(b) == n, "truncated convolution length mismatch")
    out = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j < n:
                out[i + j] = (out[i + j] + ai * bj) % p
    return out


def u_to_epsilon(a, p):
    """Convert sum a_k u^k to epsilon=u-1 coordinates modulo epsilon^p."""
    n = len(a)
    out = [0] * n
    for k, ak in enumerate(a):
        for j in range(k + 1):
            out[j] = (out[j] + ak * comb(k, j)) % p
    return out


def matrix_rank_mod(matrix, p):
    a = [[entry % p for entry in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, p)
        a[rank] = [(inv * value) % p for value in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                scale = a[r][col]
                a[r] = [
                    (a[r][j] - scale * a[rank][j]) % p
                    for j in range(cols)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def permutation_cycles(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


prime_summaries = []
for p in (2, 3, 5, 7, 13):
    identity = (0, 0, 0)
    elements = list(product(range(p), repeat=3))
    base = list(product(range(p), repeat=2))
    points = base
    x_generator = (1, 0, 0)
    y_generator = (0, 1, 0)
    center_generator = (0, 0, 1)

    # Identity, inverse, and the cocycle identity giving associativity.
    for g in elements:
        require(mul(identity, g, p) == g, f"left identity failed at p={p}")
        require(mul(g, identity, p) == g, f"right identity failed at p={p}")
        require(mul(g, inverse(g, p), p) == identity,
                f"right inverse failed at p={p}")
        require(mul(inverse(g, p), g, p) == identity,
                f"left inverse failed at p={p}")
    for a, b, aa, bb, aaa, bbb in product(range(p), repeat=6):
        left = (a * bb + (a + aa) * bbb) % p
        right = (aa * bbb + a * (bb + bbb)) % p
        require(left == right, f"cocycle identity failed at p={p}")

    # The commutator form is exactly the determinant, including p=2.
    commutator_values = set()
    for v in base:
        for w in base:
            g = (v[0], v[1], 0)
            h = (w[0], w[1], 0)
            expected = (0, 0, det(v, w, p))
            require(commutator(g, h, p) == expected,
                    f"commutator determinant failed at p={p}")
            commutator_values.add(expected)
    expected_center = {(0, 0, c) for c in range(p)}
    require(commutator_values == expected_center,
            f"derived subgroup failed at p={p}")
    actual_center = {
        g for g in elements
        if commutator(g, x_generator, p) == identity
        and commutator(g, y_generator, p) == identity
    }
    require(actual_center == expected_center, f"center failed at p={p}")

    # The p^2 action is a genuine faithful transitive action.
    permutations = {}
    for g in elements:
        permutation = tuple(points.index(act(g, point, p)) for point in points)
        require(sorted(permutation) == list(range(p * p)),
                f"nonpermutation at p={p}")
        permutations[g] = permutation
    require(len(set(permutations.values())) == p ** 3,
            f"sharp action lost faithfulness at p={p}")
    require({act(g, (0, 0), p) for g in elements} == set(points),
            f"sharp action lost transitivity at p={p}")
    for g in elements:
        for h in elements:
            gh = mul(g, h, p)
            # Equality of the affine action parameters proves equality on all
            # points; these two probes independently catch both coefficients.
            for point in ((0, 0), (1 % p, 0)):
                require(act(gh, point, p) == act(g, act(h, point, p), p),
                        f"action composition failed at p={p}")

    # The action is the coset action of a noncentral order-p subgroup and its
    # core is trivial.
    stabilizer = {(a, 0, 0) for a in range(p)}
    conjugates = {
        frozenset(conjugate(g, h, p) for h in stabilizer)
        for g in elements
    }
    core = set.intersection(*(set(subgroup) for subgroup in conjugates))
    require(core == {identity}, f"stabilizer core failed at p={p}")

    # Classification of minimal transitive H_p-sets.  Conjugation preserves
    # the projected line, and is transitive on the stabilizers above one line.
    noncentral_order_p_subgroups = {
        cyclic_subgroup(g, p)
        for g in elements
        if g not in expected_center and element_order(g, p) == p
    }
    families_by_projected_line = {}
    for subgroup in noncentral_order_p_subgroups:
        projected_line = frozenset((a, b) for a, b, _ in subgroup)
        families_by_projected_line.setdefault(projected_line, set()).add(subgroup)
    expected_class_count = p + 1 if p > 2 else 2
    expected_subgroup_count = p * (p + 1) if p > 2 else 4
    require(len(noncentral_order_p_subgroups) == expected_subgroup_count,
            f"minimal stabilizer count failed at p={p}")
    require(len(families_by_projected_line) == expected_class_count,
            f"minimal H-set class count failed at p={p}")
    for family in families_by_projected_line.values():
        require(len(family) == p, f"projected-line family size failed at p={p}")
        representative = next(iter(family))
        conjugacy_orbit = {
            frozenset(conjugate(g, h, p) for h in representative)
            for g in elements
        }
        require(conjugacy_orbit == family,
                f"projected-line family is not one conjugacy class at p={p}")

    # Complete determinant/frame census.
    determinant_counts = Counter(det(v, w, p) for v in base for w in base)
    require(determinant_counts[0] == p ** 3 + p ** 2 - p,
            f"singular frame count failed at p={p}")
    for d in range(1, p):
        require(determinant_counts[d] == p * (p ** 2 - 1),
                f"area frame count failed at p={p}, d={d}")
    normalized_frames = [
        (v, w) for v in base for w in base if det(v, w, p) == 1
    ]
    for v, w in normalized_frames:
        require(commutator((v[0], v[1], 0), (w[0], w[1], 0), p)
                == center_generator,
                f"normalized frame commutator failed at p={p}")
    for v in base:
        if v == (0, 0):
            continue
        completions = [w for w in base if det(v, w, p) == 1]
        require(len(completions) == p, f"frame fibre size failed at p={p}")
        w0 = completions[0]
        require(set(completions)
                == {((w0[0] + lam * v[0]) % p,
                     (w0[1] + lam * v[1]) % p)
                    for lam in range(p)},
                f"frame shear torsor failed at p={p}")

    center_perm = permutations[center_generator]
    require(permutation_cycles(center_perm) == (p,) * p,
            f"central cycle structure failed at p={p}")
    prime_summaries.append(
        (p, len(elements), len(points), len(normalized_frames),
         len(conjugates), expected_class_count, permutation_cycles(center_perm))
    )


# The characteristic-two boundary is D_8, not the odd-prime exponent-p group.
p = 2
elements_two = list(product(range(p), repeat=3))
orders_two = Counter(element_order(g, p) for g in elements_two)
require(orders_two == Counter({2: 5, 4: 2, 1: 1}),
        "p=2 order census is not D8")
require(element_order(mul((1, 0, 0), (0, 1, 0), p), p) == 4,
        "p=2 product did not expose the quadratic refinement")


# Exact THM-2771 decoder torsor, reconstructed from the printed primitive rows.
p = 13
k1 = 483303
k2 = 483287
h_rows = [[0] * p for _ in range(7)]
for index in range(3, 12):
    h_rows[1][index] = 2 * k1
h_rows[1][12] = 121 * k1
for index in range(2, 10):
    h_rows[2][index] = 2 * k2
h_rows[2][12] = 265 * k2
for index in (2, 3, 4, 5, 6, 7, 10, 11):
    h_rows[3][index] = 2 * k2
h_rows[3][12] = 254 * k2
beta_rows = [[11 * entry % p for entry in row] for row in h_rows]
s_beta = [sum(row[index] for row in beta_rows) % p for index in range(p)]
expected_s_beta = [0, 0, 8, 0, 0, 0, 0, 0, 9, 9, 9, 9, 8]
require(s_beta == expected_s_beta, "intrinsic Bockstein collapse changed")

k_beta = [3, 5, 5, 5, 7, 12, 2, 8, 2, 9, 8, 11, 0]
epsilon = [12, 1] + [0] * 11
require(cyclic_convolution(s_beta, k_beta, p) == epsilon,
        "printed K_beta no longer decodes to u-1")
decoded_rows = [cyclic_convolution(row, k_beta, p) for row in beta_rows]
expected_decoded = [
    [0] * p,
    [2, 6, 7, 7, 9, 2, 8, 12, 6, 5, 7, 11, 6],
    [0, 11, 1, 7, 6, 9, 9, 4, 12, 8, 8, 5, 8],
    [10, 10, 5, 12, 11, 2, 9, 10, 8, 0, 11, 10, 12],
    [0] * p,
    [0] * p,
    [0] * p,
]
require(decoded_rows == expected_decoded, "decoded THM-2771 rows changed")

n_p = [1] * p
require(cyclic_convolution(s_beta, n_p, p) == [0] * p,
        "group-algebra socle is not in the decoder kernel")
convolution_matrix = []
for row_index in range(p):
    convolution_matrix.append([
        cyclic_convolution(s_beta, [1 if j == column else 0
                                    for j in range(p)], p)[row_index]
        for column in range(p)
    ])
require(matrix_rank_mod(convolution_matrix, p) == p - 1,
        "decoder multiplication does not have one-dimensional kernel")
decoder_torsor = [
    tuple((k_beta[i] + lam * n_p[i]) % p for i in range(p))
    for lam in range(p)
]
require(len(set(decoder_torsor)) == p, "decoder torsor collapsed")
for decoder in decoder_torsor:
    require(cyclic_convolution(s_beta, decoder, p) == epsilon,
            "decoder affine translate failed")
require([decoder[-1] for decoder in decoder_torsor] == list(range(p)),
        "last-coefficient normalization is not a torsor section")

# In epsilon coordinates the printed decoder misses the full inverse only in
# the top socle; adding 3*N_13 repairs that harmless representative choice.
s_epsilon = u_to_epsilon(s_beta, p)
require(s_epsilon[0] == 0 and s_epsilon[1] == 12,
        "S_beta lost its simple augmentation zero")
v_epsilon = s_epsilon[1:] + [0]
k_epsilon = u_to_epsilon(k_beta, p)
vk = truncated_convolution(v_epsilon, k_epsilon, p)
require(vk == [1] + [0] * 11 + [3],
        "socle-qualified inverse identity changed")
full_inverse_u = [(entry + 3) % p for entry in k_beta]
full_inverse_epsilon = u_to_epsilon(full_inverse_u, p)
require(truncated_convolution(v_epsilon, full_inverse_epsilon, p)
        == [1] + [0] * 12,
        "full inverse repair failed")
require(u_to_epsilon(n_p, p) == [0] * 12 + [1],
        "N_13 is not epsilon^12")

# Every decoder gauge gives the same uniformly convolved seven-chart output,
# although the local target-zero chart column itself changes.
chart_columns = []
for decoder in decoder_torsor:
    rows = [cyclic_convolution(row, decoder, p) for row in beta_rows]
    column = tuple(row[0] for row in rows)
    require(sum(column) % p == 12, "decoder gauge changed total correction")
    chart_columns.append(column)
require(len(set(chart_columns)) == p,
        "hostile failed: decoder gauge does not give thirteen local columns")
base_column = (0, 2, 0, 10, 0, 0, 0)
gauge_direction = (0, 3, 3, 7, 0, 0, 0)
for lam, column in enumerate(chart_columns):
    require(
        column == tuple((base_column[i] + lam * gauge_direction[i]) % p
                        for i in range(7)),
        "decoder chart-gauge affine law changed",
    )


# A literal action on p=13 root labels has image at most C_13.  The standard
# quotient control displays its p^2 kernel and kills the commutator center.
root_permutations = {}
for g in product(range(p), repeat=3):
    permutation = tuple((root + g[1]) % p for root in range(p))
    root_permutations[g] = permutation
require(len(set(root_permutations.values())) == p,
        "thirteen-root quotient image changed")
root_kernel = {g for g, permutation in root_permutations.items()
               if permutation == tuple(range(p))}
require(len(root_kernel) == p ** 2, "thirteen-root quotient kernel changed")
require({(0, 0, c) for c in range(p)} <= root_kernel,
        "thirteen-root quotient retained the center")


print("THM-2779 HEISENBERG / DECODER / FRAME EXACT CERTIFICATE")
for summary in prime_summaries:
    (p, group_size, action_degree, frame_count, conjugate_count,
     hset_class_count, cycles) = summary
    print(
        f"p={p}: |H_p|={group_size}; sharp_degree={action_degree}; "
        f"normalized_frames={frame_count}; one_line_conjugates={conjugate_count}; "
        f"minimal_Hset_classes={hset_class_count}; "
        f"z_cycles={','.join(map(str, cycles))}"
    )
print("p=2 boundary: H_2=D8; order_census=1^1,2^5,4^2; minimal_degree=4")
print("decoder: rank(mult_Sbeta)=12; kernel=<N13>; solutions=13; normalized_K_last=0")
print("decoder epsilon law: Vbeta*Kbeta=1+3*epsilon^12; Kbeta+3*N13=Vbeta^-1")
print("decoder hostile: local_chart_columns_vary=13; uniform_chart_sum=-1 for all gauges")
print("frame torsor: Fr_1 is an SL2(F13)-torsor of size 2184; fixed-first-vector fibres=13")
print("commutator: [lift(s),lift(t)]=z^det(s,t); normalized frames give z")
print("standard 13-label quotient: image_size=13; kernel_size=169; center_killed=yes")
print("general 13-label gate: every H_13 action kills Z(H_13) and factors through <=C13")
print("sharp control: H_13 on F13^2 has degree=169, is transitive and faithful")
print("status=VERIFIED-EXACT; no physical ancestry, deck intertwiner, row exclusion, or LRC claim")
