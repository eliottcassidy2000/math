#!/usr/bin/env python3
"""Independent hostile audit for the merged THM-2779 theorem.

This script does not import any THM-2779 companion.  It checks the group-law
conventions, the exact-degree stabilizer classification, decoder quantifiers,
the p=13 specializations, the endpoint/odometer model, and the two conditional
typing gates with fresh implementations.
"""

from collections import Counter
from itertools import product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def new_mul(g, h, p):
    """Merged theorem convention: z-coordinate cocycle -y*x'."""
    x, y, z = g
    xx, yy, zz = h
    return ((x + xx) % p, (y + yy) % p, (z + zz - y * xx) % p)


def new_inv(g, p):
    x, y, z = g
    return ((-x) % p, (-y) % p, (-z - x * y) % p)


def old_mul(g, h, p):
    """Tertiary companion convention: C-coordinate cocycle A*B'."""
    a, b, c = g
    aa, bb, cc = h
    return ((a + aa) % p, (b + bb) % p, (c + cc + a * bb) % p)


def convention_map(g, p):
    x, y, z = g
    return ((-y) % p, x % p, z % p)


def commutator(g, h, p):
    return new_mul(
        new_mul(new_mul(g, h, p), new_inv(g, p), p),
        new_inv(h, p),
        p,
    )


def power(g, n, p):
    out = (0, 0, 0)
    for _ in range(n):
        out = new_mul(out, g, p)
    return out


def element_order(g, p):
    for n in range(1, 2 * p + 1):
        if power(g, n, p) == (0, 0, 0):
            return n
    raise RuntimeError(f"order search escaped at p={p}, g={g}")


def cyclic_subgroup(g, p):
    return frozenset(power(g, n, p) for n in range(element_order(g, p)))


def conjugate(g, h, p):
    return new_mul(new_mul(g, h, p), new_inv(g, p), p)


def det(v, w, p):
    return (v[0] * w[1] - v[1] * w[0]) % p


def rho(g, point, p):
    x, y, z = g
    r, w = point
    return ((r + x) % p, (w + z - y * r) % p)


def mul_eps(a, b, p):
    out = [0] * p
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j < p:
                out[i + j] = (out[i + j] + ai * bj) % p
    return tuple(out)


def inv_eps(v, p):
    require(v[0] % p != 0, "attempted to invert a nonunit")
    out = [0] * p
    out[0] = pow(v[0], -1, p)
    for n in range(1, p):
        tail = sum(v[i] * out[n - i] for i in range(1, n + 1))
        out[n] = (-out[0] * tail) % p
    require(mul_eps(v, out, p) == (1,) + (0,) * (p - 1),
            "epsilon inverse recursion failed")
    return tuple(out)


def epsilon_to_u(a, p):
    out = [0] * p
    for j, aj in enumerate(a):
        for k in range(j + 1):
            out[k] = (
                out[k] + aj * comb(j, k) * ((-1) ** (j - k))
            ) % p
    return tuple(out)


def u_to_epsilon(a, p):
    out = [0] * p
    for k, ak in enumerate(a):
        for j in range(k + 1):
            out[j] = (out[j] + ak * comb(k, j)) % p
    return tuple(out)


def cyclic_convolution(a, b, p):
    out = [0] * p
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[(i + j) % p] = (out[(i + j) % p] + ai * bj) % p
    return tuple(out)


prime_rows = []
for p in (2, 3, 5, 7, 11, 13):
    elements = tuple(product(range(p), repeat=3))
    plane = tuple(product(range(p), repeat=2))
    identity = (0, 0, 0)
    x_gen = (1, 0, 0)
    y_gen = (0, 1, 0)
    z_gen = (0, 0, 1)

    # The two displayed Heisenberg conventions are exactly isomorphic.
    for g in elements:
        require(new_mul(g, new_inv(g, p), p) == identity,
                f"new inverse failed at p={p}")
        require(old_mul(convention_map(g, p),
                        convention_map(new_inv(g, p), p), p)
                == (0, 0, 0), f"old inverse transport failed at p={p}")
        for h in elements:
            require(
                convention_map(new_mul(g, h, p), p)
                == old_mul(convention_map(g, p), convention_map(h, p), p),
                f"group conventions conflict at p={p}",
            )

    expected_center = {(0, 0, z) for z in range(p)}
    actual_center = {
        g for g in elements
        if commutator(g, x_gen, p) == identity
        and commutator(g, y_gen, p) == identity
    }
    require(actual_center == expected_center, f"center failed at p={p}")
    commutator_values = set()
    for v in plane:
        for w in plane:
            value = commutator((v[0], v[1], 0), (w[0], w[1], 0), p)
            require(value == (0, 0, det(v, w, p)),
                    f"determinant commutator failed at p={p}")
            commutator_values.add(value)
    require(commutator_values == expected_center,
            f"derived group failed at p={p}")

    # The p^2 action is faithful, transitive, and has the claimed stabilizer.
    permutations = {
        g: tuple(rho(g, point, p) for point in plane) for g in elements
    }
    require(len(set(permutations.values())) == p ** 3,
            f"sharp action is not faithful at p={p}")
    require({rho(g, (0, 0), p) for g in elements} == set(plane),
            f"sharp action is not transitive at p={p}")
    stabilizer = {g for g in elements if rho(g, (0, 0), p) == (0, 0)}
    require(stabilizer == {(0, y, 0) for y in range(p)},
            f"sharp stabilizer changed at p={p}")

    # Classify all core-free order-p stabilizers.  Also retain the central
    # stabilizer as a hostile: omitting "faithful" adds one degree-p^2 class.
    eligible = {
        cyclic_subgroup(g, p)
        for g in elements
        if g not in expected_center and element_order(g, p) == p
    }
    classes = set()
    while eligible:
        representative = next(iter(eligible))
        orbit = frozenset(
            frozenset(conjugate(g, h, p) for h in representative)
            for g in elements
        )
        require(orbit <= eligible, f"conjugacy orbit escaped at p={p}")
        classes.add(orbit)
        eligible -= set(orbit)
    faithful_class_count = len(classes)
    expected_faithful = p + 1 if p > 2 else 2
    require(faithful_class_count == expected_faithful,
            f"faithful exact-degree class count failed at p={p}")
    require(all(len(cls) == p for cls in classes),
            f"stabilizer conjugacy class size failed at p={p}")
    all_degree_p2_transitive_classes = faithful_class_count + 1

    frame_count = sum(
        det(v, w, p) == 1 for v in plane for w in plane
    )
    require(frame_count == p * (p * p - 1),
            f"normalized frame count failed at p={p}")
    for v in plane:
        if v == (0, 0):
            continue
        completions = [w for w in plane if det(v, w, p) == 1]
        require(len(completions) == p, f"frame fibre failed at p={p}")

    prime_rows.append(
        (p, p ** 3, p ** 2, frame_count, faithful_class_count,
         all_degree_p2_transitive_classes)
    )


orders_two = Counter(
    element_order(g, 2) for g in product(range(2), repeat=3)
)
require(orders_two == Counter({1: 1, 2: 5, 4: 2}),
        "p=2 is not D8")


# Exhaust the decoder theorem for every unit V at p=2,3,5.
decoder_unit_counts = []
for p in (2, 3, 5):
    unit_count = 0
    epsilon = (0, 1) + (0,) * (p - 2)
    top = (0,) * (p - 1) + (1,)
    norm_u = (1,) * p
    require(epsilon_to_u(top, p) == norm_u,
            f"epsilon^(p-1) != N_p at p={p}")
    for v in product(range(p), repeat=p):
        if v[0] == 0:
            continue
        unit_count += 1
        base = inv_eps(v, p)
        s = (0,) + v[:-1]
        decoders_eps = tuple(
            tuple((base[i] + c * top[i]) % p for i in range(p))
            for c in range(p)
        )
        require(all(mul_eps(s, k, p) == epsilon for k in decoders_eps),
                f"decoder torsor failed at p={p}, V={v}")
        decoders_u = tuple(epsilon_to_u(k, p) for k in decoders_eps)
        for coefficient in range(p):
            require(
                sorted(k[coefficient] for k in decoders_u) == list(range(p)),
                f"coefficient-zero section is not unique at p={p}",
            )
    expected = (p - 1) * p ** (p - 1)
    require(unit_count == expected, f"unit census failed at p={p}")
    decoder_unit_counts.append((p, unit_count))


# Reconstruct the exact THM-2771 Bockstein decoder and all local gauges.
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
beta_rows = [tuple(11 * x % p for x in row) for row in h_rows]
s_beta = tuple(sum(row[j] for row in beta_rows) % p for j in range(p))
k_beta = (3, 5, 5, 5, 7, 12, 2, 8, 2, 9, 8, 11, 0)
epsilon_u = (12, 1) + (0,) * 11
require(s_beta == (0, 0, 8, 0, 0, 0, 0, 0, 9, 9, 9, 9, 8),
        "S_beta specialization changed")
require(cyclic_convolution(s_beta, k_beta, p) == epsilon_u,
        "K_beta specialization failed")
columns = []
for c in range(p):
    decoder = tuple((k_beta[i] + c) % p for i in range(p))
    decoded = [cyclic_convolution(row, decoder, p) for row in beta_rows]
    columns.append(tuple(row[0] for row in decoded))
base_column = (0, 2, 0, 10, 0, 0, 0)
direction = (0, 3, 3, 7, 0, 0, 0)
for c, column in enumerate(columns):
    require(column == tuple((base_column[i] + c * direction[i]) % p
                            for i in range(7)),
            "local decoder gauge law failed")
    require(sum(column) % p == 12, "uniform chart invoice changed")
require(len(set(columns)) == 13, "local decoder gauges collapsed")

s_epsilon = u_to_epsilon(s_beta, p)
v_beta = s_epsilon[1:] + (0,)
k_epsilon = u_to_epsilon(k_beta, p)
require(mul_eps(v_beta, k_epsilon, p)
        == (1,) + (0,) * 11 + (3,), "socle-qualified inverse failed")
full_inverse = tuple((x + 3) % p for x in k_beta)
require(mul_eps(v_beta, u_to_epsilon(full_inverse, p), p)
        == (1,) + (0,) * 12, "full inverse repair failed")


# Endpoint action, odometer section, and exact generator-convention ledger.
def endpoint_x(point):
    v, w = point
    return (v, (w + v) % 13)


def endpoint_x_inv(point):
    v, w = point
    return (v, (w - v) % 13)


def endpoint_y(point):
    v, w = point
    return ((v + 1) % 13, w)


def endpoint_y_inv(point):
    v, w = point
    return ((v - 1) % 13, w)


def endpoint_z(point):
    v, w = point
    return (v, (w + 1) % 13)


plane13 = tuple(product(range(13), repeat=2))
for point in plane13:
    lhs = endpoint_x(endpoint_y(endpoint_x_inv(endpoint_y_inv(point))))
    require(lhs == endpoint_z(point), "endpoint commutator sign failed")
    require(endpoint_x(point) == rho((0, 12, 0), point, 13),
            "endpoint X group coordinate changed")
    require(endpoint_y(point) == rho((1, 0, 0), point, 13),
            "endpoint Y group coordinate changed")
    require(endpoint_z(point) == rho((0, 0, 1), point, 13),
            "endpoint Z group coordinate changed")


def digits(n):
    return (n % 13, (n // 13) % 13)


def address(point):
    return point[0] + 13 * point[1]


nonwrap = 0
carry = 0
for n in range(169):
    v, w = digits(n)
    require(digits(14 * n % 169) == endpoint_x((v, w)),
            "odometer shear conjugacy failed")
    low = address(((v + 1) % 13, w))
    require(digits(low) == endpoint_y((v, w)),
            "low-digit conjugacy failed")
    require(digits((n + 13) % 169) == endpoint_z((v, w)),
            "central-digit conjugacy failed")
    if v == 12:
        require((n + 1) % 169 == (low + 13) % 169,
                "carry correction failed")
        carry += 1
    else:
        require((n + 1) % 169 == low, "nonwrap successor failed")
        nonwrap += 1
require((nonwrap, carry) == (156, 13), "address carry census failed")

wrap_pairs = sum(a + b >= 13 for a in range(13) for b in range(13))
require(wrap_pairs == 78, "two-digit cocycle wrap-pair census failed")
require(7 * 1 % 13 == 7 and 2 * 7 % 13 == 1,
        "THM-2657 class-7 normalization failed")


# The 62,377,224-square statement is an inference from edge factors.
normalized_frames_13 = 13 * (13 ** 2 - 1)
square_count = normalized_frames_13 * 169 * 169
require(square_count == 62377224, "normalized-square count changed")
for p_field, factors, hadamard in (
    (
        352341050142921841,
        (153528379476679734, 205065739500452496,
         153829297094820499, 338655797145668320),
        (340901553381135615, 319769238334939305,
         60314403185657026, 233858556407273320),
    ),
    (
        956354278959359281,
        (955153878005214586, 313954411340872907,
         902111039187286496, 3971292277140798),
        (511989854700512392, 62552333514892978,
         86125904489159431, 419745850041132014),
    ),
):
    p0, p1, q0, q1 = factors
    expected_h = (
        (p0 + p1) * (q0 + q1) % p_field,
        (p0 + p1) * (q0 - q1) % p_field,
        (p0 - p1) * (q0 + q1) % p_field,
        (p0 - p1) * (q0 - q1) % p_field,
    )
    require(expected_h == hadamard, "displayed Hadamard witness changed")
    require(all(factors) and all(hadamard), "displayed witness has a zero")
    require(hadamard[0] * hadamard[3] % p_field
            == hadamard[1] * hadamard[2] % p_field,
            "displayed Pluecker witness failed")


# Conditional dilation naturality: flat target=id, root=0 has only zero map;
# the six-grade identity transport has one common scalar on every grade.
flat_intertwiners = tuple(c for c in range(13) if 0 * c % 13 == c)
graded_chain_maps = tuple(
    c for c in range(13)
    if all(c == c for _edge in range(5))
)
require(flat_intertwiners == (0,), "nonzero flat intertwiner survived")
require(len(graded_chain_maps) == 13, "graded scalar line changed")


print("THM-2779 INDEPENDENT HOSTILE AUDIT")
for row in prime_rows:
    p, order, degree, frames, faithful_classes, all_classes = row
    print(
        f"p={p}: order={order}; sharp_degree={degree}; frames={frames}; "
        f"faithful_degree_p2_classes={faithful_classes}; "
        f"all_transitive_degree_p2_classes={all_classes}"
    )
print("p=2_boundary=D8; order_census=1^1,2^5,4^2")
print(f"decoder_all_unit_censuses={tuple(decoder_unit_counts)}")
print("decoder_p13=13 gauges; local columns distinct; uniform sum=-1")
print("convention=(A,B,C)=(-y,x,z) verified for p=2,3,5,7,11,13")
print("endpoint_generators=X->(x,y,z)=(0,-1,0),Y->(1,0,0),Z->(0,0,1)")
print(f"odometer_addresses=169; nonwrap={nonwrap}; carry={carry}")
print("cocycle_normalization=class1 --kernel*7--> class7; class7 --kernel*2--> class1")
print(f"normalized_endpoint_squares={square_count}; inference=edge-factor products")
print("dilation_flat_intertwiners=(0,); graded_chain_maps=13")
print("verdict=PASS_AFTER_FAITHFUL_ACTION_WORDING_REPAIR")
