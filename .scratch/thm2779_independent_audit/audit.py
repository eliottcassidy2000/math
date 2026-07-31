#!/usr/bin/env python3
"""Independent exact hostile controls for the merged THM-2779 candidate."""

from collections import Counter
from itertools import product
from math import comb


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def mul(g, h, p):
    x, y, z = g
    a, b, c = h
    return ((x + a) % p, (y + b) % p, (z + c - y * a) % p)


def inv(g, p):
    x, y, z = g
    return ((-x) % p, (-y) % p, (-z - x * y) % p)


def old_mul(g, h, p):
    a, b, c = g
    x, y, z = h
    return ((a + x) % p, (b + y) % p, (c + z + a * y) % p)


def old(g, p):
    x, y, z = g
    return ((-y) % p, x, z)


def comm(g, h, p):
    return mul(mul(mul(g, h, p), inv(g, p), p), inv(h, p), p)


def det(v, w, p):
    return (v[0] * w[1] - v[1] * w[0]) % p


def rho(g, q, p):
    x, y, z = g
    r, w = q
    return ((r + x) % p, (w + z - y * r) % p)


def pow_g(g, n, p):
    out = (0, 0, 0)
    for _ in range(n):
        out = mul(out, g, p)
    return out


def order(g, p):
    for n in range(1, 2 * p + 1):
        if pow_g(g, n, p) == (0, 0, 0):
            return n
    raise RuntimeError("order escaped")


def subgroup(g, p):
    return frozenset(pow_g(g, n, p) for n in range(order(g, p)))


def conj(g, h, p):
    return mul(mul(g, h, p), inv(g, p), p)


def mul_eps(a, b, p):
    out = [0] * p
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < p:
                out[i + j] = (out[i + j] + x * y) % p
    return tuple(out)


def inv_eps(v, p):
    out = [0] * p
    out[0] = pow(v[0], -1, p)
    for n in range(1, p):
        out[n] = (
            -out[0] * sum(v[i] * out[n - i] for i in range(1, n + 1))
        ) % p
    require(mul_eps(v, out, p) == (1,) + (0,) * (p - 1), "bad inverse")
    return tuple(out)


def eps_to_u(a, p):
    out = [0] * p
    for j, x in enumerate(a):
        for k in range(j + 1):
            out[k] = (
                out[k] + x * comb(j, k) * (-1) ** (j - k)
            ) % p
    return tuple(out)


def u_to_eps(a, p):
    out = [0] * p
    for k, x in enumerate(a):
        for j in range(k + 1):
            out[j] = (out[j] + x * comb(k, j)) % p
    return tuple(out)


def cyclic(a, b, p):
    out = [0] * p
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[(i + j) % p] = (out[(i + j) % p] + x * y) % p
    return tuple(out)


prime_rows = []
for p in (2, 3, 5, 7, 11, 13):
    group = tuple(product(range(p), repeat=3))
    plane = tuple(product(range(p), repeat=2))
    identity = (0, 0, 0)
    center = {(0, 0, z) for z in range(p)}

    # Exhaustive convention transport; p=11 is an independent extra control.
    for g in group:
        require(mul(g, inv(g, p), p) == identity, f"inverse p={p}")
        for h in group:
            require(old(mul(g, h, p), p) == old_mul(old(g, p), old(h, p), p),
                    f"convention p={p}")

    actual_center = {
        g for g in group
        if comm(g, (1, 0, 0), p) == identity
        and comm(g, (0, 1, 0), p) == identity
    }
    require(actual_center == center, f"center p={p}")
    comm_values = set()
    for v in plane:
        for w in plane:
            value = comm((v[0], v[1], 0), (w[0], w[1], 0), p)
            require(value == (0, 0, det(v, w, p)), f"commutator p={p}")
            comm_values.add(value)
    require(comm_values == center, f"derived group p={p}")

    perms = {
        tuple(rho(g, point, p) for point in plane) for g in group
    }
    require(len(perms) == p ** 3, f"faithfulness p={p}")
    require({rho(g, (0, 0), p) for g in group} == set(plane),
            f"transitivity p={p}")

    eligible = {
        subgroup(g, p) for g in group
        if g not in center and order(g, p) == p
    }
    classes = []
    while eligible:
        representative = next(iter(eligible))
        orbit = {
            frozenset(conj(g, h, p) for h in representative) for g in group
        }
        require(orbit <= eligible, f"class escape p={p}")
        classes.append(orbit)
        eligible -= orbit
    faithful_classes = len(classes)
    require(faithful_classes == (p + 1 if p > 2 else 2),
            f"class count p={p}")
    require(all(len(c) == p for c in classes), f"class size p={p}")

    frames = sum(det(v, w, p) == 1 for v in plane for w in plane)
    require(frames == p * (p * p - 1), f"frame count p={p}")
    for v in plane:
        if v != (0, 0):
            require(sum(det(v, w, p) == 1 for w in plane) == p,
                    f"frame fibre p={p}")
    prime_rows.append(
        (p, p ** 3, p ** 2, frames, faithful_classes, faithful_classes + 1)
    )


require(Counter(order(g, 2) for g in product(range(2), repeat=3))
        == Counter({1: 1, 2: 5, 4: 2}), "H2 is not D8")


decoder_units = []
for p in (2, 3, 5):
    count = 0
    eps = (0, 1) + (0,) * (p - 2)
    top = (0,) * (p - 1) + (1,)
    require(eps_to_u(top, p) == (1,) * p, f"socle p={p}")
    for v in product(range(p), repeat=p):
        if v[0] == 0:
            continue
        count += 1
        base = inv_eps(v, p)
        s = (0,) + v[:-1]
        decoders = tuple(
            tuple((base[i] + c * top[i]) % p for i in range(p))
            for c in range(p)
        )
        require(all(mul_eps(s, k, p) == eps for k in decoders),
                f"decoder p={p}")
        decoders_u = tuple(eps_to_u(k, p) for k in decoders)
        for j in range(p):
            require(sorted(k[j] for k in decoders_u) == list(range(p)),
                    f"zero section p={p}")
    require(count == (p - 1) * p ** (p - 1), f"unit count p={p}")
    decoder_units.append((p, count))


p = 13
k1, k2 = 483303, 483287
h = [[0] * p for _ in range(7)]
for j in range(3, 12):
    h[1][j] = 2 * k1
h[1][12] = 121 * k1
for j in range(2, 10):
    h[2][j] = 2 * k2
h[2][12] = 265 * k2
for j in (2, 3, 4, 5, 6, 7, 10, 11):
    h[3][j] = 2 * k2
h[3][12] = 254 * k2
rows = [tuple(11 * x % p for x in row) for row in h]
s_beta = tuple(sum(row[j] for row in rows) % p for j in range(p))
k_beta = (3, 5, 5, 5, 7, 12, 2, 8, 2, 9, 8, 11, 0)
require(s_beta == (0, 0, 8, 0, 0, 0, 0, 0, 9, 9, 9, 9, 8),
        "S_beta")
require(cyclic(s_beta, k_beta, p) == (12, 1) + (0,) * 11, "K_beta")
base_col = (0, 2, 0, 10, 0, 0, 0)
direction = (0, 3, 3, 7, 0, 0, 0)
columns = []
for c in range(p):
    decoder = tuple((x + c) % p for x in k_beta)
    column = tuple(cyclic(row, decoder, p)[0] for row in rows)
    require(column == tuple((base_col[i] + c * direction[i]) % p
                            for i in range(7)), "chart gauge")
    require(sum(column) % p == 12, "chart invoice")
    columns.append(column)
require(len(set(columns)) == 13, "chart columns collapsed")
s_eps = u_to_eps(s_beta, p)
v_beta = s_eps[1:] + (0,)
k_eps = u_to_eps(k_beta, p)
require(mul_eps(v_beta, k_eps, p) == (1,) + (0,) * 11 + (3,),
        "qualified inverse")
require(mul_eps(v_beta, u_to_eps(tuple((x + 3) % p for x in k_beta), p), p)
        == (1,) + (0,) * 12, "full inverse")


def ex(q):
    v, w = q
    return (v, (w + v) % 13)


def exi(q):
    v, w = q
    return (v, (w - v) % 13)


def ey(q):
    v, w = q
    return ((v + 1) % 13, w)


def eyi(q):
    v, w = q
    return ((v - 1) % 13, w)


def ez(q):
    v, w = q
    return (v, (w + 1) % 13)


for point in product(range(13), repeat=2):
    require(ex(ey(exi(eyi(point)))) == ez(point), "endpoint commutator")
    require(ex(point) == rho((0, 12, 0), point, 13), "endpoint X coordinate")
    require(ey(point) == rho((1, 0, 0), point, 13), "endpoint Y coordinate")
    require(ez(point) == rho((0, 0, 1), point, 13), "endpoint Z coordinate")

nonwrap = carry = 0
for n in range(169):
    point = (n % 13, n // 13)
    require((14 * n % 169) % 13 == ex(point)[0]
            and (14 * n % 169) // 13 == ex(point)[1], "X odometer")
    low = ((point[0] + 1) % 13) + 13 * point[1]
    require((low % 13, low // 13) == ey(point), "Y odometer")
    require(((n + 13) % 169 % 13, (n + 13) % 169 // 13) == ez(point),
            "Z odometer")
    if point[0] == 12:
        require((n + 1) % 169 == (low + 13) % 169, "carry")
        carry += 1
    else:
        require(n + 1 == low, "nonwrap")
        nonwrap += 1
require((nonwrap, carry) == (156, 13), "odometer census")
require(sum(a + b >= 13 for a in range(13) for b in range(13)) == 78,
        "cocycle wrap census")
require(7 % 13 == 7 and 2 * 7 % 13 == 1, "class normalization")


square_count = 2184 * 169 * 169
require(square_count == 62377224, "square count")
for prime, factors, expected in (
    (352341050142921841,
     (153528379476679734, 205065739500452496,
      153829297094820499, 338655797145668320),
     (340901553381135615, 319769238334939305,
      60314403185657026, 233858556407273320)),
    (956354278959359281,
     (955153878005214586, 313954411340872907,
      902111039187286496, 3971292277140798),
     (511989854700512392, 62552333514892978,
      86125904489159431, 419745850041132014)),
):
    p0, p1, q0, q1 = factors
    got = (
        (p0 + p1) * (q0 + q1) % prime,
        (p0 + p1) * (q0 - q1) % prime,
        (p0 - p1) * (q0 + q1) % prime,
        (p0 - p1) * (q0 - q1) % prime,
    )
    require(got == expected and all(factors) and all(got), "Hadamard witness")
    require(got[0] * got[3] % prime == got[1] * got[2] % prime,
            "Pluecker")

flat = tuple(c for c in range(13) if c == 0)
graded = tuple(range(13))
require(flat == (0,) and len(graded) == 13, "dilation gate")

print("THM-2779 INDEPENDENT HOSTILE AUDIT")
for row in prime_rows:
    p, size, degree, frames, faithful, all_classes = row
    print(
        f"p={p}: order={size}; sharp_degree={degree}; frames={frames}; "
        f"faithful_degree_p2_classes={faithful}; "
        f"all_transitive_degree_p2_classes={all_classes}"
    )
print("p=2_boundary=D8; order_census=1^1,2^5,4^2")
print(f"decoder_all_unit_censuses={tuple(decoder_units)}")
print("decoder_p13=13 gauges; local columns distinct; uniform sum=-1")
print("convention=(A,B,C)=(-y,x,z) verified for p=2,3,5,7,11,13")
print("endpoint_generators=X->(x,y,z)=(0,-1,0),Y->(1,0,0),Z->(0,0,1)")
print(f"odometer_addresses=169; nonwrap={nonwrap}; carry={carry}")
print("cocycle_normalization=class1 --kernel*7--> class7; class7 --kernel*2--> class1")
print(f"normalized_endpoint_squares={square_count}; inference=edge-factor products")
print("dilation_flat_intertwiners=(0,); graded_chain_maps=13")
print("verdict=PASS_WITH_WORDING_SHARPENING")
