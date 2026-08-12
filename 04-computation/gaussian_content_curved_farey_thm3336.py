#!/usr/bin/env python3
"""Exact companion for THM-3336; dependency-free and -O safe."""

from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gpow(z, n):
    out = (1, 0)
    base = z
    while n:
        if n & 1:
            out = gmul(out, base)
        base = gmul(base, base)
        n //= 2
    return out


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def primitive(z):
    return z != (0, 0) and gcd(abs(z[0]), abs(z[1])) == 1


def content(z):
    return gcd(abs(z[0]), abs(z[1]))


def tau(z):
    require(primitive(z), f"tau requires primitive input: {z}")
    return int(z[0] % 2 != 0 and z[1] % 2 != 0)


def dcontent(z, w):
    require(primitive(z) and primitive(w), "dcontent requires primitive inputs")
    return content(gmul(z, w))


def star(z, w):
    zw = gmul(z, w)
    d = dcontent(z, w)
    out = (zw[0] // d, zw[1] // d)
    require(primitive(out), f"star output not primitive: {z},{w}")
    return out


def kappa(z, w):
    return dcontent(z, w) // (2 ** (tau(z) * tau(w)))


def phi(z):
    a, b = z
    return (a * a - b * b, 2 * a * b, norm(z))


def primitive_triple(z):
    scale = 2 ** tau(z)
    return tuple(x // scale for x in phi(z))


def box(x, y):
    a, b, c = x
    d, e, f = y
    return (a * d - b * e, a * e + b * d, c * f)


def hypotenuse(z):
    return norm(z) // (2 ** tau(z))


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1]


def vadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def lorentz(x, y):
    return x[2] * y[2] - x[0] * y[0] - x[1] * y[1]


def gdiv_exact(z, p):
    n = norm(p)
    q = gmul(z, gconj(p))
    if q[0] % n or q[1] % n:
        return None
    return (q[0] // n, q[1] // n)


def valuation(z, p):
    out = 0
    while True:
        q = gdiv_exact(z, p)
        if q is None:
            return out
        z = q
        out += 1


def charge(z, p):
    return valuation(z, p) - valuation(z, gconj(p))


def allocation_cube(prime_rows):
    out = {}
    for bits in product((0, 1), repeat=len(prime_rows)):
        z = (1, 0)
        for bit, (pi, _, exponent) in zip(bits, prime_rows):
            factor = gconj(pi) if bit else pi
            z = gmul(z, gpow(factor, exponent))
        out[bits] = z
    return out


def xor(x, y):
    return tuple(a ^ b for a, b in zip(x, y))


def main():
    canonical_small = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    all_small = [
        (a, b)
        for a in range(-8, 9)
        for b in range(-8, 9)
        if primitive((a, b))
    ]
    primitive_bank = canonical_small + [z for z in all_small if z not in canonical_small][:52]
    torsion = []
    torsion_state = (1, 0)
    for _ in range(8):
        torsion_state = star(torsion_state, (1, 1))
        torsion.append(torsion_state)
    require(
        torsion
        == [(1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)],
        "signed primitive C8 torsion failed",
    )
    cocycle_rows = 0
    for z in primitive_bank:
        for w in primitive_bank:
            d = dcontent(z, w)
            require((d & -d).bit_length() - 1 == tau(z) * tau(w), f"2-content law failed: {z},{w}")
            zw = star(z, w)
            require(tau(zw) == (tau(z) ^ tau(w)), f"tau XOR failed: {z},{w}")
            require(
                box(primitive_triple(z), primitive_triple(w))
                == tuple(kappa(z, w) ** 2 * x for x in primitive_triple(zw)),
                f"Pythagorean content law failed: {z},{w}",
            )
            require(
                hypotenuse(zw) * kappa(z, w) ** 2 == hypotenuse(z) * hypotenuse(w),
                f"hypotenuse coboundary failed: {z},{w}",
            )
        for w in primitive_bank[:24]:
            for t in primitive_bank[:12]:
                require(
                    dcontent(z, w) * dcontent(star(z, w), t)
                    == dcontent(w, t) * dcontent(z, star(w, t)),
                    f"d cocycle failed: {z},{w},{t}",
                )
                require(
                    kappa(z, w) * kappa(star(z, w), t)
                    == kappa(w, t) * kappa(z, star(w, t)),
                    f"kappa cocycle failed: {z},{w},{t}",
                )
                require(star(star(z, w), t) == star(z, star(w, t)), "star associativity failed")
                cocycle_rows += 1

    # Charge-lattice cancellation at two split primes.
    split_primes = (((2, 1), 5), ((3, 2), 13))
    charge_rows = 0
    charge_objects = {}
    for l5 in range(-3, 4):
        for l13 in range(-3, 4):
            z = (1, 0)
            for exponent, (pi, _) in zip((l5, l13), split_primes):
                if exponent >= 0:
                    z = gmul(z, gpow(pi, exponent))
                else:
                    z = gmul(z, gpow(gconj(pi), -exponent))
            require(primitive(z), "charge object not primitive")
            require((charge(z, split_primes[0][0]), charge(z, split_primes[1][0])) == (l5, l13), "charge readback failed")
            charge_objects[(l5, l13)] = z
    for lz, z in charge_objects.items():
        for lw, w in charge_objects.items():
            got = kappa(z, w)
            expected = 1
            for i, (_, rational_prime) in enumerate(split_primes):
                defect = (abs(lz[i]) + abs(lw[i]) - abs(lz[i] + lw[i])) // 2
                expected *= rational_prime**defect
            require(got == expected, f"charge cancellation failed: {lz},{lw}")
            sw = star(z, w)
            require(
                tuple(charge(sw, pi) for pi, _ in split_primes)
                == tuple(lz[i] + lw[i] for i in range(2)),
                "charge addition failed",
            )
            fz = 5 ** max(0, -lz[0]) * 13 ** max(0, -lz[1])
            fw = 5 ** max(0, -lw[0]) * 13 ** max(0, -lw[1])
            ls = (lz[0] + lw[0], lz[1] + lw[1])
            fs = 5 ** max(0, -ls[0]) * 13 ** max(0, -ls[1])
            require(got * fs == fz * fw, "orientation-cochain coboundary failed")
            charge_rows += 1

    # Farey edges and fixed-multiplier curved faces.
    rays = [(a, b) for a in range(0, 17) for b in range(0, 17) if primitive((a, b))]
    farey_edges = [(u, v) for u in rays for v in rays if det(u, v) == 1]
    for u, v in farey_edges:
        d = dcontent(u, v)
        require(d == gcd(norm(u), norm(v)), f"Farey content law failed: {u},{v}")
        require(d % 2 == 1, f"Farey content should be odd: {u},{v}")
        h = dot(u, v)
        require(norm(u) * norm(v) == h * h + 1, "Farey dot/determinant identity failed")
        require(norm(star(u, v)) == (h * h + 1) // (d * d), "Farey star norm failed")
        require((h * h + 1) % (d * d) == 0, "Farey congruence failed")
    require(dcontent((2, 1), (2, 1)) == 1, "nonadjacent equal-norm hostile failed")
    require(dcontent((2, 1), (7, 4)) == 5 and star((2, 1), (7, 4)) == (2, 3), "Farey positive hostile failed")

    multipliers = [(a, b) for a in range(1, 10) for b in range(0, 10) if primitive((a, b))]
    face_rows = 0
    for alpha in multipliers:
        A = norm(alpha)
        for u0, u1 in farey_edges[:700]:
            u2 = vadd(u0, u1)
            ds = (dcontent(alpha, u0), dcontent(alpha, u1), dcontent(alpha, u2))
            require(gcd(ds[0], ds[1]) == gcd(ds[1], ds[2]) == gcd(ds[2], ds[0]) == 1, "face contents not coprime")
            require(A % (ds[0] * ds[1] * ds[2]) == 0, "face curvature not integral")
            curvature = A // (ds[0] * ds[1] * ds[2])
            xs = (star(alpha, u0), star(alpha, u1), star(alpha, u2))
            require(
                (ds[2] * xs[2][0], ds[2] * xs[2][1])
                == (ds[0] * xs[0][0] + ds[1] * xs[1][0], ds[0] * xs[0][1] + ds[1] * xs[1][1]),
                "weighted face relation failed",
            )
            edge_dets = (abs(det(xs[0], xs[1])), abs(det(xs[1], xs[2])), abs(det(xs[2], xs[0])))
            require(edge_dets == (curvature * ds[2], curvature * ds[0], curvature * ds[1]), "curved determinant law failed")
            require(
                (lorentz(phi(xs[0]), phi(xs[1])), lorentz(phi(xs[1]), phi(xs[2])), lorentz(phi(xs[2]), phi(xs[0])))
                == tuple(2 * d * d for d in edge_dets),
                "raw Lorentz shell failed",
            )
            normalized_shells = (
                lorentz(primitive_triple(xs[0]), primitive_triple(xs[1])),
                lorentz(primitive_triple(xs[1]), primitive_triple(xs[2])),
                lorentz(primitive_triple(xs[2]), primitive_triple(xs[0])),
            )
            require(
                normalized_shells
                == (
                    2 * edge_dets[0] ** 2 // (2 ** (tau(xs[0]) + tau(xs[1]))),
                    2 * edge_dets[1] ** 2 // (2 ** (tau(xs[1]) + tau(xs[2]))),
                    2 * edge_dets[2] ** 2 // (2 ** (tau(xs[2]) + tau(xs[0]))),
                ),
                "normalized Lorentz shell failed",
            )
            face_rows += 1

    alpha = (8, 1)
    face = ((1, 1), (2, 3), (3, 4))
    ds = tuple(dcontent(alpha, u) for u in face)
    xs = tuple(star(alpha, u) for u in face)
    require(ds == (1, 13, 5), "composite face contents failed")
    require((abs(det(xs[0], xs[1])), abs(det(xs[1], xs[2])), abs(det(xs[2], xs[0]))) == (5, 1, 13), "composite face determinants failed")
    parity_face = ((1, 0), (0, 1), (1, 1))
    parity_ds = tuple(dcontent((1, 1), u) for u in parity_face)
    parity_xs = tuple(star((1, 1), u) for u in parity_face)
    require(parity_ds == (1, 1, 2), "ramified-two face failed")
    require(
        (abs(det(parity_xs[0], parity_xs[1])), abs(det(parity_xs[1], parity_xs[2])), abs(det(parity_xs[2], parity_xs[0])))
        == (2, 1, 1),
        "ramified-two determinants failed",
    )

    # Boolean fibres and state-dependent groupoid flips.
    rows65 = (((2, 1), 5, 1), ((3, 2), 13, 1))
    rows1105 = rows65 + (((4, 1), 17, 1),)
    folded_by_c = {}
    for rows in (rows65, rows1105):
        cube = allocation_cube(rows)
        c = 1
        for _, rational_prime, exponent in rows:
            c *= rational_prime**exponent
        require(all(norm(z) == c and primitive(z) for z in cube.values()), "allocation cube failed")
        r = len(rows)
        for x, z in cube.items():
            for y, w in cube.items():
                direction = xor(x, y)
                expected = 1
                for bit, (_, rational_prime, exponent) in zip(direction, rows):
                    if bit:
                        expected *= rational_prime**exponent
                require(kappa(z, w) == expected, f"Boolean kappa failed at c={c},{x},{y}")
            for mask in product((0, 1), repeat=r):
                multiplier = (1, 0)
                expected_content = 1
                for bit, chosen_bit, (pi, rational_prime, exponent) in zip(mask, x, rows):
                    if not bit:
                        continue
                    chosen = gconj(pi) if chosen_bit else pi
                    multiplier = gmul(multiplier, gpow(gconj(chosen), 2 * exponent))
                    expected_content *= rational_prime**exponent
                require(primitive(multiplier), "state-dependent multiplier not primitive")
                require(dcontent(multiplier, z) == expected_content, "groupoid flip content failed")
                require(star(multiplier, z) == cube[xor(x, mask)], "groupoid bit flip failed")
        folded = set()
        for direction in product((0, 1), repeat=r):
            # Zero and all-ones are the same zero direction after quotienting
            # raw allocations by global conjugation.
            if not any(direction) or all(direction):
                continue
            value = 1
            for bit, (_, rational_prime, exponent) in zip(direction, rows):
                if bit:
                    value *= rational_prime**exponent
            folded.add(tuple(sorted((value, c // value))))
        folded_by_c[c] = sorted(folded)
    require(folded_by_c[65] == [(5, 13)], "c=65 folded weight failed")
    require(folded_by_c[1105] == [(5, 221), (13, 85), (17, 65)], "c=1105 folded weights failed")

    # Multiplication does not descend through conjugation without a section.
    pi5 = (2, 1)
    require(hypotenuse(star(pi5, pi5)) == 25, "conjugation hostile same-orientation failed")
    require(hypotenuse(star(gconj(pi5), pi5)) == 1, "conjugation hostile opposite-orientation failed")

    # With an explicit orientation section, one fixed multiplier can swap the
    # two c=65 quotient points, but it exits the grade on conjugate lifts.
    pi13 = (3, 2)
    z0, z1 = gmul(pi5, pi13), gmul(pi5, gconj(pi13))
    section_multiplier = gpow(gconj(pi5), 2)
    require((z0, z1, section_multiplier) == ((4, 7), (8, -1), (3, -4)), "section setup failed")
    require(star(section_multiplier, z0) == (8, 1), "section swap first point failed")
    require(star(section_multiplier, z1) == (4, -7), "section swap second point failed")
    require(hypotenuse(star(section_multiplier, gconj(z0))) == 1625, "section conjugate hostile 0 failed")
    require(hypotenuse(star(section_multiplier, gconj(z1))) == 1625, "section conjugate hostile 1 failed")

    pi = (2, 1)
    fixed = gpow(pi, 2)
    require(fixed == (3, 4), "fixed multiplier setup failed")
    require(gmul(fixed, gconj(pi)) == (10, 5) and star(fixed, gconj(pi)) == pi, "one-sided flip failed")
    require(gmul(fixed, pi) == (2, 11) and hypotenuse(star(fixed, pi)) == 125, "opposite orientation hostile failed")

    # The universal grade-preservers have primitive hypotenuse one.
    preservers = []
    for alpha in primitive_bank:
        preserves_bank = all(hypotenuse(star(alpha, z)) == hypotenuse(z) for z in primitive_bank)
        if preserves_bank:
            preservers.append(alpha)
        require(preserves_bank == (hypotenuse(alpha) == 1), f"grade-preserver classification failed at {alpha}")

    print("THM-3336 VERIFIED-EXACT")
    print(f"signed_primitive_torsion_C8={torsion}")
    print(f"primitive_cocycle_triples={cocycle_rows}")
    print(f"charge_lattice_pairs={charge_rows}; split_primes=5,13")
    print(f"farey_edges={len(farey_edges)}; multiplier_face_rows={face_rows}")
    print("farey_positive=(2+i,7+4i): content=5; output=2+3i")
    print("composite_face=alpha(8+i); d=(1,13,5); edge_det=(5,1,13)")
    print("ramified_two_face=d=(1,1,2); edge_det=(2,1,1)")
    print(f"folded_weights_65={folded_by_c[65]}")
    print(f"folded_weights_1105={folded_by_c[1105]}")
    print("conjugation_hostile=(2+i)^2:grade25; (2-i)*(2+i):grade1")
    print("section_swap_c65=works_on_(4+7i,8-i); conjugates_exit_grade1625")
    print("fixed_multiplier_hostile=(3+4i)*(2-i)->grade5; *(2+i)->grade125")
    print(f"universal_grade_preservers_in_bank={preservers}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
