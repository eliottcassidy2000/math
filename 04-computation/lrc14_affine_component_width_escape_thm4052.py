#!/usr/bin/env python3
"""Exact certificate for THM-4052's affine component-width escape.

This primary path uses the closed formulas inherited from THM-4041 and
checks their consequences, endpoint controls, and one fully typed physical
row.  The independent audit reconstructs the spoiled set from circle walls.
No Python assertions are used, so every gate survives ``python -O``.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd, prod


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def distance(value):
    value -= value.numerator // value.denominator
    return min(value, 1 - value)


def clearance(speeds, time):
    return min(distance(speed * time) for speed in speeds)


def fully_spoiled(d, exceptions, y):
    for label in range(d):
        lift = Fraction(y + label, d)
        if not any(distance(speed * lift) < Fraction(1, 14)
                   for speed in exceptions):
            return False
    return True


def d2_component_formula(alpha, beta):
    """Return the complete sorted multiset of open component lengths."""
    require(alpha < beta and alpha % 2 == beta % 2 == 1,
            "d2 formula requires ordered positive odd exceptions")
    g = gcd(alpha, beta)
    a, b = alpha // g, beta // g
    widths = []
    r = 1
    while 7 * r < a + b:
        ell = min(Fraction(1, 7 * b),
                  Fraction(a + b - 7 * r, 14 * a * b))
        widths.extend([2 * ell / g] * (2 * g))
        r += 2
    return tuple(sorted(widths))


def d2_max_width(alpha, beta):
    require(alpha < beta and alpha % 2 == beta % 2 == 1,
            "d2 maximum requires ordered positive odd exceptions")
    g = gcd(alpha, beta)
    a, b = alpha // g, beta // g
    if a + b <= 7:
        return Fraction(0)
    return min(Fraction(2, 7 * beta),
               Fraction(alpha + beta - 7 * g, 7 * alpha * beta))


def frac_text(value):
    return f"{value.numerator}/{value.denominator}"


def semantic_rows(limit):
    rows = []
    for alpha in range(1, limit + 1, 2):
        for beta in range(alpha + 2, limit + 1, 2):
            widths = d2_component_formula(alpha, beta)
            rows.append(
                f"{alpha},{beta}:" + ",".join(frac_text(w) for w in widths)
            )
    return rows


def factor(n):
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def atlas_admissible(total):
    fs = factor(total)
    return bool(fs) and all(p % 3 == 2 and exponent <= 2
                            for p, exponent in fs.items())


def main():
    # The cited pack margins are just the exact Lipschitz losses
    # 1/12-1/14 and 1/11-1/14, doubled to obtain arc length.
    d2_margin = 2 * (Fraction(1, 12) - Fraction(1, 14))
    d34_margin = 2 * (Fraction(1, 11) - Fraction(1, 14))
    require(d2_margin == Fraction(1, 42), "d2 pack margin changed")
    require(d34_margin == Fraction(3, 77), "d3/d4 pack margin changed")

    # Check the exact d2 scalar reformulation throughout a hostile box.
    scalar_checks = 0
    for alpha in range(1, 128, 2):
        for beta in range(alpha + 2, 128, 2):
            g = gcd(alpha, beta)
            a, b = alpha // g, beta // g
            width = d2_max_width(alpha, beta)
            if a + b <= 7:
                require(width == 0, "certificate-negative pair has a component")
                continue
            expected = min(Fraction(2, 7 * beta),
                           Fraction(alpha + beta - 7 * g, 7 * alpha * beta))
            require(width == expected, "d2 maximum-width formula changed")
            for M in range(1, 65):
                topological = Fraction(1, 42 * M) >= width
                scalar = (beta >= 12 * M or
                          alpha * beta >= 6 * M * (alpha + beta - 7 * g))
                require(topological == scalar, "d2 scalar wedge mismatch")
                scalar_checks += 1
    require(scalar_checks == 126848, "d2 scalar audit universe changed")

    # Exact positive controls.  The displayed y is the pack phase and the
    # displayed x is one of its labelled lifts.
    controls = (
        (2, tuple(range(1, 12)), (1, 133), Fraction(1, 12),
         Fraction(13, 154), Fraction(167, 308), Fraction(1, 14)),
        (3, tuple(range(1, 11)), (1, 110, 23), Fraction(1, 11),
         Fraction(13, 140), Fraction(51, 140), Fraction(1, 14)),
        (4, tuple(range(1, 11)), (2, 185, 11), Fraction(1, 11),
         Fraction(137, 1540), Fraction(4757, 6160), Fraction(137, 1540)),
    )
    control_text = []
    for d, pack, exceptions, y0, y, x, expected_clearance in controls:
        require(fully_spoiled(d, exceptions, y0),
                "declared central pack phase is not fully spoiled")
        label = d * x - y
        require(label.denominator == 1 and 0 <= label < d,
                "escape time is not a labelled lift")
        speeds = tuple(d * h for h in pack) + exceptions
        actual = clearance(speeds, x)
        require(actual == expected_clearance, "control clearance changed")
        control_text.append(f"d{d}:{frac_text(x)}:{frac_text(actual)}")

    # The d2 method hostile: this argument cannot escape below its wedge.
    hostile_pack = tuple(range(1, 12))
    hostile_exceptions = (1, 11)
    hostile_y0 = Fraction(1, 12)
    hostile_radius = Fraction(1, 84 * max(hostile_pack))
    require(fully_spoiled(2, hostile_exceptions, hostile_y0),
            "method hostile lost central spoilage")
    require(fully_spoiled(2, hostile_exceptions, hostile_y0 - hostile_radius),
            "method hostile lost lower endpoint")
    require(fully_spoiled(2, hostile_exceptions, hostile_y0 + hostile_radius),
            "method hostile lost upper endpoint")

    # Fully typed THM-3818/4024 d2 control outside THM-4049's residue bank.
    omissions = (37, 43, 61, 67, 73, 79, 97, 103, 127)
    P = 3 * 5 * prod(omissions)
    require(P == 713721382004055345, "physical control product changed")
    body = tuple(2 * P // r for r in omissions) + (P // 3, P // 5)
    pair = (2, 8)
    row = body + pair
    require(gcd(*body) == 1, "physical body is not primitive")
    require(len(set(row)) == 13, "physical row lost distinctness")
    require(sum(row) == 574570283588268864, "physical row sum changed")
    Q = 91 ** 6
    B = Q ** 2
    require(sum(row) < B, "physical row left the THM-3818 box")
    require(sum(value % 2 == 0 for value in body) == 9,
            "physical body left c2=9")

    # Ten decoder-tree edges: eight among the nine even omission labels and
    # one to each odd exception.  The final pair component has ratio 1:4.
    body_by_name = {str(r): 2 * P // r for r in omissions}
    body_by_name.update({"P/3": P // 3, "P/5": P // 5})
    edge_names = (
        ("37", "73"), ("73", "43"), ("73", "97"), ("43", "67"),
        ("67", "103"), ("37", "79"), ("37", "127"), ("103", "61"),
        ("79", "P/3"), ("37", "P/5"),
    )
    tree_sums = []
    adjacency = {name: set() for name in body_by_name}
    for left, right in edge_names:
        common = gcd(body_by_name[left], body_by_name[right])
        tree_sums.append((body_by_name[left] + body_by_name[right]) // common)
        adjacency[left].add(right)
        adjacency[right].add(left)
    require(tuple(tree_sums) == (110, 116, 170, 110, 170, 116, 164, 164, 85, 47),
            "physical decoder-tree sums changed")
    reached = {next(iter(adjacency))}
    frontier = list(reached)
    while frontier:
        vertex = frontier.pop()
        for neighbour in adjacency[vertex] - reached:
            reached.add(neighbour)
            frontier.append(neighbour)
    require(reached == set(adjacency) and len(edge_names) == len(adjacency) - 1,
            "physical decoder edges are not a spanning tree")
    require(all(atlas_admissible(total) for total in tuple(tree_sums) + (5,)),
            "physical decoder tree left the admissible atlas")
    require(max(omissions) <= Q, "internal relation height exceeds Q")

    body_gcd_floor = min(gcd(x, y) for i, x in enumerate(body)
                         for y in body[i + 1:])
    one_body_floor = min(body) // 10
    two_body_floor = body_gcd_floor // 2
    require(one_body_floor == 1123970680321347,
            "one-body crossing lower bound changed")
    require(two_body_floor == 54561683510745,
            "two-body crossing lower bound changed")
    require(min(one_body_floor, two_body_floor) > Q,
            "a bounded crossing relation entered the physical control")

    H = tuple(P // r for r in omissions) + (1, 4)
    residues = tuple(h % 56 for h in H)
    forbidden = {0, 11, 14, 22, 23, 28, 33, 34, 42, 45}
    require(residues == (37, 11, 5, 43, 41, 23, 17, 47, 39, 1, 4),
            "physical residue profile changed")
    require(set(residues) & forbidden == {11, 23},
            "physical control no longer defeats the residue firewall")
    alpha, beta = P // 5, P // 3
    M = max(H)
    require(M == P // 37 and beta > 12 * M,
            "physical control left the width cone")
    require(d2_max_width(alpha, beta) == Fraction(1, 7 * P),
            "physical spoiled-component width changed")
    require(Fraction(1, 42 * M) == Fraction(37, 42 * P),
            "physical pack margin changed")
    require(clearance(row, Fraction(1, 14)) == Fraction(1, 14),
            "physical consequence control lost its boundary witness")

    rows = semantic_rows(79)
    require(len(rows) == 780, "d2 component universe changed")
    digest = sha256("\n".join(rows).encode()).hexdigest()

    print("LRC14 AFFINE COMPONENT-WIDTH ESCAPE CONES")
    print("status=PROVED_RELATIVE_TO_LRCUpTo13_AND_AFFINE_BOUNDARIES;LRC14=OPEN")
    print(f"pack_arc_lengths=d2:{frac_text(d2_margin)}/M,d3d4:{frac_text(d34_margin)}/M")
    print("coarse_cones=d2:E>=12M;d3:E>=11M;d4:3E>=44M")
    print("d2_failure_wedge=a+b>7;beta<12M;alpha*beta<6M(alpha+beta-7g)")
    print(f"d2_scalar_checks={scalar_checks}")
    print("controls=" + ",".join(control_text))
    print("method_hostile=H1..11;exceptions(1,11);central_and_both_margin_endpoints_spoiled")
    print(f"physical=P:{P};sum:{sum(row)};Q:{Q};crossing_floor:{min(one_body_floor,two_body_floor)}")
    print("physical_residues=" + ",".join(map(str, residues)) + ";firewall_hits=11,23")
    print(f"physical_width={frac_text(d2_max_width(alpha,beta))};witness=1/14")
    print(f"semantic_digest={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
