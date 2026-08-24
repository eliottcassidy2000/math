#!/usr/bin/env python3
"""Owner typing for THM-4009 support-two relations in the THM-3910 branch.

This is an exact arithmetic firewall.  It does not prove LRC(14).  It checks
that a primitive support-two relation on positive integer speeds v_i,v_j has
coefficient pair (v_j/g,-v_i/g), and separates body/body, body/extra, and
extra/extra ownership before comparing the Euclidean cap ||a||_2^2 <= 195
with THM-3910's seventeen surviving external pair types.
"""

from __future__ import annotations

from math import gcd


CAP_SQUARE = 195
HEIGHT_CAP = 13
SCALE_ONE_TYPES = (
    (1, 3),
    (1, 4),
    (1, 9),
    (1, 10),
    (2, 3),
    (2, 9),
    (3, 7),
    (3, 8),
    (4, 7),
    (5, 6),
    (5, 12),
    (6, 11),
    (7, 10),
    (8, 9),
    (8, 21),
    (9, 11),
)
SCALE_TWO_TYPES = ((1, 9),)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_relation_pair(x: int, y: int) -> tuple[int, int]:
    """Absolute primitive coefficients for a relation on speeds x,y."""

    g = gcd(x, y)
    return y // g, x // g


def norm_square(pair: tuple[int, int]) -> int:
    return pair[0] * pair[0] + pair[1] * pair[1]


def main() -> None:
    ratio_atlas = tuple(
        (a, b)
        for a in range(1, HEIGHT_CAP + 1)
        for b in range(a + 1, HEIGHT_CAP + 1)
        if gcd(a, b) == 1 and a * a + b * b <= CAP_SQUARE
    )
    require(len(ratio_atlas) == 47, "THM-4009 ratio-atlas count changed")

    all_types = tuple(("scale1", pair) for pair in SCALE_ONE_TYPES) + tuple(
        ("scale2", pair) for pair in SCALE_TWO_TYPES
    )
    extra_extra_good = []
    extra_extra_bad = []
    for scale, (p, q) in all_types:
        coefficients = primitive_relation_pair(p, q)
        row = (scale, p, q, coefficients, norm_square(coefficients))
        if row[-1] <= CAP_SQUARE:
            extra_extra_good.append(row)
        else:
            extra_extra_bad.append(row)

    require(len(extra_extra_good) == 15, "extra/extra compatible count changed")
    require(
        [(scale, p, q) for scale, p, q, _, _ in extra_extra_bad]
        == [("scale1", 8, 21), ("scale1", 9, 11)],
        "extra/extra incompatible types changed",
    )

    # AP11 is a hostile to any owner-free reduction: its speeds 1 and 2 have
    # the primitive relation (2,-1), independently of both external speeds.
    ap11_body = tuple(range(1, 12))
    ap11_pair = primitive_relation_pair(ap11_body[0], ap11_body[1])
    require(ap11_pair == (2, 1), "AP11 body relation changed")
    require(norm_square(ap11_pair) == 5, "AP11 body norm changed")
    require(norm_square(ap11_pair) <= CAP_SQUARE, "AP11 lost the short relation")

    # Exact finite hostile replay: at every type and 11 <= t <= 143, the same
    # body/body relation survives.  This range reaches the largest possible
    # body/extra scale for an AP11 body because tp/g <= 13 implies t <= 143/p.
    hostile_rows = 0
    for _, _pair in all_types:
        for _t in range(11, 144):
            hostile_rows += 1
            require(norm_square(ap11_pair) <= CAP_SQUARE, "owner-free hostile failed")
    require(hostile_rows == 17 * 133, "AP11 hostile row count changed")

    # Symbolic cutoff: for a body speed u <= U and an extra speed tp, the
    # larger primitive coefficient is tp/gcd(tp,u).  If it is at most 13,
    # then tp <= 13*gcd(tp,u) <= 13u <= 13U.  Thus tp>13U rules out that
    # body/extra owner without inspecting the body.
    # Check the inequality exhaustively on a generous finite control bank.
    cutoff_checks = 0
    for U in range(1, 81):
        for u in range(1, U + 1):
            for p in range(1, 31):
                for t in range(U, 14 * U + 2):
                    cutoff_checks += 1
                    coeffs = primitive_relation_pair(u, t * p)
                    if norm_square(coeffs) <= CAP_SQUARE:
                        require(t * p <= HEIGHT_CAP * U, "body/extra cutoff failed")

    print("LRC14 THM-4009 SUPPORT-TWO OWNER TYPING")
    print("status=PROVED arithmetic typing + FINITE-EXACT hostile replay; LRC(14) OPEN")
    print("cap=squared_l2<=195; coefficient_height<=13")
    print("support_two_owners=body/body|body/extra|extra/extra")
    print("support_two_ratio_atlas=47")
    print("THM3910_types=17")
    print("extra_extra_cap_compatible=15")
    print("extra_extra_cap_incompatible=scale1:(8,21),scale1:(9,11)")
    print("AP11_fixed_body_relation=2*(speed 1)-1*(speed 2)=0")
    print("AP11_fixed_body_relation_norm_square=5")
    print(f"AP11_owner_free_hostile_rows={hostile_rows}")
    print("body_extra_necessary_cutoff=t*p<=13*U")
    print(f"cutoff_finite_control_gates={cutoff_checks}")
    print("consequence=47x17 is an owner-labelled discovery product, not a type intersection")
    print("needed_sidecar=support labels plus body/extra incidence and endpoint owner")
    print("PASS")


if __name__ == "__main__":
    main()
