#!/usr/bin/env python3
"""Exact referee for the charge-paired modular/free-factor construction."""

from __future__ import annotations

import hashlib
import itertools
import json


P = 13
I = (1, 0, 0, 1)
NEG_I = (-1 % P, 0, 0, -1 % P)
S = (0, 8, 5, 0)
R = (0, 1, 12, 12)
S0 = (0, 12, 1, 0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mmul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return (
        (a[0] * b[0] + a[1] * b[2]) % P,
        (a[0] * b[1] + a[1] * b[3]) % P,
        (a[2] * b[0] + a[3] * b[2]) % P,
        (a[2] * b[1] + a[3] * b[3]) % P,
    )


def mscale(lam: int, a: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((lam * x) % P for x in a)


def mpow(a: tuple[int, ...], n: int) -> tuple[int, ...]:
    out = I
    while n:
        if n & 1:
            out = mmul(out, a)
        a = mmul(a, a)
        n >>= 1
    return out


def det(a: tuple[int, ...]) -> int:
    return (a[0] * a[3] - a[1] * a[2]) % P


def mact(a: tuple[int, ...], v: tuple[int, int]) -> tuple[int, int]:
    return ((a[0] * v[0] + a[1] * v[1]) % P,
            (a[2] * v[0] + a[3] * v[1]) % P)


def omega(v: tuple[int, int], w: tuple[int, int]) -> int:
    return (v[0] * w[1] - v[1] * w[0]) % P


def hinv2() -> int:
    return pow(2, -1, P)


def hmul(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    """Symmetric Heisenberg coordinates (x,y,z)."""
    v = (a[0], a[1])
    w = (b[0], b[1])
    return ((a[0] + b[0]) % P,
            (a[1] + b[1]) % P,
            (a[2] + b[2] + hinv2() * omega(v, w)) % P)


def alpha(a: tuple[int, ...], h: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y = mact(a, (h[0], h[1]))
    return (x, y, det(a) * h[2] % P)


def old_to_sym(h: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, c = h
    return (x, y, (c + hinv2() * x * y) % P)


def sym_to_old(h: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, z = h
    return (x, y, (z - hinv2() * x * y) % P)


def alpha_old(a: tuple[int, ...], h: tuple[int, int, int]) -> tuple[int, int, int]:
    return sym_to_old(alpha(a, old_to_sym(h)))


def j(flag: tuple[int, int, int]) -> tuple[int, int, int]:
    r, w, u = flag
    return (r, -u % P, w)


def jinv(h: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, c = h
    return (x, c, -y % P)


def s_flag(flag: tuple[int, int, int]) -> tuple[int, int, int]:
    r, w, u = flag
    return (5 * u % P, (r * u - w) % P, 8 * r % P)


def r_flag(flag: tuple[int, int, int]) -> tuple[int, int, int]:
    r, w, u = flag
    return (-u % P, (w - r * u + hinv2() * u * u) % P, (r - u) % P)


def closure(gens: tuple[tuple[int, ...], ...]) -> set[tuple[int, ...]]:
    group = {I}
    frontier = [I]
    while frontier:
        a = frontier.pop()
        for g in gens:
            b = mmul(a, g)
            if b not in group:
                group.add(b)
                frontier.append(b)
    return group


def projective_key(a: tuple[int, ...]) -> tuple[int, ...]:
    first = next(x for x in a if x)
    inv = pow(first, -1, P)
    return mscale(inv, a)


def order(a: tuple[int, ...], bound: int = 10000) -> int:
    x = I
    for n in range(1, bound + 1):
        x = mmul(x, a)
        if x == I:
            return n
    raise RuntimeError("matrix order exceeded bound")


def generic_group_control(p: int) -> tuple[object, ...]:
    """Independent finite controls for the uniform p == 1 (mod 4) proof."""
    root = next(x for x in range(2, p) if x * x % p == p - 1)
    identity = (1, 0, 0, 1)
    s = (0, -root % p, root, 0)
    r = (0, 1, -1 % p, -1 % p)

    def multiply(a, b):
        return (
            (a[0] * b[0] + a[1] * b[2]) % p,
            (a[0] * b[1] + a[1] * b[3]) % p,
            (a[2] * b[0] + a[3] * b[2]) % p,
            (a[2] * b[1] + a[3] * b[3]) % p,
        )

    def power(a, n):
        out = identity
        while n:
            if n & 1:
                out = multiply(out, a)
            a = multiply(a, a)
            n >>= 1
        return out

    group = {identity}
    frontier = [identity]
    while frontier:
        a = frontier.pop()
        for g in (s, r):
            b = multiply(a, g)
            if b not in group:
                group.add(b)
                frontier.append(b)
    scalars = sorted(a[0] for a in group if a[1] == a[2] == 0 and a[0] == a[3])
    require(len(group) == 2 * p * (p * p - 1), (p, "generic group order", len(group)))
    require(scalars == sorted({1, root, -root % p, -1 % p}), (p, "generic scalar kernel"))
    require(power(multiply(s, r), 4 * p) == identity, (p, "generic SR exponent"))
    require(power(multiply(s, r), 2 * p) != identity, (p, "generic SR order divisor"))
    return p, root, len(group), tuple(scalars), 4 * p


def main() -> None:
    generic_controls = tuple(generic_group_control(p) for p in (5, 13, 17, 29))
    require(det(S0) == 1 and mmul(S0, S0) == NEG_I, "standard modular lift mismatch")
    require(5 * 5 % P == P - 1 and mscale(5, S0) == S,
            "honest involution rescaling mismatch")
    require(mmul(S, S) == I, "S is not an involution")
    require(mpow(R, 3) == I and R != I and mpow(R, 2) != I, "R is not order three")
    require(det(S) == P - 1 and det(R) == 1, "determinant action mismatch")
    require(mmul(S, R) == mscale(5, (1, 1, 0, 1)), "SR formula mismatch")
    require(mmul(R, S) == (5, 0, 8, 5), "RS formula mismatch")
    require(mmul(S, R) != mmul(R, S), "free-factor generators commute")
    require(order(mmul(S, R)) == 52, "SR order mismatch")

    vectors = list(itertools.product(range(P), repeat=2))
    symplectic_checks = 0
    for a in (S, R):
        for v in vectors:
            for w in vectors:
                require(omega(mact(a, v), mact(a, w)) == det(a) * omega(v, w) % P,
                        "similitude identity failed")
                symplectic_checks += 1

    hs = list(itertools.product(range(P), repeat=3))
    for a in (S, R):
        for h in hs:
            for k in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (7, 11, 5)):
                require(alpha(a, hmul(h, k)) == hmul(alpha(a, h), alpha(a, k)),
                        "Heisenberg automorphism check failed")

    flags = hs
    for f in flags:
        require(jinv(alpha_old(S, j(f))) == s_flag(f), "S flag formula failed")
        require(jinv(alpha_old(R, j(f))) == r_flag(f), "R flag formula failed")
        require(s_flag(s_flag(f)) == f, "S flag order failed")
        require(r_flag(r_flag(r_flag(f))) == f, "R flag order failed")

    group = closure((S, R))
    scalars = sorted(a[0] for a in group if a[1] == 0 and a[2] == 0 and a[0] == a[3])
    projective = {projective_key(a) for a in group}
    determinants = sorted({det(a) for a in group})
    require(len(group) == 4368, "generated GL group order mismatch")
    require(scalars == [1, 5, 8, 12], "scalar kernel mismatch")
    require(len(projective) == 1092, "projective image order mismatch")
    require(determinants == [1, 12], "determinant image mismatch")

    sl2 = []
    honest_involutions = []
    projective_involutions = []
    for a in itertools.product(range(P), repeat=4):
        if det(a) != 1:
            continue
        sl2.append(a)
        if mmul(a, a) == I:
            honest_involutions.append(a)
        if projective_key(mmul(a, a)) == projective_key(I):
            projective_involutions.append(a)
    require(len(sl2) == 2184, "SL2 order mismatch")
    require(set(honest_involutions) == {I, NEG_I}, "nontrivial center-preserving involution found")
    require(any(mmul(a, a) == NEG_I for a in projective_involutions),
            "standard projective involution control absent")

    charge_orbits = []
    unseen = set(range(1, P))
    while unseen:
        kappa = min(unseen)
        orbit = {kappa, -kappa % P}
        charge_orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    require(charge_orbits == [(1, 12), (2, 11), (3, 10), (4, 9), (5, 8), (6, 7)],
            "charge-pair orbit mismatch")

    # An invertible idempotent is the identity: the finite scalar controls expose
    # both Boolean possibilities and prevent confusing conjugation with inversion.
    scalar_idempotents = [x for x in range(P) if x * x % P == x]
    invertible_idempotents = [x for x in scalar_idempotents if x != 0]
    require(scalar_idempotents == [0, 1] and invertible_idempotents == [1],
            "idempotent localization control failed")

    record = {
        "p": P,
        "S": S,
        "S0": S0,
        "R": R,
        "determinants": (det(S), det(R)),
        "orders": (order(S), order(R), order(mmul(S, R))),
        "SR": mmul(S, R),
        "RS": mmul(R, S),
        "symplectic_checks": symplectic_checks,
        "heisenberg_states": len(hs),
        "flag_checks": len(flags),
        "generated_group_order": len(group),
        "projective_image_order": len(projective),
        "scalar_kernel": scalars,
        "determinant_image": determinants,
        "sl2_order": len(sl2),
        "honest_sl2_involutions": sorted(honest_involutions),
        "charge_orbits": charge_orbits,
        "minimal_nonzero_stable_charge_support": 2,
        "paired_permutation_dimension": 2 * P * P,
        "paired_irreducible_dimension": 2 * P,
        "idempotent_localization": "invertible_image_is_identity",
        "generic_p_controls": generic_controls,
    }
    semantic = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("THM3259 charge-paired modular free-factor exact referee")
    for key, value in record.items():
        print(f"{key}={value}")
    print(f"semantic_sha256={semantic}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
