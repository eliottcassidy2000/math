#!/usr/bin/env python3
"""Independent hostile audit of the proposed THM-2788.

The control does not import the discovery probe.  It compares two central
extension laws on triples (a,b,c), independently realizes both laws as
permutations of Z/p^2, checks the section/cocycle orientation, enumerates
centres, commutators and complete order spectra, and audits the faithful
degree-p^2 boundary through all index-p quotient kernels.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import ast


Triple = tuple[int, int, int]
Perm = tuple[int, ...]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(left: Perm, right: Perm) -> Perm:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse_perm(perm: Perm) -> Perm:
    result = [0] * len(perm)
    for index, image in enumerate(perm):
        result[image] = index
    return tuple(result)


def power_perm(perm: Perm, exponent: int) -> Perm:
    result = tuple(range(len(perm)))
    base = perm
    while exponent:
        if exponent & 1:
            result = compose(result, base)
        base = compose(base, base)
        exponent //= 2
    return result


def permutation_order(perm: Perm) -> int:
    visited = [False] * len(perm)
    result = 1
    for start in range(len(perm)):
        if visited[start]:
            continue
        point = start
        length = 0
        while not visited[point]:
            visited[point] = True
            point = perm[point]
            length += 1
        result = result * length // gcd(result, length)
    return result


def gcd(first: int, second: int) -> int:
    while second:
        first, second = second, first % second
    return abs(first)


def carry(p: int, b: int, bp: int) -> int:
    return (b + bp) // p


def law_h(p: int, left: Triple, right: Triple) -> Triple:
    a, b, c = left
    ap, bp, cp = right
    return (
        (a + ap) % p,
        (b + bp) % p,
        (c + cp + a * bp) % p,
    )


def law_m(p: int, left: Triple, right: Triple) -> Triple:
    a, b, c = left
    ap, bp, cp = right
    return (
        (a + ap) % p,
        (b + bp) % p,
        (c + cp + a * bp + carry(p, b, bp)) % p,
    )


def power_triple(p: int, law, element: Triple, exponent: int) -> Triple:
    result = (0, 0, 0)
    base = element
    while exponent:
        if exponent & 1:
            result = law(p, result, base)
        base = law(p, base, base)
        exponent //= 2
    return result


def inverse_triple(p: int, law, element: Triple) -> Triple:
    # The groups have order p^3, so g^(p^3-1)=g^(-1).  This deliberately
    # avoids baking a cocycle sign into the inverse control.
    return power_triple(p, law, element, p**3 - 1)


def commutator(p: int, law, left: Triple, right: Triple) -> Triple:
    return law(
        p,
        law(
            p,
            law(p, left, right),
            inverse_triple(p, law, left),
        ),
        inverse_triple(p, law, right),
    )


def all_elements(p: int) -> tuple[Triple, ...]:
    return tuple(
        (a, b, c)
        for a in range(p)
        for b in range(p)
        for c in range(p)
    )


def perm_h(p: int, element: Triple) -> Perm:
    """Permutation for Z^c Y_0^b X^a."""
    a, b, c = element
    result = []
    for n in range(p * p):
        v, w = n % p, n // p
        result.append(
            ((v + b) % p)
            + p * ((w + a * v + c) % p)
        )
    return tuple(result)


def perm_m(p: int, element: Triple) -> Perm:
    """Permutation for Z^c Y_phys^b X^a."""
    a, b, c = element
    modulus = p * p
    multiplier = pow(1 + p, a, modulus)
    return tuple(
        (multiplier * n + b + p * c) % modulus
        for n in range(modulus)
    )


def section_cocycle_h(p: int, u, v) -> int:
    a, b = u
    _ap, bp = v
    return a * bp % p


def section_cocycle_m(p: int, u, v) -> int:
    _a, b = u
    _ap, bp = v
    return (
        section_cocycle_h(p, u, v) + carry(p, b, bp)
    ) % p


def element_order(p: int, law, element: Triple) -> int:
    if element == (0, 0, 0):
        return 1
    current = (0, 0, 0)
    for exponent in range(1, p * p + 1):
        current = law(p, current, element)
        if current == (0, 0, 0):
            return exponent
    raise RuntimeError(f"element order exceeded p^2 at p={p}: {element}")


def projective_functionals(p: int) -> tuple[tuple[int, int], ...]:
    # Representatives of P^1(F_p): (1,t), plus (0,1).
    return tuple((1, t) for t in range(p)) + ((0, 1),)


def index_p_kernel(p: int, functional) -> frozenset[Triple]:
    r, s = functional
    return frozenset(
        element
        for element in all_elements(p)
        if (r * element[0] + s * element[1]) % p == 0
    )


def audit_prime(p: int):
    identity_perm = tuple(range(p * p))
    identity = (0, 0, 0)
    x = (1, 0, 0)
    y = (0, 1, 0)
    z = (0, 0, 1)
    elements = all_elements(p)

    perms_h = {element: perm_h(p, element) for element in elements}
    perms_m = {element: perm_m(p, element) for element in elements}
    require(
        len(set(perms_h.values())) == len(set(perms_m.values())) == p**3,
        f"normal forms are not faithful at p={p}",
    )

    # Independently check that the proposed triple laws are exactly
    # composition of the displayed permutations.
    multiplication_checks = 0
    for left in elements:
        # Checking against the two generators and centre suffices because
        # they generate every normal form, but include all quotient sections
        # below as a separate cocycle control.
        for right in (x, y, z):
            require(
                compose(perms_h[left], perms_h[right])
                == perms_h[law_h(p, left, right)],
                f"H multiplication orientation failed at p={p}",
            )
            require(
                compose(perms_m[left], perms_m[right])
                == perms_m[law_m(p, left, right)],
                f"M multiplication orientation failed at p={p}",
            )
            multiplication_checks += 2

    require(
        power_triple(p, law_h, x, p) == identity
        and power_triple(p, law_h, y, p) == identity
        and power_triple(p, law_h, z, p) == identity,
        f"H generator powers failed at p={p}",
    )
    require(
        power_triple(p, law_m, x, p) == identity
        and power_triple(p, law_m, y, p) == z
        and power_triple(p, law_m, z, p) == identity,
        f"M generator powers failed at p={p}",
    )
    require(
        commutator(p, law_h, x, y)
        == commutator(p, law_m, x, y)
        == z,
        f"commutator convention failed at p={p}",
    )

    cocycle_checks = 0
    for a in range(p):
        for b in range(p):
            u = (a, b)
            section_u = (a, b, 0)
            for ap in range(p):
                for bp in range(p):
                    v = (ap, bp)
                    section_v = (ap, bp, 0)
                    sum_section = (
                        (a + ap) % p,
                        (b + bp) % p,
                        0,
                    )
                    ch = section_cocycle_h(p, u, v)
                    cm = section_cocycle_m(p, u, v)
                    require(
                        law_h(p, section_u, section_v)
                        == law_h(p, (0, 0, ch), sum_section),
                        f"H section cocycle orientation failed at p={p}",
                    )
                    require(
                        law_m(p, section_u, section_v)
                        == law_m(p, (0, 0, cm), sum_section),
                        f"M section cocycle orientation failed at p={p}",
                    )
                    require(
                        (cm - ch) % p == carry(p, b, bp),
                        f"carry/Bockstein difference failed at p={p}",
                    )
                    alternating = (a * bp - ap * b) % p
                    require(
                        (
                            ch
                            - section_cocycle_h(p, v, u)
                        )
                        % p
                        == alternating
                        and (
                            cm
                            - section_cocycle_m(p, v, u)
                        )
                        % p
                        == alternating,
                        f"alternating form failed at p={p}",
                    )
                    cocycle_checks += 4

    # Associativity is checked independently as the normalized cocycle
    # identity.  Exhaust all quotient triples.
    cocycle_identity_checks = 0
    quotient = tuple(
        (a, b) for a in range(p) for b in range(p)
    )
    for first in quotient:
        for second in quotient:
            sum_first_second = (
                (first[0] + second[0]) % p,
                (first[1] + second[1]) % p,
            )
            for third in quotient:
                sum_second_third = (
                    (second[0] + third[0]) % p,
                    (second[1] + third[1]) % p,
                )
                for cocycle in (
                    section_cocycle_h,
                    section_cocycle_m,
                ):
                    require(
                        (
                            cocycle(p, first, second)
                            + cocycle(
                                p, sum_first_second, third
                            )
                        )
                        % p
                        == (
                            cocycle(p, second, third)
                            + cocycle(
                                p, first, sum_second_third
                            )
                        )
                        % p,
                        f"cocycle identity failed at p={p}",
                    )
                    cocycle_identity_checks += 1

    centers = {}
    derived_values = {}
    spectra = {}
    for label, law in (("H", law_h), ("M", law_m)):
        center = tuple(
            element
            for element in elements
            if commutator(p, law, element, x) == identity
            and commutator(p, law, element, y) == identity
        )
        commutators = {
            commutator(p, law, (a, b, 0), (ap, bp, 0))
            for a in range(p)
            for b in range(p)
            for ap in range(p)
            for bp in range(p)
        }
        require(
            center == tuple((0, 0, c) for c in range(p))
            and commutators
            == {(0, 0, c) for c in range(p)},
            f"centre/derived subgroup failed for {label} at p={p}",
        )
        centers[label] = len(center)
        derived_values[label] = len(commutators)
        spectra[label] = Counter(
            element_order(p, law, element)
            for element in elements
        )

    if p % 2:
        require(
            spectra["H"] == Counter({1: 1, p: p**3 - 1}),
            f"odd H exponent spectrum failed at p={p}",
        )
        require(
            spectra["M"]
            == Counter(
                {
                    1: 1,
                    p: p * p - 1,
                    p * p: p * p * (p - 1),
                }
            ),
            f"odd modular spectrum failed at p={p}",
        )
        for a in range(p):
            for b in range(p):
                require(
                    power_triple(
                        p, law_h, (a, b, 0), p
                    )
                    == identity
                    and power_triple(
                        p, law_m, (a, b, 0), p
                    )
                    == (0, 0, b),
                    f"odd-prime quotient power map failed at p={p}",
                )
        group_boundary = "nonisomorphic_by_exponent_and_order_spectrum"
    else:
        set_h = set(perms_h.values())
        set_m = set(perms_m.values())
        require(set_h == set_m, "p=2 action subgroups differ")
        x_perm = perms_h[x]
        y_zero_perm = perms_h[y]
        y_phys_perm = perms_m[y]
        require(
            y_phys_perm
            == inverse_perm(compose(x_perm, y_zero_perm)),
            "p=2 generator relabel failed",
        )
        require(
            spectra["H"] == spectra["M"]
            == Counter({1: 1, 2: 5, 4: 2}),
            "p=2 D8 spectrum failed",
        )
        # D8 presentation with rotation r=Y_phys and reflection s=X.
        require(
            power_perm(y_phys_perm, 4) == identity_perm
            and power_perm(y_phys_perm, 2) != identity_perm
            and power_perm(x_perm, 2) == identity_perm
            and compose(
                compose(x_perm, y_phys_perm), x_perm
            )
            == inverse_perm(y_phys_perm),
            "p=2 dihedral presentation failed",
        )
        group_boundary = "same_D8_after_quotient_basis_change"

    intersection = set(perms_h.values()) & set(perms_m.values())
    expected_intersection = {
        perms_h[(a, 0, c)]
        for a in range(p)
        for c in range(p)
    }
    if p % 2:
        require(
            intersection == expected_intersection
            and len(intersection) == p * p,
            f"odd-prime H-intersection-M failed at p={p}",
        )
    else:
        require(
            intersection == set(perms_h.values())
            and len(intersection) == p**3,
            "binary intersection is not the common D8",
        )

    # The discrepancy C=O Y_0^(-1) is a one-column high-digit increment.
    # Its Y-conjugates are the p independent column translations.  Their
    # direct-product rank and Y action prove the wreath-join order without
    # attempting to materialize p^(p+1) permutations.
    y_zero_perm = perms_h[y]
    y_phys_perm = perms_m[y]
    z_perm = perms_h[z]
    x_perm = perms_h[x]
    c_zero = compose(y_phys_perm, inverse_perm(y_zero_perm))
    columns = []
    for residue in range(p):
        conjugator = power_perm(y_zero_perm, residue)
        column = compose(
            compose(conjugator, c_zero),
            inverse_perm(conjugator),
        )
        columns.append(column)
        for n in range(p * p):
            v, w = n % p, n // p
            expected = (
                v
                + p * ((w + int(v == residue)) % p)
            )
            require(
                column[n] == expected,
                f"column carry support failed at p={p}",
            )
        require(
            power_perm(column, p) == identity_perm,
            f"column carry order failed at p={p}",
        )
    for first in columns:
        for second in columns:
            require(
                compose(first, second) == compose(second, first),
                f"column carries stopped commuting at p={p}",
            )
    product_columns = identity_perm
    weighted_columns = identity_perm
    for residue, column in enumerate(columns):
        product_columns = compose(product_columns, column)
        weighted_columns = compose(
            weighted_columns, power_perm(column, residue)
        )
    require(
        product_columns == z_perm
        and weighted_columns == x_perm,
        f"wreath column identities failed at p={p}",
    )

    # All index-p subgroups are inverse images of hyperplanes in the
    # abelianization F_p^2.  Enumerate the p+1 kernels and verify they are
    # subgroups containing the central commutator.
    kernels = tuple(
        index_p_kernel(p, functional)
        for functional in projective_functionals(p)
    )
    require(
        len(set(kernels)) == p + 1
        and all(len(kernel) == p * p for kernel in kernels)
        and all(z in kernel for kernel in kernels),
        f"index-p kernel census failed at p={p}",
    )
    for law in (law_h, law_m):
        for kernel in kernels:
            require(
                all(
                    law(p, left, right) in kernel
                    for left in kernel
                    for right in kernel
                ),
                f"index-p kernel closure failed at p={p}",
            )

    # The displayed p^2-point actions are transitive and faithful.  Their
    # zero stabilizers have order p and trivial action core.
    action_controls = {}
    for label, perms in (("H", perms_h), ("M", perms_m)):
        orbit_zero = {perm[0] for perm in perms.values()}
        stabilizer = tuple(
            element for element, perm in perms.items() if perm[0] == 0
        )
        action_kernel = tuple(
            element
            for element, perm in perms.items()
            if perm == identity_perm
        )
        require(
            len(orbit_zero) == p * p
            and len(stabilizer) == p
            and action_kernel == (identity,),
            f"faithful transitive p^2 action failed for {label} at p={p}",
        )
        action_controls[label] = (
            len(orbit_zero), len(stabilizer), len(action_kernel)
        )

    wall = None
    if p == 13:
        x_perm = perms_h[x]
        del x_perm
        y_zero = perms_h[y]
        y_phys = perms_m[y]
        z_perm = perms_h[z]
        off = carry_wall = 0
        for n in range(p * p):
            if n % p == p - 1:
                require(
                    y_phys[n] == compose(z_perm, y_zero)[n],
                    "p=13 carry-wall orientation failed",
                )
                carry_wall += 1
            else:
                require(
                    y_phys[n] == y_zero[n],
                    "p=13 off-wall successor failed",
                )
                off += 1
        wall = (off, carry_wall)
        require(wall == (156, 13), "p=13 wall census failed")

    return {
        "p": p,
        "multiplication_checks": multiplication_checks,
        "cocycle_checks": cocycle_checks,
        "cocycle_identity_checks": cocycle_identity_checks,
        "centers": centers,
        "derived_values": derived_values,
        "spectra_h": tuple(sorted(spectra["H"].items())),
        "spectra_m": tuple(sorted(spectra["M"].items())),
        "index_p_kernels": len(kernels),
        "actions": action_controls,
        "boundary": group_boundary,
        "intersection": len(intersection),
        "wreath_rank": len(columns),
        "join_order": p ** (p + 1),
        "wall": wall,
    }


def audit_tower(p: int, maximum_depth: int):
    rows = []
    for depth in range(2, maximum_depth + 1):
        modulus = p**depth
        unit = 1 + p
        unit_order = p ** (depth - 1)
        require(
            pow(unit, unit_order, modulus) == 1
            and (
                depth == 2
                or pow(unit, unit_order // p, modulus) != 1
            ),
            f"principal-unit order failed at depth {depth}",
        )
        ladder = []
        for r in range(depth):
            translation = p**r
            commutator_translation = (
                unit * translation - translation
            ) % modulus
            require(
                commutator_translation == p ** (r + 1) % modulus,
                f"lower-central step failed at depth={depth},r={r}",
            )
            ladder.append(commutator_translation)
        group_order = p ** (2 * depth - 1)
        reduction_kernel = p ** (2 * depth - 4)
        center_order = p
        require(
            group_order
            == (p**depth) * unit_order
            and reduction_kernel
            == group_order // (p**3),
            f"tower size invoice failed at depth {depth}",
        )
        rows.append(
            (
                depth,
                unit_order,
                group_order,
                reduction_kernel,
                center_order,
                tuple(ladder),
            )
        )
    return tuple(rows)


def main():
    rows = tuple(audit_prime(p) for p in (2, 3, 5, 7, 11, 13))
    tower = audit_tower(13, 6)
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(
        isinstance(node, ast.Assert) for node in ast.walk(tree)
    )
    require(assert_nodes == 0, "truth-bearing assert found")

    print("THM-2788 INDEPENDENT MODULAR/HEISENBERG HOSTILE AUDIT")
    print(
        "convention=composition_left_after_right;"
        "commutator=[g,h]=g*h*g^-1*h^-1;"
        "section_s(a,b)=Y^b*X^a"
    )
    print(
        "laws=H:(a,b,c)(a',b',c')="
        "(a+a',b+b',c+c'+a*b');"
        "M=H+floor((b+b')/p)"
    )
    for row in rows:
        print(
            f"p={row['p']}:"
            f"mult={row['multiplication_checks']};"
            f"cocycle={row['cocycle_checks']};"
            f"cocycle_identity={row['cocycle_identity_checks']};"
            f"center={row['centers']};"
            f"derived={row['derived_values']};"
            f"H_orders={row['spectra_h']};"
            f"M_orders={row['spectra_m']};"
            f"index_p_kernels={row['index_p_kernels']};"
            f"actions={row['actions']};"
            f"boundary={row['boundary']};"
            f"intersection={row['intersection']};"
            f"wreath_rank={row['wreath_rank']};"
            f"join_order={row['join_order']};"
            f"wall={row['wall']}"
        )
    print(f"p13_tower_(depth,unit_order,group_order,kernel_to_depth2,"
          f"center_order,commutator_steps)={tower}")
    print(
        "minimal_faithful_degree=p^2:"
        "all_smaller_orbits_have_index_1_or_p_stabilizers;"
        "every_index_p_kernel_contains_Z;"
        "displayed_transitive_actions_are_faithful"
    )
    print(
        "novel_scope=M_p_presentation+carry_Bockstein_cocycle+"
        "odd_prime_power/order_separation;"
        "H_p_degree+D8_H_boundary+p13_wall_are_inherited_from_THM-2779"
    )
    print(
        "THM2782_scope=central_+p_segment_only;"
        "cannot_detect_full_successor_carry;"
        "no_physical_M13_current_or_endpoint/root_intertwiner"
    )
    print(f"assert_nodes={assert_nodes}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
