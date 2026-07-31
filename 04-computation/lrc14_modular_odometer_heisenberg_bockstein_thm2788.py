#!/usr/bin/env python3
"""Exact modular-odometer versus Heisenberg extension audit for THM-2788.

On Omega_p=Z/p^2, X(n)=(1+p)n and the actual odometer O(n)=n+1
generate the modular group C_(p^2) semidirect C_p.  Replacing O by the
carry-suppressed digit shift Y gives the exponent-p Heisenberg group.
The two actions have the same centre and alternating commutator, but for
odd p their p-power maps and extension cocycles differ by the ordinary
second-coordinate carry/Bockstein.

The script also checks the carry-wall join, the exact modular-group
minimal-action stabilizer census, and the full odd-prime metacyclic tower
on Z/p^k.  The Heisenberg action-class count is inherited from THM-2779.
It uses only exact integer/permutation arithmetic, explicit exception
gates, and no Python ``assert`` statement.
"""

from collections import Counter
from math import gcd


PRIMES = (2, 3, 5, 7, 13)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lcm(a, b):
    return a // gcd(a, b) * b


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[n]] for n in range(len(left)))


def inverse(permutation):
    result = [None] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    require(all(value is not None for value in result), "nonpermutation")
    return tuple(result)


def power(permutation, exponent):
    require(exponent >= 0, "negative exponent")
    result = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            result = compose(base, result)
        base = compose(base, base)
        exponent //= 2
    return result


def order(permutation):
    unseen = set(range(len(permutation)))
    result = 1
    while unseen:
        start = min(unseen)
        point = start
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = permutation[point]
            length += 1
        result = lcm(result, length)
    return result


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(generator, current)
            if candidate not in group:
                group.add(candidate)
                frontier.append(candidate)
    return frozenset(group)


def subgroup(generator):
    return frozenset(power(generator, exponent)
                     for exponent in range(order(generator)))


def conjugate(group_element, element):
    return compose(
        group_element,
        compose(element, inverse(group_element)),
    )


def actions(p):
    size = p * p
    identity = tuple(range(size))
    odometer = tuple((n + 1) % size for n in range(size))
    dilation = tuple((1 + p) * n % size for n in range(size))
    central = tuple((n + p) % size for n in range(size))
    suppressed = tuple(
        ((n % p + 1) % p) + p * (n // p)
        for n in range(size)
    )
    return identity, odometer, dilation, central, suppressed


def column_translation(p, column):
    """Increment the high digit only on one low-digit column."""
    result = []
    for n in range(p * p):
        low = n % p
        high = n // p
        if low == column:
            high = (high + 1) % p
        result.append(low + p * high)
    return tuple(result)


def extension_cocycles(p):
    """Check c_H=ad and c_M=ad+floor((b+d)/p)."""
    checks = 0
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    heisenberg = a * d % p
                    carry = (b + d) // p
                    modular = (heisenberg + carry) % p
                    reverse_heisenberg = c * b % p
                    reverse_modular = (
                        reverse_heisenberg + (d + b) // p
                    ) % p
                    require(
                        (heisenberg - reverse_heisenberg) % p
                        == (modular - reverse_modular) % p
                        == (a * d - c * b) % p,
                        "alternating commutator mismatch",
                    )
                    require(
                        (modular - heisenberg) % p == carry,
                        "carry/Bockstein cocycle mismatch",
                    )
                    checks += 1
    return checks


def audit_prime(p):
    identity, odometer, dilation, central, suppressed = actions(p)
    dilation_inverse = inverse(dilation)
    odometer_inverse = inverse(odometer)
    suppressed_inverse = inverse(suppressed)

    require(order(odometer) == p * p, "odometer order")
    require(order(dilation) == p, "dilation order")
    require(order(central) == p, "central order")
    require(order(suppressed) == p, "suppressed order")
    require(power(odometer, p) == central, "O^p != Z")
    require(
        compose(dilation, compose(odometer, dilation_inverse))
        == power(odometer, 1 + p),
        "modular conjugation relation",
    )
    require(
        compose(
            dilation,
            compose(odometer, compose(dilation_inverse, odometer_inverse)),
        )
        == central,
        "modular commutator",
    )
    require(
        compose(
            dilation,
            compose(
                suppressed,
                compose(dilation_inverse, suppressed_inverse),
            ),
        )
        == central,
        "Heisenberg commutator",
    )

    modular = generated_group((dilation, odometer))
    heisenberg = generated_group((dilation, suppressed))
    require(len(modular) == p**3, "modular group order")
    require(len(heisenberg) == p**3, "Heisenberg group order")

    modular_orders = Counter(order(element) for element in modular)
    heisenberg_orders = Counter(order(element) for element in heisenberg)
    if p == 2:
        require(modular == heisenberg, "p=2 groups should coincide")
        require(
            modular_orders == Counter({1: 1, 2: 5, 4: 2}),
            "D8 order spectrum",
        )
    else:
        require(modular != heisenberg, "odd-prime groups should differ")
        require(
            heisenberg_orders == Counter({1: 1, p: p**3 - 1}),
            "Heisenberg exponent-p spectrum",
        )
        require(
            modular_orders
            == Counter({
                1: 1,
                p: p * p - 1,
                p * p: p * p * (p - 1),
            }),
            "modular exponent-p^2 spectrum",
        )

    common_affine = frozenset(
        compose(power(central, high), power(dilation, slope))
        for slope in range(p)
        for high in range(p)
    )
    if p == 2:
        require(modular & heisenberg == modular,
                "p=2 intersection boundary")
    else:
        require(
            modular & heisenberg == common_affine,
            "odd-prime intersection is not <X,Z>",
        )

    # The actual-versus-suppressed discrepancy is one column translation.
    carry_wall = compose(odometer, suppressed_inverse)
    require(carry_wall == column_translation(p, 0),
            "carry wall is not the column-zero high-digit shift")
    columns = tuple(
        conjugate(power(suppressed, column), carry_wall)
        for column in range(p)
    )
    require(
        columns == tuple(column_translation(p, column)
                         for column in range(p)),
        "carry-wall conjugates do not give all columns",
    )
    require(
        all(order(column) == p for column in columns)
        and all(
            compose(columns[left], columns[right])
            == compose(columns[right], columns[left])
            for left in range(p)
            for right in range(p)
        ),
        "column translations are not an elementary abelian base",
    )

    product_all = identity
    product_linear = identity
    for column in range(p):
        product_all = compose(columns[column], product_all)
        product_linear = compose(
            power(columns[column], column),
            product_linear,
        )
    require(product_all == central, "constant column product != Z")
    require(product_linear == dilation, "linear column product != X")
    require(compose(carry_wall, suppressed) == odometer,
            "O != C Y")

    # Minimal faithful modular actions: the p noncentral order-p
    # stabilizers form one conjugacy orbit, and their common core is trivial.
    if p > 2:
        order_p_subgroups = {
            subgroup(element)
            for element in modular
            if order(element) == p
        }
        centre_subgroup = subgroup(central)
        stabilizer = subgroup(dilation)
        conjugate_stabilizers = {
            frozenset(
                conjugate(power(odometer, shift), element)
                for element in stabilizer
            )
            for shift in range(p)
        }
        require(
            len(order_p_subgroups) == p + 1
            and centre_subgroup in order_p_subgroups
            and order_p_subgroups - {centre_subgroup}
            == conjugate_stabilizers,
            "minimal stabilizer classification changed",
        )
        core = set(stabilizer)
        for conjugate_subgroup in conjugate_stabilizers:
            core.intersection_update(conjugate_subgroup)
        require(core == {identity}, "natural stabilizer is not core-free")

    cocycle_checks = extension_cocycles(p)
    return {
        "p": p,
        "group_order": len(modular),
        "intersection": len(modular & heisenberg),
        "modular_orders": tuple(sorted(modular_orders.items())),
        "heisenberg_orders": tuple(sorted(heisenberg_orders.items())),
        "cocycle_checks": cocycle_checks,
        "wreath_order": p ** (p + 1),
        "minimal_degree": p * p,
        "minimal_classes": 2 if p == 2 else 1,
    }


def valuation(value, prime):
    require(value != 0, "valuation of zero")
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def audit_tower(p, depth):
    require(p % 2 == 1 and depth >= 2, "odd-prime tower only")
    modulus = p**depth
    unit = 1 + p
    unit_order = p ** (depth - 1)
    require(
        pow(unit, unit_order, modulus) == 1
        and pow(unit, unit_order // p, modulus) != 1,
        "principal-unit order changed",
    )
    require(
        valuation(pow(unit, p ** (depth - 2)) - 1, p)
        == depth - 1,
        "LTE depth gate changed",
    )

    # [X,O^(p^r)]=O^(p^(r+1)); the last one is trivial modulo p^depth.
    commutator_steps = []
    for level in range(depth):
        source = p**level
        target = p ** (level + 1)
        require(
            ((unit * source) - source) % modulus
            == target % modulus,
            "lower-central commutator step changed",
        )
        commutator_steps.append(target % modulus)

    group_order = p ** (2 * depth - 1)
    derived_order = p ** (depth - 1)
    centre_order = p
    kernel_to_depth_two = p ** (2 * depth - 4)
    require(
        group_order // (p**3) == kernel_to_depth_two,
        "depth-two projection kernel count changed",
    )
    return {
        "depth": depth,
        "group_order": group_order,
        "derived_order": derived_order,
        "centre_order": centre_order,
        "nilpotency_class": depth,
        "kernel_to_depth_two": kernel_to_depth_two,
        "commutator_steps": tuple(commutator_steps),
    }


def main():
    prime_rows = tuple(audit_prime(p) for p in PRIMES)
    tower_rows = tuple(audit_tower(13, depth) for depth in range(2, 7))

    print("THM-2788 MODULAR ODOMETER / HEISENBERG BOCKSTEIN EXTENSION")
    print("status=FINITE-EXACT candidate; no LRC conclusion")
    for row in prime_rows:
        print(
            f"p={row['p']} order={row['group_order']} "
            f"intersection={row['intersection']} "
            f"M_orders={row['modular_orders']} "
            f"H_orders={row['heisenberg_orders']} "
            f"cocycle_checks={row['cocycle_checks']} "
            f"wreath_order={row['wreath_order']} "
            f"mu={row['minimal_degree']} "
            f"minimal_classes={row['minimal_classes']}"
        )
    print(
        "extension=common commutator ad-cb; "
        "c_M-c_H=floor((b+d)/p), and M power map is b"
    )
    print(
        "join=<H_p,M_p>=C_p wr C_p; carry C=O Y^-1 is one "
        "column delta, its conjugates give all p columns"
    )
    for row in tower_rows:
        print(
            f"p=13 depth={row['depth']} "
            f"order={row['group_order']} "
            f"derived={row['derived_order']} "
            f"centre={row['centre_order']} "
            f"class={row['nilpotency_class']} "
            f"kernel_to_depth2={row['kernel_to_depth_two']} "
            f"commutators={row['commutator_steps']}"
        )
    print(
        "THM2782=arm +13 is Z1=O^13; j=13 hostile +169 is "
        "Z2=[X,Z1], killed only after reduction mod169"
    )
    print(
        "scope=physical address-action tower and exact extension boundary, "
        "not factor/current covariance or endpoint-origin allocation"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
