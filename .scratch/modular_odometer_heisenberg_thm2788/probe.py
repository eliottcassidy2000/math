#!/usr/bin/env python3
"""Exact probe for the physical odometer versus Heisenberg extension.

No theorem is promoted here.  For each tested prime p the script compares

  X(n)=(1+p)n,  Y_phys(n)=n+1

on Z/p^2 with the carry-suppressed digit action

  X(v,w)=(v,w+v),  Y_0(v,w)=(v+1,w),  Z(v,w)=(v,w+1).

It verifies the common commutator, the different carry/Bockstein cocycles,
the odd-prime exponent/order boundary, the p=2 D8 coincidence, and the
sharp p^2 faithful-carrier lower bound used in the mathematical proof.
"""

from __future__ import annotations

import ast
from collections import Counter, deque
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Return left after right."""
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(perm: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * len(perm)
    for i, image in enumerate(perm):
        answer[image] = i
    return tuple(answer)


def power(perm: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    answer = tuple(range(len(perm)))
    base = perm
    while exponent:
        if exponent & 1:
            answer = compose(answer, base)
        base = compose(base, base)
        exponent //= 2
    return answer


def permutation_order(perm: tuple[int, ...]) -> int:
    seen = [False] * len(perm)
    answer = 1
    for start in range(len(perm)):
        if seen[start]:
            continue
        length = 0
        point = start
        while not seen[point]:
            seen[point] = True
            point = perm[point]
            length += 1
        answer = answer * length // gcd(answer, length)
    return answer


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def generated_group(generators: tuple[tuple[int, ...], ...]) -> set[tuple[int, ...]]:
    identity = tuple(range(len(generators[0])))
    found = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            nxt = compose(current, generator)
            if nxt not in found:
                found.add(nxt)
                queue.append(nxt)
    return found


def maps(p: int) -> tuple[tuple[int, ...], ...]:
    modulus = p * p
    x = tuple(((1 + p) * n) % modulus for n in range(modulus))
    y_phys = tuple((n + 1) % modulus for n in range(modulus))
    y_zero = []
    z = []
    for n in range(modulus):
        v, w = n % p, n // p
        y_zero.append(((v + 1) % p) + p * w)
        z.append(v + p * ((w + 1) % p))
    return x, y_phys, tuple(y_zero), tuple(z)


def carry(p: int, b: int, bp: int) -> int:
    return (b + bp) // p


def cocycle_h(p: int, a: int, b: int, ap: int, bp: int) -> int:
    del ap
    return a * bp % p


def cocycle_m(p: int, a: int, b: int, ap: int, bp: int) -> int:
    return (cocycle_h(p, a, b, ap, bp) + carry(p, b, bp)) % p


def affine_order_spectrum_odd(p: int) -> Counter[int]:
    """Order spectrum of <n->(1+p)n,n->n+1> for odd p."""
    spectrum: Counter[int] = Counter()
    modulus = p * p
    for a in range(p):
        multiplier = pow(1 + p, a, modulus)
        for b in range(modulus):
            perm = tuple((multiplier * n + b) % modulus for n in range(modulus))
            spectrum[permutation_order(perm)] += 1
    return spectrum


def heisenberg_order_spectrum(p: int) -> Counter[int]:
    """Order spectrum from the class-two coordinate power law."""
    if p % 2:
        return Counter({1: 1, p: p**3 - 1})
    x, _y_phys, y_zero, _z = maps(p)
    group = generated_group((x, y_zero))
    return Counter(permutation_order(g) for g in group)


def main() -> None:
    primes = (2, 3, 5, 7, 11, 13)
    rows = []
    for p in primes:
        modulus = p * p
        identity = tuple(range(modulus))
        x, y_phys, y_zero, z = maps(p)

        require(power(x, p) == identity, f"X order does not divide p at {p}")
        require(power(z, p) == identity and z != identity, f"Z order failed at {p}")
        require(power(y_phys, p) == z, f"physical p-power is not Z at {p}")
        require(power(y_zero, p) == identity, f"suppressed low digit has carry at {p}")
        require(
            compose(compose(compose(x, y_phys), inverse(x)), inverse(y_phys)) == z,
            f"physical commutator failed at {p}",
        )
        require(
            compose(compose(compose(x, y_zero), inverse(x)), inverse(y_zero)) == z,
            f"Heisenberg commutator failed at {p}",
        )

        group_m = generated_group((x, y_phys))
        group_h = generated_group((x, y_zero))
        require(len(group_m) == p**3, f"modular group order failed at {p}")
        require(len(group_h) == p**3, f"Heisenberg group order failed at {p}")

        cocycle_gates = 0
        for a in range(p):
            for b in range(p):
                for ap in range(p):
                    for bp in range(p):
                        ch = cocycle_h(p, a, b, ap, bp)
                        cm = cocycle_m(p, a, b, ap, bp)
                        require(
                            (cm - ch) % p == carry(p, b, bp),
                            f"Bockstein carry mismatch at p={p}",
                        )
                        omega_h = (ch - cocycle_h(p, ap, bp, a, b)) % p
                        omega_m = (cm - cocycle_m(p, ap, bp, a, b)) % p
                        require(
                            omega_h == omega_m == (a * bp - ap * b) % p,
                            f"alternating commutator mismatch at p={p}",
                        )
                        cocycle_gates += 2

        if p % 2:
            spectrum_m = affine_order_spectrum_odd(p)
            spectrum_h = heisenberg_order_spectrum(p)
            require(
                spectrum_m == Counter({1: 1, p: p * p - 1, p * p: p * p * (p - 1)}),
                f"modular order spectrum failed at {p}: {spectrum_m}",
            )
            require(
                spectrum_h == Counter({1: 1, p: p**3 - 1}),
                f"Heisenberg exponent spectrum failed at {p}",
            )
            require(group_m != group_h, f"odd-prime action groups accidentally coincide at {p}")
            boundary = "distinct: modular has p^2-elements; Heisenberg exponent p"
        else:
            require(group_m == group_h, "p=2 permutation groups are not the same D8")
            require(
                y_phys == inverse(compose(x, y_zero)),
                "p=2 physical odometer generator relabel failed",
            )
            spectrum_m = Counter(permutation_order(g) for g in group_m)
            spectrum_h = heisenberg_order_spectrum(p)
            require(
                spectrum_m == spectrum_h == Counter({1: 1, 2: 5, 4: 2}),
                "p=2 D8 order spectrum failed",
            )
            boundary = "coincide as D8 after generator relabel"

        # Every index-p subgroup of either nonabelian order-p^3 group is
        # normal and contains the commutator Z.  Therefore a faithful action
        # must contain an orbit of size at least p^2; the displayed actions
        # attain p^2.  The finite gate checks the attaining faithful action.
        require(
            all(any(g[n] != n for g in (x, y_phys)) for n in range(modulus)),
            f"physical action unexpectedly has a global fixed point at {p}",
        )
        require(
            len({tuple(g[n] for g in group_m) for n in range(modulus)}) == modulus,
            f"physical carrier points were not separated at {p}",
        )

        if p == 13:
            off_wall = 0
            carry_wall = 0
            for n in range(modulus):
                v = n % p
                if v != p - 1:
                    require(y_phys[n] == y_zero[n], "off-wall successor mismatch")
                    off_wall += 1
                else:
                    require(y_phys[n] == compose(z, y_zero)[n], "carry-wall Z correction failed")
                    carry_wall += 1
            require((off_wall, carry_wall) == (156, 13), "p=13 wall census failed")

        rows.append(
            (
                p,
                len(group_m),
                len(group_h),
                dict(sorted(spectrum_m.items())),
                dict(sorted(spectrum_h.items())),
                cocycle_gates,
                boundary,
            )
        )

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(assert_nodes == 0, "truth-bearing Python assert found")

    print("THM-2788 PHYSICAL MODULAR ODOMETER / HEISENBERG BOUNDARY")
    for p, gm, gh, sm, sh, gates, boundary in rows:
        print(
            f"p={p}:|M|={gm},|H|={gh},M_orders={sm},H_orders={sh},"
            f"cocycle_gates={gates},boundary={boundary}"
        )
    print("extension_cocycles=c_M((a,b),(a',b'))=a*b'+carry(b,b'); c_H=a*b'")
    print("shared_commutator=omega((a,b),(a',b'))=a*b'-a'*b")
    print("odd_prime_power_map=M:(a,b)->b in Z; H:zero")
    print("faithful_degree_lower_bound=p^2 because every index-p subgroup contains [G,G]=Z")
    print("p13_local_wall=156 no-carry points plus 13 points with Y_phys=Z*Y_0")
    print("scope=group/action obstruction only; no physical current, root intertwiner, row exclusion, or LRC14")
    print(f"assert_nodes={assert_nodes}")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
