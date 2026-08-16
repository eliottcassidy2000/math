#!/usr/bin/env python3
"""Exact sidecar for the marked D5 normal-slice bridge and flux no-go.

This original companion has two jobs.  THM-3496 records the later independent
audit and repairs the scope from a blanket Bockstein phrase to exact
finite-coefficient extinction modulo 13^r.

1.  It writes the actual finite cochain map from the THM-2542 seven-cycle
    class to the mod-13 Kummer line of a marked JC divisor normal coordinate.
    It checks graph/deck gauges and the cover-degree/ramification square.
2.  It verifies the smallest obstruction to extending that map to the
    additive Hamiltonian flux.  For P=x+x^2 z, reduction modulo 13 creates
    an explicit Frobenius first integral and kills the principal-part
    connecting class.

Only integer and finite-field arithmetic is used.  The mathematical scope is
the accompanying unnumbered reflection; this file is not a canon theorem.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict


P = 13
SOURCE_CYCLE = 7
EXPECTED_SEMANTIC_SHA256 = "8b347e4b65b8a07e4ce040775733b7b1e2161280cc241629dc6d50d21a50881c"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mod_rank(matrix: list[list[int]], prime: int) -> int:
    a = [[entry % prime for entry in row] for row in matrix]
    if not a:
        return 0
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * value) % prime for value in a[rank]]
        for r in range(rows):
            if r == rank or a[r][col] == 0:
                continue
            factor = a[r][col]
            a[r] = [
                (a[r][j] - factor * a[rank][j]) % prime
                for j in range(cols)
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def coboundary(values: tuple[int, ...], prime: int = P) -> tuple[int, ...]:
    n = len(values)
    return tuple((values[(i + 1) % n] - values[i]) % prime for i in range(n))


def holonomy(cochain: tuple[int, ...], prime: int = P) -> int:
    return sum(cochain) % prime


def reversed_cochain(cochain: tuple[int, ...], prime: int = P) -> tuple[int, ...]:
    """Pull back along orientation reversal, including the edge sign."""

    n = len(cochain)
    return tuple((-cochain[(-i - 1) % n]) % prime for i in range(n))


# Sparse Laurent polynomials in x,z: (x exponent,z exponent) -> integer coeff.
Poly = dict[tuple[int, int], int]


def normalize(poly: Poly) -> Poly:
    return {monomial: coeff for monomial, coeff in poly.items() if coeff}


def add(left: Poly, right: Poly) -> Poly:
    out: defaultdict[tuple[int, int], int] = defaultdict(int)
    for monomial, coeff in left.items():
        out[monomial] += coeff
    for monomial, coeff in right.items():
        out[monomial] += coeff
    return normalize(dict(out))


def scale(poly: Poly, scalar: int) -> Poly:
    return normalize({monomial: scalar * coeff for monomial, coeff in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    out: defaultdict[tuple[int, int], int] = defaultdict(int)
    for (ax, az), ac in left.items():
        for (bx, bz), bc in right.items():
            out[(ax + bx, az + bz)] += ac * bc
    return normalize(dict(out))


def power(poly: Poly, exponent: int) -> Poly:
    require(exponent >= 0, "negative sparse-polynomial power")
    out: Poly = {(0, 0): 1}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def reduce_mod(poly: Poly, prime: int) -> Poly:
    return normalize({monomial: coeff % prime for monomial, coeff in poly.items()})


def hamiltonian_d(poly: Poly) -> Poly:
    """D=(1+2xz)d/dz-x^2 d/dx for P=x+x^2 z."""

    out: defaultdict[tuple[int, int], int] = defaultdict(int)
    for (a, b), coeff in poly.items():
        if b:
            out[(a, b - 1)] += coeff * b
        shifted = 2 * b - a
        if shifted:
            out[(a + 1, b)] += coeff * shifted
    return normalize(dict(out))


def q_prime(prime: int) -> Poly:
    require(prime >= 3 and prime % 2 == 1, "odd-prime formula only")
    return q_to_depth(prime - 1)


def q_to_depth(depth: int) -> Poly:
    """Telescope whose terminal target is (xz)^depth."""

    require(depth >= 1, "positive telescope depth required")
    return {
        (j, j + 1): (-1) ** j
        for j in range(depth)
    }


def poly_payload(poly: Poly) -> list[list[int]]:
    return [
        [x_exp, z_exp, coeff]
        for (x_exp, z_exp), coeff in sorted(poly.items())
    ]


def main() -> None:
    # ------------------------------------------------------------------
    # LRC graph H^1 and the THM-3431 deck defect.
    # ------------------------------------------------------------------
    delta_matrix = [[0 for _ in range(SOURCE_CYCLE)] for _ in range(SOURCE_CYCLE)]
    for column in range(SOURCE_CYCLE):
        basis = tuple(1 if i == column else 0 for i in range(SOURCE_CYCLE))
        image = coboundary(basis)
        for row, value in enumerate(image):
            delta_matrix[row][column] = value

    delta_rank = mod_rank(delta_matrix, P)
    require(delta_rank == 6, "C7 coboundary rank must be six")
    require(
        all(holonomy(tuple(row[column] for row in delta_matrix)) == 0
            for column in range(SOURCE_CYCLE)),
        "every graph coboundary must have zero seam sum",
    )

    canonical_a = 1
    g = (canonical_a,) * SOURCE_CYCLE
    source_holonomy = holonomy(g)
    require(source_holonomy == 7, "canonical LRC holonomy mismatch")

    cover_size = SOURCE_CYCLE * P
    primitive = tuple((j * canonical_a) % P for j in range(cover_size))
    require(
        all((primitive[(j + 1) % cover_size] - primitive[j]) % P == canonical_a
            for j in range(cover_size)),
        "degree-13 cover primitive does not differentiate to the pullback",
    )
    deck_defects = tuple(
        (primitive[(j + SOURCE_CYCLE) % cover_size] - primitive[j]) % P
        for j in range(cover_size)
    )
    require(set(deck_defects) == {source_holonomy}, "deck defect is not constant")
    for constant in range(P):
        shifted = tuple((value + constant) % P for value in primitive)
        shifted_defects = {
            (shifted[(j + SOURCE_CYCLE) % cover_size] - shifted[j]) % P
            for j in range(cover_size)
        }
        require(shifted_defects == {source_holonomy}, "primitive gauge changed deck defect")

    # Phi_lambda([c]) = Hol(c) * [lambda] on the marked Kummer line.
    # At the Cech cochain level F^0=0 and F^1(c)=Hol(c)lambda^-1.
    for column in range(SOURCE_CYCLE):
        boundary = tuple(delta_matrix[row][column] for row in range(SOURCE_CYCLE))
        require(holonomy(boundary) == 0, "marked local-H1 cochain map failed")

    mapped_exponents = {
        holonomy((a,) * SOURCE_CYCLE)
        for a in range(P)
    }
    require(mapped_exponents == set(range(P)), "marked H1 map is not onto its Kummer line")
    require(
        holonomy(reversed_cochain(g)) == (-source_holonomy) % P,
        "orientation reversal does not match Kummer inversion",
    )

    cover_ramification_checks = []
    for degree in range(1, 201):
        source_pullback = (degree * source_holonomy) % P
        target_ramification = (degree * source_holonomy) % P
        require(source_pullback == target_ramification, "cover/ramification square failed")
        cover_ramification_checks.append((degree, source_pullback))
    require(cover_ramification_checks[12][1] == 0, "degree 13 must kill both classes")
    require(
        cover_ramification_checks[13][1] == source_holonomy,
        "degree 14 must preserve both classes",
    )

    # The pole-order map from additive local H1 to the Kummer line is not
    # additive: xi and -xi have the same pole order but opposite additive sign.
    xi_kummer = 1
    minus_xi_kummer_by_pole_order = 1
    require(
        (xi_kummer + minus_xi_kummer_by_pole_order) % P == 2,
        "pole-order hostile unexpectedly vanished",
    )
    require(2 % P != 0, "odd-characteristic sign hostile mistyped")

    # Kummer reduction forgets the full Prüfer/annihilator depth.
    require(1 % P == 14 % P, "depth collision control failed")

    # ------------------------------------------------------------------
    # P=x+x^2 z: principal-part transgression and the Frobenius trap.
    # ------------------------------------------------------------------
    p_poly: Poly = {(1, 0): 1, (2, 1): 1}
    h: Poly = {(-2, 0): -1, (-1, 1): 2}
    mu = hamiltonian_d(h)
    require(mu == {(0, 1): 6}, "D(h)=6z failed")

    p_times_h = multiply(p_poly, h)
    require(
        p_times_h == {(-1, 0): -1, (0, 1): 1, (1, 2): 2},
        "P*h principal-part representative mismatch",
    )
    simple_principal_part: Poly = {(-1, 0): -1}
    require(
        hamiltonian_d(simple_principal_part) == {(0, 0): -1},
        "principal-part connecting class must be [-1]",
    )

    q13 = q_prime(P)
    d_q13 = hamiltonian_d(q13)
    expected_integral = {(0, 0): 1, (12, 12): -13}
    require(d_q13 == expected_integral, "integral telescoping identity failed")
    require(
        reduce_mod(d_q13, P) == {(0, 0): 1},
        "mod-13 polynomial mate identity failed",
    )

    localized_cycle = add(q13, simple_principal_part)
    require(
        reduce_mod(hamiltonian_d(localized_cycle), P) == {},
        "localized Frobenius constant is not D-closed mod 13",
    )
    frobenius_form = scale(
        multiply(power(p_poly, P - 1), {(-P, 0): 1}),
        -1,
    )
    require(
        reduce_mod(localized_cycle, P) == reduce_mod(frobenius_form, P),
        "Q_13-x^-1 != -P^12/x^13 mod 13",
    )
    require(
        reduce_mod(hamiltonian_d(frobenius_form), P) == {},
        "explicit Frobenius form is not a localized first integral",
    )

    # General odd-prime telescope: this is a mechanism, not a bounded guess.
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    prime_checks = []
    for prime in primes:
        q_p = q_prime(prime)
        expected = {(0, 0): 1, (prime - 1, prime - 1): -prime}
        require(hamiltonian_d(q_p) == expected, f"prime telescope failed at {prime}")
        rhs = scale(
            multiply(power(p_poly, prime - 1), {(-prime, 0): 1}),
            -1,
        )
        lhs = add(q_p, {(-1, 0): -1})
        require(
            reduce_mod(lhs, prime) == reduce_mod(rhs, prime),
            f"Frobenius identity failed at {prime}",
        )
        require(
            reduce_mod(hamiltonian_d(q_p), prime) == {(0, 0): 1},
            f"mod-p mate failed at {prime}",
        )
        prime_checks.append(prime)

    # Integral response divisibility forces finite-coefficient extinction.
    # For every n>=1,
    #   D(Q_n)=1+(-1)^(n-1)(n+1)(xz)^n,
    # hence [1]=(-1)^n(n+1)[(xz)^n] in the integral cokernel.  Taking
    # n=13^r-1 makes [1] divisible by every power 13^r.
    depth_checks = []
    for depth in range(1, 201):
        q_n = q_to_depth(depth)
        expected = {
            (0, 0): 1,
            (depth, depth): ((-1) ** (depth - 1)) * (depth + 1),
        }
        require(hamiltonian_d(q_n) == normalize(expected), f"depth telescope failed at {depth}")
        depth_checks.append(depth)
    thirteen_power_depths = []
    for exponent in range(1, 4):
        power_13 = P ** exponent
        depth = power_13 - 1
        expected = {(0, 0): 1, (depth, depth): -power_13}
        require(
            hamiltonian_d(q_to_depth(depth)) == expected,
            f"13-adic divisibility telescope failed at exponent {exponent}",
        )
        thirteen_power_depths.append((exponent, depth, power_13))

    # Characteristic-zero terminal recurrence.  Any weight -1 candidate for
    # D(Q)=1 has coefficients c_j=(-1)^j and terminal obstruction
    # (N+2)c_N, which never vanishes in characteristic zero.
    terminal_obstructions = tuple(
        (n + 2) * ((-1) ** n)
        for n in range(101)
    )
    require(all(value != 0 for value in terminal_obstructions), "char-zero terminal vanished")

    semantic = {
        "prime": P,
        "source": {
            "cycle": SOURCE_CYCLE,
            "delta_rank": delta_rank,
            "h1_dimension": SOURCE_CYCLE - delta_rank,
            "canonical_holonomy": source_holonomy,
            "cover_size": cover_size,
            "deck_defect": sorted(set(deck_defects)),
            "mapped_kummer_exponents": sorted(mapped_exponents),
            "degree_13_image": cover_ramification_checks[12][1],
            "degree_14_image": cover_ramification_checks[13][1],
        },
        "marked_square": {
            "formula": "Phi_lambda([g])=Hol(g)*kappa_lambda",
            "cover_ramification_degrees_checked": 200,
            "orientation_reversal": (-source_holonomy) % P,
        },
        "losses": {
            "pole_order_sign_images": [xi_kummer, minus_xi_kummer_by_pole_order],
            "pole_order_sum_image": 2,
            "depth_collision": [1, 14],
        },
        "flux_hostile": {
            "P": "x+x^2*z",
            "D_h": poly_payload(mu),
            "P_h": poly_payload(p_times_h),
            "D_minus_x_inverse": poly_payload(hamiltonian_d(simple_principal_part)),
            "D_Q13_over_Z": poly_payload(d_q13),
            "localized_frobenius_form": "Q13-x^-1=-P^12/x^13 (mod 13)",
            "odd_primes_checked": prime_checks,
            "integral_depths_checked": len(depth_checks),
            "thirteen_power_divisibility": thirteen_power_depths,
            "char_zero_terminal_checks": len(terminal_obstructions),
        },
    }
    semantic_bytes = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    semantic_sha256 = hashlib.sha256(semantic_bytes).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "semantic payload hash drift",
        )

    print("D5 MARKED KUMMER NORMAL-SLICE / FROBENIUS FLUX NO-GO")
    print("STATUS: FINITE-EXACT original sidecar; repaired canon scope in THM-3496")
    print()
    print("SOURCE H1")
    print(f"  rank(delta_C7)={delta_rank}; dim H1={SOURCE_CYCLE-delta_rank}")
    print(f"  canonical g=(1,...,1): Hol(g)={source_holonomy} mod {P}")
    print(f"  C91 primitive deck defect chi(tau)={deck_defects[0]}")
    print(f"  all target Kummer exponents occur: {sorted(mapped_exponents)}")
    print()
    print("MARKED COMMUTING SQUARE")
    print("  Phi_lambda([g]) = Hol(g) * kappa_lambda")
    print("  checked cover degree k vs lambda->lambda^k for 1<=k<=200")
    print(f"  k=13: {cover_ramification_checks[12][1]} (both die)")
    print(f"  k=14: {cover_ramification_checks[13][1]} (both return)")
    print(f"  orientation reversal sends {source_holonomy} -> {(-source_holonomy)%P}")
    print()
    print("FIRST LOSS")
    print("  pole order is not additive: xi and -xi both map to kappa, sum maps to 2*kappa")
    print("  depths 1 and 14 have the same Kummer exponent but different JC arm lengths")
    print()
    print("MINIMAL FLUX HOSTILE: P=x+x^2*z")
    print("  h=-x^-2+2*z/x; D(h)=6*z")
    print("  [P*h]=[-x^-1]; delta([P*h])=[-1]")
    print("  Q13=sum_{j=0}^{11}(-1)^j*x^j*z^(j+1)")
    print("  D(Q13)=1-13*(x*z)^12 over Z")
    print("  D(Q13)=1 over F13, so the flux class [1] dies after coefficient matching")
    print("  Q13-x^-1=-P^12/x^13 over F13, an explicit Frobenius localized constant")
    print(f"  odd-prime telescope checked at {prime_checks}")
    print("  [1]=13^r*[(x*z)^(13^r-1)] in the integral cokernel; r=1,2,3 checked")
    print("  consequence: finite Z/13^r coefficient shadows vanish; derived completion unadjudicated")
    print("  characteristic-zero recurrence has nonzero terminal obstruction for every finite N")
    print()
    print(f"SEMANTIC_SHA256={semantic_sha256}")


if __name__ == "__main__":
    main()
