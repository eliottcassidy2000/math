#!/usr/bin/env python3
"""Exact companion for the proved global literal half-twist cap-seven theorem.

The global classification is an analytic strong-induction argument built from
proved predecessor theorems.  This standard-library companion freezes its
elementary interfaces and sharp boundary controls:

* direct strict-mask reconstruction and period pullback;
* the exact odd-cofactor fixed fibre for ``v_2(Q)=1,2``;
* every proposed support atom, its witness, and its exact rank through seven;
* the mixed-prime capacity list and the two dyadic exclusion ledgers;
* the residual ``2*17^b`` quotient-density and two-sided-mass arithmetic;
* direct no-cap negative searches at ``Q=6,34,578``.

The literal masks are evaluated at the fixed common centre ``1/(2Q)``.  The
companion does not classify arbitrary centres or construct a nonzero current.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3405",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3434",
        ROOT / "01-canon/theorems/THM-3434-seventeen-fibre-two-sided-mass-closure.md",
        "52a5a5caed75ad48ab35ce287c15f6cece074c88dddda2b87329792a3df70af7",
    ),
    (
        "THM-3435",
        ROOT / "01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md",
        "e2d202f97b1f39af032597767633c55f1abc4da261662824058b05b70724bb52",
    ),
    (
        "THM-3445",
        ROOT / "01-canon/theorems/THM-3445-prime-even-half-twist-cap-seven-classification.md",
        "06b28dd920ab525da087a6eba7bc3ecd25b51a878ce7e81f07725fb00e40b6bc",
    ),
    (
        "THM-3451",
        ROOT / "01-canon/theorems/THM-3451-prime-quarter-half-twist-cap-seven-classification.md",
        "87b7c2cf42c63f80d7be8e9f64ba8aff2dc72039026682ad912163d8d171c2a1",
    ),
)

SUPPORT_ATOMS = (8, 9, 10, 11, 12, 13, 14, 15, 23, 25, 29, 38, 51, 68, 148)
CAP_SIX_ATOMS = (8, 9, 10, 11, 12, 15, 23, 25)
EXACT_RANK = {
    8: 4, 9: 4, 10: 5, 11: 6, 12: 5, 13: 7, 14: 7, 15: 6,
    23: 6, 25: 6, 29: 7, 38: 7, 51: 7, 68: 7, 148: 7,
}
WITNESSES = {
    8: (1, 3, 5, 7),
    9: (1, 5, 6, 7),
    10: (1, 3, 4, 7, 9),
    11: (1, 2, 3, 5, 7, 9),
    12: (1, 5, 7, 8, 11),
    13: (1, 2, 3, 5, 7, 9, 11),
    14: (1, 3, 4, 5, 9, 11, 13),
    15: (1, 4, 6, 7, 8, 10),
    23: (1, 4, 5, 7, 9, 11),
    25: (1, 9, 10, 11, 19, 21),
    29: (1, 5, 7, 8, 12, 13, 22),
    38: (1, 9, 17, 20, 21, 29, 37),
    51: (1, 11, 12, 18, 23, 34, 35),
    68: (8, 11, 23, 24, 45, 56, 57),
    148: (8, 33, 41, 100, 107, 115, 140),
}
MIXED_PRIMES = (3, 5, 11, 17, 23, 29)
EXPECTED_SEMANTIC_SHA256 = "5f5d2c29d47a15abc8cf74d2b0f40769c220c3a646b4ba4c6bb21ca5d49cc1fa"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def in_arc(word: int, modulus: int) -> bool:
    residue = word % modulus
    return 14 * min(residue, modulus - residue) < modulus


def half_mask(modulus: int, coefficient: int) -> int:
    answer = 0
    for sheet in range(modulus):
        if in_arc(coefficient * (2 * sheet + 1), 2 * modulus):
            answer |= 1 << sheet
    return answer


def fraction_hit(modulus: int, coefficient: int, sheet: int) -> bool:
    phase = Fraction(coefficient * (2 * sheet + 1), 2 * modulus)
    fractional = phase - phase.numerator // phase.denominator
    return min(fractional, 1 - fractional) < Fraction(1, 14)


def reference_mask_certificate(limit: int = 32):
    rows = cells = endpoints = 0
    for modulus in range(2, limit + 1):
        for coefficient in range(1, 2 * modulus):
            direct = half_mask(modulus, coefficient)
            reference = 0
            for sheet in range(modulus):
                word = coefficient * (2 * sheet + 1) % (2 * modulus)
                if 14 * min(word, 2 * modulus - word) == 2 * modulus:
                    endpoints += 1
                if fraction_hit(modulus, coefficient, sheet):
                    reference |= 1 << sheet
                cells += 1
            require(direct == reference,
                    ("fraction reference", modulus, coefficient))
            rows += 1
    return limit, rows, cells, endpoints


def direct_union(modulus: int, coefficients: tuple[int, ...]) -> int:
    answer = 0
    for coefficient in coefficients:
        answer |= half_mask(modulus, coefficient)
    return answer


def odd_prime_factors(value: int) -> tuple[int, ...]:
    answer = []
    residual = value
    divisor = 3
    while divisor * divisor <= residual:
        if residual % divisor == 0:
            answer.append(divisor)
            while residual % divisor == 0:
                residual //= divisor
        divisor += 2
    if residual > 1:
        answer.append(residual)
    return tuple(answer)


def period_descent_certificate(limit: int = 96):
    rows = cells = 0
    for modulus in range(2, limit + 1):
        for lower in range(2, modulus + 1):
            if modulus % lower:
                continue
            scale = modulus // lower
            for coefficient in range(1, 2 * lower):
                base = half_mask(lower, coefficient)
                pulled = sum(
                    1 << sheet
                    for sheet in range(modulus)
                    if base >> (sheet % lower) & 1
                )
                direct = half_mask(modulus, scale * coefficient)
                require(direct == pulled,
                        ("period descent", modulus, lower, coefficient))
                rows += 1
                cells += modulus
    return limit, rows, cells


def fixed_odd_cofactor_certificate(odd_limit: int = 149):
    factor_rows = owner_rows = active_rows = sheet_cells = 0
    for exponent in (1, 2):
        dyadic = 1 << exponent
        for odd_core in range(3, odd_limit + 1, 2):
            modulus = dyadic * odd_core
            for prime in odd_prime_factors(odd_core):
                local_modulus = dyadic * prime
                odd_cofactor = odd_core // prime
                fixed_base = (odd_cofactor - 1) // 2
                require(2 * fixed_base + 1 == odd_cofactor,
                        ("fixed base", odd_core, prime))
                factor_rows += 1
                for coefficient in range(1, 2 * modulus):
                    if coefficient % prime:
                        active_rows += 1
                        require(coefficient % local_modulus != 0,
                                ("active local transversality", modulus,
                                 prime, coefficient))
                    for local_sheet in range(local_modulus):
                        global_sheet = fixed_base + odd_cofactor * local_sheet
                        direct = in_arc(
                            coefficient * (2 * global_sheet + 1),
                            2 * modulus,
                        )
                        local = in_arc(
                            coefficient * (2 * local_sheet + 1),
                            2 * local_modulus,
                        )
                        require(direct == local,
                                ("fixed fibre", modulus, prime,
                                 coefficient, local_sheet))
                        sheet_cells += 1
                    owner_rows += 1
    expected = (216, 103488, 89952, 11582552)
    result = (factor_rows, owner_rows, active_rows, sheet_cells)
    require(result == expected, ("fixed-fibre universe", result, expected))
    return odd_limit, result


def canonical_bank(modulus: int):
    grouped = {}
    for coefficient in range(1, modulus):
        mask = half_mask(modulus, coefficient)
        if not mask:
            continue
        old = grouped.get(mask)
        if old is None or coefficient < old:
            grouped[mask] = coefficient
    unique = tuple(sorted((mask, coefficient) for mask, coefficient in grouped.items()))
    maximal = tuple(
        item
        for item in unique
        if not any(
            item[0] != other[0] and item[0] | other[0] == other[0]
            for other in unique
        )
    )
    return unique, maximal


def exact_cover(full: int, items: tuple[tuple[int, int], ...], cap: int):
    masks = tuple(mask for mask, _coefficient in items)
    coefficients = tuple(coefficient for _mask, coefficient in items)
    width = full.bit_length()
    coverers = tuple(
        tuple(index for index, mask in enumerate(masks) if mask >> bit & 1)
        for bit in range(width)
    )
    joined = 0
    for mask in masks:
        joined |= mask
    if joined != full:
        return (), 1, 0, 0

    nodes = branches = 0

    @lru_cache(maxsize=None)
    def solve(state: int, slots: int):
        nonlocal nodes, branches
        nodes += 1
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted(
            ((mask & missing).bit_count() for mask in masks), reverse=True
        )
        if sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing >> bit & 1)
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(masks[index] | state != state for index in coverers[bit]),
                bit,
            ),
        )
        options = sorted(
            (
                index
                for index in coverers[pivot]
                if masks[index] | state != state
            ),
            key=lambda index: (
                -(masks[index] & missing).bit_count(), coefficients[index]
            ),
        )
        for index in options:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    chosen = solve(0, cap)
    hits = solve.cache_info().hits
    if chosen is None:
        return (), nodes, branches, hits
    witness = tuple(sorted(coefficients[index] for index in chosen))
    return witness, nodes, branches, hits


def atom_certificate():
    reports = []
    for modulus in SUPPORT_ATOMS:
        supplied = WITNESSES[modulus]
        full = (1 << modulus) - 1
        require(len(set(supplied)) == len(supplied) <= 7,
                ("atom owner count", modulus, supplied))
        require(all(0 < coefficient < modulus for coefficient in supplied),
                ("atom sign representatives", modulus, supplied))
        require(direct_union(modulus, supplied) == full,
                ("atom witness", modulus, supplied))
        require(gcd(modulus, *supplied) == 1,
                ("atom joint period", modulus, supplied))

        unique, maximal = canonical_bank(modulus)
        rank = EXACT_RANK[modulus]
        # Search the complete unique-mask bank.  The maximal-bank count is
        # reported as a hostile control, but no dominance pruning is used for
        # the exact-rank verdict.
        below = exact_cover(full, unique, rank - 1)
        exact = exact_cover(full, unique, rank)
        require(not below[0] and exact[0],
                ("exact atom rank", modulus, rank, below, exact))
        require(len(exact[0]) == rank,
                ("exact witness length", modulus, exact))
        sizes = tuple(half_mask(modulus, coefficient).bit_count()
                      for coefficient in supplied)
        multiplicities = [0] * modulus
        for coefficient in supplied:
            mask = half_mask(modulus, coefficient)
            for sheet in range(modulus):
                multiplicities[sheet] += (mask >> sheet) & 1
        distribution = tuple(
            (value, multiplicities.count(value))
            for value in sorted(set(multiplicities))
        )
        reports.append((
            modulus, rank, len(unique), len(maximal), supplied, sizes,
            distribution, below[1:], exact,
        ))
    require(tuple(EXACT_RANK) == SUPPORT_ATOMS, "rank/support ordering")
    require(
        all(
            not any(modulus % base == 0 for base in CAP_SIX_ATOMS)
            for modulus in (13, 14, 29, 38, 51, 68, 148)
        ),
        "rank-seven atom contains cap-six base",
    )
    return tuple(reports)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor = 3 if divisor == 2 else divisor + 2
    return True


def mixed_prime_certificate():
    survivors = tuple(
        prime
        for prime in range(3, 1000, 2)
        if is_prime(prime) and prime <= 6 * ((prime + 6) // 7)
    )
    require(survivors == MIXED_PRIMES, ("mixed primes", survivors))
    require(6 * ((31 + 6) // 7) < 31, "p=31 stopping boundary")
    require(
        all(6 * (prime + 6) < 7 * prime for prime in (37, 101, 997)),
        "large-prime capacity cutoff",
    )

    residual = {}
    for exponent in (1, 2):
        dyadic = 1 << exponent
        residual[exponent] = tuple(
            prime
            for prime in MIXED_PRIMES
            if not any((dyadic * prime) % atom == 0 for atom in SUPPORT_ATOMS)
        )
    require(residual == {1: (3, 17), 2: ()},
            ("dyadic mixed residual", residual))

    all_active_support = {1: (7, 19), 2: (17, 37)}
    local_atoms = {
        exponent: tuple((1 << exponent) * prime for prime in primes)
        for exponent, primes in all_active_support.items()
    }
    require(local_atoms == {1: (14, 38), 2: (68, 148)}, local_atoms)
    return survivors, tuple(sorted(residual.items())), tuple(sorted(local_atoms.items()))


def tower_certificate(max_exponent: int = 8):
    density_rows = []
    invoice_rows = []
    for exponent in range(1, max_exponent + 1):
        odd_order = 17 ** exponent
        even_order = 2 * odd_order
        even_ceiling = (even_order + 6) // 7
        require(17 * even_ceiling <= 3 * even_order,
                ("even quotient density", exponent, even_order, even_ceiling))
        require((2 * odd_order + 6) // 7 <= 2 * ((odd_order + 6) // 7),
                ("active size", exponent, odd_order))
        density_rows.append((
            exponent, even_order, even_ceiling,
            3 * even_order - 17 * even_ceiling,
        ))

    require(half_mask(2, 1) == 0, "order-two quotient is not empty")
    for exponent in range(1, max_exponent + 1):
        base_count = 2 * 17 ** (exponent - 1)
        active_order = 17 ** exponent
        left_numerator = 274 * base_count
        active_upper = 12 * ((active_order + 6) // 7)
        right_numerator = 17 * active_upper
        invoice_rows.append((
            exponent, base_count, left_numerator, right_numerator,
            left_numerator > right_numerator,
        ))
    require(not invoice_rows[0][-1], ("b=1 boundary", invoice_rows[0]))
    require(all(row[-1] for row in invoice_rows[1:]),
            ("tower tail", invoice_rows))
    require(184 * 34 > 1224, "uniform b>=2 contradiction")
    return (
        tuple(density_rows), tuple(invoice_rows),
        (184, 1224, 34),
    )


def negative_boundary_certificate():
    reports = []
    expected = {
        6: (1, 1, 1, 0, 0),
        34: (32, 24, 154, 164, 11),
        578: (576, 568, 10413, 10412, 0),
    }
    for modulus in (6, 34, 578):
        unique, maximal = canonical_bank(modulus)
        result = exact_cover((1 << modulus) - 1, unique, 7)
        require(not result[0], ("negative boundary cover", modulus, result[0]))
        row = (len(unique), len(maximal), result[1], result[2], result[3])
        require(row == expected[modulus],
                ("negative boundary search", modulus, row, expected[modulus]))
        union = 0
        for mask, _coefficient in unique:
            union |= mask
        reports.append((modulus, row, union.bit_count(), modulus - union.bit_count()))

    q6_union = direct_union(6, tuple(range(1, 6)))
    q6_support = tuple(sheet for sheet in range(6) if q6_union >> sheet & 1)
    require(q6_support == (1, 4), ("Q6 support", q6_support))
    tiny = tuple(
        (modulus, direct_union(modulus, tuple(range(1, modulus))).bit_count())
        for modulus in (2, 4)
    )
    require(tiny == ((2, 0), (4, 0)), ("tiny periods", tiny))
    return tuple(reports), q6_support, tiny


def security_certificate(source: Path):
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    require(not (called & forbidden_calls),
            ("forbidden calls", called & forbidden_calls))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "fractions", "functools", "hashlib", "math",
        "pathlib",
    }
    require(imported <= allowed, ("unexpected imports", imported - allowed))
    return tuple(sorted(imported)), tuple(sorted(forbidden_calls))


def main() -> None:
    source = Path(__file__)
    dependency_hashes = tuple(
        (label, lf_sha256(path)) for label, path, _expected in DEPENDENCIES
    )
    expected_hashes = tuple(
        (label, expected) for label, _path, expected in DEPENDENCIES
    )
    require(dependency_hashes == expected_hashes,
            ("dependency hash drift", dependency_hashes, expected_hashes))
    security = security_certificate(source)
    reference = reference_mask_certificate()
    period_descent = period_descent_certificate()
    fixed_fibre = fixed_odd_cofactor_certificate()
    atoms = atom_certificate()
    mixed = mixed_prime_certificate()
    tower = tower_certificate()
    negatives = negative_boundary_certificate()

    semantic_payload = (
        dependency_hashes, security, reference, period_descent, fixed_fibre,
        atoms, mixed, tower, negatives,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("THM-3453 global literal half-twist cap-seven support companion")
    print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"support_atoms={SUPPORT_ATOMS}")
    print(f"dependency_hashes={dependency_hashes}")
    print(f"security_certificate={security}")
    print(f"strict_fraction_reference=(limit,rows,cells,endpoints)={reference}")
    print(f"period_descent=(limit,rows,cells)={period_descent}")
    print(f"fixed_odd_cofactor_fibre=(odd_limit,(factor_rows,owner_rows,active_rows,sheet_cells))={fixed_fibre}")
    print("atom_exact_ranks=" + repr(tuple(
        (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[8][0])
        for row in atoms
    )))
    print(f"mixed_prime_sieve=(capacity_primes,residual_by_v2,all_active_local_atoms)={mixed}")
    print("tower_density=(b,m_even,ceil(m_even/7),density_slack)="
          f"{tower[0]}")
    print("tower_invoice=(b,N,274N,17_active_upper,contradiction)="
          f"{tower[1]}")
    print(f"tower_uniform=(left_coefficient,right_constant,N_min)={tower[2]}")
    print(f"negative_boundaries=(searches,Q6_support,Q2_Q4)={negatives}")
    print("classification_candidate=rho_H(Q)<=7 iff some support atom divides Q")
    print("scope=fixed half-twist common centre and zero cochain; no arbitrary-centre classification, nonzero current, decrement, or LRC(14) consequence")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
