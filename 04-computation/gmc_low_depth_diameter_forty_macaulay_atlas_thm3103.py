#!/usr/bin/env python3
"""Finite exact companion for THM-3103.

For every translated four-slot factorial support

    (n,n+a,n+b,n+M),  0 <= n <= 3,  0 < a < b < M <= 40,

the script constructs the normalized quadratic, cubic, and quartic forms
after eliminating the first moment, and evaluates one fixed 36-by-36
degree-seven Macaulay minor over F_1000003.  A nonzero modular determinant
is an exact certificate that the corresponding rational minor is nonzero.

The fast constructor expands the original four coefficient variables before
the linear elimination.  A pinned, independently organized inclusion-
exclusion/rising-factorial constructor from THM-2921 checks representative
forms and minors.  Repeated-exponent controls make the fixed minor vanish.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from itertools import combinations
from math import factorial, isqrt
from pathlib import Path

try:
    from flint import nmod_mat
except ModuleNotFoundError as error:
    raise RuntimeError("THM-3103 exact replay requires python-flint") from error


PRIME = 1_000_003
DEPTHS = tuple(range(4))
MAX_DIAMETER = 40
ORDERS = (2, 3, 4)
TARGET_DEGREE = 7
SELECTED_ROWS = (
    tuple(range(20))
    + tuple(range(21, 30))
    + (35,)
    + tuple(range(36, 42))
)

HERE = Path(__file__).resolve().parent
CONSTRUCTOR = HERE / "gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py"
CONSTRUCTOR_SHA256 = (
    "42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_digest(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value) + 1, 2))


def compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first,) + tail
        for first in range(total + 1)
        for tail in compositions(total - first, parts - 1)
    )


MONOMIALS3 = {
    degree: compositions(degree, 3)
    for degree in range(TARGET_DEGREE + 1)
}
TARGET_MONOMIALS = MONOMIALS3[TARGET_DEGREE]
TARGET_INDEX = {
    monomial: index
    for index, monomial in enumerate(TARGET_MONOMIALS)
}


def multinomial(exponents: tuple[int, ...]) -> int:
    answer = factorial(sum(exponents))
    for exponent in exponents:
        answer //= factorial(exponent)
    return answer


# Each source count vector contributes a scalar depending on the support and
# a fixed signed expansion of c_3^d=(-x-y-z)^d.
COUNT_BANK = {}
for order in ORDERS:
    entries = []
    for counts in compositions(order, 4):
        expansions = []
        for tail in compositions(counts[3], 3):
            target = tuple(counts[index] + tail[index] for index in range(3))
            signed_copy = multinomial(tail)
            if counts[3] % 2:
                signed_copy = -signed_copy
            expansions.append((target, signed_copy % PRIME))
        entries.append((counts, multinomial(counts) % PRIME, tuple(expansions)))
    COUNT_BANK[order] = tuple(entries)


MAX_EXPONENT = max(DEPTHS) + MAX_DIAMETER
FACTORIAL_MOD = [1]
for value in range(1, 4 * MAX_EXPONENT + 1):
    FACTORIAL_MOD.append(FACTORIAL_MOD[-1] * value % PRIME)
INV_FACTORIAL_MOD = [pow(value, PRIME - 2, PRIME) for value in FACTORIAL_MOD]


def direct_form_mod(
    depth: int,
    order: int,
    offsets: tuple[int, int, int, int],
) -> dict[tuple[int, int, int], int]:
    """Expand four variables, then substitute c_3=-(x+y+z)."""
    support = tuple(depth + offset for offset in offsets)
    answer = {monomial: 0 for monomial in MONOMIALS3[order]}
    base_inverse = (
        pow(FACTORIAL_MOD[depth], order, PRIME)
        * INV_FACTORIAL_MOD[order * depth]
    ) % PRIME
    support_inverse = tuple(INV_FACTORIAL_MOD[value] for value in support)
    for counts, copies, expansions in COUNT_BANK[order]:
        total_exponent = sum(
            counts[index] * support[index]
            for index in range(4)
        )
        coefficient = copies * FACTORIAL_MOD[total_exponent] % PRIME
        coefficient = coefficient * base_inverse % PRIME
        for index in range(4):
            coefficient = (
                coefficient
                * pow(support_inverse[index], counts[index], PRIME)
            ) % PRIME
        for target, signed_copy in expansions:
            answer[target] = (
                answer[target] + coefficient * signed_copy
            ) % PRIME
    return answer


ROW_DESCRIPTORS = []
for order in ORDERS:
    for multiplier in MONOMIALS3[TARGET_DEGREE - order]:
        ROW_DESCRIPTORS.append((order, multiplier))
require(len(ROW_DESCRIPTORS) == 46, "Macaulay row count changed")
SELECTED_DESCRIPTORS = tuple(ROW_DESCRIPTORS[index] for index in SELECTED_ROWS)


def selected_matrix(
    depth: int,
    offsets: tuple[int, int, int, int],
) -> list[list[int]]:
    forms = {
        order: direct_form_mod(depth, order, offsets)
        for order in ORDERS
    }
    rows = []
    for order, multiplier in SELECTED_DESCRIPTORS:
        row = [0] * len(TARGET_MONOMIALS)
        for monomial, coefficient in forms[order].items():
            target = tuple(
                multiplier[index] + monomial[index]
                for index in range(3)
            )
            row[TARGET_INDEX[target]] = coefficient
        rows.append(row)
    require(len(rows) == len(TARGET_MONOMIALS) == 36, "selected chart changed")
    return rows


def selected_determinant(
    depth: int,
    offsets: tuple[int, int, int, int],
) -> int:
    return int(nmod_mat(selected_matrix(depth, offsets), PRIME).det())


def load_independent_constructor():
    require(CONSTRUCTOR.is_file(), "missing THM-2921 constructor")
    require(
        lf_digest(CONSTRUCTOR) == CONSTRUCTOR_SHA256,
        "THM-2921 constructor hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm2921_for_3103", CONSTRUCTOR)
    require(spec is not None and spec.loader is not None, "import spec failed")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def independent_constructor_audit() -> int:
    source = load_independent_constructor()
    checks = 0
    controls = (
        (0, (0, 1, 2, 3)),
        (0, (0, 7, 19, 40)),
        (1, (0, 3, 27, 40)),
        (2, (0, 13, 37, 40)),
        (3, (0, 1, 39, 40)),
    )
    for depth, offsets in controls:
        local_forms = tuple(
            direct_form_mod(depth, order, offsets)
            for order in ORDERS
        )
        independent_forms = []
        for order in ORDERS:
            rational_form = source.moment_form(depth, order, offsets)
            reduced_form = {
                monomial: (
                    value.numerator
                    * pow(value.denominator % PRIME, PRIME - 2, PRIME)
                ) % PRIME
                for monomial, value in rational_form.items()
            }
            require(
                reduced_form == local_forms[order - 2],
                f"constructor mismatch: depth={depth}, offsets={offsets}, order={order}",
            )
            independent_forms.append(reduced_form)
            checks += 1
        independent_rows = source.macaulay_rows_mod(tuple(independent_forms))
        independent_minor = source.determinant_mod(
            [independent_rows[index] for index in SELECTED_ROWS]
        )
        require(
            independent_minor == selected_determinant(depth, offsets),
            f"minor mismatch: depth={depth}, offsets={offsets}",
        )
        checks += 1
    return checks


def digest_update(digest, record: str) -> None:
    digest.update(record.encode("ascii"))
    digest.update(b"\n")


def main() -> None:
    require(is_prime(PRIME), "audit modulus is not prime")
    require(4 * MAX_EXPONENT < PRIME, "factorial denominator meets modulus")
    require(SELECTED_ROWS == tuple(range(20)) + tuple(range(21, 30)) + (35,) + tuple(range(36, 42)), "selected rows changed")

    constructor_checks = independent_constructor_audit()

    # Repeated exponents make the eliminated system effectively binary, so a
    # full ternary degree-seven Macaulay chart must vanish.  These are exact
    # negative controls outside the distinct-support universe.
    repeated_controls = (
        (0, (0, 1, 1, 3)),
        (3, (0, 7, 7, 40)),
        (2, (0, 0, 19, 40)),
    )
    for depth, offsets in repeated_controls:
        require(
            selected_determinant(depth, offsets) == 0,
            f"repeated-exponent hostile did not vanish: {depth}, {offsets}",
        )

    total = 0
    global_digest = sha256()
    depth_records = []
    diameter_counts = {diameter: 0 for diameter in range(3, MAX_DIAMETER + 1)}
    for depth in DEPTHS:
        count = 0
        residue_digest = sha256()
        for diameter in range(3, MAX_DIAMETER + 1):
            for first, second in combinations(range(1, diameter), 2):
                offsets = (0, first, second, diameter)
                residue = selected_determinant(depth, offsets)
                require(
                    residue != 0,
                    f"zero fixed minor: depth={depth}, offsets={offsets}",
                )
                record = f"{depth}:{first}:{second}:{diameter}:{residue}"
                digest_update(residue_digest, record)
                digest_update(global_digest, record)
                count += 1
                total += 1
                diameter_counts[diameter] += 1
        require(count == 9880, f"depth-{depth} census changed")
        depth_records.append((depth, count, residue_digest.hexdigest()))

    require(total == 39_520, "total atlas census changed")
    for diameter, count in diameter_counts.items():
        require(
            count == 4 * ((diameter - 1) * (diameter - 2) // 2),
            f"diameter-{diameter} census changed",
        )

    print("THM-3103 LOW-DEPTH DIAMETER-FORTY MACAULAY ATLAS")
    print(
        f"prime={PRIME} selected_rows=36 constructor_checks={constructor_checks} "
        f"repeated_zero_controls={len(repeated_controls)}"
    )
    print(
        f"universe=0<=n<=3;0<a<b<M<=40 total={total} "
        "zero_minors=0 unresolved=0"
    )
    for depth, count, digest in depth_records:
        print(f"depth={depth} supports={count} residue_sha256={digest}")
    print(
        f"diameter3_count={diameter_counts[3]} "
        f"diameter40_count={diameter_counts[40]} "
        f"global_residue_sha256={global_digest.hexdigest()}"
    )
    print(
        "scope=finite_exact_modular_nonvanishing;fixed_macaulay_chart;"
        "no_sign;no_all_depth;no_width5"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
