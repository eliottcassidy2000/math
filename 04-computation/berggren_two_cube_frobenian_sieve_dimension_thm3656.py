#!/usr/bin/env python3
"""Exact local companion for THM-3656's dimension-3/2 sieve.

The analytic upper-bound sieve is applied in the theorem.  This script pins
the local bad-line geometry, the joint mod-24 character table, and CRT
multiplicativity on exact finite hostile controls.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT = (
    ROOT / "01-canon/theorems/"
    "THM-3645-berggren-two-cube-local-prime-support-and-mod24-cones.md"
)
EXPECTED_PARENT_SHA256 = "436cc88dd18150e625b92f6e04f6c07538b5f056bdd0ed29902e9e24f12240db"
EXPECTED_SEMANTIC_SHA256 = "2ea01388cebc907078d58844ff4799e8433d509f342195c6ef0240829476ae82"
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def prime_list(limit: int):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for value in range(2, int(limit**0.5) + 1):
        if sieve[value]:
            sieve[value * value:limit + 1:value] = b"\x00" * (
                (limit - value * value) // value + 1
            )
    return tuple(index for index, flag in enumerate(sieve) if flag)


def legendre(value: int, prime: int) -> int:
    residue = pow(value % prime, (prime - 1) // 2, prime)
    require(residue in (1, prime - 1), ("Legendre zero", value, prime))
    return 1 if residue == 1 else -1


def character_pair(prime: int):
    return legendre(6, prime), legendre(-2, prime)


def rho_formula(prime: int) -> int:
    chi24, chi8 = character_pair(prime)
    e24 = int(chi24 == -1)
    e8 = int(chi8 == -1)
    if (e24, e8) == (0, 0):
        return 0
    if (e24, e8) == (0, 1):
        return prime
    if (e24, e8) == (1, 0):
        return 2 * prime - 1
    return 3 * prime - 2


def is_bad(n: int, v: int, prime: int) -> bool:
    chi24, chi8 = character_pair(prime)
    return ((chi24 == -1 and (n % prime == 0 or v % prime == 0))
            or (chi8 == -1 and (2 * n + v) % prime == 0))


def first_prime_in_class(residue: int) -> int:
    for prime in prime_list(10_000):
        if prime >= 5 and prime % 24 == residue:
            return prime
    raise RuntimeError(("missing residue prime", residue))


def main() -> None:
    require(raw_sha256(PARENT) == EXPECTED_PARENT_SHA256, "THM-3645 drift")

    residue_table = []
    coefficient_sum = Fraction(0)
    for residue in range(1, 24):
        if gcd(residue, 24) != 1:
            continue
        representative = first_prime_in_class(residue)
        chi24, chi8 = character_pair(representative)
        coefficient = 2 * int(chi24 == -1) + int(chi8 == -1)
        identity = Fraction(3, 2) - chi24 - Fraction(chi8, 2)
        require(Fraction(coefficient) == identity,
                ("character identity", residue))
        coefficient_sum += coefficient
        residue_table.append((residue, chi24, chi8, coefficient))
    require(len(residue_table) == 8, "unit classes mod 24")
    mean_coefficient = coefficient_sum / len(residue_table)
    require(mean_coefficient == Fraction(3, 2), "sieve dimension mean")

    primes = tuple(prime for prime in prime_list(251) if prime >= 5)
    local_records = []
    for prime in primes:
        chi24, chi8 = character_pair(prime)
        observed = sum(
            ((chi24 == -1 and (n == 0 or v == 0))
             or (chi8 == -1 and (2 * n + v) % prime == 0))
            for n in range(prime) for v in range(prime)
        )
        expected = rho_formula(prime)
        require(observed == expected, ("rho", prime, observed, expected))
        local_records.append((prime, *character_pair(prime), observed))

    crt_records = []
    for factors in ((5, 7), (5, 7, 11)):
        modulus = 1
        expected = 1
        for prime in factors:
            modulus *= prime
            expected *= rho_formula(prime)
        observed = sum(
            all(is_bad(n, v, prime) for prime in factors)
            for n in range(modulus) for v in range(modulus)
        )
        require(observed == expected,
                ("CRT multiplicativity", factors, observed, expected))
        crt_records.append((factors, modulus, observed))

    # The universal remainder majorant used in the theorem follows from
    # rho(p)<=3p and multiplicativity.  Pin it on every squarefree product of
    # at most three tested primes as a hostile sign/order control.
    test_primes = (5, 7, 11, 13, 17, 19)
    remainder_controls = []
    for mask in range(1, 1 << len(test_primes)):
        factors = tuple(test_primes[index] for index in range(len(test_primes))
                        if mask >> index & 1)
        if len(factors) > 3:
            continue
        modulus = 1
        rho = 1
        for prime in factors:
            modulus *= prime
            rho *= rho_formula(prime)
        bound = (3 ** len(factors)) * modulus
        require(rho <= bound, ("rho majorant", factors, rho, bound))
        remainder_controls.append((factors, rho, bound))

    semantic_record = {
        "parent": EXPECTED_PARENT_SHA256,
        "residue_table": tuple(residue_table),
        "dimension": (mean_coefficient.numerator, mean_coefficient.denominator),
        "local_prime_limit": 251,
        "local_records_digest": digest(tuple(local_records)),
        "crt_records": tuple(crt_records),
        "remainder_controls_digest": digest(tuple(remainder_controls)),
        "scope": "local geometry and CRT only; analytic fundamental lemma is cited in theorem",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3656 two-cube Frobenian sieve dimension ==")
    print("parent_sha256_raw=" + EXPECTED_PARENT_SHA256)
    print("residue_table=(class,chi24,chi-8,2e24+e-8)=" + repr(tuple(residue_table)))
    print("mean_bad_line_coefficient=3/2")
    print(f"rho_prime_checks={len(primes)};limit=251;digest={semantic_record['local_records_digest']}")
    print("crt_records=" + repr(tuple(crt_records)))
    print(f"remainder_majorant_controls={len(remainder_controls)};digest={semantic_record['remainder_controls_digest']}")
    print("semantic_sha256=" + semantic)
    print(f"CHECKS={CHECKS}")
    print("status=FINITE-EXACT LOCAL/CRT COMPANION;ANALYTIC UPPER SIEVE IN THEOREM")
    print("scope=no lower bound/asymptotic/constant/primitive tail/Pell/compiler/two-cube existence")


if __name__ == "__main__":
    main()
