#!/usr/bin/env python3
"""Numerical Euler-product probe for the two-cube coupled constant.

This evaluates the convergent products in kps-S190 through p<=5,000,000.
It is a numerical diagnostic, not a proof of the correlation asymptotic.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np


LIMIT = 5_000_000
CHECKPOINTS = (1_000, 10_000, 100_000, 1_000_000, 5_000_000)
EXPECTED = (
    (1_000, 0.308466264898, 0.405052002764, 0.958770858860, 0.009238064020),
    (10_000, 0.308457000284, 0.405039867778, 0.958714629054, 0.009236690618),
    (100_000, 0.308456299582, 0.405038951963, 0.958710359602, 0.009236586635),
    (1_000_000, 0.308456242687, 0.405038877337, 0.958710008338, 0.009236578142),
    (5_000_000, 0.308456238399, 0.405038871708, 0.958709981736, 0.009236577500),
)
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


def primes_up_to(limit: int) -> np.ndarray:
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime:limit + 1:prime] = False
    return np.nonzero(sieve)[0]


def character(value: int, prime: int) -> int:
    return 1 if pow(value % prime, (prime - 1) // 2, prime) == 1 else -1


def constants(log_p24: float, log_pm8: float, log_corr: float):
    l24 = math.log(5 + 2 * math.sqrt(6)) / math.sqrt(6)
    lm8 = math.pi / (2 * math.sqrt(2))
    c24 = math.sqrt(l24 * math.exp(log_p24) * (1 - 1 / 2) * (1 - 1 / 3)
                    / math.pi)
    cm8 = math.sqrt(lm8 * math.exp(log_pm8) * (1 - 1 / 2) / math.pi)
    corr = math.exp(log_corr)
    kappa = c24 * c24 * cm8 * corr / 4
    return c24, cm8, corr, kappa


def main() -> None:
    log_p24 = 0.0
    log_pm8 = 0.0
    log_corr = 0.0
    records = []
    checkpoint_index = 0

    primes = primes_up_to(LIMIT)
    require(int(primes[-1]) == 4_999_999, "last prime")
    for raw_prime in primes:
        prime = int(raw_prime)
        if prime < 5:
            continue
        while (checkpoint_index < len(CHECKPOINTS)
               and prime > CHECKPOINTS[checkpoint_index]):
            values = constants(log_p24, log_pm8, log_corr)
            records.append((CHECKPOINTS[checkpoint_index],)
                           + tuple(round(value, 12) for value in values))
            checkpoint_index += 1

        chi24 = character(6, prime)
        chim8 = character(-2, prime)
        if chi24 == -1:
            log_p24 += math.log1p(-1 / (prime * prime))
        if chim8 == -1:
            log_pm8 += math.log1p(-1 / (prime * prime))
        if chi24 == chim8 == 1:
            log_corr += math.log1p(-1 / (prime * prime))
        elif chi24 == chim8 == -1:
            log_corr += math.log(1 - 2 / prime) - 2 * math.log1p(-1 / prime)

    while checkpoint_index < len(CHECKPOINTS):
        values = constants(log_p24, log_pm8, log_corr)
        records.append((CHECKPOINTS[checkpoint_index],)
                       + tuple(round(value, 12) for value in values))
        checkpoint_index += 1

    records_tuple = tuple(records)
    if EXPECTED != "TO_BE_PINNED":
        require(records_tuple == EXPECTED, ("Euler records", records_tuple))
    require(all(0.30 < row[1] < 0.32 and 0.40 < row[2] < 0.41
                and 0.95 < row[3] < 0.97 and 0.009 < row[4] < 0.0095
                for row in records_tuple), "range controls")

    finite_normalized = 633_416 * math.log(49_999) ** 1.5 / 49_999**2
    ratio = finite_normalized / records_tuple[-1][4]
    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== two-cube coupled Euler-product constant probe ==")
    print("records=(prime_cap,C24,Cminus8,correction,kappa):" + repr(records_tuple))
    print(f"finite_N49999_normalized={finite_normalized:.12f};ratio_to_kappa={ratio:.12f}")
    print("semantic_sha256=" + digest(records_tuple))
    print(f"CHECKS={CHECKS}")
    print("status=VERIFIED-NUMERICAL EULER-PRODUCT TRUNCATION")
    print("scope=p<=5000000;no proof of infinite product error or correlation asymptotic")


if __name__ == "__main__":
    main()
