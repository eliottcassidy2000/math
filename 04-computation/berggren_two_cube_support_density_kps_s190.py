#!/usr/bin/env python3
"""Exact finite census for the THM-3645 mod-prime support gate.

The output is a bounded structural count, not a proof of the conjectural
N^2/(log N)^(3/2) asymptotic discussed in the accompanying reflection.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


MAX_N = 50_000
CHECKPOINTS = (401, 997, 1_999, 4_999, 9_999, 19_999, 49_999)
EXPECTED = (
    (401, 8_195, 2_053, 104),
    (997, 50_507, 12_654, 481),
    (1_999, 202_713, 50_741, 1_668),
    (4_999, 1_266_692, 316_812, 8_941),
    (9_999, 5_065_838, 1_266_211, 31_859),
    (19_999, 20_264_681, 5_066_584, 114_230),
    (49_999, 126_652_918, 31_664_297, 633_416),
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


def primes_up_to(limit: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p:limit + 1:p] = (
                b"\x00" * (((limit - p * p) // p) + 1)
            )
    return tuple(index for index in range(2, limit + 1) if sieve[index])


def legendre(value: int, prime: int) -> int:
    if prime == 2:
        return 1
    residue = pow(value % prime, (prime - 1) // 2, prime)
    return 0 if residue == 0 else (1 if residue == 1 else -1)


def semigroup_mask(limit: int, character) -> np.ndarray:
    good = np.ones(limit + 1, dtype=np.bool_)
    good[0] = False
    for prime in primes_up_to(limit):
        if character(prime) != 1:
            good[prime::prime] = False
    return good


def main() -> None:
    good_6 = semigroup_mask(MAX_N, lambda p: legendre(6, p))
    good_minus2 = semigroup_mask(MAX_N, lambda p: legendre(-2, p))
    candidate_total = 0
    three_eligible_total = 0
    survivor_total = 0
    records = []

    for n in range(3, MAX_N + 1, 2):
        v = np.arange(1, n, 2, dtype=np.int64)
        candidate = ((v + n) % 4 == 0) & (np.gcd(v, n) == 1)
        candidate_total += int(np.count_nonzero(candidate))

        three_eligible = candidate & (v % 3 == n % 3) & (n % 3 != 0)
        three_eligible_total += int(np.count_nonzero(three_eligible))

        if good_6[n]:
            quotient = (2 * n + v) // 3
            survive = three_eligible & good_6[v] & good_minus2[quotient]
            survivor_total += int(np.count_nonzero(survive))

        if n in CHECKPOINTS:
            records.append((n, candidate_total, three_eligible_total,
                            survivor_total))

    records_tuple = tuple(records)
    require(records_tuple == EXPECTED, ("checkpoint census", records_tuple))
    require(records_tuple[0][1:] == (8_195, 2_053, 104), "THM-3640 control")
    require(records_tuple[1][1:] == (50_507, 12_654, 481), "THM-3645 control")
    require(all(0 < survivors < eligible < candidates
                for _n, candidates, eligible, survivors in records_tuple),
            "strict count nesting")

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== two-cube mod-prime support density probe ==")
    print("coordinates=U:n-m;V:2m-n;T:2m+n;support=n*V*T")
    print("gate=p|nV:(6/p)=1;p|T,p!=3:(-2/p)=1;p3:m=n!=0")
    print(f"records=(N,candidates,p3_eligible,survivors):{records_tuple}")
    print(f"semantic_sha256={digest(records_tuple)}")
    print(f"CHECKS={CHECKS}")
    print("status=FINITE-EXACT STRUCTURAL CENSUS")
    print("scope=N<=49999;mod-prime gate only;no p-adic/Pell/global/asymptotic claim")


if __name__ == "__main__":
    main()
