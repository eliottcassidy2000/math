#!/usr/bin/env python3
"""Exact replay of the historical mac-mini-S105 8,260-row band bank.

This audits only the bank that S105 actually generated: interval cores
``{1,...,k}`` for ``k = 11, 10, 9``, its stated outlier pool and ordering,
and its cap of 4,000 accepted rows per k.  All witness decisions use integer
arithmetic.  In particular, this does not promote the capped bank to a
complete LRC(14) residual classification.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations


N = 14
Q25 = 25
Q_AUDIT = 38
W0 = 475
EXPECTED_BANK_HASH = "711bddb07352201565df55f1a916631a61a23b87a6d026269593ab45b8dc632e"
EXPECTED_CERT_HASH = "d6cae4b7ba7426dd189e93a7fa3b04146b191a4ad68d00b5b7f7c41a6d75bf0b"
EXPECTED_FAIL_HASH = "a9606cae75abe2926079126b684b104c6f6dc5f0c1f163460453b4a7c0f3bbcf"


def covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % d == 0 for v in speeds) for d in range(2, N + 1))


def is_witness(speeds: tuple[int, ...], q: int, a: int) -> bool:
    """Exact test for min_v ||a v / q|| >= 1/14."""

    return all(
        N * min((a * v) % q, q - ((a * v) % q)) >= q
        for v in speeds
    )


def first_witness(speeds: tuple[int, ...]) -> tuple[int, int] | None:
    for q in range(2, Q_AUDIT + 1):
        for a in range(1, q):
            if is_witness(speeds, q, a):
                return q, a
    return None


def candidate_pool() -> tuple[int, ...]:
    return tuple(sorted({
        d * j
        for d in range(2, N + 1)
        for j in range(1, W0 // d + 1)
        if 15 <= d * j <= W0
    }))


def replay_bank(pool: tuple[int, ...]):
    for k in (11, 10, 9):
        accepted = 0
        core = tuple(range(1, k + 1))
        for outliers in combinations(pool, 13 - k):
            if max(outliers) <= 220:
                continue
            speeds = core + outliers
            if not covering(speeds):
                continue
            witness = first_witness(speeds)
            assert witness is not None, ("no witness through q=38", speeds)
            yield k, speeds, witness
            accepted += 1
            if accepted >= 4000:
                break


def hexdigest(lines: list[str]) -> str:
    return sha256("".join(lines).encode("ascii")).hexdigest()


def main() -> None:
    pool = candidate_pool()
    rows = list(replay_bank(pool))

    by_k: dict[int, Counter[int]] = {k: Counter() for k in (11, 10, 9)}
    count_by_k = Counter()
    q_hist = Counter()
    failures = []
    for k, speeds, (q, _a) in rows:
        count_by_k[k] += 1
        by_k[k][q] += 1
        q_hist[q] += 1
        if q > Q25:
            failures.append(speeds)

    bank_lines = [
        f"{k};{','.join(map(str, speeds))}\n"
        for k, speeds, _witness in rows
    ]
    cert_lines = [
        f"{k};{','.join(map(str, speeds))};{q};{a}\n"
        for k, speeds, (q, a) in rows
    ]
    fail_lines = [f"{','.join(map(str, speeds))}\n" for speeds in failures]
    bank_hash = hexdigest(bank_lines)
    cert_hash = hexdigest(cert_lines)
    fail_hash = hexdigest(fail_lines)

    expected_by_k = {
        11: Counter({25: 207, 27: 43, 28: 1, 37: 8, 38: 1}),
        10: Counter({23: 3304, 24: 219, 25: 439, 26: 11, 27: 26, 34: 1}),
        9: Counter({21: 1588, 23: 2296, 24: 56, 25: 60}),
    }
    expected_hist = Counter({
        21: 1588, 23: 5600, 24: 275, 25: 706, 26: 11,
        27: 69, 28: 1, 34: 1, 37: 8, 38: 1,
    })
    max_rows = [(speeds, witness) for _k, speeds, witness in rows if witness[0] == 38]

    assert len(pool) == 371
    assert count_by_k == Counter({10: 4000, 9: 4000, 11: 260})
    assert by_k == expected_by_k
    assert q_hist == expected_hist
    assert len(failures) == 91
    assert max_rows == [
        ((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 338, 420), (38, 3))
    ]
    assert bank_hash == EXPECTED_BANK_HASH
    assert cert_hash == EXPECTED_CERT_HASH
    assert fail_hash == EXPECTED_FAIL_HASH

    print("S105 q<=25 bank audit (exact integer replay)")
    print(f"candidate pool: {len(pool)}")
    print(f"accepted rows: {len(rows)} (k=11: {count_by_k[11]}, k=10: {count_by_k[10]}, k=9: {count_by_k[9]})")
    for k in (11, 10, 9):
        dist = " ".join(f"q{q}={count}" for q, count in sorted(by_k[k].items()))
        above = sum(count for q, count in by_k[k].items() if q > Q25)
        print(f"k={k}: {dist}; q>25={above}")
    print("overall qmin: " + " ".join(f"q{q}={count}" for q, count in sorted(q_hist.items())))
    print(f"q>25 failures: {len(failures)}/8260")
    print("maximum qmin: 38")
    print("unique qmin=38 row: (1,2,3,4,5,6,7,8,9,10,11,338,420); witness=3/38")
    print(f"bank sha256: {bank_hash}")
    print(f"full qmin certificate sha256: {cert_hash}")
    print(f"q>25 failure-list sha256: {fail_hash}")
    print("scope: historical capped/interval-core S105 bank only; not a uniform residual theorem")


if __name__ == "__main__":
    main()
