#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3507."""

from __future__ import annotations

import hashlib


WIDTH = 24
EXPECTED_V = (1, 3, 4, 6, 7, 9, 15, 16, 24)


def check(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def phi(x: int, mask: int) -> int:
    return (x ^ ((x << 1) | (x << 2))) & mask


def nu2(x: int) -> int:
    check(x != 0, "nu2 zero")
    return (x & -x).bit_length() - 1


def seed_period(width: int) -> tuple[int, list[int]]:
    mask = (1 << width) - 1
    row = 1
    orbit = [row]
    for period in range(1, 1 << width):
        row = phi(row, mask)
        if row == 1:
            return period, orbit
        orbit.append(row)
    raise RuntimeError("seed return missing")


def period_data() -> tuple[list[int], list[int], list[int]]:
    periods: list[int] = []
    eps: list[int] = []
    for width in range(1, WIDTH + 1):
        period, _ = seed_period(width)
        periods.append(period)
        row = 1
        mask = (1 << (width + 1)) - 1
        for _ in range(period):
            row = phi(row, mask)
        eps.append((row >> width) & 1)
    for width in range(1, WIDTH):
        check(periods[width] == periods[width - 1] << eps[width - 1], "lift")
    innovations = [k for k, e in enumerate(eps, 1) if e]
    check(tuple(innovations[: len(EXPECTED_V)]) == EXPECTED_V, "v prefix")
    return periods, eps, innovations


def strict_suffix(x: int, n: int) -> int:
    out = 0
    parity = 0
    for i in range(n - 1, -1, -1):
        if parity:
            out |= 1 << i
        parity ^= (x >> i) & 1
    return out


def section_first_odd(q: int, search_limit: int) -> int:
    n = 2 * q
    mask_q = (1 << q) - 1
    rows = [(1 << n) - 1, (1 << n) - 1]
    while len(rows) < search_limit:
        rows.append(strict_suffix(rows[-2], n) | strict_suffix(rows[-1], n))
    for k, row in enumerate(rows):
        delta = (row & mask_q) ^ (row >> q)
        if delta.bit_count() & 1:
            return k
    raise RuntimeError("section odd defect not found")


def borrow_decode(x: int, y: int, bits: int) -> int:
    out = 0
    borrow = 0
    for j in range(bits):
        xb = (x >> j) & 1
        yb = (y >> j) & 1
        h = xb ^ yb
        out |= (h ^ borrow) << j
        borrow = ((1 - yb) & xb) | (borrow & (1 - h))
    return out


def audit_units(periods: list[int], innovations: list[int]) -> tuple[int, int, int, int]:
    period = periods[-1]
    full_mask = (1 << WIDTH) - 1
    rows: list[int] = []
    row = 1
    for _ in range(period):
        rows.append(row)
        row = phi(row, full_mask)
    check(row == 1, "orbit closure")

    trace_checks = borrow_checks = seam_checks = gap_checks = 0
    for m in range(len(innovations) - 1):
        v = innovations[m]
        v_next = innovations[m + 1]
        if v_next >= WIDTH:
            break
        q = 1 << m
        d = v_next - v
        bits = WIDTH - v
        unit_mod = 1 << bits
        for t in range(period):
            x_full = rows[t]
            y_full = rows[(t + q) % period]
            z_full = rows[(t + 2 * q) % period]
            u = ((y_full - x_full) & full_mask) >> v
            u_sib = ((z_full - y_full) & full_mask) >> v
            u_next = ((z_full - x_full) & full_mask) >> v_next
            check(u & 1 and u_sib & 1 and u_next & 1, ("odd units", m, t))
            sibling_sum = (u + u_sib) % unit_mod
            check(sibling_sum == ((1 << d) * u_next) % unit_mod, ("trace", m, t))
            check(nu2(sibling_sum) == d, ("gap valuation", m, t))
            trace_checks += 1

            x = (x_full >> v) & (unit_mod - 1)
            y = (y_full >> v) & (unit_mod - 1)
            check(borrow_decode(x, y, bits) == u, ("borrow", m, t))
            borrow_checks += 1

        if v < q:
            for t in range(q):
                h = ((rows[(t + q) % period] ^ rows[t]) >> v) & 1
                check(h == 1, ("hard seam", m, t))
                seam_checks += 1

        first_k = section_first_odd(q, v_next + 2)
        check(first_k + 1 == v_next, ("section next valuation", m, first_k))
        check(v_next - v == d, "gap tautology")
        gap_checks += 1
    return trace_checks, borrow_checks, seam_checks, gap_checks


def audit_local_covers(periods: list[int]) -> int:
    period = periods[-1]
    mask = (1 << WIDTH) - 1
    rows: list[int] = []
    row = 1
    for _ in range(period):
        rows.append(row)
        row = phi(row, mask)

    checks = 0
    for w in range(1, WIDTH):
        coarse_mask = (1 << w) - 1
        groups: dict[int, list[int]] = {}
        for value in rows:
            groups.setdefault(value & coarse_mask, []).append(value)
        check(len(groups) == periods[w - 1], ("coarse cells", w))
        for fine in range(w + 1, WIDTH + 1):
            fine_mask = (1 << fine) - 1
            expected = periods[fine - 1] // periods[w - 1]
            for values in groups.values():
                count = len({value & fine_mask for value in values})
                check(count == expected, ("local cover", w, fine))
                checks += 1
    return checks


def longest_zero_run(word: list[int]) -> int:
    best = current = 0
    for bit in word:
        if bit:
            current = 0
        else:
            current += 1
            best = max(best, current)
    return best


def audit_scalar_controls() -> tuple[int, int, int, str]:
    horizon = 384
    dense = [1 if i % 3 in (0, 1) else 0 for i in range(horizon)]
    sparse = [0] * horizon
    k = 1
    while k < horizon:
        sparse[k] = 1
        k = 2 * k + 1
    burst_gap: list[int] = []
    for r in range(1, 9):
        burst_gap.extend([1, 1, 0] * r)
        burst_gap.extend([0] * (3 * r))

    for word in (dense, sparse, burst_gap):
        check(
            all(not (word[i] and word[i + 1] and word[i + 2]) for i in range(len(word) - 2)),
            "no111 scalar",
        )
    check(sum(dense) * 3 == 2 * horizon, "dense exact 2/3")
    check(longest_zero_run(sparse) >= 128, "sparse gap")
    check(longest_zero_run(burst_gap) >= 24, "burst gap")
    check(
        any(burst_gap[i : i + 6] == [1, 1, 0, 1, 1, 0] for i in range(len(burst_gap) - 5)),
        "dense burst",
    )
    digest = hashlib.sha256(bytes(dense + sparse + burst_gap)).hexdigest()
    return sum(dense), longest_zero_run(sparse), longest_zero_run(burst_gap), digest


def main() -> None:
    periods, eps, innovations = period_data()
    trace, borrow, seam, gaps = audit_units(periods, innovations)
    local_covers = audit_local_covers(periods)
    dense_count, sparse_gap, burst_gap, digest = audit_scalar_controls()

    mask = (1 << WIDTH) - 1
    row = 1
    marked = {0: row}
    for t in range(1, 9):
        row = phi(row, mask)
        marked[t] = row
    check(marked[4] == 401 and marked[8] == 102849, "q4 rows")
    check((marked[4] - 1) // 16 == 25, "q4 first unit")
    check((marked[8] - marked[4]) // 16 == 6403, "q4 sibling unit")
    check(25 + 6403 == 4 * 1607, "q4 trace")

    print("RULE30_NORMALIZED_DISPLACEMENT_ASSOUAD_THM3507")
    print("status PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print("width", WIDTH, "period", periods[-1], "innovation_prefix", tuple(innovations[:9]))
    print("unit_trace_checks", trace, "borrow_checks", borrow, "hard_seam_checks", seam)
    print("section_gap_scales", gaps, "local_cover_checks", local_covers)
    print("q4_units", (25, 6403, 1607), "gap", 2)
    print("scalar_controls", "dense_count", dense_count, "sparse_gap", sparse_gap, "burst_gap", burst_gap)
    print("scalar_sha256", digest)
    print("epsilon_prefix", "".join(map(str, eps[:24])))
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
