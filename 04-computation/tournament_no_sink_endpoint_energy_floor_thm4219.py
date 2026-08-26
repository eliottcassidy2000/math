#!/usr/bin/env python3
"""Exact audit for the provisional THM-4219 endpoint-energy packet.

The script uses only literal tournament adjacency and subset Hamilton-path
DP.  It checks the all-order proof identities on every labelled tournament
through order five, evaluates the independent one-bad-permutation model, and
checks the two closed near-ordinal tower formulae against literal counts.
"""

from __future__ import annotations

from itertools import permutations


def tournament_from_word(n: int, word: int) -> tuple[int, ...]:
    out = [0] * n
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (word >> p) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            p += 1
    return tuple(out)


def tower(n: int, r: int) -> tuple[int, ...]:
    """T(n,r): transitive base and z beating 0,...,r-1,n-3."""
    assert n >= r + 4
    z = n - 1
    special = set(range(r)) | {n - 3}
    out = [0] * n
    for i in range(z):
        for j in range(i + 1, z):
            out[i] |= 1 << j
    for i in range(z):
        if i in special:
            out[z] |= 1 << i
        else:
            out[i] |= 1 << z
    return tuple(out)


def stats(out: tuple[int, ...]) -> tuple[int, tuple[int, ...], int, int, int, int, int]:
    n = len(out)
    size = 1 << n
    full = size - 1
    finish = [[0] * n for _ in range(size)]
    h = [0] * size
    h[0] = 1
    for i in range(n):
        finish[1 << i][i] = 1
    for mask in range(1, size):
        for last in range(n):
            if not ((mask >> last) & 1):
                continue
            rest = mask ^ (1 << last)
            if rest:
                finish[mask][last] = sum(
                    finish[rest][prev]
                    for prev in range(n)
                    if ((rest >> prev) & 1) and ((out[prev] >> last) & 1)
                )
        h[mask] = sum(finish[mask])
    H = h[full]
    v = tuple(
        sum(
            finish[mask][i] * h[full ^ mask]
            for mask in range(1, size)
            if (mask >> i) & 1
        )
        for i in range(n)
    )
    a = sum(v)
    W = a - H
    b = W + 2 * H
    m = sum(x * x for x in v)
    delta = b * b - H * H - 3 * m
    core_gap = b * b - H * H - 5 * m
    return H, v, W, b, m, delta, core_gap


def one_bad_counts(out: tuple[int, ...]) -> tuple[int, ...]:
    n = len(out)
    counts = [0] * n
    for order in permutations(range(n)):
        bad = [
            k
            for k in range(n - 1)
            if not ((out[order[k]] >> order[k + 1]) & 1)
        ]
        if len(bad) == 1:
            counts[order[bad[0]]] += 1
    return tuple(counts)


def sigma(out: tuple[int, ...], H: int, v: tuple[int, ...]) -> tuple[int, ...]:
    n = len(out)
    return tuple(
        H
        + sum(v[j] for j in range(n) if (out[j] >> i) & 1)
        - v[i]
        for i in range(n)
    )


def partition_profile(n: int, r: int) -> tuple[int, tuple[int, ...]]:
    """Evaluate the exact L^up,z,R^up sums in THM-4219 equation (24)."""
    m = n - 1
    full = (1 << m) - 1
    special = sum(1 << i for i in set(range(r)) | {m - 2})

    h = [0] * (1 << m)
    for U in range(1 << m):
        L = U
        while True:
            R = U ^ L
            left_ok = L == 0 or not ((special >> (L.bit_length() - 1)) & 1)
            right_ok = R == 0 or ((special >> ((R & -R).bit_length() - 1)) & 1)
            if left_ok and right_ok:
                h[U] += 1
            if L == 0:
                break
            L = (L - 1) & U

    def F(C: int) -> int:
        value = 1
        for q in range(m):
            if ((C >> q) & 1) and not ((special >> q) & 1):
                value += 1 << ((C & ((1 << q) - 1)).bit_count())
        return value

    v = []
    for i in range(m):
        upper = full ^ ((1 << (i + 1)) - 1)
        first = sum(h[upper | lower] for lower in range(1 << i))
        second = 0
        for q in range(i + 1):
            if not ((special >> q) & 1):
                continue
            if q == i:
                second += F(full ^ (1 << i))
                continue
            middle_width = i - q - 1
            for middle in range(1 << middle_width):
                R = (1 << q) | (1 << i) | (middle << (q + 1))
                second += F(full ^ R)
        v.append(first + second)
    v.append(F(full))
    return h[full], tuple(v)


def expected_t1(n: int) -> tuple[int, tuple[int, ...], int, int, int, int, int]:
    x = 2 ** (n - 4)
    y = 3 ** (n - 4)
    H = 3 * x + 3
    v = [H]
    v.extend(3**i * 2 ** (n - i - 3) + 3 * 2**i for i in range(1, n - 4))
    v.extend((2 * y + 3 * x - 1, 4 * y + 6 * x - 1, 4 * y + 6 * x - 3, 6 * x - 1))
    W = 14 * y + 18 * x - 12
    b = 14 * y + 24 * x - 6
    m_num = 609 * x * x + 570 * x * y - 330 * x + 196 * y * y - 180 * y + 45
    assert m_num % 5 == 0
    m = m_num // 5
    delta_num = 2 * (504 * x * x + 825 * x * y - 270 * x + 196 * y * y - 150 * y)
    assert delta_num % 5 == 0
    delta = delta_num // 5
    gap = 6 * (17 * x * y + 2 * y - 7 * x * x + 4 * x - 3)
    return H, tuple(v), W, b, m, delta, gap


def expected_t2(n: int) -> tuple[int, tuple[int, ...], int, int, int, int, int]:
    x = 2 ** (n - 5)
    y = 3 ** (n - 5)
    H = 9 * x + 1
    v = [9 * x + 2, 12 * x + 5]
    v.extend(
        15 * 3 ** (i - 2) * 2 ** (n - i - 3) + 2 ** (i + 1)
        for i in range(2, n - 4)
    )
    v.extend((10 * y + 4 * x - 3, 20 * y + 8 * x - 3, 20 * y + 4 * x - 9, 12 * x - 3))
    W = 70 * y + 14 * x - 20
    b = 70 * y + 32 * x - 18
    m_num = 871 * x * x + 1800 * x * y - 540 * x + 2940 * y * y - 1620 * y + 347
    assert m_num % 3 == 0
    m = m_num // 3
    delta = 2 * (36 * x * x + 1340 * x * y - 315 * x + 980 * y * y - 450 * y - 12)
    gap_num = 2 * (2220 * x * y + 270 * y - 763 * x * x - 405 * x - 383)
    assert gap_num % 3 == 0
    gap = gap_num // 3
    return H, tuple(v), W, b, m, delta, gap


def audit_labelled() -> None:
    total = 0
    no_sink = 0
    weak_equal = 0
    by_order = []
    for n in range(1, 6):
        rows = 1 << (n * (n - 1) // 2)
        order_no_sink = 0
        for word in range(rows):
            out = tournament_from_word(n, word)
            H, v, _W, _b, _m, delta, _gap = stats(out)
            D = one_bad_counts(out)
            assert all(v[i] == H + D[i] for i in range(n))
            sig = sigma(out, H, v)
            assert all(x >= 0 for x in sig)
            assert delta == 2 * sum(v[i] * sig[i] for i in range(n))
            total += 1
            if all(out[i] != 0 for i in range(n)):
                no_sink += 1
                order_no_sink += 1
                assert all(x >= H for x in v)
                assert delta >= n * (n - 1) * H * H
                if delta == 2 * n * H * H:
                    weak_equal += 1
                    assert n == 3
        by_order.append(order_no_sink)
    assert total == 1099
    assert no_sink == 738
    assert weak_equal == 2  # the two labelled orientations of C3
    print(f"labelled_orders_1_to_5={total} no_sink={no_sink} no_sink_by_order={by_order}")
    print(f"weak_floor_equal_labelled={weak_equal} object_identity=PASS canonical_sigma=PASS")


def audit_towers() -> None:
    for n in range(5, 13):
        got = stats(tower(n, 1))
        want = expected_t1(n)
        assert got == want
        assert partition_profile(n, 1) == got[:2]
        assert got[-1] > 0
        print(
            f"T1 n={n} H={got[0]} W={got[2]} b={got[3]} "
            f"m={got[4]} Delta={got[5]} G={got[6]}"
        )
    for n in range(6, 13):
        got = stats(tower(n, 2))
        want = expected_t2(n)
        assert got == want
        assert partition_profile(n, 2) == got[:2]
        assert got[-1] > 0
        print(
            f"T2 n={n} H={got[0]} W={got[2]} b={got[3]} "
            f"m={got[4]} Delta={got[5]} G={got[6]}"
        )
    for n in range(5, 65):
        assert expected_t1(n)[-1] > 0
    for n in range(6, 65):
        assert expected_t2(n)[-1] > 0
    print("tower_closed_form_literal_DP_and_partition_sums=T1[5,12],T2[6,12] positivity_integer_control_through_64=PASS")


def hostile_controls() -> None:
    transitive3 = tournament_from_word(3, 0b111)
    transitive4 = tournament_from_word(4, 0b111111)
    s3 = stats(transitive3)
    s4 = stats(transitive4)
    assert (s3[0], s3[2], s3[4], s3[6]) == (1, 6, 21, -42)
    assert (s4[0], s4[2], s4[4], s4[6]) == (1, 14, 85, -170)
    c3 = tournament_from_word(3, 0b101)
    sc = stats(c3)
    assert (sc[0], sc[1], sc[2], sc[4], sc[6]) == (3, (3, 3, 3), 6, 27, 0)
    print("controls=C3_equal;transitive3_G=-42;transitive4_G=-170")


def main() -> None:
    audit_labelled()
    hostile_controls()
    audit_towers()
    print("THM4219_PRIMARY_AUDIT=PASS")


if __name__ == "__main__":
    main()
