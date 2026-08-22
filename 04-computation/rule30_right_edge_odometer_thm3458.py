#!/usr/bin/env python3
"""Exact audit for THM-3458's Rule-30 edge/odometer compiler.

Finite computations audit two independent row implementations, every state of
the first twelve triangular quotients, seed periods through width thirty,
fixed-offset rational words, periodic-ring prefixes, de Bruijn preimage counts,
and the characteristic/lift boundaries.  Universal statements are proved in
the theorem text; no fitted period prefix is extrapolated.
"""

from __future__ import annotations

from itertools import product


def rule30(l: int, c: int, r: int) -> int:
    return l ^ c ^ r ^ (c & r)


def phi(x: int, width: int | None = None) -> int:
    y = x ^ ((x << 1) | (x << 2))
    if width is not None:
        y &= (1 << width) - 1
    return y


def phi_inverse(y: int, width: int) -> int:
    x = 0
    for k in range(width):
        lower_one = ((x >> (k - 1)) & 1) if k >= 1 else 0
        lower_two = ((x >> (k - 2)) & 1) if k >= 2 else 0
        bit = ((y >> k) & 1) ^ (lower_one | lower_two)
        x |= bit << k
    return x


def direct_rows(horizon: int) -> list[int]:
    row = [1]
    packed = []
    for _ in range(horizon + 1):
        packed.append(sum(bit << k for k, bit in enumerate(reversed(row))))
        padded = [0, 0] + row + [0, 0]
        row = [rule30(padded[i], padded[i + 1], padded[i + 2]) for i in range(len(padded) - 2)]
    return packed


def seed_period(width: int) -> int:
    x = 1
    period = 0
    while True:
        x = phi(x, width)
        period += 1
        if x == 1:
            return period


def is_power_two(n: int) -> bool:
    return n > 0 and n & (n - 1) == 0


def bit_word(offset: int, period: int) -> str:
    x = 1
    bits = []
    width = offset + 1
    for _ in range(period):
        bits.append(str((x >> offset) & 1))
        x = phi(x, width)
    return "".join(bits)


def ring_step(state: tuple[int, ...]) -> tuple[int, ...]:
    n = len(state)
    return tuple(rule30(state[(i - 1) % n], state[i], state[(i + 1) % n]) for i in range(n))


def ring_orbit(n: int) -> tuple[int, int, list[tuple[int, ...]]]:
    state = (1,) + (0,) * (n - 1)
    seen: dict[tuple[int, ...], int] = {}
    orbit = []
    while state not in seen:
        seen[state] = len(orbit)
        orbit.append(state)
        state = ring_step(state)
    mu = seen[state]
    return mu, len(orbit) - mu, orbit


def mat_mul(a: tuple[tuple[int, ...], ...], b: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    n = len(a)
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)) for i in range(n))


def mat_identity(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(1 if i == j else 0 for j in range(n)) for i in range(n))


def mat_trace(a: tuple[tuple[int, ...], ...]) -> int:
    return sum(a[i][i] for i in range(len(a)))


def debruijn_matrices() -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    mats = []
    for output in (0, 1):
        m = [[0] * 4 for _ in range(4)]
        for a, b, c in product((0, 1), repeat=3):
            if rule30(a, b, c) == output:
                m[2 * a + b][2 * b + c] += 1
        mats.append(tuple(tuple(row) for row in m))
    return mats[0], mats[1]


def cyclic_preimages(word: tuple[int, ...]) -> int:
    a0, a1 = debruijn_matrices()
    product_matrix = mat_identity(4)
    for bit in word:
        product_matrix = mat_mul(product_matrix, (a0, a1)[bit])
    return mat_trace(product_matrix)


def brute_cyclic_preimages(word: tuple[int, ...]) -> int:
    n = len(word)
    return sum(1 for state in product((0, 1), repeat=n) if ring_step(state) == word)


def determinant_bareiss(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((r for r in range(k + 1, n) if a[r][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot_value - a[i][k] * a[k][j]) // previous
        previous = pivot_value
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[-1][-1]


def cyclic_jacobian(n: int, value: int) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for j in range(n):
        matrix[j][(j - 1) % n] += 1
        matrix[j][j] += 1 + value
        matrix[j][(j + 1) % n] += 1 + value
    return matrix


def audit_rows() -> list[int]:
    direct = direct_rows(512)
    x = 1
    packed = []
    for _ in range(513):
        packed.append(x)
        x = phi(x)
    assert packed == direct
    assert packed[:9] == [1, 7, 25, 111, 401, 1783, 6409, 28479, 102849]
    assert all(value.bit_length() == 2 * t + 1 for t, value in enumerate(packed))
    return packed


def audit_quotients() -> tuple[list[int], str]:
    for width in range(1, 13):
        for x in range(1 << width):
            assert phi_inverse(phi(x, width), width) == x
            assert phi(phi_inverse(x, width), width) == x

    periods = [seed_period(w) for w in range(1, 31)]
    for w, period in enumerate(periods, 1):
        assert is_power_two(period)
        assert period <= 1 << (w - 1)
        assert 2 * period >= w
        if w < len(periods):
            assert periods[w] in (period, 2 * period)

    eps = []
    for w, period in enumerate(periods, 1):
        lifted = 1
        for _ in range(period):
            lifted = phi(lifted, w + 1)
        epsilon = (lifted >> w) & 1
        eps.append(str(epsilon))
        if w < len(periods):
            assert periods[w] == period * (1 << epsilon)
    return periods, "".join(eps)


def audit_fixed_offsets(periods: list[int]) -> dict[int, tuple[int, str, int]]:
    data = {}
    for offset in range(13):
        period = periods[offset]
        word = bit_word(offset, period)
        data[offset] = (period, word, word.count("1"))
    assert data[5] == (8, "00010101", 3)
    return data


def audit_moving_observer(rows: list[int]) -> tuple[int, int, int, int]:
    c7 = (rows[7] >> 7) & 1
    c15 = (rows[15] >> 15) & 1
    assert rows[7] & 63 == rows[15] & 63 == 63
    assert (c7, c15) == (0, 1)
    return rows[7] & 63, rows[15] & 63, c7, c15


def audit_rings(rows: list[int]) -> dict[int, tuple[int, int, int]]:
    table = {}
    for n in range(2, 13):
        mu, lam, orbit = ring_orbit(n)
        for t in range(n):
            index = t if t < len(orbit) else mu + (t - mu) % lam
            assert orbit[index][0] == ((rows[t] >> t) & 1)
        cycle_ones = sum(orbit[mu + r][0] for r in range(lam))
        table[n] = (mu, lam, cycle_ones)
    assert table[5] == (1, 5, 3)
    return table


def audit_debruijn() -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    a0, a1 = debruijn_matrices()
    for n in range(1, 9):
        for word in product((0, 1), repeat=n):
            assert cyclic_preimages(word) == brute_cyclic_preimages(word)
        assert cyclic_preimages((0,) * n) == 2
        assert cyclic_preimages((1,) * n) == (3 if n % 3 == 0 else 0)
    return a0, a1


def audit_jacobian_boundary() -> dict[int, tuple[int, int]]:
    dets = {}
    for n in range(3, 15):
        det_zero = determinant_bareiss(cyclic_jacobian(n, 0))
        det_minus_one = determinant_bareiss(cyclic_jacobian(n, -1))
        assert det_zero == 0 if n % 3 == 0 else abs(det_zero) == 3
        assert abs(det_minus_one) == 1
        assert det_zero != det_minus_one
        dets[n] = (det_zero, det_minus_one)

    # Same Boolean truth table, two incompatible characteristic-zero lifts.
    for x, a, b in product((0, 1), repeat=3):
        anf = x + a + b + a * b
        or_real = a + b - a * b
        faithful = x + or_real - 2 * x * or_real
        assert (anf & 1) == faithful == (x ^ (a | b))
    assert 0 + 1 + 1 + 1 * 1 == 3
    return dets


def main() -> None:
    rows = audit_rows()
    periods, eps = audit_quotients()
    offsets = audit_fixed_offsets(periods)
    hostile = audit_moving_observer(rows)
    rings = audit_rings(rows)
    a0, a1 = audit_debruijn()
    dets = audit_jacobian_boundary()

    print("THM-3458 EXACT AUDIT")
    print("packed_rows_0..8=", rows[:9], sep="")
    print("direct_vs_packed_rows=0..512")
    print("triangular_inverse_all_states_width=1..12")
    print("seed_periods_w1..30=", periods, sep="")
    print("period_lift_bits_eps1..30=", eps, sep="")
    print("fixed_offsets_k0..12=(period,word,ones)=", offsets, sep="")
    print("moving_observer_hostile=(low6_t7,low6_t15,c7,c15)=", hostile, sep="")
    print("single_seed_ring_(mu,lambda,cycle_ones)_n2..12=", rings, sep="")
    print("debruijn_A0=", a0, sep="")
    print("debruijn_A1=", a1, sep="")
    print("cyclic_jacobian_(det_at_0,det_at_-1)_n3..14=", dets, sep="")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
