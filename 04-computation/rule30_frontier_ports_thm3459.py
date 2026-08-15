#!/usr/bin/env python3
"""Exact audit for THM-3459's Rule 30 frontier ports.

Universe:
* every Boolean local input;
* every Rule-30 initial cone through horizon 8;
* every one-step open block through output width 10 and selected iterates;
* every mask family on a three-point ground set through six owners;
* exact sparse polynomial factorial moments through power 8;
* periodic Rule-30 rings through size 12 for the all-one-preimage law.

Only the finite checks are computational.  The all-horizon trace codec,
preimage count, ternary mask compiler, and stopping statements are proved in
the theorem text.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
from math import ceil, comb, factorial, log


def check(condition: bool, message: str = "exact audit failed") -> None:
    """An optimization-stable audit gate."""
    if not condition:
        raise RuntimeError(message)


def rule30(l: int, c: int, r: int) -> int:
    return l ^ c ^ r ^ (c & r)


def rule60(l: int, c: int, r: int) -> int:
    del r
    return l ^ c


def rule90(l: int, c: int, r: int) -> int:
    del c
    return l ^ r


def rule150(l: int, c: int, r: int) -> int:
    return l ^ c ^ r


def truth_lift(l: int, c: int, r: int) -> int:
    return l + c + r - c * r - 2 * l * c - 2 * l * r + 2 * l * c * r


def trace_from_cone(bits: tuple[int, ...], rule=rule30) -> tuple[int, ...]:
    """The full center trace from a minimal odd initial cone."""
    row = list(bits)
    out = []
    while row:
        out.append(row[len(row) // 2])
        if len(row) == 1:
            break
        row = [rule(row[i], row[i + 1], row[i + 2]) for i in range(len(row) - 2)]
    return tuple(out)


def reconstruct_cone(trace: tuple[int, ...], right: tuple[int, ...], rule=rule30) -> tuple[int, ...]:
    """Invert the trace/right-half codec by solving fresh left bits."""
    horizon = len(trace) - 1
    check(len(right) == horizon)
    x = {0: trace[0]}
    for j, bit in enumerate(right, 1):
        x[j] = bit
    for s in range(1, horizon + 1):
        x[-s] = 0
        trial = tuple(x[j] for j in range(-s, s + 1))
        discrepancy = trace_from_cone(trial, rule)[-1] ^ trace[s]
        x[-s] = discrepancy
        check(trace_from_cone(tuple(x[j] for j in range(-s, s + 1)), rule)[-1] == trace[s])
    return tuple(x[j] for j in range(-horizon, horizon + 1))


def block_step(bits: tuple[int, ...], rule=rule30) -> tuple[int, ...]:
    return tuple(rule(bits[i], bits[i + 1], bits[i + 2]) for i in range(len(bits) - 2))


def block_iter(bits: tuple[int, ...], depth: int, rule=rule30) -> tuple[int, ...]:
    row = bits
    for _ in range(depth):
        row = block_step(row, rule)
    return row


def single_seed_trace(rule, horizon: int) -> tuple[int, ...]:
    row = {0: 1}
    out = []
    for t in range(horizon + 1):
        out.append(row.get(0, 0))
        row = {
            j: rule(row.get(j - 1, 0), row.get(j, 0), row.get(j + 1, 0))
            for j in range(-t - 1, t + 2)
        }
    return tuple(out)


def base3_code(i: int, depth: int) -> tuple[int, ...]:
    digits = []
    for _ in range(depth):
        digits.append(i % 3)
        i //= 3
    return tuple(digits)


def rule30_overlap_defect(masks: tuple[int, ...], ground_mask: int) -> int:
    """Union of all overlaps recovered by the rotated ternary compiler."""
    count = len(masks)
    if count <= 1:
        return 0
    depth = ceil(log(count, 3))
    codes = [base3_code(i, depth) for i in range(count)]
    defect = 0
    for d in range(depth):
        groups = []
        for colour in range(3):
            union = 0
            for i, mask in enumerate(masks):
                if codes[i][d] == colour:
                    union |= mask
            groups.append(union)
        a, b, c = groups
        parity = a ^ b ^ c
        wa = a ^ (b | c)
        wb = b ^ (c | a)
        wc = c ^ (a | b)
        recovered = (wa ^ parity) | (wb ^ parity) | (wc ^ parity)
        check(recovered == ((a & b) | (b & c) | (c & a)))
        defect |= recovered
    return defect & ground_mask


def actual_overlap(masks: tuple[int, ...]) -> int:
    defect = 0
    for i in range(len(masks)):
        for j in range(i + 1, len(masks)):
            defect |= masks[i] & masks[j]
    return defect


Poly = dict[tuple[int, int, int], int]


def poly_mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for a, ca in p.items():
        for b, cb in q.items():
            e = (a[0] + b[0], a[1] + b[1], a[2] + b[2])
            out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c}


def factorial_readout(p: Poly) -> int:
    return sum(c * factorial(a) * factorial(b) * factorial(d) for (a, b, d), c in p.items())


def truth_lift_moments(limit: int) -> list[int]:
    g: Poly = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 1, 1): -1,
        (1, 1, 0): -2,
        (1, 0, 1): -2,
        (1, 1, 1): 2,
    }
    power: Poly = {(0, 0, 0): 1}
    values = []
    for _ in range(limit + 1):
        values.append(factorial_readout(power))
        power = poly_mul(power, g)
    return values


def derangements(limit: int) -> list[int]:
    values = [1]
    if limit >= 1:
        values.append(0)
    for n in range(2, limit + 1):
        values.append((n - 1) * (values[n - 1] + values[n - 2]))
    return values


def truth_lift_compiler(limit: int) -> list[int]:
    d = derangements(2 * limit)
    m = [sum((-1) ** j * comb(k, j) * d[j] ** 2 for j in range(k + 1)) for k in range(2 * limit + 1)]
    return [
        sum(
            comb(power, q)
            * factorial(q)
            * sum(comb(q, s) * (-2) ** s * m[power - q + s] for s in range(q + 1))
            for q in range(power + 1)
        )
        for power in range(limit + 1)
    ]


def anf_lift_moments(limit: int) -> tuple[list[int], list[int]]:
    p: Poly = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1, (0, 1, 1): 1}
    power: Poly = {(0, 0, 0): 1}
    mu = []
    for _ in range(limit + 1):
        mu.append(factorial_readout(power))
        power = poly_mul(power, p)
    u = [mu[m] // factorial(m) for m in range(limit + 1)]
    check(all(mu[m] == factorial(m) * u[m] for m in range(limit + 1)))
    for n in range(4, limit + 1):
        check(u[n] == (n + 2) * u[n - 1] - (2 * n - 1) * u[n - 2] + (n - 3) * u[n - 3] + u[n - 4])
    return mu, u


def ring_step(state: tuple[int, ...]) -> tuple[int, ...]:
    n = len(state)
    return tuple(rule30(state[(i - 1) % n], state[i], state[(i + 1) % n]) for i in range(n))


def audit_truth_table() -> None:
    outputs = []
    for l, c, r in product((0, 1), repeat=3):
        expected = l ^ (c | r)
        check(rule30(l, c, r) == expected == truth_lift(l, c, r))
        det_f2 = ((c + 1) * (r + 1) - (l + 1)) & 1
        check(det_f2 == expected)
        outputs.append(expected)
    check(outputs == [0, 1, 1, 1, 1, 0, 0, 0])


def audit_trace_codec() -> None:
    for horizon in range(9):
        counts: Counter[tuple[int, ...]] = Counter()
        for bits in product((0, 1), repeat=2 * horizon + 1):
            trace = trace_from_cone(bits)
            counts[trace] += 1
            if horizon <= 6:
                rebuilt = reconstruct_cone(trace, bits[horizon + 1 :])
                check(rebuilt == bits)
        check(len(counts) == 1 << (horizon + 1))
        check(set(counts.values()) == {1 << horizon})

    for width in range(1, 11):
        counts = Counter(block_step(bits) for bits in product((0, 1), repeat=width + 2))
        check(len(counts) == 1 << width)
        check(set(counts.values()) == {4})

    for width in range(1, 6):
        for depth in range(1, 4):
            counts = Counter(
                block_iter(bits, depth)
                for bits in product((0, 1), repeat=width + 2 * depth)
            )
            check(len(counts) == 1 << width)
            check(set(counts.values()) == {4**depth})


def audit_single_seed_hostiles() -> None:
    check(single_seed_trace(rule90, 128) == (1,) + (0,) * 128)
    check(single_seed_trace(rule60, 128) == (1,) * 129)
    check(single_seed_trace(rule150, 128) == (1,) * 129)


def audit_mask_compiler() -> None:
    ground_size = 3
    ground_mask = (1 << ground_size) - 1
    families = 0
    for count in range(1, 7):
        for masks in product(range(1 << ground_size), repeat=count):
            families += 1
            check(rule30_overlap_defect(masks, ground_mask) == actual_overlap(masks))
    check(families == sum((1 << ground_size) ** r for r in range(1, 7)))

    partition = (0b001001, 0b010010, 0b100100)
    check(actual_overlap(partition) == 0)
    check(partition[0] | partition[1] | partition[2] == 0b111111)


def audit_factorial_ports() -> tuple[list[int], list[int], list[int]]:
    truth = truth_lift_moments(8)
    check(truth == truth_lift_compiler(8))
    check(truth[:6] == [1, 0, 6, 78, 17520, 4101240])
    anf_mu, anf_u = anf_lift_moments(10)
    check(anf_u == [sum(factorial(d) * comb(n + d + 2, n - d) for d in range(n + 1)) for n in range(11)])
    check(anf_mu[:6] == [1, 4, 26, 270, 4416, 109560])
    check(anf_u[:6] == [1, 4, 13, 45, 184, 913])
    # x^2 and x are the same Boolean function, but factorial readout differs.
    check(factorial(2) - factorial(1) == 1)
    return truth, anf_mu, anf_u


def audit_periodic_gluing() -> dict[int, int]:
    counts = {}
    for n in range(1, 13):
        all_one = (1,) * n
        preimages = [state for state in product((0, 1), repeat=n) if ring_step(state) == all_one]
        expected = 3 if n % 3 == 0 else 0
        check(len(preimages) == expected)
        if preimages:
            base = tuple((1 if i % 3 == 0 else 0) for i in range(n))
            rotations = {base[k:] + base[:k] for k in range(3)}
            check(set(preimages) == rotations)
        counts[n] = len(preimages)
    return counts


def main() -> None:
    audit_truth_table()
    audit_trace_codec()
    audit_single_seed_hostiles()
    audit_mask_compiler()
    truth, anf_mu, anf_u = audit_factorial_ports()
    ring_counts = audit_periodic_gluing()

    print("THM-3459 EXACT AUDIT")
    print("truth_table=00011110")
    print("trace_horizons=0..8 uniform_fibre=2^T codec_roundtrip=0..6")
    print("open_blocks=width1..10 fibre4; iterates=width1..5 depth1..3 fibre4^depth")
    print("single_seed_hostiles: rule90=1,0^infinity rule60=rule150=1^infinity")
    print("mask_families_ground3_owners1..6=", sum(8**r for r in range(1, 7)), sep="")
    print("q6_partition_masks=({0,3},{1,4},{2,5})")
    print("truth_lift_factorial_moments_0..8=", truth, sep="")
    print("anf_lift_factorial_moments_0..10=", anf_mu, sep="")
    print("anf_normalized_u_0..10=", anf_u, sep="")
    print("periodic_all_one_preimages_n1..12=", ring_counts, sep="")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
