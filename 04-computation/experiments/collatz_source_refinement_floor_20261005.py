#!/usr/bin/env python3
"""Exact controls for source floors, signed bills, and refinement energy.

The complete inherited cost-eight certificate is rechecked, not rediscovered.
Finite orbit checks have a hard cap; they are not a universal convergence test.
All checks remain enabled under python -O. Standard library only.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(value, witness=None):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(witness)


def step(n):
    if n < 1 or not n & 1:
        raise ValueError(n)
    x = 3 * n + 1
    a = (x & -x).bit_length() - 1
    return x >> a, a


@lru_cache(None)
def weight(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


@lru_cache(None)
def labels(n):
    """Finite discovery with a resource guard and a strict first visit to 1."""
    if n == 1:
        return 0, 0
    word = []
    x = n
    for _ in range(20000):
        x, a = step(x)
        word.append(a)
        if x == 1:
            return len(word) - 1, sum((a - 1) // 2 for a in word)
    raise RuntimeError(("orbit cap; no convergence conclusion", n))


def W(n):
    return weight(*labels(n))


def unsibling(n):
    j = 0
    while n % 8 == 5:
        n = (n - 1) // 4
        j += 1
    return n, j


def inverse(n, a):
    q, r = divmod((1 << a) * n - 1, 3)
    if r or q < 1 or not q & 1:
        raise ValueError((n, a))
    return q


def rho(n):
    return inverse(n, {1: 6, 2: 5, 4: 4, 5: 1, 7: 2, 8: 3}[n % 9])


def verify_kernel(rows, q):
    bank = {b: (L, K, p, j) for b, L, K, p, j in rows if K <= q}
    check(bank.get(1) == (0, 0, None, None), "root record")
    for b, (L, K, p, j) in bank.items():
        check(unsibling(b) == (b, 0))
        check(labels(b) == (L, K), ("independent forward labels", b))
        if b != 1:
            check(p in bank and unsibling(step(b)[0]) == (p, j))
            check((L, K) == (bank[p][0] + 1, bank[p][1] + j))
        for depth in range(q - K + 1):
            n = 4 ** depth * b + (4 ** depth - 1) // 3
            if n % 3 == 0 or n == 1:
                continue
            child = inverse(n, 1 if n % 3 == 2 else 2)
            check(bank.get(child) == (L + 1, K + depth, b, depth),
                  ("incomplete cost boundary", q, child))
    return bank


def rational_floor(eta, a, b):
    """Quote one aggregate signed bill, using a rational endpoint cutoff."""
    if not F(0) < eta <= F(1, 3):
        raise ValueError(eta)
    m = max(0, a, b)
    if m == 0:
        return eta
    s = F(1, 2)
    while 4 * s * s > eta:
        s /= 2
    return (eta - 2 * s * s) * (s * (1 - s)) ** m


def bill(n, h, common):
    def prefix(x):
        length = cost = 0
        while x != common:
            check(x != 1 and length < 20000, ("bad common endpoint", n, h, common))
            x, a = step(x)
            length += 1
            cost += (a - 1) // 2
        return cost, length
    x, y = prefix(n), prefix(h)
    return x[0] - y[0], x[1] - y[1]


def children(n):
    start = 2 if n % 3 == 1 else 3
    return tuple(inverse(n, a) for a in range(start, start + 5, 2)
                 if inverse(n, a) % 3)


def tree_member(n):
    depth = 0
    while n != 5:
        if n < 5 or not n & 1 or n % 3 == 0:
            return None
        p, _ = step(n)
        if p < 5 or p >= n or n not in children(p):
            return None
        n = p
        depth += 1
    return depth


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    inherited_path = ROOT / "05-knowledge/results/collatz_adaptive_mixture_flow_20261005.json"
    inherited = json.loads(inherited_path.read_text())
    rows = inherited["cost_kernel"]
    bank = verify_kernel(rows, 8)
    # Hostile: a forward-correct subset is not a complete kernel.
    rejected = False
    try:
        verify_kernel([row for row in rows if row[0] != 9], 2)
    except RuntimeError:
        rejected = True
    check(rejected, "must reject a missing zero-cost continuation")
    exceptional = []
    constants = {J: 8 * (F(4, J) + 4 ** J) for J in range(1, 6)}
    for J in range(1, 6):
        q = 2 * J - 2
        low_bank = verify_kernel(rows, q)
        # A bijection: each nonroot base maps to its unit target with same K.
        units = {step(b)[0] for b in low_bank if b != 1}
        by_siblings = set()
        for b, (L, K, _, _) in low_bank.items():
            for j in range(q - K + 1):
                n = 4 ** j * b + (4 ** j - 1) // 3
                if n > 1 and n % 3:
                    by_siblings.add(n)
        check(units == by_siblings and len(units) == len(low_bank) - 1)
        worst, source = F(0), None
        for n in units:
            L, K = labels(n)
            ratio = weight(L, K) ** 3 / weight(L + 1, K + J) ** 2
            check(ratio <= constants[J], ("exception fails", J, n, ratio))
            if ratio > worst:
                worst, source = ratio, n
        exceptional.append(dict(J=J, cost_bound=q, unit_count=len(units),
                                squared_bill_bound=str(constants[J]),
                                largest_formal_bill=str(worst), source=source))
    check(sorted(step(b)[0] for b, _, K, _, _ in rows if b != 1 and K <= 2)
          == [5, 7, 11, 13, 17])

    # The designated section is much stronger than an arbitrary bounded edge.
    section_exceptions = {step(b)[0] for b, _, K, _, _ in rows if b != 1 and K <= 6}
    section_max = max((W(n) ** 3 / W(rho(n)) ** 2, n) for n in section_exceptions)
    check(section_max == (F(441, 160), 85) and rho(85) == 453)
    for n in section_exceptions:
        check(441 * W(rho(n)) ** 2 >= 160 * W(n) ** 3)
    check(F(3025, 1296) < F(441, 160))
    for K in range(3, 80):
        for L in range(80):
            x, k = L + 2, K + 1
            gap = ((x + k) ** 3 * (x + 1) ** 2
                   - x ** 3 * (x + k + 3) ** 2)
            expanded = ((k - 4) * x ** 4 + 2 * (k * k - 4) * x ** 3
                        + (3 * k + 6 * k * k + k ** 3) * x ** 2
                        + (3 * k * k + 2 * k ** 3) * x + k ** 3)
            check(gap == expanded and gap > 0)
        check((K + 3) ** 5 - (K + 6) ** 2 * (K + 1) ** 3
              == 15 * K ** 3 + 125 * K ** 2 + 285 * K + 207)

    max_section_ratio = (F(0), None)
    unit_count = edge_count = 0
    for n in range(5, 1 << 15, 2):
        if n % 3 == 0:
            continue
        unit_count += 1
        L, K = labels(n)
        check(K >= 1 and W(n) <= F(1, 3))
        z = rho(n)
        check(z % 3 == 0 and step(z)[0] == n and 3 * z < 64 * n)
        ratio = W(n) ** 3 / W(z) ** 2
        check(ratio <= F(441, 160), (n, z, ratio))
        if ratio > max_section_ratio[0]:
            max_section_ratio = ratio, n
        for a in range(1, 13):
            if ((1 << a) * n - 1) % 3:
                continue
            j = (a - 1) // 2
            y = inverse(n, a)
            check(labels(y) == (L + 1, K + j))
            if j <= 5:
                J = max(1, j)
                check(constants[J] * W(y) ** 2 >= W(n) ** 3)
            edge_count += 1

    sharp = []
    for t in range(12):
        k = 9 * t + 1
        n = (4 ** (k + 1) - 1) // 3
        z = (2 * n - 1) // 3
        check(n % 9 == 5 and rho(n) == z)
        check(W(n) == F(2, (k + 1) * (k + 2)))
        check(W(z) == F(4, (k + 1) * (k + 2) * (k + 3)))
        check(W(z) / W(n) == F(2, k + 3))
        sharp.append(dict(k=k, W=str(W(n)), leaf=str(W(z))))

    check(bill(13, 3, 5) == (1, 0) and bill(3, 13, 5) == (-1, 0))
    eta = W(3)
    degraded = eta
    for _ in range(5):
        degraded = rational_floor(degraded, 1, 0)
        degraded = rational_floor(degraded, -1, 0)
    check(degraded < eta and rational_floor(eta, 0, 0) == eta)
    for t in range(32):
        n, h = 155 + 2048 * t, 111 + 1458 * t
        common = step(h)[0]
        a, b = bill(n, h, common)
        check((a, b) == (1, 6))
        check(W(n) >= W(h) ** 4 / 8192)
        check(W(n) >= rational_floor(W(h), a, b))

    # Fixed-source selector energy on nested residue partitions, positive and
    # deleted atom controls. Infinite geometric masses are evaluated exactly.
    energies = []
    previous = None
    for q in (2, 4, 8, 16, 32, 64):
        positive = F(1, 4) / (1 - F(1, 2 ** q))
        deleted = F(1, 3 * (2 ** q - 1))
        check(positive >= F(1, 4) and deleted > 0)
        if previous:
            oldq, oldpos, olddel = previous
            check(q % oldq == 0)
            for c0, c1 in ((oldpos, positive), (olddel, deleted)):
                # Direct two-region norm squared, versus telescoping energy.
                direct = c1 * (1 / c1 - 1 / c0) ** 2 + (c0 - c1) / c0 ** 2
                check(direct == 1 / c1 - 1 / c0 and direct > 0)
        energies.append(dict(modulus=q, positive_norm2=str(1 / positive),
                             deleted_norm2=str(1 / deleted)))
        previous = q, positive, deleted

    # An intrinsically recognized binary inverse tree, independent of finding
    # arbitrary ROOT paths. Strict parent decrease makes recognition total.
    level = {5}
    seen = set()
    tree_count = 0
    for depth in range(12):
        check(len(level) == 2 ** depth and not level & seen)
        seen |= level
        next_level = set()
        for n in level:
            check(tree_member(n) == depth)
            L, K = labels(n)
            bits = n.bit_length()
            check(L == depth and K <= 3 * depth + 1 and depth <= 3 * bits)
            floor = weight(3 * bits, 9 * bits + 1)
            check(W(n) >= floor)
            z = rho(n)
            # A rational weakening of the sharp Holder floor, no square root.
            check(W(z) >= floor ** 2 / 12)
            for c in children(n):
                check(c > n and step(c)[0] == n and c % 3)
                check(3 * (c - 1) >= 4 * (n - 1))
                next_level.add(c)
            check(len(children(n)) == 2)
            tree_count += 1
        level = next_level
    check(tree_member(7) is None and W(7) > 0,
          "the binary tree is a proper rooted subfamily")
    check(tree_member(3) is None and tree_member(1) is None)

    # Incoming Pascal-note audit: discount changes the split; zero unit defect
    # is compatible with a positive unit weight.
    D = lambda L, K: weight(L, K) * F(2, (L + 1) * (L + 2))
    check(D(0, 0) == 1 and D(1, 0) + D(0, 1) == F(5, 9))
    for L in range(10):
        for K in range(10):
            check(D(L, K) == F(L + 3, L + 1) * D(L + 1, K) + D(L, K + 1))
    check(W(5) == F(1, 3))  # its unit incoming deficit is zero by the fibre sum

    result = dict(
        status="PROVED scoped plus FINITE-EXACT; universal positivity OPEN",
        checks=CHECKS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        inherited_kernel_sha256=sha256(inherited_path.read_bytes()).hexdigest(),
        kernel_bases=len(bank), exceptional_kernels=exceptional,
        sharp_section=dict(squared_coefficient="160/441", equality_source=85,
                           equality_leaf=453, exceptional_unit_count=len(section_exceptions),
                           large_cost_bill_bound="3025/1296"),
        finite_unit_targets=unit_count, finite_inverse_edges=edge_count,
        largest_observed_section_bill=dict(value=str(max_section_ratio[0]),
                                           source=max_section_ratio[1]),
        sharp_exponent_controls=sharp,
        canceling_bill=dict(aggregate_floor=str(eta),
                           subdivided_floor_less_than_original=True),
        selector_energy_controls=energies,
        binary_tree=dict(levels="0..11", sources=tree_count,
                         rejected_rooted_control=7),
        repaired_discount_split="D00=1; D10+D01=5/9",
        limitations=["No arbitrary-source positive anchor was proved.",
                     "General Holder constants J=1..5 use the verified cost-eight boundary.",
                     "Finite exact checks supplement, not replace, the analytic proofs.",
                     "No fixed-source uniform selector-energy bound was proved."])
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
