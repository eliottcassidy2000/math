#!/usr/bin/env python3
"""Exact hostile audit for THM-3344.

C1 checks the odd-prefix repair formula through d=256.
C2 exhaustively constructs and audits every word through M=16.
C3 checks the signed polynomial cancellation through d=256.
C4 exhausts every orientation-sign word for d=2,4,8 and confirms that the
   donor enumerator has (1+x)-adic valuation at most one.
"""

from itertools import product
from math import comb


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def critical(word):
    for i, bit in enumerate(word[1:], 1):
        if bit != word[0]:
            return i
    return None


def layer_count(length, n, w):
    tail = length - n - 1
    return choose(tail, w - 1) + choose(tail, w - n)


def prefix_data(d, w):
    """Data at odd prefix H=2d-1."""
    H = 2 * d - 1
    donor = layer_count(H, d, w)
    heads = sum(
        layer_count(H, n, w)
        for n in range(d + 1, H)
        if (n - d) % 2 == 1
    )
    tails = sum(
        layer_count(H, n, w)
        for n in range(d + 1, H)
        if (n - d) % 2 == 0
    )
    defect = heads - tails
    total = donor + heads + tails
    tau = -1 if w % 2 else 1
    repair = (donor - defect + tau) // 2
    return donor, heads, tails, defect, total, tau, repair


def check_c1():
    layers = 0
    for d in [2, 4, 8, 16, 32, 64, 128, 256]:
        H = 2 * d - 1
        for w in range(1, H):
            A, E, O, D, B, tau, q = prefix_data(d, w)
            assert B == choose(H - d, w) + choose(H - d, w - d)
            assert B % 2 == 1
            assert 0 <= D <= A
            assert (A - D) % 2 == 1
            assert 0 <= q <= A
            assert 2 * q - A + D == tau
            layers += 1
    print(f"[PASS] C1 odd-prefix repair: {layers} layers through d=256")


def polynomial_add(a, b):
    n = max(len(a), len(b))
    return [
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(n)
    ]


def polynomial_mul_1x(a):
    return polynomial_add(a, [0] + a)


def check_c3():
    cases = 0
    for d in [2, 4, 8, 16, 32, 64, 128, 256]:
        M, H = 2 * d, 2 * d - 1
        signed_prefix = [0] * (M + 1)
        for w in range(1, H):
            A, E, O, D, B, tau, q = prefix_data(d, w)
            signed_prefix[w] = 2 * q - A + D
            assert signed_prefix[w] == (-1 if w % 2 else 1)
        expanded = polynomial_mul_1x(signed_prefix)
        expanded += [0] * (M + 1 - len(expanded))
        expanded[1] += 1       # 0^(M-1)1 is heads
        expanded[H] -= 1       # 1^(M-1)0 is tails
        assert all(c == 0 for c in expanded)
        cases += 1
    print(f"[PASS] C3 exact signed cancellation in {cases} dyadic annuli")


def check_words(d):
    M, H = 2 * d, 2 * d - 1
    prefixes = [
        y for y in product((0, 1), repeat=H)
        if critical(y) is not None and critical(y) >= d
    ]
    donor_prefixes = [y for y in prefixes if critical(y) == d]
    labels = {}

    for y in prefixes:
        n = critical(y)
        if n > d:
            labels[y] = int((n - d) % 2 == 1)

    for w in range(1, H):
        candidates = sorted(y for y in donor_prefixes if sum(y) == w)
        q = prefix_data(d, w)[-1]
        labels.update({y: int(i < q) for i, y in enumerate(candidates)})

    assert len(labels) == len(prefixes)
    full_shell = [
        x for x in product((0, 1), repeat=M)
        if critical(x) is not None and critical(x) >= d
    ]
    full_labels = {}
    for x in full_shell:
        n = critical(x)
        if n == H:
            full_labels[x] = int(x[0] == 0)
        else:
            full_labels[x] = labels[x[:H]]

    for w in range(1, M):
        layer = [x for x in full_shell if sum(x) == w]
        assert 2 * sum(full_labels[x] for x in layer) == len(layer)

    for x in full_shell:
        n = critical(x)
        if n == d:
            assert x[:H] in labels
        elif n < H:
            assert full_labels[x] == int((n - d) % 2 == 1)
        else:
            assert n == H
    return len(prefixes), len(full_shell)


def check_c2():
    prefixes = words = 0
    for d in [2, 4, 8]:
        a, b = check_words(d)
        prefixes += a
        words += b
    print(f"[PASS] C2 exhaustive literal extractors: {prefixes} prefixes, {words} full words through M=16")


def orient_poly(M, n, orientation):
    out = [0] * (M + 1)
    shift = 1 if orientation == 0 else n
    for j in range(M - n):
        out[shift + j] += choose(M - n - 1, j)
    return out


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def divide_1x(a):
    a = trim(a[:])
    q = [0] * (len(a) - 1)
    q[-1] = a[-1]
    for k in range(len(a) - 2, 0, -1):
        q[k - 1] = a[k] - q[k]
    assert a[0] == q[0]
    return trim(q)


def valuation_1x(a):
    a = trim(a[:])
    value = 0
    while any(a) and sum(c * (-1) ** i for i, c in enumerate(a)) == 0:
        a = divide_1x(a)
        value += 1
    return value


def check_c4():
    configurations = 0
    for d in [2, 4, 8]:
        M = 2 * d
        donor = polynomial_add(orient_poly(M, d, 0), orient_poly(M, d, 1))
        rows = [
            orient_poly(M, n, orientation)
            for n in range(d + 1, M)
            for orientation in (0, 1)
        ]
        max_value = 0
        attained_one = False
        for signs in product((-1, 1), repeat=len(rows)):
            signed = [sum(s * row[i] for s, row in zip(signs, rows)) for i in range(M + 1)]
            difference = [donor[i] - signed[i] for i in range(M + 1)]
            if any(c % 2 for c in difference):
                continue
            repair = [c // 2 for c in difference]
            if any(not 0 <= repair[i] <= donor[i] for i in range(M + 1)):
                continue
            value = valuation_1x(repair)
            assert value <= 1
            max_value = max(max_value, value)
            attained_one |= value == 1
            configurations += 1
        assert max_value == 1 and attained_one
    print(f"[PASS] C4 orientation-sign hostile census: {configurations} feasible configurations, max v_(1+x)=1")


def main():
    check_c1()
    check_c2()
    check_c3()
    check_c4()
    print("THM-3344 ORIENTATION-SPLIT AUDIT: ALL CHECKS PASS")


if __name__ == "__main__":
    main()
