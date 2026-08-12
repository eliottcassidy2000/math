#!/usr/bin/env python3
"""Exact hostile audit for THM-3343 (shifted dyadic donors).

The proof is combinatorial.  This companion independently checks:

  C1  the binomial layer formulas and repair bounds for every dyadic pair
      d<M<=512;
  C2  the parity characterization for every 1<=d<M<=64;
  C3  both rotation bijections/injections by exhaustive word enumeration for
      every dyadic pair M<=16;
  C4  the literal lexicographic donor coloring bisects every layer for the
      same exhaustive cubes;
  C5  the donor enumerator is nonzero at x=-1, so no rule with the prescribed
      higher-run labels can decide every donor branch before bit M.

Only Python standard-library integer arithmetic is used.
"""

from itertools import product
from math import comb


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def is_power_two(n):
    return n > 0 and n & (n - 1) == 0


def critical(word):
    a = word[0]
    for i, bit in enumerate(word[1:], 1):
        if bit != a:
            return i
    return None


def rotate_left(word, k=1):
    return word[k:] + word[:k]


def rotate_right(word, k=1):
    return word[-k:] + word[:-k]


def layer_count(M, n, w):
    """Words of length M, weight w, and critical value n."""
    tail = M - n - 1
    return choose(tail, w - 1) + choose(tail, w - n)


def layer_data(d, M, w):
    donor = layer_count(M, d, w)
    heads = sum(
        layer_count(M, n, w)
        for n in range(d + 1, M)
        if (n - d) % 2 == 1
    )
    tails = sum(
        layer_count(M, n, w)
        for n in range(d + 1, M)
        if (n - d) % 2 == 0
    )
    total = donor + heads + tails
    defect = heads - tails
    return donor, heads, tails, total, defect


def check_c1():
    pairs = 0
    layers = 0
    for d in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        for M in [2, 4, 8, 16, 32, 64, 128, 256, 512]:
            if not d < M:
                continue
            pairs += 1
            alt = 0
            for w in range(1, M):
                A, E, O, B, D = layer_data(d, M, w)
                direct = sum(layer_count(M, n, w) for n in range(d, M))
                closed = choose(M - d, w) + choose(M - d, w - d)
                assert B == direct == closed
                assert B % 2 == 0
                assert 0 <= D <= A
                assert (A - D) % 2 == 0
                r = (A - D) // 2
                assert 0 <= r <= A
                assert E + r == B // 2
                alt += (-1) ** w * r
                layers += 1
            assert alt == (-1 if d == 1 else 1)
    print(f"[PASS] C1 dyadic layer formulas/repair bounds: {pairs} pairs, {layers} layers through M=512")


def check_c2():
    tested = 0
    positives = []
    for d in range(1, 64):
        for M in range(d + 1, 65):
            even = all(
                (choose(M - d, w) + choose(M - d, w - d)) % 2 == 0
                for w in range(1, M)
            )
            expected = is_power_two(d) and is_power_two(M)
            assert even == expected, (d, M, even, expected)
            if even:
                positives.append((d, M))
            tested += 1
    print(f"[PASS] C2 parity iff d,M are powers of two: {tested} pairs; {len(positives)} positives")


def shell_words(d, M):
    for word in product((0, 1), repeat=M):
        n = critical(word)
        if n is not None and n >= d:
            yield word


def check_pair_maps_and_coloring(d, M):
    words = list(shell_words(d, M))
    donors = [x for x in words if critical(x) == d]
    heads = [x for x in words if critical(x) > d and (critical(x) - d) % 2 == 1]
    tails = [x for x in words if critical(x) > d and (critical(x) - d) % 2 == 0]

    # Rotation by one: tails biject to the heads ending in their initial bit.
    target = [x for x in heads if x[-1] == x[0]]
    images = [rotate_left(x) for x in tails]
    assert len(images) == len(set(images))
    assert set(images) == set(target)
    for x in tails:
        y = rotate_left(x)
        assert rotate_right(y) == x
        assert critical(y) == critical(x) - 1
        assert y[-1] == y[0]

    # Unmatched heads inject into the donor branch.  The terminal run recovers
    # the rotation distance and therefore the source word.
    unmatched = [x for x in heads if x[-1] != x[0]]
    donor_images = []
    for x in unmatched:
        n = critical(x)
        k = n - d
        y = rotate_left(x, k)
        assert critical(y) == d
        terminal = 0
        for bit in reversed(y):
            if bit == y[0]:
                terminal += 1
            else:
                break
        assert terminal == k
        assert rotate_right(y, terminal) == x
        donor_images.append(y)
    assert len(donor_images) == len(set(donor_images))
    assert set(donor_images).issubset(set(donors))

    # Literal deterministic repair, independently in each Hamming layer.
    labels = {x: 1 for x in heads}
    labels.update({x: 0 for x in tails})
    for w in range(1, M):
        layer = [x for x in words if sum(x) == w]
        donor_layer = sorted(x for x in donors if sum(x) == w)
        A, E, O, B, D = layer_data(d, M, w)
        assert len(layer) == B and len(donor_layer) == A
        r = (A - D) // 2
        labels.update({x: int(i < r) for i, x in enumerate(donor_layer)})
        assert sum(labels[x] for x in layer) * 2 == len(layer)
    assert len(labels) == len(words)
    return len(words), len(unmatched)


def check_c3_c4():
    pairs = 0
    words = 0
    unmatched = 0
    for M in [2, 4, 8, 16]:
        for d in [1, 2, 4, 8]:
            if d >= M:
                continue
            count, miss = check_pair_maps_and_coloring(d, M)
            pairs += 1
            words += count
            unmatched += miss
    print(f"[PASS] C3 exhaustive rotation maps: {pairs} dyadic pairs, {words} shell words, {unmatched} injected defects")
    print("[PASS] C4 lexicographic donor repair bisects every exhaustive Hamming layer")


def check_c5():
    cases = 0
    for d in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        for M in [2, 4, 8, 16, 32, 64, 128, 256, 512]:
            if d >= M:
                continue
            value = sum(
                (-1) ** w * (layer_data(d, M, w)[0] - layer_data(d, M, w)[4]) // 2
                for w in range(1, M)
            )
            assert value == (-1 if d == 1 else 1)
            cases += 1
    print(f"[PASS] C5 donor enumerator R(-1) is nonzero in all {cases} dyadic pairs")


def main():
    check_c1()
    check_c2()
    check_c3_c4()
    check_c5()
    print("THM-3343 SHIFTED-DONOR AUDIT: ALL CHECKS PASS")


if __name__ == "__main__":
    main()
