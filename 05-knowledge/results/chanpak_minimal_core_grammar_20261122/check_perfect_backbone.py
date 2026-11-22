"""Exact grammar for every orientation of the five-vertex lower-rail backbone."""

from collections import Counter, defaultdict
from itertools import combinations, product
from math import comb, prod


ROWS = (2, 8, 3, 37, 15, 151, 63, 95)
PORTS = (6, 4, 1)
BACKBONE = (0, 2, 3, 5, 7)
EDGES = tuple(combinations(BACKBONE, 2))
FULL = sum(1 << j for j in BACKBONE)
LOWER = {16: 32, 4: -4, 1: -1}


def source(mask):
    rows = list(ROWS)
    for i, (a, b) in enumerate(EDGES):
        if mask >> i & 1:
            rows[a] ^= 1 << b
            rows[b] ^= 1 << a
    return rows


def compressed_paths(rows):
    """Each backbone vertex is used once. Port runs lie in transitive intervals."""
    layer = {(0, (0, 0, 0), -1): 1}
    signature = Counter()
    while layer:
        advance = defaultdict(int)
        for (seen, counts, end), ways in layer.items():
            if seen == FULL and all(counts):
                signature[counts] += ways
            available = (rows[end] if end >= 0 else 255) & ~seen
            while available:
                bit = available & -available
                available ^= bit
                vertex = bit.bit_length() - 1
                if vertex in PORTS:
                    i = PORTS.index(vertex)
                    new = list(counts)
                    new[i] += 1
                    advance[seen, tuple(new), vertex] += ways
                else:
                    advance[seen | bit, counts, vertex] += ways
        layer = advance
    return signature


def difference(signature):
    """Collect coefficients of base**t for H(mask) minus the lower rail."""
    coeff = Counter()
    for runs, ways in signature.items():
        for slots in product(*(range(1, r + 1) for r in runs)):
            i, j, k = slots
            sign_and_choices = prod((-1) ** (r - s) * comb(r, s)
                                    for r, s in zip(runs, slots))
            coeff[i * i * j * k] += ways * i**3 * sign_and_choices
    coeff.subtract(LOWER)
    return tuple(sorted((base, c) for base, c in coeff.items() if c))


def evaluate(coeff, t):
    return sum(c * b**t for b, c in coeff)


def subtract(coeff, other):
    shifted = Counter(dict(coeff))
    shifted.subtract(dict(other))
    return tuple(sorted((b, c) for b, c in shifted.items() if c))


def tail_sign(coeff):
    if not coeff:
        return 0, 1
    last, lead = coeff[-1]
    t = 1
    while abs(lead) * last**t <= sum(abs(c) * b**t for b, c in coeff[:-1]):
        t += 1
    return (1 if lead > 0 else -1), t


def all_exponents(readouts, target):
    identities = []
    exceptions = defaultdict(list)
    max_cutoff = 0
    for coeff, masks in readouts.items():
        terms = subtract(coeff, ((1, target),))
        if not terms:
            identities.extend(masks)
            continue
        _, cutoff = tail_sign(terms)
        for t in range(1, cutoff):
            if evaluate(terms, t) == 0:
                exceptions[t].extend(masks)
        max_cutoff = max(max_cutoff, cutoff)
    return sorted(identities), {t: sorted(ms) for t, ms in exceptions.items()}, max_cutoff


def width_inverse(signatures, target):
    # Any wider live slot with two or three runs would already exceed target.
    limit = (target + 2).bit_length() - 1
    hits = []
    for signature in signatures:
        assert all(i <= 2 and j <= 3 and k <= 2 for (i, j, k), _ in signature)
        for b in range(1, limit + 1):
            fill = (0, 1, 2**b - 2, 3**b - 3 * 2**b + 3)
            coeff = [[0, 0], [0, 0]]
            for (i, j, k), ways in signature:
                coeff[i - 1][k - 1] += ways * fill[j]
            (base, third), (first, mixed) = coeff
            if base > target:
                break
            for c in range(1, limit + 1):
                y = 2**c - 2
                remainder = target - base - third * y
                if remainder < 0:
                    break
                slope = first + mixed * y
                if slope:
                    q, r = divmod(remainder, slope)
                    power = q + 2
                    if r == 0 and power & (power - 1) == 0:
                        hits.append((signature, (power.bit_length() - 1, b, c)))
                elif remainder == 0:
                    hits.append((signature, (1, b, c)))
    return hits


def check():
    readouts = defaultdict(list)
    signatures = Counter()
    for mask in range(1 << len(EDGES)):
        signature = compressed_paths(source(mask))
        signatures[tuple(sorted(signature.items()))] += 1
        readouts[difference(signature)].append(mask)
    assert Counter(signatures.values()) == {2: 499, 4: 5, 6: 1}
    assert Counter(map(len, readouts.values())) == {2: 471, 4: 19, 6: 1}
    assert len(readouts) == 491
    assert readouts[()] == [0, 1]
    assert compressed_paths(source(0)) == compressed_paths(source(1))
    results = {target: all_exponents(readouts, target) for target in (0, 2)}
    assert results[0] == ([0, 1], {1: [3, 16, 148, 149, 162, 163]}, 4)
    assert results[2] == ([], {}, 4)
    print("All-exponent solutions for changes 0 and 2:", results)
    lo = ((1, 2), (2, -1), (4, 16), (8, -8))
    hi = ((1, 2), (4, 16))
    assert readouts[lo] == [3, 16] and readouts[hi] == [18, 19]
    for coeff in readouts:
        sign, start = tail_sign(coeff)
        above_hi, start_hi = tail_sign(subtract(coeff, hi))
        above_lo, start_lo = tail_sign(subtract(coeff, lo))
        assert max(start, start_hi, start_lo) <= 4
        assert sign <= 0 or above_hi >= 0
        assert sign >= 0 or above_lo <= 0
    for t in (1, 2, 3):
        values = [evaluate(c, t) for c in readouts]
        adjacent = (max(d for d in values if d < 0), min(d for d in values if d > 0))
        assert adjacent == {1: (-16, 18), 2: (-258, 252),
                            3: (evaluate(lo, t), evaluate(hi, t))}[t]
    for a, pair in enumerate(EDGES):
        coeff = difference(compressed_paths(source(1 << a)))
        print("single flip", pair, coeff, "at 61:", evaluate(coeff, 29))
    at_61 = {evaluate(c, 29): (c, masks) for c, masks in readouts.items()}
    assert len(at_61) == len(readouts)
    value = min(key for key in at_61 if key > 0)
    assert value == 2 + 2**62
    assert at_61[value] == (((1, 2), (4, 16)), [18, 19])
    print("Smallest strictly positive change at 61:", value, "masks:", at_61[value][1])
    perfect = 2**60 * (2**61 - 1)
    assert width_inverse(signatures, perfect + 1) == []
    hits = width_inverse(signatures, perfect - 1)
    expected = tuple(sorted(compressed_paths(source(0)).items()))
    assert hits == [(expected, (61, 29, 29))]
    assert signatures[expected] == 2
    print("All independent positive widths at 61: no upper hit; only the original lower profile.")


if __name__ == "__main__":
    check()
