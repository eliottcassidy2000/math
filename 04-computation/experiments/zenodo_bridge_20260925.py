"""Two-channel arithmetic sewing and a terminating affine-word decoder.

Pure exact controls for the accompanying result note. No external data/code
are needed to reproduce the arithmetic experiment.
"""

from fractions import Fraction
from itertools import product
from pathlib import Path
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def encode(word, sigma=1):
    a, b, d = 1, 0, 1
    for digit in word:
        if digit:
            a, b = 3 * a, 3 * b + sigma * d
        d *= 2
    return a, b, d


def decode(matrix, sigma=1):
    a, b, d = matrix
    if sigma not in (-1, 1) or any(type(x) is not int for x in matrix):
        raise ValueError("integer domain/sign")
    if a <= 0 or a % 2 == 0 or d <= 0 or d & (d - 1):
        raise ValueError("positive odd slope/power-of-two denominator")
    word = []
    while d > 1:
        digit = b % 2
        if digit:
            if a % 3:
                raise ValueError("odd digit needs a factor three")
            a //= 3
            b = (b - sigma * a) // 2
        else:
            b //= 2
        d //= 2
        word.append(digit)
    if (a, b, d) != (1, 0, 1):
        raise ValueError("nonidentity terminal matrix")
    return tuple(word)


def compose(prefix, suffix):
    a, b, d = prefix
    e, f, g = suffix
    return e * a, e * b + f * d, g * d


def step(n, sigma):
    return (3 * n + sigma) // 2 if n % 2 else n // 2


def follow_word(n, word, sigma):
    for digit in word:
        if n % 2 != digit:
            return None
        n = step(n, sigma)
    return n


def boundaries(n, target, prefix, suffix):
    a, b, d = prefix
    e, f, g = suffix
    return Fraction(a * n + b, d), Fraction(g * target - f, e)


def residual(n, target, prefix, suffix):
    a, b, d = prefix
    e, f, g = suffix
    return e * (a * n + b) + d * (f - g * target)


def main():
    decoder_count = guard_count = sewing_count = 0
    for sigma in (1, -1):
        seen = set()
        for length in range(13):
            for word in product((0, 1), repeat=length):
                matrix = encode(word, sigma)
                check(decode(matrix, sigma) == word, "decode/encode")
                check(matrix not in seen, "word matrix collision")
                seen.add(matrix)
                decoder_count += 1
                if length <= 8:
                    a, b, d = matrix
                    for n in range(1, 65):
                        direct = follow_word(n, word, sigma)
                        candidate = Fraction(a * n + b, d)
                        check((candidate.denominator == 1) == (direct is not None),
                              "endpoint versus intermediate parity")
                        if direct is not None:
                            check(candidate == direct, "word endpoint")
                        guard_count += 1
                if 2 <= length <= 8:
                    for split in range(1, length):
                        prefix = encode(word[:split], sigma)
                        suffix = encode(word[split:], sigma)
                        check(compose(prefix, suffix) == matrix, "sewing composition")
                        for n in (1, 9, 17):
                            target = n
                            for _ in word:
                                target = step(target, sigma)
                            for r in (target, target + 1):
                                left, right = boundaries(n, r, prefix, suffix)
                                zero = residual(n, r, prefix, suffix) == 0
                                check(zero == (left == right), "boundary equality")
                                check(zero == (follow_word(n, word, sigma) == r), "root guard")
                                if zero:
                                    check(left.denominator == right.denominator == 1,
                                          "dyadic/triadic intersection")
                                sewing_count += 1

    rank_controls = []
    for sigma, h, k, n, r in product((1, -1), range(1, 7), range(1, 7), (1, 9), (4, 7)):
        prefixes = [encode((0,) * h, sigma), encode((1,) + (0,) * (h - 1), sigma)]
        suffixes = [encode((0,) * k, sigma), encode((1,) + (0,) * (k - 1), sigma)]
        m = [[residual(n, r, u, v) for v in suffixes] for u in prefixes]
        determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0]
        expected = -(1 << h) * (2 * n + sigma) * (sigma + (1 << (k + 1)) * r)
        check(determinant == expected and determinant != 0, "exact rank-two minor")
        rank_controls.append(determinant)

    rejected = []
    for matrix in ((2, 0, 2), (3, 0, 2), (1, 1, 2), (1, 0, 3), (1, 2, 1)):
        try:
            decode(matrix)
        except ValueError:
            rejected.append(matrix)
        else:
            raise RuntimeError("invalid matrix accepted")

    n, word, states = 9, [], [9]
    while n != 4:
        word.append(n % 2)
        n = step(n, 1)
        states.append(n)
    prefix, suffix = encode(word[:2]), encode(word[2:])
    left, right = boundaries(9, 4, prefix, suffix)
    check(prefix == (3, 1, 4) and suffix == (243, 347, 512), "nine certificate")
    check(left == right == 7, "nine seam")
    check(encode((0, 1)) == (3, 2, 4) and encode((1, 0)) == (3, 1, 4), "order hostile")

    data = {
        "status": "FINITE-EXACT controls; all-depth statements proved in note",
        "decoder_universe": "both signs, every binary word of length0..12",
        "decoder_checks": decoder_count,
        "guard_universe": "both signs, lengths0..8, n1..64",
        "guard_checks": guard_count,
        "sewing_universe": "both signs, lengths2..8, every split, n1/9/17, actual target and target+1",
        "sewing_checks": sewing_count,
        "rank_two_nonzero_minors": len(rank_controls),
        "invalid_matrices_rejected": rejected,
        "nine_certificate": {"states": states, "word": word, "prefix": prefix,
                             "suffix": suffix, "seam": str(left), "endpoint": 4},
        "ordered_carry_hostile": {"01": encode((0, 1)), "10": encode((1, 0)),
                                  "n9_endpoints": [str(Fraction(29, 4)), "7"]},
    }
    out = Path(__file__).resolve().parents[2] / "05-knowledge/results/zenodo_bridge_20260925.out"
    out.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("PASS", decoder_count, "word decodes;", guard_count, "guards;", sewing_count,
          "sewing cases;", len(rank_controls), "rank-two minors; five rejection controls")


if __name__ == "__main__":
    main()
